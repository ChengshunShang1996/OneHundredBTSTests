from cmath import pi
import KratosMultiphysics
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import KratosMultiphysics.DEMApplication as DEM
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage
from KratosMultiphysics.DEMApplication import DEM_procedures as DEM_procedures
import math
import datetime

class BrazilianSplitTest(DEMAnalysisStage):

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

        self.parameters = parameters

    def Initialize(self):
        super().Initialize()
        self.InitializeMaterialTest()
        self.PrepareDataForGraph()    
        self._GetSolver().cplusplus_strategy.HealAllBonds()
        ParallelBondUtilities().SetCurrentIndentationAsAReferenceInParallelBondsForPBM(self.spheres_model_part)
        self.spheres_model_part.ProcessInfo[DELTA_TIME] = self.parameters["MaxTimeStep"].GetDouble()
        self.dt = self.spheres_model_part.ProcessInfo[DELTA_TIME]
        self.end_sim = 0
             
    def RunSolutionLoop(self):

        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        self.PrintGraph(self.time)
        self.CheckSimulationEnd()
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.MeasureForcesAndPressure()

    def Finalize(self):
        super().Finalize()
        # TODO: After self.CleanUpOperations() in base class!!
        self.FinalizeGraphs()

    def KeepAdvancingSolutionLoop(self):
        
        return (self.time < self.end_time and self.end_sim < 1)

    def InitializeMaterialTest(self):

        self.top_mesh_nodes = [] 
        self.bot_mesh_nodes = []

        self.graph_counter = 0  
        self.length_correction_factor = 1.0
        self.graph_frequency = int(self.parameters["GraphExportFreq"].GetDouble()/self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME))
        self.strain = 0.0
        self.total_stress_top = 0.0
        self.total_stress_bot = 0.0
        self.total_stress_mean = 0.0
        self.total_tensile_stress = 0.0
        self.LoadingVelocity = 0.0
        self.MeasuringSurface = 0.0
        self.total_tensile_stress_list = []
        self.total_tensile_stress_max = 0.0
        self.total_tensile_stress_max_time = 0.0

        if "material_test_settings" in self.parameters.keys():
            self.thickness = self.parameters["material_test_settings"]["SpecimenThickness"].GetDouble()
            self.diameter = self.parameters["material_test_settings"]["SpecimenDiameter"].GetDouble()
            self.test_type = self.parameters["material_test_settings"]["TestType"].GetString()
            self.width = self.parameters["material_test_settings"]["WallWidth"].GetDouble()
            self.length = self.parameters["material_test_settings"]["WallLength"].GetDouble()
        else:
            self.thickness = self.parameters["SpecimenThickness"].GetDouble()
            self.diameter = self.parameters["SpecimenDiameter"].GetDouble()
            self.test_type = self.parameters["TestType"].GetString()
            self.width = self.parameters["material_test_settings"]["WallWidth"].GetDouble()
            self.length = self.parameters["material_test_settings"]["WallLength"].GetDouble()

        self.ComputeLoadingVelocity()
        self.ComputeMeasuringSurface()
        self.problem_name = self.parameters["problem_name"].GetString()
        self.initial_time = datetime.datetime.now()
        absolute_path_to_file = os.path.join(self.graphs_path, self.problem_name + "_Parameter_chart.grf")
        self.chart = open(absolute_path_to_file, 'w')
        self.aux = AuxiliaryUtilities()
        self.PreUtilities = PreUtilities()
        self.PrepareTests()
        domain_volume = math.pi * 0.5 * 0.5 * self.diameter * self.diameter * self.thickness
        DEM_procedures.GranulometryUtils(domain_volume, self.spheres_model_part)

    def ComputeLoadingVelocity(self):
        top_vel = bot_vel = 0.0
        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                top_vel = smp[LINEAR_VELOCITY_Y]
            if smp[IDENTIFIER] == 'BOTTOM':
                bot_vel = smp[LINEAR_VELOCITY_Y]
        self.LoadingVelocity = top_vel - bot_vel

    def MeasureForcesAndPressure(self):

        dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.strain += -100 * self.length_correction_factor * self.LoadingVelocity * dt / self.diameter

        total_force_top = 0.0
        for node in self.top_mesh_nodes:
            force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_top += force_node_y
        self.total_stress_top = total_force_top / self.MeasuringSurface

        total_force_bot = 0.0
        for node in self.bot_mesh_nodes:
            force_node_y = -node.GetSolutionStepValue(ELASTIC_FORCES)[1]
            total_force_bot += force_node_y
        self.total_stress_bot = total_force_bot / self.MeasuringSurface
        
        self.total_stress_mean = 0.5 * (self.total_stress_bot + self.total_stress_top)

        self.total_tensile_stress = 2.0 * 0.5 * (total_force_top + total_force_bot) / (math.pi * self.diameter * self.thickness)

        self.total_tensile_stress_list.append(self.total_tensile_stress)

        if self.total_tensile_stress_max < self.total_tensile_stress:
            self.total_tensile_stress_max = self.total_tensile_stress
            self.total_tensile_stress_max_time = self.time

    def CheckSimulationEnd(self):

        if self.total_tensile_stress_max_time < 0.5 * self.time:
            self.end_sim = 2   # means end the simulation

    def ComputeMeasuringSurface(self):
        self.MeasuringSurface = self.width * self.length

    def PrepareTests(self):

        absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph_mean.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_graph_top.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_bot.grf")
        absolute_path_to_file4 = os.path.join(self.graphs_path, self.problem_name + "_graph_tensile.grf")
        self.graph_export_1 = open(absolute_path_to_file1, 'w')
        self.graph_export_2 = open(absolute_path_to_file2, 'w')
        self.graph_export_3 = open(absolute_path_to_file3, 'w')
        self.graph_export_4 = open(absolute_path_to_file4, 'w')

    def PrepareDataForGraph(self):

        prepare_check = [0,0,0,0]
        self.total_check = 0

        for smp in self.rigid_face_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                self.top_mesh_nodes = smp.Nodes
                prepare_check[0] = 1
            if smp[IDENTIFIER] == 'BOTTOM':
                self.bot_mesh_nodes = smp.Nodes
                prepare_check[1] = 1

        for smp in self.spheres_model_part.SubModelParts:
            if smp[IDENTIFIER] == 'TOP':
                self.top_mesh_nodes = smp.Nodes
                prepare_check[2] = -1

            if smp[IDENTIFIER] == 'BOTTOM':
                self.bot_mesh_nodes = smp.Nodes
                prepare_check[3] = -1

        for it in range(len(prepare_check)):
            self.total_check += prepare_check[it]

        if math.fabs(self.total_check) != 2:
            self.Procedures.KratosPrintWarning(" ERROR in the definition of TOP BOT groups. Both groups are required to be defined, they have to be either on FEM groups or in DEM groups")

    def PrintGraph(self, time):

        if self.graph_counter == self.graph_frequency:
            self.graph_counter = 0
            self.graph_export_1.write(str("%.6g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_mean)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_2.write(str("%.8g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_top)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_3.write(str("%.8g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_stress_bot)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_4.write(str("%.8g"%self.strain).rjust(13) + "  " + str("%.6g"%(self.total_tensile_stress)).rjust(13) + "  " + str("%.8g"%time).rjust(12) + '\n')
            self.graph_export_1.flush()
            self.graph_export_2.flush()
            self.graph_export_3.flush()
            self.graph_export_4.flush()
        self.graph_counter += 1
    
    def FinalizeGraphs(self):
        # Create a copy and renaming
        absolute_path_to_file1 = os.path.join(self.graphs_path, self.problem_name + "_graph.grf")
        absolute_path_to_file2 = os.path.join(self.graphs_path, self.problem_name + "_bts.grf")
        absolute_path_to_file3 = os.path.join(self.graphs_path, self.problem_name + "_graph_VOL.grf")
        for filename in os.listdir("."):
            if filename.startswith(absolute_path_to_file1):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", absolute_path_to_file1 + str(self.initial_time).replace(":", "") + ".grf")
            if filename.startswith(absolute_path_to_file2):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", absolute_path_to_file2 + str(self.initial_time).replace(":", "") + ".grf")
            if filename.startswith(absolute_path_to_file3):
                shutil.copy(filename, filename + "COPY")
                os.rename(filename+"COPY", absolute_path_to_file3 + str(self.initial_time).replace(":", "") + ".grf")
        self.graph_export_1.close()
        self.graph_export_2.close()
        self.graph_export_3.close()
        self.graph_export_4.close()

if __name__ == "__main__":

    with open("ProjectParametersDEM.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    BrazilianSplitTest(model, parameters).Run()
