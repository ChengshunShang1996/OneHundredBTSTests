import os
import shutil
#import numpy as np

#sigma_limit_list = np.arange(100, 100000, 0.5).tolist()
#tension_limit_list = np.arange(100, 100000, 0.5).tolist()
#sigma_limit_list = [10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000]
#tension_limit_list = [10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000]
sigma_limit_list = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000]
shear_limit_list = [10000, 50000, 100000, 200000, 300000, 400000, 500000, 600000]
#sigma_limit_list = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900]
#shear_limit_list = [512000, 514000, 516000, 518000, 520000, 522000, 524000, 526000, 528000, 530000]

# creat the cases_run.sh
cases_run_path_and_name = os.path.join(os.getcwd(),'cases_run.sh')
with open(cases_run_path_and_name, "w") as f_w_cases_run:
    f_w_cases_run.write('#!/bin/bash'+'\n')

    for sigma_limit in sigma_limit_list:

        for shear_limit in shear_limit_list:

            #creat new folder
            new_folder_name = 'BTS_Sigma' + str(sigma_limit) + '_Shear' + str(shear_limit)
            aim_path = os.path.join(os.getcwd(),'Generated_BTS_cases', new_folder_name)
            if os.path.exists(aim_path):
                shutil.rmtree(aim_path)
            os.mkdir(aim_path)

            #copy source file
            seed_file_name_list = ['BTS_PBM_220903.py', 'BTStestDEM_FEM_boundary.mdpa', 'BTStestDEM.mdpa', 'ProjectParametersDEM.json', 'MaterialsDEM.json', 'run_omp.sh']
            for seed_file_name in seed_file_name_list:
                seed_file_path_and_name = os.path.join(os.getcwd(),'BTS_seed_files',seed_file_name)
                aim_file_path_and_name = os.path.join(aim_path, seed_file_name)

                if seed_file_name == 'MaterialsDEM.json':
                    with open(seed_file_path_and_name, "r") as f_material:
                        with open(aim_file_path_and_name, "w") as f_material_w:
                            for line in f_material.readlines():
                                if "BOND_SIGMA_MAX" in line:
                                    line = line.replace("1e3", str(sigma_limit))
                                if "BOND_TAU_ZERO" in line:
                                    line = line.replace("2.6e6", str(shear_limit))
                                f_material_w.write(line)
                elif seed_file_name == 'run_omp.sh':
                    with open(seed_file_path_and_name, "r") as f_run_omp:
                        with open(aim_file_path_and_name, "w") as f_run_omp_w:
                            for line in f_run_omp.readlines():
                                if "BTS-Q-Ep6.2e10-T1e3-f0.1" in line:
                                    line = line.replace("BTS-Q-Ep6.2e10-T1e3-f0.1", new_folder_name)
                                f_run_omp_w.write(line)
                else:
                    shutil.copyfile(seed_file_path_and_name, aim_file_path_and_name) 

            # write the cases_run.sh
            f_w_cases_run.write('cd '+ aim_path + '\n')
            f_w_cases_run.write('sbatch run_omp.sh' + '\n')


