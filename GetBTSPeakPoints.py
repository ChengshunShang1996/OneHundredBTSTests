import os

sigma_limit_list = [10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000]
tension_limit_list = [10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000]

# creat the BTS_peak_points.dat
BTS_peak_points_path_and_name = os.path.join(os.getcwd(),'BTS_peak_points.dat')
with open(BTS_peak_points_path_and_name, "w") as f_w_peak_points:

    for sigma_limit in sigma_limit_list:

        for tension_limit in tension_limit_list:

            aim_folder_name = 'BTS_Sigma' + str(sigma_limit) + '_Tension' + str(tension_limit)
            aim_path_and_name = os.path.join(os.getcwd(),'Generated_BTS_cases', aim_folder_name, 'BTStest_Graphs', 'BTStest_graph_tensile.grf')

            tensile_data_list = []
            with open(aim_path_and_name, 'r') as tensile_data:
                for line in tensile_data:
                    values = [float(s) for s in line.split()]
                    tensile_data_list.append(values[1]) 

            # write BTS_peak_points.dat
            f_w_peak_points.write(str(sigma_limit) + ' ' + str(tension_limit) + ' ' + str(max(tensile_data_list)) + '\n')
