import os

sigma_limit_list = [1,2]
tension_limit_list = [3,5]

for sigma in sigma_limit_list:
    for tension in tension_limit_list:
        new_folder_name = 'BTS_Sigma' + str(sigma) + '_Tension' + str(tension)
        #aim_folder_path = os.path.join(os.getcwd(),Generated_BTS_cases,new_folder_name)
        aim_path = os.path.join(os.getcwd(),'Generated_BTS_cases', new_folder_name)
        os.mkdir(aim_path)
        