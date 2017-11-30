import os
import shutil

parent_dir = os.getcwd()
parent_lis = os.listdir(parent_dir)

child1_dir = []
for files in parent_lis:
    if files.startswith('Alp'):
        child1_dir.append(os.path.join(parent_dir, files))

title_list = ['Veloc_VS_t.png', 'Thet_VS_t.png', 'XVelo_VS_t.png', 'YVelo_VS_t.png', 'XVel_Var_VS_t.png', 'YVel_Var_VS_t.png', 'Spin_VS_t.png', 'Spin_Var_VS_t.png']
for folders_loc in child1_dir:
    os.chdir(folders_loc)
    child_lis = os.listdir(folders_loc)
    name_arr = folders_loc.split('/')
    name_video = str(name_arr[-1]) + ".mp4"
    name_text = str(name_arr[-1]) + ".txt"
    for files in child_lis:
        if files.endswith('.mp4'):
            files_path = os.path.join(folders_loc, files)
            updated_files_path = os.path.join(folders_loc, name_video)
            os.rename(files_path, updated_files_path)
        if files.endswith('.txt'):
            files_path = os.path.join(folders_loc, files)
            updated_files_path = os.path.join(folders_loc, name_text)
            os.rename(files_path, updated_files_path)

        if files.endswith('.png'):
            for i in range(2,10):
                name_image = title_list[i-2]
                if files.startswith(str(i)):
                    files_path = os.path.join(folders_loc, files)
                    updated_files_path = os.path.join(folders_loc, name_image)
                    os.rename(files_path, updated_files_path)


