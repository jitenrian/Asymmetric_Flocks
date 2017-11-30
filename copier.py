#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu  Nov 21:55:15 2017

@author: jitendrian
"""

import os
import shutil

parent_dir = os.getcwd()
parent_lis = os.listdir(parent_dir)
naam = parent_dir.split('/')
os.mkdir(str(naam[-1]))
os.mkdir(str(naam[-1]) + "_all")
destination0 = os.path.join(parent_dir, naam[-1])
destination1 = os.path.join(parent_dir, naam[-1] + "_all")
for folders in parent_lis:
    copying_from = []
    copy_to = []
    file_name = []
    if folders.startswith('Alp'):
        child1_dir = os.path.join(parent_dir, folders)
        os.chdir(child1_dir)
        child1_lis = os.listdir(child1_dir)
        for files in child1_lis:
            if files.endswith('.png'):
                copy_file = os.path.join(child1_dir, files)
                copying_from.append(copy_file)
                file_name.append(files)
            if files.endswith('.mp4'):
                copy_file = os.path.join(child1_dir, files)
                copying_from.append(copy_file)
                file_name.append(files)
        
        os.chdir(destination0)
        os.mkdir(folders)
        cousin1_dir = os.path.join(destination0, folders)
        for i in range(len(file_name)):
            final_copy_dir = os.path.join(cousin1_dir, file_name[i])
            shutil.copy(copying_from[i], final_copy_dir)
        
        moved_to = os.path.join(destination1, folders)
        shutil.move(child1_dir, moved_to)
        
