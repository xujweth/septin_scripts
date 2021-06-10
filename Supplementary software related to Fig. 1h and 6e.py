#!/bin/env/python
#This version is orginally used to analyze the distance between bacterial membrane and the septin filaments
#This script is written by Jingwei Xu, ETH Zurich

import math, time, statistics
import numpy as np
from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	usage="""python model_dist.py [model1] [model2]"""
	parser=OptionParser(usage=usage)
	parser.add_option("--interpolate", dest="interpolate", action="store_true", help="To interpolate the curve based on current points. Only available now in U shape membrane.", default=False)
	parser.add_option("--resampling", dest="resampling", type=int, help="The resampling number in interpolation of membrane curve. The default is 1.", default=1)
	parser.add_option("--by_x_max", dest="by_x_max", action="store_true", help="To seperate the memrbane curve based on the X maximum value instead of X minimum value. The default is False.", default=False)
	parser.add_option("--seperate", dest="seperate", action="store_true", help="To seperate the membrane curve into two parts for interpolation. Just in case when the long axis of membrane is not in the y axis, which will result one X value has two Y values. The default is False.", default=False)
	parser.add_option("--write_tmp", dest="write_tmp", action="store_true", help="To write the interpolate file? The default is False.", default=False)
	parser.add_option("--restrictZ", dest="restrictZ", action="store_true", help="To restrict the calculation only in the cell membrane Z sections. The default is false.", default=False)
	(options, args)=parser.parse_args()


	file1=args[0]
	file2=args[1]
	membrane_file=file1
	septin_file=file2

#	To define the files for processing
	if (len(open(file1,"r").readline().split())) > 3:
		membrane_file=file2
		septin_file=file1
#	print membrane_file, septin_file

	membrane_xyz=open(membrane_file,"r").readlines()
	septin_xyz=open(septin_file,"r").readlines()

	min_z, max_z = z_range(membrane_xyz)
	print "The z range of model is from %.1f to %.1f."%(min_z, max_z)

	membrane_xyz_mod = ""
	if options.interpolate:
		membrane_xyz_mod=reorganize_lst(membrane_xyz, options.resampling, options.write_tmp, options.seperate, options.by_x_max)
	else:
		membrane_xyz_mod=membrane_xyz
	print "There are total %d point in membrane model file."%(len(membrane_xyz_mod))

	out_file = "%s_mod.txt"%(septin_file.split('.')[0])
	out = open(out_file, "w")
	combine_lst = []
	temp_lst = []
	fil_num_cur = 0
	for i in septin_xyz:
		fil_num, xsep, ysep, zsep = i.split()[:]
		dist = []
		if options.restrictZ:
			if (float(zsep) >= min_z) and (float(zsep) <= max_z):
				for j in membrane_xyz_mod:
					xm, ym, zm = j.split()[:]
					dist_temp=math.sqrt((float(xsep)-float(xm))**2+(float(ysep)-float(ym))**2+(float(zsep)-float(zm))**2)
					dist.append(dist_temp)
		else:
			for j in membrane_xyz_mod:
				xm, ym, zm = j.split()[:]
				dist_temp=math.sqrt((float(xsep)-float(xm))**2+(float(ysep)-float(ym))**2+(float(zsep)-float(zm))**2)
				dist.append(dist_temp)
		mins = np.nanmin(dist)
		out_info = "%.2f\t%d\t%.1f\t%.1f\t%.1f\n"%(mins,float(i.split()[0]),float(i.split()[1]),float(i.split()[2]),float(i.split()[3]))
		out.write(out_info)
		if fil_num != fil_num_cur:
			if len(temp_lst) != 0:
				if len(temp_lst) != 0:
					combine_lst.append(temp_lst)
					temp_lst = []
				temp_lst.append([mins])
				fil_num_cur = fil_num
			else:
				temp_lst.append([mins])

#	for i in combine_lst:
#		avg = np.mean(i)
#		std = np.std(i)
#		print avg, std
						
	out.close()
	
	
def z_range(xyz_file):
	z_lst = []
	for i in xyz_file:
		z = float(i.split()[-1])
		z_lst.append(z)	
	z_max = max(z_lst)
	z_min = min(z_lst)
	return z_min, z_max

def reorganize_lst(xyz_lst, resampling, write_tmp, seperate, by_x_max):
	total_lst = []
	z_value = 0
	temp_lst = []
	for i in xyz_lst:
		x = float(i.split()[0])
		y = float(i.split()[1])
		z = float(i.split()[2])
		if z != z_value:
			if len(temp_lst) != 0:
				total_lst.append(temp_lst)
				temp_lst = []
			temp_lst.append([x,y,z])
			z_value = z
		else:
			temp_lst.append([x,y,z])
	total_lst.append(temp_lst)

	total_interp_lst = []
	for i in total_lst:
		interp_lst = interpolation(i, resampling, seperate, by_x_max)
		total_interp_lst.append(interp_lst)
	
	out = ""
	if write_tmp:
		out = open('interp_tmp.txt', "w")
	out_format_lst = []
	for i in total_interp_lst:
		for j in i:
			x_interp = j[0]
			y_interp = j[1]
			z_interp = j[2]
			interp_xyz = "%.1f\t%.1f\t%.1f\n"%(x_interp, y_interp, z_interp)
			out_format_lst.append(interp_xyz)
			if write_tmp:
				out.write(interp_xyz)
	if write_tmp:
		out.close()
	return out_format_lst

def interpolation(lst, resampling, seperate, by_x_max):
	x_lst = []
	y_lst = []
	z_value = 0
	total_lst = []
	for point in lst:
		x_lst.append(point[0])
		y_lst.append(point[1])
		z_value = point[2]
	x_np = np.array(x_lst, dtype=np.float32)
	y_np = np.array(y_lst, dtype=np.float32)
	x_min_index = np.where(x_np == min(x_np))[0][0] + 1
	x_max_index = np.where(x_np == max(x_np))[0][0] + 1

	if seperate:
		x_index = x_min_index	
		if by_x_max:
			x_index = x_max_index
		
		x_part1 = x_np[:x_index]
		y_part1 = y_np[:x_index]
		f1 = interp1d(x_part1,y_part1)
		x_new_part1 = np.linspace(min(x_part1),max(x_part1),num=(max(x_part1)-min(x_part1))*resampling,endpoint=True)
		y_new_part1 = f1(x_new_part1)
	
		x_part2 = x_np[x_index:]
		y_part2 = y_np[x_index:]
		f2 = interp1d(x_part2,y_part2)
		x_new_part2 = np.linspace(min(x_part2),max(x_part2),num=(max(x_part2)-min(x_part2))*resampling,endpoint=True)
		y_new_part2 = f2(x_new_part2)
	
		x_comb = np.concatenate((x_new_part1,x_new_part2))
		y_comb = np.concatenate((y_new_part1,y_new_part2))
		z_comb = np.full((len(x_comb)),z_value)

		for num in range(len(x_comb)):
			interp_lst = [x_comb[num], y_comb[num], z_comb[num]]
			total_lst.append(interp_lst)
	else:
		f = interp1d(x_np, y_np)
		x_comb = np.linspace(min(x_np),max(x_np),num=(max(x_np)-min(x_np))*resampling,endpoint=True)
		y_comb = f(x_comb)
		z_comb = np.full((len(x_comb)),z_value)

		for num in range(len(x_comb)):
			interp_lst = [x_comb[num], y_comb[num], z_comb[num]]
			total_lst.append(interp_lst)

	return total_lst
	

if __name__ == "__main__":
	main()	
