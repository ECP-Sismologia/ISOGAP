#                           ISOGAP 0.1
#This script calculates the primary azimuthal gap for a certain network geometry
#Input:
#- A csv file with station name, longitude and latitude columns
#- Corners of a rectangle to create the grided area
#- Grid step
## Author: Nelson David Perez e-mail:ndperezg@gmail.com 

import os
import pandas as pd
import numpy as np
import nvector as nv
import sys
from itertools import combinations
import concurrent.futures as cf

#Calculates gaps for a list of azimuths 
def gaps(Azz):
	azz = sorted(Azz)
	print(azz)
	gaps_ = []
	for i in range(len(azz)):
		if i != 0:
			alpha = azz[i]-azz[i-1]
		else:
			alpha = azz[0]+ 360 - azz[-1]
		gaps_.append(alpha)
	return gaps_

#Calculates azimuths between two points
def azimuth(Lon1,Lat1,Lon2,Lat2):
	wgs84 = nv.FrameE(name='WGS84')
	pointA = wgs84.GeoPoint(latitude=Lat1, longitude=Lon1, z=0, degrees=True)
	pointB = wgs84.GeoPoint(latitude=Lat2, longitude=Lon2, z=0, degrees=True)
	p_AB_N = pointA.delta_to(pointB)
	print(p_AB_N.azimuth_deg )
	azim = p_AB_N.azimuth_deg
	if azim < 0:
		azim += 360
	else:
		pass
	return azim

#calculates isogap for each point
def each_gap(lon,lat,net):
	azz=[]
	for sta in net:
		print(lon,lat,net[sta][0], net[sta][1])
		azim = azimuth(lon,lat,net[sta][0], net[sta][1])
		azz.append(azim)
	GAP = max(gaps(azz))
	return GAP

def iter_read_stations(arc,sta_number=[]):
	"""
	Parameters:
	-----------
	arc : str
		csv path
	sta_number: list of int
		Each element  of the list is the number of available stations 
	out_folder = str
		Path of the folder to export the results

	Results:
	--------
		dict:
			key : int (number of stations)
			value: list (list of dicts)
				dict: 
					key:str (name of station)
					value: list ([lon,lat])
	"""
	def get_combinations(df,sta_number):
		## intento paralelizacion
		# def _get_combinations(index):
		# 	my_df = df.loc[index,:]
		# 	NET = {}
		# 	for _,row in my_df.iterrows():
		# 		NET[row["STA"]] = [row["LON"],row["LAT"]]
		# 	return NET

		# index_combinations = list(combinations(df.index,sta_number))
		# with cf.ThreadPoolExecutor(max_workers=10) as executor:
		# 	combinations_net = list(executor.map(_get_combinations,
		# 							index_combinations))
		########

		combinations_net = []
		for index in list(combinations(df.index,sta_number)):
			my_df = df.loc[index,:]
			NET = {}
			for _,row in my_df.iterrows():
				NET[row["STA"]] = [row["LON"],row["LAT"]]

			combinations_net.append(NET )
		return combinations_net


	df = pd.read_csv(arc)
	if sta_number:
		all_df = {}
		for number in sta_number:
			combinations_df = get_combinations(df,number)
			all_df[number] = combinations_df
	else:
		l = len(df)
		all_df = {l:get_combinations(df,len(df))}
	
	print(all_df)
	return all_df

def export_gap(NET, minlon, maxlon, minlat, maxlat, 
				step,out_path):
	if os.path.isdir(os.path.dirname(out_path)) == False:
		os.makedirs(os.path.dirname(out_path))

	Lons = np.arange(minlon,maxlon,step)
	Lats = np.arange(minlat,maxlat,step)

	out = open(out_path, 'w')
	out.write('LON,LAT,GAP\n')

	for i in Lons:
		for j in Lats:
			az_gap = each_gap(i,j,NET)
			print(i,j,az_gap)
			out.write('%s,%s,%4.2f\n'%(round(i,5),round(j,5),az_gap))

	out.close()

def run_itergap(sta_csv,out_folder,grid,step,
				sta_number=[]):
	"""
	out_folder = str
		Path of the folder to export the results
	"""
	minlon, maxlon,	minlat, maxlat = grid


	all_df = iter_read_stations(arc=sta_csv,sta_number=sta_number)

	for number_sta,combinations in all_df.items():
		path = os.path.join(out_folder,str(number_sta)+"_stations")
		for i,network in enumerate(combinations,1):
			out_path = os.path.join(path,f"net_{i}_S{step}.csv")
			export_gap(network, minlon, maxlon,
						minlat, maxlat, step,
						out_path)

#Ask for longitudes and latitudes for the study area
def input_area():
	lons= input("Enter min and max longitudes separated by a comma:\n")
	lats= input("Enter min and max latitudes separated by a comma:\n")
	if len(lons.split(','))!=2 or len(lats.split(','))!=2:
		print("Bad input, try again\n")
		sys.exit()
	minlon = float(lons.split(',')[0])
	maxlon = float(lons.split(',')[1])
	minlat = float(lats.split(',')[0])
	maxlat = float(lats.split(',')[1])
	if (minlon>=maxlon) or (minlat>=maxlat):
		print("Wrong values, try again\n")
		sys.exit()
	return minlon, maxlon, minlat, maxlat

if __name__ == "__main__":
	step = 0.005
	grid = (-73.6,-73.2,3.8,4.3)
	sta_csv = "/home/emmanuel/Ecopetrol/ISOGAP/data/stations.csv"
	out_folder = "/home/emmanuel/Ecopetrol/ISOGAP/outs"
	run_itergap(sta_csv,out_folder,grid,
				step,sta_number=[6])

