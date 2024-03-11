# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 18:29:47 2023

@author: KhorsandiKuhanestP
"""

import numpy as np 
import pyvista as pv 
import matplotlib.pyplot as plt 
import netCDF4 
from datetime import datetime, timedelta 
import pandas as pd  
import triangle as tr 
import math
import shapely 
from shapely import Polygon 
from shapely.geometry import Point, Polygon
import matplotlib.tri as mtri 
import math
import pickle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
import time 
t_0 = time.time() 

ncfile_currentmesh =  r'C:\\\\Users\\\\KhorsandiKuhanestP\\\\OneDrive - University of Twente\\\\D_Flow_FM\\\\221117-Question one\\\\COMBINED_V2\\\\VERSION 4\\\\CurrentMesh_Z.nc' 
# Show all variables
ncfid_1 = netCDF4.Dataset(ncfile_currentmesh)
projected_coordinate_system_1=np.array(ncfid_1.variables['projected_coordinate_system'])
mesh2d_1=np.array(ncfid_1.variables['mesh2d'])
mesh2d_node_x_1=np.array(ncfid_1.variables['mesh2d_node_x'])
mesh2d_node_y_1=np.array(ncfid_1.variables['mesh2d_node_y'])
mesh2d_node_z_1=np.array(ncfid_1.variables['mesh2d_node_z'])

x_Cmesh=list()
y_Cmesh=list()
Z_Cmesh=list()
for i in range (0, (len(mesh2d_node_x_1))):
    x_1=mesh2d_node_x_1[i]
    y_1=mesh2d_node_y_1[i]
    z_1=mesh2d_node_z_1[i]
    x_Cmesh.append(x_1)
    y_Cmesh.append(y_1)
    Z_Cmesh.append(z_1)
mesh_points=list(zip(x_Cmesh, y_Cmesh, Z_Cmesh))

#sortting the mesh points based on high x 
mesh_points.sort(key=lambda row: (row[0]))

#finding the number of the cross sections
n=1
i=0
resolution=10
for u in range(i,len(mesh_points)-1):
    Diff_x=float(mesh_points[u+1][0])-float(mesh_points[u][0])
    if Diff_x>float(resolution):
        n=n+1
        i=u

#finding the coss section based on the distance in the x direction
Cross_sections = []
threshold=1
i=0

for f in range(0,n):
    innerlist_Cross_section = []
    innerlist_Cross_section.append(mesh_points[i])
    for u in range(i+1,len(mesh_points)):
        Diff_x=float(mesh_points[u][0])-float(mesh_points[u-1][0])
        if Diff_x<=float(resolution):
            innerlist_Cross_section.append(mesh_points[u])
          
        if Diff_x>float(resolution):
            i=u
            break
    Cross_sections.append(innerlist_Cross_section)

#Sorting each cross section
CROSS_section_sorted=list()
for f in range(0,n):
    Cross_sections_sorted=Cross_sections[f]
    Cross_sections_sorted.sort(key=lambda row: (row[1]), reverse=True)
    CROSS_section_sorted.append(Cross_sections_sorted)

def closest_item(lst, target):
    return min(lst, key=lambda x: abs(x - target))

#mesh resolution input and water level steps
mesh_reso=40
WLR=0.5 #water level resolution

#Main channel information (based on the most down stream cross section)
X_start_mainchannel=1100
X_end_mainchannel=600
Width_X_S =[]
Width_X_E =[]
Width_L_S =[]
Width_L_E =[]
LOWEST_min =[]
Highest_Y=[]
X_CROSS =[]
Y_CROSS =[]
L_CROSS =[]
WATER_LEVEL =[]
WATER_DEPTH =[]
for y in range (len(CROSS_section_sorted)):
    #defining series:
    WATER_level=[] #the related water level, every WLR m
    Width_x_S=[]
    Width_x_E=[]
    Width_l_S=[]
    Width_l_E=[]
    X_cross=list()
    Y_cross=list()
    L_cross=list()
    highest_y=min(CROSS_section_sorted[y][0][2], CROSS_section_sorted[y][len(CROSS_section_sorted[y])-1][2])
    for u in range (0, len(CROSS_section_sorted[y])):
        X_cross.append(CROSS_section_sorted[y][u][1])
        Y_cross.append(CROSS_section_sorted[y][u][2])
        L_cross.append(CROSS_section_sorted[y][u][0])

    X_start_mainchannel=closest_item(X_cross, X_start_mainchannel)
    X_end_mainchannel=closest_item(X_cross, X_end_mainchannel)
    for t in range (len(X_cross)):
        if X_cross[t]==X_start_mainchannel:
            index_start_mainchannel=t
        if X_cross[t]==X_end_mainchannel:
            index_end_mainchannel=t
    lowest_2= min(Y_cross[index_start_mainchannel:index_end_mainchannel]) 
    #print(\"lowest_2:\", lowest_2)
    w=(highest_y-lowest_2)/WLR
    w=int(w)-2
    for i in range (len(X_cross)):
        if Y_cross[i]==lowest_2:
            num_lowest=i
    star_right=num_lowest
    star_left=num_lowest    

    WATER_depth=[]
    for wl in range (w):
        Water_level=float(lowest_2+wl*WLR+3*WLR)
        Water_depth=wl*WLR+3*WLR
        #for the right part of the river
        for i in range(star_right, len(X_cross)):
            if Y_cross[i] >= Water_level:
                break
            if i<len(X_cross)-1:
                j=i+1
            if Y_cross[i] < Water_level <= Y_cross[j]:
                Water_level_x_e=(((X_cross[j]-X_cross[i])/(Y_cross[j]-Y_cross[i]))*(Water_level-Y_cross[j]))+(X_cross[j])
                Water_level_l_e=(((L_cross[j]-L_cross[i])/(Y_cross[j]-Y_cross[i]))*(Water_level-Y_cross[j]))+(L_cross[j])
        #for the left part of the river
        for k in range(0, star_left):
            i=star_left-k-1
            if Y_cross[i+1] > Water_level:
                break
            if i<star_left:
                j=i+1
            if Y_cross[j] < Water_level <= Y_cross[i]:
                Water_level_x_s=(((X_cross[j]-X_cross[i])/(Y_cross[j]-Y_cross[i]))*(Water_level-Y_cross[i]))+(X_cross[i])
                Water_level_l_s=(((L_cross[j]-L_cross[i])/(Y_cross[j]-Y_cross[i]))*(Water_level-Y_cross[i]))+(L_cross[i])
        width=float(Water_level_x_s)-float(Water_level_x_e)
        WATER_level.append(Water_level) 
        WATER_depth.append(Water_depth)
        Width_x_S.append(float(Water_level_x_s))
        Width_x_E.append(float(Water_level_x_e))
        Width_l_S.append(float(Water_level_l_s))
        Width_l_E.append(float(Water_level_l_e))
    Width_X_S.append(Width_x_S)
    Width_X_E.append(Width_x_E)
    Width_L_S.append(Width_l_S)
    Width_L_E.append(Width_l_E)
    LOWEST_min.append(lowest_2)
    Highest_Y.append(highest_y)
    X_CROSS.append(X_cross)
    Y_CROSS.append(Y_cross)
    L_CROSS.append(L_cross)
    WATER_LEVEL.append(WATER_level)
    WATER_DEPTH.append(WATER_depth)

DR=0.1 #water level resolution for calculating the area and wetted perimiter 
X_start_mainchannel=1100 
X_end_mainchannel=600 
Depth_Wetted_Area=[] 
for y in range (0, len(CROSS_section_sorted)): 
    D_AR=[] 
    X_start_mainchannel=closest_item(X_CROSS[y], X_start_mainchannel) 
    X_end_mainchannel=closest_item(X_CROSS[y], X_end_mainchannel) 
    for t in range (len(X_CROSS[y])): 
        if X_CROSS[y][t]==X_start_mainchannel: 
            index_start_mainchannel=t 
        if X_CROSS[y][t]==X_end_mainchannel: 
            index_end_mainchannel=t 
    lowest_2= min(Y_CROSS[y][index_start_mainchannel:index_end_mainchannel]) 
    for i in range (len(X_CROSS[y])): 
        if Y_CROSS[y][i]==lowest_2: 
            num_lowest=i 
    star_right=num_lowest 
    star_left=num_lowest 
    highest_y=min(CROSS_section_sorted[y][0][2], CROSS_section_sorted[y][len(CROSS_section_sorted[y])-1][2])
    D=int((highest_y-lowest_2)/DR) 
    for wl in range (D): 
        Water_depth=wl*DR+DR 
        Water_level=float(lowest_2+Water_depth)
        #print(Water_depth, Water_level)
        #print(star_right, len(X_CROSS[y]))
        #finding the nodes of collision of water level and border  & calculatin the area for each water level 
        area=0 
        wet=0 
        #for the right part of the river 
        for i in range(star_right, len(X_CROSS[y])): 

            if Y_CROSS[y][i] >= Water_level: 
                break 
            if i<len(X_CROSS[y])-1: 
                j=i+1 
            if Y_CROSS[y][i] < Water_level <= Y_CROSS[y][j]: 
                Water_level_x_e=(((X_CROSS[y][j]-X_CROSS[y][i])/(Y_CROSS[y][j]-Y_CROSS[y][i]))*(Water_level-Y_CROSS[y][j]))+(X_CROSS[y][j]) 
                Water_level_l_e=(((L_CROSS[y][j]-L_CROSS[y][i])/(Y_CROSS[y][j]-Y_CROSS[y][i]))*(Water_level-Y_CROSS[y][j]))+(L_CROSS[y][j]) 
                #print(\"yes\", X_CROSS[y][j], Water_level_x_e, X_CROSS[y][i])
                area=area+(((Water_level-Y_CROSS[y][i]))*(X_CROSS[y][i]-Water_level_x_e)/2) 
                wet = wet + math.sqrt((Water_level - Y_CROSS[y][i])**2 + (Water_level_x_e - X_CROSS[y][i])**2) 
            if Y_CROSS[y][j] < Water_level and Y_CROSS[y][i] < Water_level: 
                area= area+(((Water_level-Y_CROSS[y][i])+(Water_level-Y_CROSS[y][j]))*(X_CROSS[y][i]-X_CROSS[y][j])/2) 
                wet=wet+(((Y_CROSS[y][j]-Y_CROSS[y][i])**2+(X_CROSS[y][j]-X_CROSS[y][i])**2)**(1/2)) 
        for k in range(0, star_left): 
            i=star_left-k-1 
            #print(\"star_left-k-1 \", i)
            if Y_CROSS[y][i+1] > Water_level: 
                break 
            if i<star_left: 
                j=i+1 
            if Y_CROSS[y][j] < Water_level <= Y_CROSS[y][i]: 
                Water_level_x_s=(((X_CROSS[y][j]-X_CROSS[y][i])/(Y_CROSS[y][j]-Y_CROSS[y][i]))*(Water_level-Y_CROSS[y][i]))+(X_CROSS[y][i])
                Water_level_l_s=(((L_CROSS[y][j]-L_CROSS[y][i])/(Y_CROSS[y][j]-Y_CROSS[y][i]))*(Water_level-Y_CROSS[y][i]))+(L_CROSS[y][i]) 
                #print(\"yayyy\", i, Water_level_l_s)
                area=area+((Water_level-Y_CROSS[y][j]))*(Water_level_x_s-X_CROSS[y][j])/2 
                wet=wet+(((Water_level-Y_CROSS[y][j])**2+(X_CROSS[y][j]-Water_level_x_s)**2)**(1/2)) 
            if Y_CROSS[y][j] < Water_level and Y_CROSS[y][i] < Water_level: 
                area= area+(((Water_level-Y_CROSS[y][i])+(Water_level-Y_CROSS[y][j]))*(X_CROSS[y][i]-X_CROSS[y][j])/2) 
                wet=wet+(((Y_CROSS[y][j]-Y_CROSS[y][i])**2+(X_CROSS[y][j]-X_CROSS[y][i])**2)**(1/2)) 
        AR=area*((area/wet)**(2/3)) 
        D_AR.append([Water_depth, Water_level, AR, Water_level_x_e,Water_level_l_e, Water_level_x_s, Water_level_l_s]) 
    Depth_Wetted_Area.append(D_AR)

PLYGONS_points=list()

for r in range(0,n-1):
    points_Polygon_x=list()
    points_Polygon_y=list()
    for g in range (len(CROSS_section_sorted[r])):
        points_Polygon_x.append(CROSS_section_sorted[r][g][0])
        points_Polygon_y.append(CROSS_section_sorted[r][g][1])
    for e in range (len(CROSS_section_sorted[r+1])):
        points_Polygon_x.append(CROSS_section_sorted[r+1][len(CROSS_section_sorted[r+1])-1-e][0])
        points_Polygon_y.append(CROSS_section_sorted[r+1][len(CROSS_section_sorted[r+1])-1-e][1])
    points_Polygon_x.append(CROSS_section_sorted[r][0][0])    
    points_Polygon_y.append(CROSS_section_sorted[r][0][1])  

    points_Polygon=list(zip(points_Polygon_x, points_Polygon_y))
    PLYGONS_points.append(points_Polygon)
p1 = Polygon(PLYGONS_points[0])

x,y = p1.exterior.xy
plt.plot(x,y)

p2 = Polygon(PLYGONS_points[1])
x,y = p2.exterior.xy
plt.plot(x,y)

p3 = Polygon(PLYGONS_points[2])
x,y = p3.exterior.xy
plt.plot(x,y)


#reading the bed level file and creatinf surface from it:
    
input_path=r"C:\\Users\\KhorsandiKuhanestP\\OneDrive - University of Twente\\D_Flow_FM\\221117-Question one\\COMBINED_V2\\VERSION 4/"  
output_path = r"C:\\Users\\KhorsandiKuhanestP\\OneDrive - University of Twente\\D_Flow_FM\\221117-Question one\\COMBINED_V2\\VERSION 4/" 
dataname="Bathymetry_Combined_V4.xyz"
point_cloud= np.loadtxt(input_path+dataname,skiprows=1)
BED_points=list()
X_P=list()
Y_P=list()
Z_P=list()
for o in range (len(point_cloud)):
    x=point_cloud[o][0]
    y=point_cloud[o][1]
    z=point_cloud[o][2]
    X_P.append(x)
    Y_P.append(y)
    Z_P.append(z)
    
    
BED_points=list(zip(X_P, Y_P, Z_P))
BED_points.sort(key=lambda row: (row[0]))

# import time 
# t_0 = time.time() 
# POINTS_inside = [] 
 
# bed_points_length = len(BED_points) 
# for f in range(len(PLYGONS_points)): 
#     points_inside = [] 
#     p1 = Polygon(PLYGONS_points[f]) 
#     print(f) 
#     points_inside.extend( 
#         [Point(point[0], point[1], point[2]) for point in point_cloud[:bed_points_length] if 
#          Point(point[0], point[1], point[2]).within(p1) or Point(point[0], point[1], point[2]).touches(p1)] 
#     ) 
#     POINTS_inside.append(points_inside) 
 
# T = time.time() - t_0 
# print(T) 
 
# import pickle 
# # Save the data to a pickle file 
# with open('0706-COMBINED-INpolygones_V4.pkl', 'wb') as picklefile: 
#     pickle.dump(POINTS_inside, picklefile) 
#print(POINTS_inside[8]) 
 
with open('C:\\Users\\KhorsandiKuhanestP\\OneDrive - University of Twente\\D_Flow_FM\\Jupyter scripts\\Volume Modifiying\\0706-COMBINED-INpolygones_V4.pkl', 'rb') as picklefile: 
    POINTS_inside = pickle.load(picklefile)
    
    
    
df = pd.DataFrame(POINTS_inside[0])


x_coords = []
y_coords = []
z_coords = []

for point in POINTS_inside[0]:
    x_coords.append(point.x)
    y_coords.append(point.y)
    z_coords.append(point.z)
for g in range (len(CROSS_section_sorted[0])):
    x_coords.append(CROSS_section_sorted[0][g][0])
    y_coords.append(CROSS_section_sorted[0][g][1])
    z_coords.append(CROSS_section_sorted[0][g][2])
for e in range (len(CROSS_section_sorted[1])):
    x_coords.append(CROSS_section_sorted[1][len(CROSS_section_sorted[1])-1-e][0])
    y_coords.append(CROSS_section_sorted[1][len(CROSS_section_sorted[1])-1-e][1])
    z_coords.append(CROSS_section_sorted[1][len(CROSS_section_sorted[1])-1-e][2])
fig = plt.figure(figsize=(15,12))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_coords, y_coords, z_coords)

ax.set_xlabel('X-Coordinate')
ax.set_ylabel('Y-Coordinate')
ax.set_zlabel('Z-Coordinate')

plt.show()


#Ploting the cross section
fig= plt.figure(figsize= (14,7))
plt.plot(X_CROSS[0], Y_CROSS[0], label = "meshed cross section")
plt.scatter(X_CROSS[0], Y_CROSS[0],c='y', label='mesh points')
# naming the x axis
plt.xlabel('X')
# naming the y axis
plt.ylabel('Y')
# giving a title to my graph
plt.title('Cross section')

# add scatter plot for Width_X_S
plt.scatter(Width_X_S[0], WATER_LEVEL[0], c='b', label='Width_X_S')
# add scatter plot for Width_X_E
plt.scatter(Width_X_E[0], WATER_LEVEL[0], c='r', label='Width_X_E')

# show a legend on the plot
plt.legend()
# function to show the plot
plt.show()


VOLUME_MESH=[]
VOLUME_BED=[]
V_Dif=[]
for o in range (1,len(CROSS_section_sorted)):
    V_Mesh=[]
    V_Bedlevel=[]
    for wl in range (min(len(Width_L_S[o-1]), len(WATER_LEVEL[o]))):
        BND_Depth = float(wl * WLR + 3 * WLR) 

        closest_sublist = min(Depth_Wetted_Area[0], key=lambda sublist: abs(sublist[0] - BND_Depth)) 
        BND_R23A = closest_sublist[2] 
  
        
        closest_R23Ao = min(Depth_Wetted_Area[o], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
        Water_Eqo = closest_R23Ao[1] 
        sublisto=closest_R23Ao
        closest_R23Ao_1 = min(Depth_Wetted_Area[o-1], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
        Water_Eqo_1 = closest_R23Ao_1[1] 
        sublisto_1=closest_R23Ao_1
         
       
        #Creating the points for bathymetry
        x_coords_WL_POLY = []
        y_coords_WL_POLY = []
        z_coords_WL_POLY = []
        x_coords_WL_POLY.append(sublisto_1[6])
        y_coords_WL_POLY.append(sublisto_1[5])
        z_coords_WL_POLY.append(Water_Eqo_1)
        for t in range (len(X_CROSS[int(o-1)])):
            if X_CROSS[int(o-1)][t] < sublisto_1[5] and X_CROSS[int(o-1)][t] > sublisto_1[3]:
                x_coords_WL_POLY.append(L_CROSS[int(o-1)][t])
                y_coords_WL_POLY.append(X_CROSS[int(o-1)][t])
                z_coords_WL_POLY.append(Y_CROSS[int(o-1)][t])
        x_coords_WL_POLY.append(sublisto_1[4])
        y_coords_WL_POLY.append(sublisto_1[3])
        z_coords_WL_POLY.append(Water_Eqo_1)
        x_coords_WL_POLY.append(sublisto[4])
        y_coords_WL_POLY.append(sublisto[3])
        z_coords_WL_POLY.append(Water_Eqo)           
        for t in range (len(X_CROSS[int(o)])):
            if X_CROSS[int(o)][len(X_CROSS[int(o)])-t-1] < sublisto[5] and X_CROSS[int(o)][len(X_CROSS[int(o)])-t-1] > sublisto[3]:
                x_coords_WL_POLY.append(L_CROSS[int(o)][len(X_CROSS[int(o)])-t-1])
                y_coords_WL_POLY.append(X_CROSS[int(o)][len(X_CROSS[int(o)])-t-1])
                z_coords_WL_POLY.append(Y_CROSS[int(o)][len(X_CROSS[int(o)])-t-1]) 
        x_coords_WL_POLY.append(sublisto[6])
        y_coords_WL_POLY.append(sublisto[5])
        z_coords_WL_POLY.append(Water_Eqo) 
        x_coords_WL_POLY.append(sublisto_1[6])
        y_coords_WL_POLY.append(sublisto_1[5])
        z_coords_WL_POLY.append(Water_Eqo_1)    
        points_Polygon_wl=list(zip(x_coords_WL_POLY, y_coords_WL_POLY))    
        points_wl_3D = list(zip(x_coords_WL_POLY, y_coords_WL_POLY,z_coords_WL_POLY ))
        p_wl = Polygon(points_Polygon_wl)
        points_inside_WL=list()
        x_coords_inside=list()
        y_coords_inside=list()
        z_coords_inside=list()
        for point in POINTS_inside[int(o-1)]:
            x_coords_inside.append(point.x)
            y_coords_inside.append(point.y)
            z_coords_inside.append(point.z)
        for i in range (len(x_coords_inside)-1):
            x=x_coords_inside[i]
            y=y_coords_inside[i]
            z=z_coords_inside[i]
            bedpoints_wl=Point(x, y, z)
            if bedpoints_wl.within(p_wl) or bedpoints_wl.touches(p_wl):
                points_inside_WL.append(bedpoints_wl)
        x_points_WL = []
        y_points_WL = []
        z_points_WL = []
        for point in points_inside_WL:
            x_points_WL.append(point.x)
            y_points_WL.append(point.y)
            z_points_WL.append(point.z)
        z_points_WL_New=[]
        for i in range (len(z_points_WL)):
            if z_points_WL[i] > Water_Eqo:
                z_points_WL[i] = Water_Eqo
            z_points_WL_New.append(z_points_WL[i])
        COM_point_x=[]
        COM_point_y=[]
        COM_point_z=[]
        for i in range (len(x_coords_WL_POLY)-1):
            COM_point_x.append(x_coords_WL_POLY[i])
            COM_point_y.append(y_coords_WL_POLY[i])
            COM_point_z.append(z_coords_WL_POLY[i])
        for i in range (len(x_points_WL)-1):
            COM_point_x.append(x_points_WL[i])
            COM_point_y.append(y_points_WL[i])
            COM_point_z.append(z_points_WL_New[i])
       # Create the triangulation object
        triang = mtri.Triangulation(COM_point_x, COM_point_y)
        # Get the list of lists of points for each triangle
      
        triangle_points = [triang.triangles[i] for i in range(triang.triangles.shape[0])]
        Volume_bed=0
        for i in range(len(triangle_points)):
        # Calculate the area
            area = 0.5 * abs(
                COM_point_x[triangle_points[i][0]] * (COM_point_y[triangle_points[i][1]] - COM_point_y[triangle_points[i][2]]) +
                COM_point_x[triangle_points[i][1]] * (COM_point_y[triangle_points[i][2]] - COM_point_y[triangle_points[i][0]]) +
                COM_point_x[triangle_points[i][2]] * (COM_point_y[triangle_points[i][0]] - COM_point_y[triangle_points[i][1]])
            )
            # Calculate the height component
            height_component = ((Water_Eqo+Water_Eqo_1)/2)-((COM_point_z[triangle_points[i][0]] + COM_point_z[triangle_points[i][1]] + COM_point_z[triangle_points[i][2]]) / 3)
            # Calculate the volume of the current triangle and add it to the total volume
            Volume_bed += area * height_component
        #_________________________________________________________________________
        #Creating surface from mesh points
        x_coords_mesh = []
        y_coords_mesh = []
        z_coords_mesh = []
        x_coords_mesh.append(sublisto_1[6])
        y_coords_mesh.append(sublisto_1[5])
        z_coords_mesh.append(Water_Eqo)
        for t in range (len(X_CROSS[o-1])):
            if X_CROSS[int(o-1)][t]< sublisto_1[5] and X_CROSS[int(o-1)][t]> sublisto_1[3]:
                x_coords_mesh.append(sublisto_1[4])
                y_coords_mesh.append(X_CROSS[int(o-1)][t])
                z_coords_mesh.append(Y_CROSS[int(o-1)][t])
        x_coords_mesh.append(sublisto_1[4])
        y_coords_mesh.append(sublisto_1[3])
        z_coords_mesh.append(Water_Eqo_1)
        x_coords_mesh.append(sublisto[6])
        y_coords_mesh.append(sublisto[5])
        z_coords_mesh.append(Water_Eqo)            
        for t in range (len(X_CROSS[int(o)])):
            if X_CROSS[int(o)][t] < sublisto[5] and X_CROSS[int(o)][t] > sublisto[3]:
                x_coords_mesh.append(L_CROSS[int(o)][t])
                y_coords_mesh.append(X_CROSS[int(o)][t])
                z_coords_mesh.append(Y_CROSS[int(o)][t]) 
        x_coords_mesh.append(sublisto[4])
        y_coords_mesh.append(sublisto[3])
        z_coords_mesh.append(Water_Eqo)
       # Create the triangulation object
        triang_mesh = mtri.Triangulation(x_coords_mesh, y_coords_mesh)
        # Plot the surface using the given points
        #surf_mesh = ax.plot_trisurf(triang_mesh, z_coords_mesh, cmap='viridis', edgecolor='k')
        # Get the list of lists of points for each triangle
        triangle_points_mesh = [triang_mesh.triangles[i] for i in range(triang_mesh.triangles.shape[0])]
        Volume_mesh=0
        for i in range(len(triangle_points_mesh)):
        # Calculate the area
            area_mesh = 0.5 * abs(
                x_coords_mesh[triangle_points_mesh[i][0]] * (y_coords_mesh[triangle_points_mesh[i][1]] - y_coords_mesh[triangle_points_mesh[i][2]]) +
                x_coords_mesh[triangle_points_mesh[i][1]] * (y_coords_mesh[triangle_points_mesh[i][2]] - y_coords_mesh[triangle_points_mesh[i][0]]) +
                x_coords_mesh[triangle_points_mesh[i][2]] * (y_coords_mesh[triangle_points_mesh[i][0]] - y_coords_mesh[triangle_points_mesh[i][1]])
            )
            # Calculate the height component
            height_component_mesh = ((Water_Eqo+Water_Eqo_1)/2)-((z_coords_mesh[triangle_points_mesh[i][0]] + z_coords_mesh[triangle_points_mesh[i][1]] + z_coords_mesh[triangle_points_mesh[i][2]]) / 3)
            #print(\"height_component_mesh\",BND_Depth,  height_component_mesh)
            # Calculate the volume of the current triangle and add it to the total volume
            Volume_mesh += area_mesh * height_component_mesh 
 
        Dif=Volume_bed-Volume_mesh
        V_Mesh.append(Volume_mesh)
        V_Bedlevel.append(Volume_bed)
        V_Dif.append(Dif)
    VOLUME_MESH.append(V_Mesh)
    VOLUME_BED.append(V_Bedlevel)

def remove_repeats(lst):
    Active_points_new = []
    previous = set()
    for sub_list in lst:
        sub_list = list(set(sub_list))
        new_numbers = [num for num in sub_list if num not in previous]
        Active_points_new.append(new_numbers)
        previous.update(sub_list)
    return Active_points_new
# finding active points for all of the cross sections
ACTIVE_POINTS=[]
NON_Zero=[]
SUB_LIST=[]
VOLUME_DIF=[]
for y in range (1, len(CROSS_section_sorted)):
    Volume_WL_DIF=[]
    CN_All=list()
    for i in range (len(Y_CROSS[y])):
        if Y_CROSS[y][i]==LOWEST_min[y]:
            num_lowest=i
    star_right=num_lowest
    star_left=num_lowest
    Active_points=list()
    for wl in range (min(len(Width_L_S[y-1]), len(WATER_LEVEL[y]))):
        Change_N=list()
        #for the right part of the river
        for i in range(int(star_right), len(Y_CROSS[y])-1):
            star_right_l=i-1
            if float(Y_CROSS[y][i])<float(WATER_LEVEL[y][wl]):
                Change_N.append(i)
            if Y_CROSS[y][i] > WATER_LEVEL[y][wl]:
                    break
        star_right=star_right_l
        #for the left part of the river
        for k in range(0, int(star_left)+1):
            j=star_left-k
            star_left_l=j
            if float(Y_CROSS[y][j])<float(WATER_LEVEL[y][wl]):
                Change_N.append(j)
            if Y_CROSS[y][j] > WATER_LEVEL[y][wl]:
                    break
        star_left=star_left_l+1
        CN= len(Change_N)
        if wl==0:
            CN=CN-1
        else:
                CN=CN-2
        CN_All.append(CN)
        Active_points.append(Change_N)
    Active_points_new = remove_repeats(Active_points)
    ACTIVE_POINTS.append(Active_points_new)
    non_zero=[item for item in range(len(CN_All)) if CN_All[item]>0]
    NON_Zero.append(non_zero)
    
    SUB_lists=[]
    for u in range (len(non_zero)):
        Dif=0
        sub_list=list()
        if u<len(non_zero)-1:
            sub_list=(CN_All[non_zero[u]:non_zero[u+1]])
        if u ==len(non_zero)-1:
            sub_list=CN_All[non_zero[-1]:len(CN_All)]
        for j in range (len(sub_list)):
            wl=non_zero[u]+j
            Dif=Dif+VOLUME_BED[y-1][wl]-VOLUME_MESH[y-1][wl]
        Volume_WL_DIF.append(Dif/len(sub_list))  
        SUB_lists.append (sub_list)
    VOLUME_DIF.append(Volume_WL_DIF)
    SUB_LIST.append(SUB_lists)
    
OFFSET_V=[]
for o in range  (1,len(CROSS_section_sorted)):
    Offset_value=[]
    for i in range (len(NON_Zero[o-1])):
        Offset_value.append(0)
    OFFSET_V.append(Offset_value)


Y_CROSS_NEW=[]
for o in range (0, len(CROSS_section_sorted)):   
    Y_CROSS_new=[]
    for t in range (len(X_CROSS[o])):
        Y_CROSS_new.append(Y_CROSS[int(o)][t])
    Y_CROSS_NEW.append(Y_CROSS_new)
    
    
    
#the cross section y is the cross section meant to change
for o in range  (1,  len(CROSS_section_sorted)):
    # Loop until the function starts to increase
    print("o=", o)
    for i in range (len(NON_Zero[o-1])):
        DF_sublist=0
        for j in range (len(SUB_LIST[o-1][i])):
            wl=NON_Zero[o-1][i]+j
            BND_Depth = float(wl * WLR + 3 * WLR) 
            print(BND_Depth) 

            closest_sublist = min(Depth_Wetted_Area[0], key=lambda sublist: abs(sublist[0] - BND_Depth)) 
            BND_R23A = closest_sublist[2] 
            closest_R23Ao = min(Depth_Wetted_Area[o], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
            Water_Eqo = closest_R23Ao[1] 
            sublisto=closest_R23Ao
            closest_R23Ao_1 = min(Depth_Wetted_Area[o-1], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
            Water_Eqo_1 = closest_R23Ao_1[1] 
            sublisto_1=closest_R23Ao_1
#-----------
            #Creating surface from mesh points
            x_coords_mesh = []
            y_coords_mesh = []
            z_coords_mesh = []
            x_coords_mesh.append(sublisto_1[6])
            y_coords_mesh.append(sublisto_1[5])
            z_coords_mesh.append(Water_Eqo)
            for t in range (len(X_CROSS[o-1])):
                if X_CROSS[int(o-1)][t]< sublisto_1[5] and X_CROSS[int(o-1)][t]> sublisto_1[3]:
                    x_coords_mesh.append(sublisto_1[4])
                    y_coords_mesh.append(X_CROSS[int(o-1)][t])
                    z_coords_mesh.append(Y_CROSS_NEW[int(o-1)][t])
            x_coords_mesh.append(sublisto_1[4])
            y_coords_mesh.append(sublisto_1[3])
            z_coords_mesh.append(Water_Eqo_1)
            x_coords_mesh.append(sublisto[6])
            y_coords_mesh.append(sublisto[5])
            z_coords_mesh.append(Water_Eqo)            
            for t in range (len(X_CROSS[int(o)])):
                if X_CROSS[int(o)][t] < sublisto[5] and X_CROSS[int(o)][t] > sublisto[3]:
                    x_coords_mesh.append(L_CROSS[int(o)][t])
                    y_coords_mesh.append(X_CROSS[int(o)][t])
                    z_coords_mesh.append(Y_CROSS_NEW[int(o)][t]) 
            x_coords_mesh.append(sublisto[4])
            y_coords_mesh.append(sublisto[3])
            z_coords_mesh.append(Water_Eqo)
            #Create the triangulation object
            triang_mesh = mtri.Triangulation(x_coords_mesh, y_coords_mesh)
            # Plot the surface using the given points
            #surf_mesh = ax.plot_trisurf(triang_mesh, z_coords_mesh, cmap='viridis', edgecolor='k')
            # Get the list of lists of points for each triangle
            triangle_points_mesh = [triang_mesh.triangles[i] for i in range(triang_mesh.triangles.shape[0])]
            Volume_mesh=0
            for p in range(len(triangle_points_mesh)):
            # Calculate the area
                area_mesh = 0.5 * abs(
                    x_coords_mesh[triangle_points_mesh[p][0]] * (y_coords_mesh[triangle_points_mesh[p][1]] - y_coords_mesh[triangle_points_mesh[p][2]]) +
                    x_coords_mesh[triangle_points_mesh[p][1]] * (y_coords_mesh[triangle_points_mesh[p][2]] - y_coords_mesh[triangle_points_mesh[p][0]]) +
                    x_coords_mesh[triangle_points_mesh[p][2]] * (y_coords_mesh[triangle_points_mesh[p][0]] - y_coords_mesh[triangle_points_mesh[p][1]])
                )
                # Calculate the height component
                height_component_mesh = ((Water_Eqo+Water_Eqo_1)/2)-((z_coords_mesh[triangle_points_mesh[p][0]] + z_coords_mesh[triangle_points_mesh[p][1]] + z_coords_mesh[triangle_points_mesh[p][2]]) / 3)
                #print(\"height_component_mesh\",BND_Depth,  height_component_mesh)
                # Calculate the volume of the current triangle and add it to the total volume
                Volume_mesh += area_mesh * height_component_mesh 
            #__________
#             points_mesh = np.array([x_coords_mesh, y_coords_mesh, z_coords_mesh]).T
#             vertices_mesh = np.column_stack((x_coords_mesh, y_coords_mesh))
#             A_mesh = dict(vertices=vertices_mesh)
#             B_mesh = tr.triangulate(A_mesh)
#             groups_mesh=B_mesh['triangles']
#             Volume_mesh=0
#             for r in range(len(groups_mesh)):
#                 # Define the points
#                 p1 = points_mesh[groups_mesh[r][0]]
#                 p2 = points_mesh[groups_mesh[r][1]]
#                 p3 = points_mesh[groups_mesh[r][2]]
#                 # Calculate the area
#                 area = 0.5 * abs(p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1]-p2[1]))
#                 Volume_mesh=Volume_mesh+(((WATER_LEVEL[int(o)][wl]-p1[2])+(WATER_LEVEL[int(o)][wl]-p2[2])+(WATER_LEVEL[int(o)][wl]-p3[2]))*area)/3
            Dif=VOLUME_BED[o-1][wl]-Volume_mesh
            DF_sublist=DF_sublist+Dif
        DF_sublist_Ave=(DF_sublist/len(SUB_LIST[o-1][i]))
        if DF_sublist_Ave < 0:
                s = 1
        elif DF_sublist_Ave > 0:
                s = -1
        Initial_Target=abs(DF_sublist_Ave)
        dx=0
        # Loop until the function starts to increase
        min_dx=0
        while True:
            DF_sublist=0
            Volume_total_mesh_changed = 0
            for j in range (len(SUB_LIST[o-1][i])):
                wl=NON_Zero[o-1][i]+j
                BND_Depth = float(wl * WLR + 3 * WLR) 
                closest_sublist = min(Depth_Wetted_Area[0], key=lambda sublist: abs(sublist[0] - BND_Depth)) 
                BND_R23A = closest_sublist[2] 
                closest_R23Ao = min(Depth_Wetted_Area[o], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
                Water_Eqo = closest_R23Ao[1] 
                sublisto=closest_R23Ao
                closest_R23Ao_1 = min(Depth_Wetted_Area[o-1], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
                Water_Eqo_1 = closest_R23Ao_1[1] 
                sublisto_1=closest_R23Ao_1
                #Creating surface from mesh points
                x_coords_mesh = []
                y_coords_mesh = []
                z_coords_mesh = []
                x_coords_mesh.append(sublisto_1[6])
                y_coords_mesh.append(sublisto_1[5])
                z_coords_mesh.append(Water_Eqo)
                for t in range (len(X_CROSS[o-1])):
                    if X_CROSS[int(o-1)][t]< sublisto_1[5] and X_CROSS[int(o-1)][t]> sublisto_1[3]:
                        x_coords_mesh.append(sublisto_1[4])
                        y_coords_mesh.append(X_CROSS[int(o-1)][t])
                        z_coords_mesh.append(Y_CROSS_NEW[int(o-1)][t])
                x_coords_mesh.append(sublisto_1[4])
                y_coords_mesh.append(sublisto_1[3])
                z_coords_mesh.append(Water_Eqo_1)
                x_coords_mesh.append(sublisto[6])
                y_coords_mesh.append(sublisto[5])
                z_coords_mesh.append(Water_Eqo)            
                for t in range (len(X_CROSS[int(o)])):
                    if X_CROSS[int(o)][t] < sublisto[5] and X_CROSS[int(o)][t] > sublisto[3]:
                        x_coords_mesh.append(L_CROSS[int(o)][t])
                        y_coords_mesh.append(X_CROSS[int(o)][t])
                        if t in ACTIVE_POINTS[o-1][NON_Zero[o-1][i]]: 
                            z_coords_mesh.append(Y_CROSS_NEW[o][t]+(dx*s)) 
                        else: 
                            z_coords_mesh.append(Y_CROSS_NEW[o][t]) 
                x_coords_mesh.append(sublisto[4])
                y_coords_mesh.append(sublisto[3])
                z_coords_mesh.append(Water_Eqo)
               # Create the triangulation object
                triang_mesh = mtri.Triangulation(x_coords_mesh, y_coords_mesh)
                # Plot the surface using the given points
                #surf_mesh = ax.plot_trisurf(triang_mesh, z_coords_mesh, cmap='viridis', edgecolor='k')
                # Get the list of lists of points for each triangle
                triangle_points_mesh = [triang_mesh.triangles[i] for i in range(triang_mesh.triangles.shape[0])]
                Volume_mesh_changed=0
                for p in range(len(triangle_points_mesh)):
                # Calculate the area
                    area_mesh = 0.5 * abs(
                        x_coords_mesh[triangle_points_mesh[p][0]] * (y_coords_mesh[triangle_points_mesh[p][1]] - y_coords_mesh[triangle_points_mesh[p][2]]) +
                        x_coords_mesh[triangle_points_mesh[p][1]] * (y_coords_mesh[triangle_points_mesh[p][2]] - y_coords_mesh[triangle_points_mesh[p][0]]) +
                        x_coords_mesh[triangle_points_mesh[p][2]] * (y_coords_mesh[triangle_points_mesh[p][0]] - y_coords_mesh[triangle_points_mesh[p][1]])
                    )
                    # Calculate the height component
                    height_component_mesh = ((Water_Eqo+Water_Eqo_1)/2)-((z_coords_mesh[triangle_points_mesh[p][0]] + z_coords_mesh[triangle_points_mesh[p][1]] + z_coords_mesh[triangle_points_mesh[p][2]]) / 3)
                    #print(\"height_component_mesh\",BND_Depth,  height_component_mesh)
                    # Calculate the volume of the current triangle and add it to the total volume
                    Volume_mesh_changed += area_mesh * height_component_mesh 
                DF=VOLUME_BED[o-1][wl]-Volume_mesh_changed#if you want to use it for different cross sections this part should be changed
                DF_sublist=DF_sublist+DF
            DF_sublist_Ave=DF_sublist/len(SUB_LIST[o-1][i])
          
            if abs(DF_sublist_Ave)<=Initial_Target and abs(dx)<0.3:
                min_dx=(dx*s)
                Initial_Target=abs(DF_sublist_Ave)
                dx += 0.02
            else:
                break
        OFFSET_V[o-1][i]=(dx-0.05)*s
        for g in range (0, len(ACTIVE_POINTS[o-1][NON_Zero[o-1][i]])):
            Y_CROSS_NEW[o][ACTIVE_POINTS[o-1][NON_Zero[o-1][i]][g]]=Y_CROSS[int(o)][ACTIVE_POINTS[o-1][NON_Zero[o-1][i]][g]]+(dx-0.05)*s
            
            
            
VOLUME_MESH_NEW=[]
for o in range (1,len(CROSS_section_sorted)):
    Volume_Mesh=[]
    for wl in range (min(len(Width_L_S[o-1]), len(WATER_LEVEL[o]))):
        BND_Depth = float(wl * WLR + 3 * WLR) 

        closest_sublist = min(Depth_Wetted_Area[0], key=lambda sublist: abs(sublist[0] - BND_Depth)) 
        BND_R23A = closest_sublist[2] 
  
        
        closest_R23Ao = min(Depth_Wetted_Area[o], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
        Water_Eqo = closest_R23Ao[1] 
        sublisto=closest_R23Ao
        closest_R23Ao_1 = min(Depth_Wetted_Area[o-1], key=lambda sublist: abs(sublist[2] - BND_R23A)) 
        Water_Eqo_1 = closest_R23Ao_1[1] 
        sublisto_1=closest_R23Ao_1
        #Creating surface from mesh points
        x_coords_mesh = []
        y_coords_mesh = []
        z_coords_mesh = []
        x_coords_mesh.append(sublisto_1[6])
        y_coords_mesh.append(sublisto_1[5])
        z_coords_mesh.append(Water_Eqo)
        for t in range (len(X_CROSS[o-1])):
            if X_CROSS[int(o-1)][t]< sublisto_1[5] and X_CROSS[int(o-1)][t]> sublisto_1[3]:
                x_coords_mesh.append(sublisto_1[4])
                y_coords_mesh.append(X_CROSS[int(o-1)][t])
                z_coords_mesh.append(Y_CROSS_NEW[int(o-1)][t])
        x_coords_mesh.append(sublisto_1[4])
        y_coords_mesh.append(sublisto_1[3])
        z_coords_mesh.append(Water_Eqo_1)
        x_coords_mesh.append(sublisto[6])
        y_coords_mesh.append(sublisto[5])
        z_coords_mesh.append(Water_Eqo)            
        for t in range (len(X_CROSS[int(o)])):
            if X_CROSS[int(o)][t] < sublisto[5] and X_CROSS[int(o)][t] > sublisto[3]:
                x_coords_mesh.append(L_CROSS[int(o)][t])
                y_coords_mesh.append(X_CROSS[int(o)][t])
                z_coords_mesh.append(Y_CROSS_NEW[int(o)][t]) 
        x_coords_mesh.append(sublisto[4])
        y_coords_mesh.append(sublisto[3])
        z_coords_mesh.append(Water_Eqo)
       # Create the triangulation object
        triang_mesh = mtri.Triangulation(x_coords_mesh, y_coords_mesh)
        # Plot the surface using the given points
        #surf_mesh = ax.plot_trisurf(triang_mesh, z_coords_mesh, cmap='viridis', edgecolor='k')
        # Get the list of lists of points for each triangle
        triangle_points_mesh = [triang_mesh.triangles[i] for i in range(triang_mesh.triangles.shape[0])]
        Volume_mesh=0
        for i in range(len(triangle_points_mesh)):
        # Calculate the area
            area_mesh = 0.5 * abs(
                x_coords_mesh[triangle_points_mesh[i][0]] * (y_coords_mesh[triangle_points_mesh[i][1]] - y_coords_mesh[triangle_points_mesh[i][2]]) +
                x_coords_mesh[triangle_points_mesh[i][1]] * (y_coords_mesh[triangle_points_mesh[i][2]] - y_coords_mesh[triangle_points_mesh[i][0]]) +
                x_coords_mesh[triangle_points_mesh[i][2]] * (y_coords_mesh[triangle_points_mesh[i][0]] - y_coords_mesh[triangle_points_mesh[i][1]])
            )
            # Calculate the height component
            height_component_mesh = ((Water_Eqo+Water_Eqo_1)/2)-((z_coords_mesh[triangle_points_mesh[i][0]] + z_coords_mesh[triangle_points_mesh[i][1]] + z_coords_mesh[triangle_points_mesh[i][2]]) / 3)
            #print(\"height_component_mesh\",BND_Depth,  height_component_mesh)
            # Calculate the volume of the current triangle and add it to the total volume
            Volume_mesh += area_mesh * height_component_mesh 
        Volume_Mesh.append(Volume_mesh)
    VOLUME_MESH_NEW.append(Volume_Mesh)
    
    
fig, ax = plt.subplots(figsize=(15,5))

# Make sure both lists have the same length
water_level = WATER_LEVEL[1][:(len(VOLUME_BED[0]))]
volume_bed = VOLUME_BED[0]
volume_mesh_NEW = VOLUME_MESH_NEW[0]

# Calculate the difference between volume_bed and volume_mesh
diff = np.array(volume_bed) - np.array(volume_mesh_NEW)
# Make sure both lists have the same length
water_level = WATER_LEVEL[1][:(len(VOLUME_BED[0]))]
volume_bed = VOLUME_BED[0]
volume_mesh = VOLUME_MESH[0]

# Calculate the difference between volume_bed and volume_mesh
diff_2 = np.array(volume_bed) - np.array(volume_mesh)

plt.plot(water_level, diff_2, label='Bedlevel Volume - VOLUME_MESH')
plt.plot(water_level, diff, label='Bedlevel Volume - VOLUME_MESH_NEW')

plt.xlabel('Water level')
plt.ylabel('Volume Difference')
plt.legend()
plt.show()



from sklearn.metrics import mean_absolute_error
mean_absolute_error_current=[]
mean_absolute_error_new=[]
for t in range (len (VOLUME_BED)):
    mean_absolute_error_current.append(mean_absolute_error(VOLUME_BED[t], VOLUME_MESH[t]))
    mean_absolute_error_new.append(mean_absolute_error(VOLUME_BED[t], VOLUME_MESH_NEW[t]))
# plot the data
fig= plt.figure(figsize= (14,7))
plt.plot(range(len(VOLUME_BED)), mean_absolute_error_current, label='current')
plt.plot(range(len(VOLUME_BED)), mean_absolute_error_new, label='new')
plt.xlabel('Num. mesh cross section')
plt.ylabel('Mean Absolute Error')
plt.legend()
plt.show()



x_Cmesh_NEW =[] 
y_Cmesh_NEW =[] 
Z_Cmesh_NEW =[] 
Z_Cmesh_Mod =[]

for y in range (len(CROSS_section_sorted)):
        for u in range (len(CROSS_section_sorted[y])):
            x_Cmesh_NEW.append(L_CROSS[y][u])
            y_Cmesh_NEW.append(X_CROSS[y][u])
            Z_Cmesh_NEW.append(Y_CROSS[y][u])
            Z_Cmesh_Mod.append(Y_CROSS_NEW[y][u])


raw_data_2 = {'x':x_Cmesh_NEW ,
                'y': y_Cmesh_NEW,
                'z': Z_Cmesh_Mod}
df_2 = pd.DataFrame(raw_data_2, columns = ['x', 'y', 'z'])
T = time.time() - t_0 
print(T) 
df_2.to_csv('20231018- AR product for slices.csv', index=False)

print(VOLUME_MESH[0])

AVG_H_LM=list()
AVG_H_L=list()  
MIN_H_LM=list()
MAX_H_LM=list()
MIN_H_L=list()
MAX_H_L=list()
for c in range (len(VOLUME_MESH)):
    for v in range (len(VOLUME_MESH[c])):
            H_LM=list()
            H_L=list()
            High_Lowmodified=VOLUME_BED[c][v]-VOLUME_MESH_NEW[c][v]
            High_Low=VOLUME_BED[c][v]-VOLUME_MESH[c][v]
            H_LM.append(High_Lowmodified)
            H_L.append(High_Low)  
    AVG_H_LM.append(np.mean(H_LM))
    AVG_H_L.append(np.mean(H_L))   
    MIN_H_LM.append(np.min(H_LM))
    MAX_H_LM.append(np.max(H_LM))
    MIN_H_L.append(np.min(H_L))
    MAX_H_L.append(np.max(H_L))
    
print("AVG_H_LM", AVG_H_LM)    
print("AVG_H_L", AVG_H_L)
print("MIN_H_LM", MIN_H_LM)  
print("MAX_H_LM", MAX_H_LM)
print("MIN_H_L", MIN_H_L)
print("MAX_H_L", MAX_H_L)
  
averages = [sum(x) / len(VOLUME_MESH_NEW) for x in zip(*VOLUME_MESH_NEW)]
averages2 = [sum(x) / len(VOLUME_BED) for x in zip(*VOLUME_BED)]
averages3 = [sum(x) / len(VOLUME_MESH) for x in zip(*VOLUME_MESH)]
std_devs = [np.std(x) for x in zip(*VOLUME_MESH_NEW)]
std_devs2 = [np.std(x) for x in zip(*VOLUME_BED)]
std_devs3 = [np.std(x) for x in zip(*VOLUME_MESH)]

print(averages)
print(averages2)
print(averages3)
print(std_devs)
print(std_devs2)
print(std_devs3)

Mesh_volume_data2 = {'VOLUME_MESH_NEW': averages,
                'VOLUME_BED': averages2,           
                'VOLUME_MESH': averages3}
df_modified_mesh_area2 = pd.DataFrame(Mesh_volume_data2, columns = ['VOLUME_MESH_NEW', 'VOLUME_BED','VOLUME_MESH'])
print(df_modified_mesh_area2)

Mesh_volume_data3 = {'std_devs_MESH_NEW': std_devs,  
                'std_devs_BED': std_devs2, 'std_devs_MESH':std_devs3}
df_modified_mesh_area3 = pd.DataFrame(Mesh_volume_data3, columns = ['std_devs_MESH_NEW', 'std_devs_BED', 'std_devs_MESH'])
print(df_modified_mesh_area3)