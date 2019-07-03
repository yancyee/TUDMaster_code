import numpy as np
import laspy as lp
import json
import math
import os
import matplotlib.pyplot as plt
import time
import datetime
import sys
from shapely.geometry import Polygon

# Determine tile containing point,
# based on whether point lies within tile's extent
def find_tile():
    for tile in tilelist:
        # tile name contains coordinates of bounding box

        new_tile = tile.strip(".las").split(",")


        tile_min_x = float(new_tile[2])
        tile_min_y = float(new_tile[3])
        tile_max_x = float(new_tile[4])
        tile_max_y = float(new_tile[5])
        if tile_min_x <= x <= tile_max_x and tile_min_y <= y <= tile_max_y:

            return (int(new_tile[0]), int(new_tile[1]))

# Search for the tiles that are adjacent to the initially selected tile
def find_tile_grid(row, col):
    iteration_list = [-1, 0, 1]
    tile_grid = []
    
    for i in iteration_list:
        for j in iteration_list:
            tilename_start = "{},{}".format((row + i), (col + j))
            for tile in tilelist:
                if tile.startswith(tilename_start):
                    tile_grid.append(tile)
                    break
    # tile_grid: list of 9 tiles
    return tile_grid

# determine height of viewpoint by sampling the ground points
# center: location of viewpoint
def getheight(tile_grid):
    center = np.array([x, y])
    pointheight = 0
    points_number = 0
    
    for tile in tile_grid:
        # read the .las file
        file_input = lp.file.File("{}{}{}".format(path, "/", tile), mode='r')

        # keep groundpoints satisfying ground_rules:
        # classification 2 for ground, inside las file
        # keep points within radius of 5 metres
        ground_rules = np.logical_and(
            file_input.raw_classification == 2,
            np.sqrt(np.sum((np.vstack((file_input.x, file_input.y)).transpose() - center) ** 2, axis=1)) <= 1)
        build_rules = np.logical_and(
            file_input.raw_classification == 6,
            np.sqrt(np.sum((np.vstack((file_input.x, file_input.y)).transpose() - center) ** 2, axis=1)) <= 1)
            
        ground_points = file_input.points[ground_rules]
        build_points = file_input.points[build_rules]

        # make array with heights of each point
        if ground_points.size > build_points.size:
            ground_point_heights = np.array((ground_points['point']['Z'])).transpose()
        else:
            ground_point_heights = np.array((build_points['point']['Z'])).transpose()
            
        if ground_point_heights.size > 0:
            pointheight += float(np.sum(ground_point_heights))
            points_number += ground_point_heights.size
    
    # get mean value of points' heights
    if points_number > 0:
        height = pointheight / points_number
        return (height+500)
    else:
        return 0


# function to get all points lying within range of the defined radius from the viewpoint
def getPoints(tile_grid, radius, view_height):
    # Viewpoint
    center = np.array([x, y])
    
    # Gather points
    arraysX, arraysY, arraysZ = [], [], []  # list of arrays of X,Y,Z coords
    arrayDistances = []  # Horizontal distances
    arrayClasses = []  # Classifications
    toBeAdded = []
    for tile in tile_grid:
        inFile = lp.file.File("{}{}{}".format(path, "/", tile), mode='r')
        coords = np.vstack((inFile.x, inFile.y)).transpose()
        elevation = inFile.z
        distances = np.sqrt(np.sum((coords - center)**2, axis=1))

        keep_points = np.logical_and(np.logical_and(np.logical_or(
            inFile.raw_classification == 1,
            inFile.raw_classification == 6),
            distances < radius),
            elevation >= view_height/1000)

        # Get coordinates
        arraysX.append(inFile.x[keep_points])
        arraysY.append(inFile.y[keep_points])
        arraysZ.append(inFile.z[keep_points])
        # Get distances
        arrayDistances.append(distances[keep_points])
        # Get classifications
        arrayClasses.append(inFile.raw_classification[keep_points])
  
    # Concatenate all information
    X, Y, Z = arraysX[0], arraysY[0], arraysZ[0]
    distances = arrayDistances[0]
    classes = arrayClasses[0]
    for arrayX, arrayY in zip(arraysX[1:], arraysY[1:]):
        X = np.hstack([X, arrayX])
        Y = np.hstack([Y, arrayY])
    for arrayZ in arraysZ[1:]:
        Z = np.hstack([Z, arrayZ])
    for arDist in arrayDistances[1:]:
        distances = np.hstack([distances, arDist])
    for arClass in arrayClasses[1:]:
        classes = np.hstack([classes, arClass])


    return X, Y, Z, distances, classes


# Create dome
def createDome(X, Y, Z, dists, classes, view_height):
    # Initialize dome
    # Indices = (Azimuth, Elevation)
    dome = np.zeros((180, 90), dtype=int)
    domeDists = np.zeros((180, 90), dtype=int)

    if X.size > 0:
        # Azimuths
        dX, dY = X - x, Y - y
        azimuths = np.arctan2(dY, dX) * 180 / math.pi - 90
        azimuths[azimuths < 0] += 360

        # Elevations
        dZ = Z - view_height / 1000
        elevations = np.arctan2(dZ, dists) * 180 / math.pi

        # Shade sectors
        # Array with dome indices, distances & classifications
        data = np.stack((azimuths // 2, elevations // 1, dists, classes), axis=-1)
        # Sort according to indices & classifications
        sortData = data[np.lexsort([data[:, 2], data[:, 1], data[:, 0]])]

        # Spot where azimuth & elevation values change
        azimuth_change = sortData[:, 0][:-1] != sortData[:, 0][1:]
        elevation_change = sortData[:, 1][:-1] != sortData[:, 1][1:]
        keep = np.where(np.logical_or(azimuth_change, elevation_change))
        # Take position of next element, plus add first row
        shortestDistance = sortData[
            np.insert(keep[0] + 1, 0, 0)]  # (inserts second element of change, first position, index of first point)
        # Define indices & classifications
        hor = shortestDistance[:, 0].astype(int)
        ver = shortestDistance[:, 1].astype(int)
        classif = shortestDistance[:, 3].astype(int)
        dists = shortestDistance[:, 2]

        # Update dome
        dome[hor, ver] = classif

        domeDists[hor, ver] = dists

        # Buildings as solids

        # Find building positions in dome
        # print dome[dome == 6].size
        if dome[dome == 6].size > 0:
            bhor, bver = np.where(dome == 6)
            # Create an array out of them
            builds = np.stack((bhor, bver), axis=-1)
            shape = (builds.shape[0] + 1, builds.shape[1])
            builds = np.append(builds, (bhor[0], bver[0])).reshape(shape)

            # Spot azimuth changes
            azimuth_change = builds[:, 0][:-1] != builds[:, 0][1:]
            keep = np.where(azimuth_change)
            # keep = np.insert(np.where(azimuth_change==True), 0, 0)
            # Change to building up to roof for each row
            roof_rows, roof_cols = builds[keep][:, 0], builds[keep][:, 1]
            for roof_row, roof_col in zip(roof_rows, roof_cols):
                condition = np.where(np.logical_or(domeDists[roof_row, :roof_col] > domeDists[roof_row, roof_col],
                                                   dome[roof_row, :roof_col] == 0))
                dome[roof_row, :roof_col][condition] = 6
    predict=plot(dome)
    #print(plot(dome))
    return dome,predict

# Plot dome
def plot(dome):
    # Create circular grid 
    theta, radius = np.mgrid[0:(2*np.pi+2*np.pi/180):2*np.pi/180, 0:90:1]
    Z = dome.copy().astype(float)
    Z = Z[0:, ::-1]  # Reverse array rows
    # assign colors depending on class
    Z[Z == 0] = 0
    Z[Z == 1] = 0.5
    Z[Z == 6] = 1

    s_list=[]
    predict={}
    for i in  np.arange(5.50,22,0.25):
        # zen=altitude angle azi=horizental angle
        Zen,Azi=sun_position([2018, 5, 27], i, 52, 4)
        #print((Z[1]))

        s_a = round(180 - Azi / 2)
        if s_a==180:
            s_a=179

        s_z = round(90 * math.cos(Zen / 180 * math.pi))
        if s_z==90:
            s_z=89

        s_list.append([s_a, s_z])
        if Z[s_a, s_z] != 0:
            predict[i]="not influenced"
            #print(i,"\t""not influenced")
        else:
            predict[i] = "influenced"
            #print(i,"\t""influenced")

    for n in s_list:
        Z[n[0],n[1]]=1.2


    if Z[Z == 6].size == 0:
        Z[0,0] = 1
    axes = plt.subplot(111, projection='polar')

    cmap = plt.get_cmap('tab20c')


    axes.pcolormesh(theta, radius, Z, cmap=cmap)
    axes.set_ylim([0, 90])
    axes.tick_params(labelleft=False)
    axes.set_theta_zero_location("N")
    #plt.savefig("1"),bbox_inches='tight')
    #plt.show()
    return predict

#calculate sun parameters
def sun_position(d,t,la,lon):
    # date difference from first day of the year
    d1=datetime.datetime(2018,1,1)
    d2=datetime.datetime(d[0],d[1],d[2])
    n_s=(str(d2-d1))
    n=float(n_s.split(" ")[0])

    # time difference between solar time and local  （summer -1）
    time_diff=(lon-15)/15-1
    t=t+time_diff
    if n<106:
        eqt = -14.2 * math.sin(math.pi * (n + 7) / 111)
    elif n>=106 and n<166:
        eqt = 4 * math.sin(math.pi * (n - 107) / 59)
    elif n >= 166 and n < 246:
        eqt = -6.5 * math.sin(math.pi * (n - 166) / 80)
    else:
        eqt = 16.4 * math.sin(math.pi * (n - 247) / 113)
    t=t+eqt/60

    # sun angle calculation(a_d=Azimuthal, h_d=Zenith, omega_rad=hour angle, sigma_rad=declination)
    sigma_rad=23.45*math.pi/180*math.sin(2*math.pi*(284+n)/365)
    omega_rad=math.pi*(12-t)/12
    la_rad=la*math.pi/180
    h=math.asin(math.sin(la_rad)*math.sin(sigma_rad)+math.cos(omega_rad)*math.cos(sigma_rad)*math.cos(la_rad))
    a_arc=(math.sin(h)*math.sin(la_rad)-math.sin(sigma_rad))/(math.cos(h)*math.cos(la_rad))
    a=math.acos(a_arc)
    h_d=h/math.pi*180
    a_d=a/math.pi*180

    # angle correction (Azimuthal angle start from north)
    if t >12:
        a_d=a_d+180
    else:
        a_d=180-a_d
    return ( h_d,a_d)

# calculate SVF, and percentage of building/vegetation obstructions
def calculate_SVF(radius, dome):
    obstructedArea = 0
    treeObstruction = 0
    buildObstruction = 0
    for i in range(0, 180):
        for j in range(0, 90):
            if dome[i, j] != 0:
                v = 90 - (j + 1)
                R = math.cos(v * math.pi / 180) * radius
                r = math.cos((v + 1) * math.pi / 180) * radius
                # calculate area of each obstructed sector (circular sector area calculation)
                cell_area = (math.pi / 180.0) * (R ** 2 - r ** 2)
                obstructedArea += cell_area

                if dome[i, j] == 1:
                    treeObstruction += cell_area
                elif dome[i, j] == 6:
                    buildObstruction += cell_area

    circleArea = math.pi * (radius ** 2)

    # SVF: proportion of open area to total area
    SVF = (circleArea - obstructedArea) / circleArea
    treeObstructionPercentage = treeObstruction / circleArea
    buildObstructionPercentage = buildObstruction / circleArea
    return SVF, treeObstructionPercentage, buildObstructionPercentage

def integer(geom):
    geometry = []
    append = geometry.append
    for point in geom:
        x, y = point[0], point[1]
        append([x, y])
    geometry = np.array(geometry) #Closed polygon
    return geometry


def inside_polygon(pt, minY, maxY, maxX, geom, fraction):
    condition1 = np.logical_and(pt[1] > minY, pt[1] <= maxY)
    condition2 = pt[0] <= maxX
    condition = np.logical_and(condition1, condition2)

    intersX = geom[:, 0][:-1][condition] + (pt[1] - geom[:, 1][:-1][condition]) * fraction[condition]
    truth = np.logical_or(geom[:, 0][:-1][condition] == geom[:, 0][1:][condition],
                          pt[0] <= intersX)

    intersections = truth[truth == True].size

    if intersections % 2 == 1:
        return pt


def run():

    row, col = find_tile()
    tile_grid = find_tile_grid(row, col)
    view_height = getheight(tile_grid)
    X, Y, Z, distances, classes = getPoints(tile_grid, radius, view_height)
    dome,predict = createDome(X, Y, Z, distances, classes, view_height)

    #SVF, tree_percentage, build_percentage = calculate_SVF(radius, dome)
    #SVF, tree_percentage, build_percentage = round(SVF*100), round(tree_percentage*100), round(build_percentage*100)
    #print ('{}%'.format(int(SVF)) + "\n" + '{}%'.format(int(tree_percentage)) + "\n" + '{}%'.format(int(build_percentage)))
    return predict

def xy(coordinate):
    global x
    x=coordinate[0]
    global y
    y=coordinate[1]
    predict=run()
    return predict


"""GLOBAL VARIABLES"""
# path for tile directory and list of tilenames
#path for tile directory and list of tilenames
path = "/Users/xuyixin/Desktop/Tiles_Complete"
tilelist = os.listdir("/Users/xuyixin/Desktop/Tiles_Complete")
# define radius
radius = 100
bufferRadius=1.5
"""END GLOBAL VARIABLES"""

#print(xy([78838.805, 457194.938]))
#x = 78613.134
#y = 454075.110
#print(xy([78834.805, 457190.938]))




