from shapely.geometry import Point, Polygon
import math
import fiona


# create potential location use two loop
def sctter_point(pt,radius,angle):
    pt_new_lis = []
    pt_new_lis.append(pt)
    ang_rad=angle/180*math.pi
    for i in range(int(2*math.pi/ang_rad)):
        for k in range(1,radius+1):

            pt_new_y=k*math.cos(i*ang_rad)
            pt_new_x=k*math.sin(i*ang_rad)
            pt_new=[pt_new_x,pt_new_y]
            pt_new_lis.append(pt_new)

    return pt_new_lis


# create potential location use bounding box
# pt= original point, radius=buffer radius, density= point distance(default=1 meter)
def bbox(pt,radius,density):
    pt_new_lis=[]
    # start from left-bottom corner
    start_point=[pt[0]-radius,pt[1]-radius]
    # create points within a rectangle
    for c in range(2*radius):
        for l in range(2*radius):
            pt_new_x=start_point[0]+l*density
            pt_new_y=start_point[1]+c*density
            pt_new = [pt_new_x, pt_new_y]
            pt_new_lis.append(pt_new)

    # remove point outside circle
    pt_new_lis2=[]
    for i in pt_new_lis:
        i=Point(i)
        if i.distance(Point(pt)) <= radius:
            pt_new_lis2.append(i.coords[0])
    #print(pt_new_lis2)
    return pt_new_lis2

def findtile(pt):
    x,y=pt[0],pt[1]
    if 76702.79<x<80201.65 and 455789.08<y<459614.32:
        return "1"
    elif 79782.22<x<84455.64 and 455766.15<y<460716.51:
        return "2"
    elif 83990.42<x<88071.27 and 455903.22<y<460007.70:
        return "3"
    elif 74086.05<x<76060.23 and 451983.07<y<455051.58:
        return "4"
    elif 75965.74<x<80078.49 and 451992.19<y<456048.60:
        return "5"
    elif 79914.23<x<84010.36 and 451925.74<y<456032.24:
        return "6"
    elif 83987.27<x<88003.98 and 451996.32<y<456048.42:
        return "7"
    elif 87996.72<x<92119.69 and 451992.42<y<456035.79:
        return "8"
    elif 71878.66<x<76168.23 and 447728.54<y<452016.89:
        return "9"
    elif 75722.74<x<80071.32 and 447643.33<y<452082.66:
        return "10"
    elif 79938.59<x<84022.68 and 447980.12<y<452043.65:
        return "11"
    else:
        return "12"

# Find nearest 10 building/traffic area around the select point, pt=point
def findNearestObject(pt,object):
    tileid=findtile(pt)
    nearestbuilding_lis=[]
    # open .shp file and find the coordinates of polygons
    with fiona.open('/Users/xuyixin/Desktop/msc thesis/svf-3d/bgt_tile/'+str(object)+'_'+tileid+'.shp', 'r') as src:
        dict_distance={}
        for i in range(len(src)):
            buildings = src[i]['geometry']['coordinates']
            id=src[i]['id']
            buildings_coordinates=(buildings[0])
            # calculate average(all vertex) distance of one polygon to the pt
            x_sum,y_sum=0,0
            for k in buildings_coordinates:
                # if data format is wrong then omit the building/road and continue
                if len(k)>3:
                    continue
                x,y=k[0],k[1]
                x_sum=x_sum+x
                y_sum=y_sum+y

            x_avg=x_sum/len(buildings_coordinates)
            y_avg=y_sum/len(buildings_coordinates)
            buildings_avg=[id,Point(x_avg,y_avg)]
            # use dictionary to store building/road id(id could be differ when select BGT layer) and distance
            distance=Point(pt).distance(buildings_avg[1])
            dict_distance[id]=distance

        # sort all buildings by avg distance and select first 30
        nearestbuilding_id=sorted(dict_distance.items(),key=lambda item:item[1])[0:30]

        # use id trace back to the polygon
        for j in nearestbuilding_id:
            for p in range(len(src)):
                if src[p].get('id') == j[0]:
                    nearestbuilding_lis.append(src[p])
                    continue

    return nearestbuilding_lis

# remove point outside the polygon
def CleanPointInsidePolygon(pt):
    # create scatter points list and buildings list
    ptlis=bbox(pt,15,1)
    buildinglis=findNearestObject(pt,"building")
    roadlis=findNearestObject(pt,"road")

    # point within building
    new_ptlis=[]

    # point outside building
    new_ptlis2=[]

    # point within traffic area
    new_ptlis3=[]

    # point outside building and traffic area
    new_ptlis4=[]

    # find points within building
    for p in ptlis:
        p=Point(p)
        for b in buildinglis:
            poly=Polygon(b['geometry']['coordinates'][0])
            # function "within" is much quicker than "disjoint"&"contain"
            if p.within(poly):
                new_ptlis.append(p.coords[0])

    # remove points within buildings
    for i in ptlis:
        if not i in new_ptlis:
            new_ptlis2.append(i)

    # find points within traffic area
    for h in new_ptlis2:
        h=Point(h)
        for t in roadlis:
            poly2=Polygon(t['geometry']['coordinates'][0])
            if h.within(poly2):
                new_ptlis3.append(h.coords[0])

    # remove points within traffic area
    for m in new_ptlis2:
        if not m in new_ptlis3:
            new_ptlis4.append(m)
    # expand search area
    if len(new_ptlis4)>0:
        return new_ptlis4

    else:
        ptlis = bbox(pt, 30, 1)
        buildinglis = findNearestObject(pt, "building")
        roadlis = findNearestObject(pt, "road")

        # point within building
        new_ptlis = []

        # point outside building
        new_ptlis2 = []

        # point within traffic area
        new_ptlis3 = []

        # point outside building and traffic area
        new_ptlis4 = []

        # find points within building
        for p in ptlis:
            p = Point(p)
            for b in buildinglis:
                poly = Polygon(b['geometry']['coordinates'][0])
                # function "within" is much quicker than "disjoint"&"contain"
                if p.within(poly):
                    new_ptlis.append(p.coords[0])

        # remove points within buildings
        for i in ptlis:
            if not i in new_ptlis:
                new_ptlis2.append(i)

        # find points within traffic area
        for h in new_ptlis2:
            h = Point(h)
            for t in roadlis:
                poly2 = Polygon(t['geometry']['coordinates'][0])
                if h.within(poly2):
                    new_ptlis3.append(h.coords[0])

        # remove points within traffic area
        for m in new_ptlis2:
            if not m in new_ptlis3:
                new_ptlis4.append(m)
        # expand search area
        #print("next new ptlist2",new_ptlis4)
        return new_ptlis4



# write points to txt file
def WritePoint2txt():
    a = CleanPointInsidePolygon([78608, 452805])
    fo = open("foo.txt", "w")
    fo.write("x,y \n")
    for i in a:
        t=str(i[0])+","+str(i[1])+"\n"
        fo.write(t)
    fo.close()



#WritePoint2txt()
