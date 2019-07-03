import polygon_test
import svfscript
import svfscript_50
import svfscript_100
import svfscript_150
import svfscript_200
#import matplotlib.pyplot as plt
import os
from osgeo.osr import SpatialReference, CoordinateTransformation


def wgs84To28992(lo,la):
    # Define the Rijksdriehoek projection system (EPSG 28992)
    epsg28992 = SpatialReference()
    epsg28992.ImportFromEPSG(28992)

    # correct the towgs84
    epsg28992.SetTOWGS84(565.237, 50.0087, 465.658, -0.406857, 0.350733, -1.87035, 4.0812)

    # Define the wgs84 system (EPSG 4326)
    epsg4326 = SpatialReference()
    epsg4326.ImportFromEPSG(4326)

    rd2latlon = CoordinateTransformation(epsg28992, epsg4326)
    latlon2rd = CoordinateTransformation(epsg4326, epsg28992)

    # Check the transformation for a point close to the centre of the projected grid
    xy = latlon2rd.TransformPoint(lo,la)
    return ([xy[0],xy[1]])

def readfile(filename):
    os.chdir('/Users/xuyixin/Desktop/msc thesis/svf-3d/each_sensor')
    fo=open(filename,"r")
    metadata=fo.readline()
    mac=metadata.split(",")[0]
    la=float(metadata.split(",")[1])
    lo=float(metadata.split(",")[2])
    coor=wgs84To28992(lo,la)

    all_incre_index_lis=[]
    for line in fo.readlines()[1:]:
        temp_lis=[]

        flat_lis=[]
        flat_index_lis=[]

        incre_lis=[]
        incre_index_lis=[]

        drop_lis=[]
        drop_index_lis=[]

        unclassify_index_lis=[]

        date=line.split(",")[0]
        for i in line.split(",")[1:]:
            temp=float(i)
            temp_lis.append(temp)

        #only select temp data between sen rise and sun set, sunrise: 5:30（22）; sun set: 9:45（87）
        for k in range(22,88):
            """
            if abs(temp_lis[k-1] - temp_lis[k + 1]) <= 0.07 :
                flat_lis.append(temp_lis[k])
                flat_index_lis.append(k)
            """
            if (temp_lis[k+1] - temp_lis[k-1]) >= 0.1 :
                incre_lis.append(temp_lis[k])
                incre_index_lis.append(k)
            elif (- temp_lis[k+1] + temp_lis[k-1]) >= 0.1 or (- temp_lis[k+2] + temp_lis[k+1]) >= 0.05 or (- temp_lis[k] + temp_lis[k-2]) >= 0.1:
                drop_index_lis.append(k)
                drop_lis.append(temp_lis[k])
            else:
                flat_lis.append(temp_lis[k])
                flat_index_lis.append(k)


        #plt.scatter(incre_index_lis,incre_lis,color="r")
        #plt.scatter(drop_index_lis,drop_lis,color="b")
        #plt.scatter(flat_index_lis,flat_lis,color='g')
        # plt.show()

        all_incre_index_lis.append(incre_index_lis)
    #print("next all incre index lis")
    #print(all_incre_index_lis)
    return [all_incre_index_lis,mac,coor]




# only return time has higher possibility influenced by sun
def incre_possibility(incre_lis):
    #print(123123123)
    poss={}
    k=5.25
    for i in range(22,88):
        k = k + 0.25
        p=0

        for t in incre_lis:
            #print(222222)
            #print(incre_lis,t)
            #print(i)
            if i in t:
                p = p + 1
        ps = p / len(incre_lis)
        #print("next ps")
        #print(ps)
        if ps > 0.5:
            poss[k] = ps
    #print("next poss")
    #print(poss)
    return poss


#compare result from sky view factor and sensor then calculate possibility
def possibility(dic_skyview,dic_sensor):
    #print(dic_skyview,dic_sensor)
    # calculate possibility
    poss_all=0
    for key in dic_sensor:
        if dic_skyview[key]== "not influenced" :
            poss=0
            poss_all=poss_all+poss
        elif dic_skyview[key]== "influenced" :
            poss=dic_sensor[key]
            poss_all = poss_all + poss
    # max radio: sensor can detect higher temperature during one period
    #sensor_max = (dic_sensor[max(dic_sensor, key=dic_sensor.get)])
    #print("next poss_all")
    #print(poss_all)
    if "influenced" not in dic_skyview.values() or len(dic_sensor)==0:
        a=0
    elif len(dic_sensor)!=0:
        # normolize possibility
        a=(poss_all/(len(dic_sensor)))

    return a

# a = sorted scatter points with possibility for one sensor b = increase possibility
def ptHeight(a,b):
    pt_3d = {}
    poss_max=(a[-5:])
    for i in poss_max:
        coor=i[0]
        simi_0=i[1]

        # svf at different height
        svf_50=svfscript_50.xy(coor)
        svf_100=svfscript_100.xy(coor)
        svf_150 = svfscript_150.xy(coor)
        svf_200 = svfscript_200.xy(coor)
        simi_50=possibility(svf_50,b)
        simi_100 = possibility(svf_100, b)
        simi_150 = possibility(svf_150, b)
        simi_200 = possibility(svf_200, b)

        simi_lis=[simi_0,simi_50,simi_100,simi_150,simi_200]
        #print(simi_lis,max(simi_lis),simi_lis.index(max(simi_lis)))
        height=simi_lis.index(max(simi_lis))*0.5
        pt_3d_key=(coor[0],coor[1],height)
        pt_3d[pt_3d_key]=max(simi_lis)
    #print(pt_3d)
    pt_3d_poss=(sorted(pt_3d.items(), key=lambda item:item[1])[-1:])
    print(pt_3d_poss)
    return pt_3d_poss

#poss=possibility(test3,incre_possibility(all_incre_index_lis))
#print(poss)


def writefile(filename):
    metadata=readfile(filename)
    #print(metadata)
    mac=metadata[1]
    sensor_loca=metadata[2]
    dict_poss={}
    pts=polygon_test.CleanPointInsidePolygon(sensor_loca)
    #print("next pts")
    #print(pts)

    all_incre_index_lis=incre_possibility(metadata[0])

    for pt in pts:

        poss=possibility(svfscript.xy(pt),all_incre_index_lis)
        #print(poss)
        dict_poss[pt]=poss
    #print("next dic_poss")
    #print(dict_poss)
    if len(dict_poss)>0:
        #poss_max=max(dict_poss, key=dict_poss.get)
        #point_max=(dict_poss[max(dict_poss, key=dict_poss.get)])
        #print(poss_max,point_max)
        #print(sorted(dict_poss.items(), key=lambda item:item[1]))

        pt_3d_simi=ptHeight(sorted(dict_poss.items(), key=lambda item:item[1]),all_incre_index_lis)

        os.chdir('/Users/xuyixin/Desktop/msc thesis/svf-3d/each_sensor_result')
        f2=open(mac+"_interpolation.csv",'w')
        f2.write("x,"+"y,"+"possibility"+"\n")
        for key in dict_poss:
            f2.write(str(key[0])+","+str(key[1])+","+str(dict_poss[key])+'\n')

        return pt_3d_simi

    else:
        return 0
#writefile("dba4.csv")