import os
import lsa_5days
import datetime
starttime = datetime.datetime.now()

os.chdir('/Users/xuyixin/Desktop/msc thesis/svf-3d/each_sensor')
f=open("all_sensor_1_2.csv",'w')
f.write("id,x,y,z,similarity\n")

dirs = os.listdir()
#print(dirs)
dirs_new=[]
for filename in dirs:
    filename_new=filename.split(".")[0]
    dirs_new.append(filename_new)


for new_name in dirs_new[79:]:
    if new_name != "all_sensor" and new_name != "all_sensor_1" and new_name != "all_sensor_1_2"and len(new_name)>3 :
        print("next new name",new_name)
        starttime_sensor = datetime.datetime.now()
        sensor_3d_simi=lsa_5days.writefile(new_name+".csv")
        endtime_sensor = datetime.datetime.now()
        print("sensor time",endtime_sensor - starttime_sensor)
        if sensor_3d_simi != 0:

            x=str(sensor_3d_simi[0][0][0])
            y=str(sensor_3d_simi[0][0][1])
            z=str(sensor_3d_simi[0][0][2])
            simi=str(sensor_3d_simi[0][1])
            f.write(new_name+","+x+","+y+","+z+","+simi+"\n")

endtime = datetime.datetime.now()
print (endtime - starttime)