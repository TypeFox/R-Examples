cat("############################################################
######################### testKmL3d ########################
############################################################")

source("../R/kml3d.r")

detectGlobal(calculTrajMean3d,1) # tapply
cent2A <- calculTrajMean3d(CLD2["traj"],p2a['clusters'])
cent2B <- calculTrajMean3d(CLD2["traj"],p2b['clusters'])
cent2C <- calculTrajMean3d(CLD2["traj"],p2c['clusters'],medianNA)

cent3A <- calculTrajMean3d(CLD3["traj"],p3a['clusters'])
cent3B <- calculTrajMean3d(CLD3["traj"],p3b['clusters'])
cent3C <- calculTrajMean3d(CLD3["traj"],p3c['clusters'])
cent3D <- calculTrajMean3d(CLD3["traj"],p3d['clusters'])
cent3E <- calculTrajMean3d(CLD3["traj"],p3e['clusters'])
cent3F <- calculTrajMean3d(CLD3["traj"],p3f['clusters'])

detectGlobal(affectIndiv3d,1) # dist3d (dans les arguments)
aC <- affectIndiv3d(CLD2["traj"],cent2A)
bC <- affectIndiv3d(CLD2["traj"],cent2B)
cC <- affectIndiv3d(CLD2["traj"],cent2C)


detectGlobal(kml3dSlow)
partInit <- initializePartition(3,8,method="randomAll")
system.time(kml3dSlow(CLD3['traj'],partInit))
partInit <- initializePartition(6,200,method="randomK")
system.time(kml3dSlow(CLD4['traj'],partInit))
#partInit <- initializePartition(6,2000,method="randomK")
#system.time(kml3dSlow(CLD5['traj'],partInit))


detectGlobal(kml3dFast)
partInit <- initializePartition(3,8,method="randomK")
system.time(kml3dFast(CLD3['traj'],partInit))

partInit <- initializePartition(6,200,method="randomK")
system.time(kml3dFast(CLD4['traj'],partInit))

partInit <- initializePartition(3,8,method="randomK")
system.time(a <- kml3dFast(CLD3['traj'],partInit))
system.time(b <- kml3dSlow(CLD3['traj'],partInit))

detectGlobal(kml3d)
kml3d(CLD3,toPlot="both")
choice(CLD3)

cat("------------------------------------------------------------
-------------------------- fin test ------------------------
---------------------------- KmL3d -------------------------
------------------------------------------------------------")

