source("../R/kml.r")
dyn.load("../src-i386/kml")
cleanProg(.partitionInitialise,,,1) # min
point <- matrix(c(0,0, 0,1, -1,0, 0,-1, 1,0),5,byrow=TRUE)
points <- rbind(point,t(t(point)+c(10,0)),t(t(point)+c(5,6)))
points <- rbind(points,t(t(points)+c(30,0)),t(t(points)+c(15,20)),t(-t(point)+c(20,10)))

paInit <- partitionInitialise(4,nrow(points),method="randomK")
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=na.omit(paInit["clustersAsInteger"])+1,type="p",cex=1,lwd=2,pch=16)

paInit <- partitionInitialise(4,nrow(points),method="randomAll")
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=paInit["clustersAsInteger"]+1,type="p",cex=1,lwd=2,pch=16)

paInit <- partitionInitialise(4,nrow(points),method="maxDist",points)
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=na.omit(paInit["clustersAsInteger"])+1,type="p",cex=1,lwd=2,pch=16)

paInit <- partitionInitialise(4,nrow(points),method="kmeans++",points)
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=na.omit(paInit["clustersAsInteger"])+1,type="p",cex=1,lwd=2,pch=16)

paInit <- partitionInitialise(4,nrow(points),method="kmeans--",points)
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=na.omit(paInit["clustersAsInteger"])+1,type="p",cex=1,lwd=2,pch=16)

paInit <- partitionInitialise(4,nrow(points),method="kmeans+",points)
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=na.omit(paInit["clustersAsInteger"])+1,type="p",cex=1,lwd=2,pch=16)

paInit <- partitionInitialise(4,nrow(points),method="kmeans-",points)
plot(points)
lines(points[!is.na(paInit["clusters"]),],col=na.omit(paInit["clustersAsInteger"])+1,type="p",cex=1,lwd=2,pch=16)




pa <- partitionInitialise(3,200,method="randomK")
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="randomAll")
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="maxDist",ld4["traj"])
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="kmeans++",ld4["traj"])
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="kmeans+",ld4["traj"])
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="kmeans--",ld4["traj"])
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="kmeans-",ld4["traj"])
plot(ld4,pa)

pa <- partitionInitialise(3,200,method="maxDist",LD4["traj"])
plot(LD4,pa)

pa <- partitionInitialise(3,200,method="kmeans++",LD4["traj"])
plot(LD4,pa)

pa <- partitionInitialise(3,200,method="kmeans+",LD4["traj"])
plot(LD4,pa)

pa <- partitionInitialise(3,200,method="kmeans--",LD4["traj"])
plot(LD4,pa)

pa <- partitionInitialise(3,200,method="kmeans-",LD4["traj"])
plot(LD4,pa)





cleanProg(calculTrajMean,,,1) # tapply
cent2a <- calculTrajMean(ld2["traj"],p2a['clusters'])
cent2b <- calculTrajMean(ld2["traj"],p2b['clusters'])
cent2c <- calculTrajMean(ld2["traj"],p2c['clusters'],medianNA)

cent3a <- calculTrajMean(ld3["traj"],p3a['clusters'])
cent3b <- calculTrajMean(ld3["traj"],p3b['clusters'])
cent3c <- calculTrajMean(ld3["traj"],p3c['clusters'])
cent3d <- calculTrajMean(ld3["traj"],p3d['clusters'])
cent3e <- calculTrajMean(ld3["traj"],p3e['clusters'])
cent3f <- calculTrajMean(ld3["traj"],p3f['clusters'])

cleanProg(calculTrajMean,,,1) # tapply
cent2A <- calculTrajMean(LD2["traj"],p2a['clusters'])
cent2B <- calculTrajMean(LD2["traj"],p2b['clusters'])
cent2C <- calculTrajMean(LD2["traj"],p2c['clusters'],medianNA)

cent3A <- calculTrajMean(LD3["traj"],p3a['clusters'])
cent3B <- calculTrajMean(LD3["traj"],p3b['clusters'])
cent3C <- calculTrajMean(LD3["traj"],p3c['clusters'])
cent3D <- calculTrajMean(LD3["traj"],p3d['clusters'])
cent3E <- calculTrajMean(LD3["traj"],p3e['clusters'])
cent3F <- calculTrajMean(LD3["traj"],p3f['clusters'])

cleanProg(affectIndiv,,,1) # dist3d (dans les arguments)
aC <- affectIndiv(ld2["traj"],cent2a)
bC <- affectIndiv(ld2["traj"],cent2b)
cC <- affectIndiv(ld2["traj"],cent2c)


cleanProg(kmlSlow)
partInit <- partitionInitialise(3,8,method="randomK")
system.time(kmlSlow(ld3['traj'],partInit))
partInit <- partitionInitialise(6,200,method="randomK")
system.time(kmlSlow(ld4['traj'],partInit))
#partInit <- partitionInitialise(6,2000,method="randomK")
#system.time(kmlSlow(ld5['traj'],partInit))

partInit <- partitionInitialise(3,8,method="randomK")
system.time(kmlSlow(LD3['traj'],partInit))
partInit <- partitionInitialise(6,200,method="randomK")
system.time(kmlSlow(LD4['traj'],partInit))
#partInit <- partitionInitialise(6,2000,method="randomK")
#system.time(kmlSlow(LD5['traj'],partInit))


cleanProg(kmlFast)
partInit <- partitionInitialise(3,8,method="randomK")
system.time(kmlFast(ld3['traj'],partInit))

partInit <- partitionInitialise(6,200,method="randomK")
system.time(kmlFast(ld4['traj'],partInit))

partInit <- partitionInitialise(6,2000,method="randomK")
system.time(a <- kmlFast(ld5['traj'],partInit))
system.time(b <- kmlSlow(ld5['traj'],partInit))
identical(a,b)

partInit <- partitionInitialise(3,8,method="randomK")
system.time(a <- kmlFast(LD3['traj'],partInit))
system.time(b <- kmlSlow(LD3['traj'],partInit))

partInit <- partitionInitialise(6,200,method="randomK")
system.time(kmlFast(LD4['traj'],partInit))
#partInit <- partitionInitialise(6,2000,method="randomK")
#system.time(kmlFast(LD5['traj'],partInit))


cleanProg(expandStartingCond)
expandStartingCond(startingCond="all",10,"")
expandStartingCond(startingCond="all",10,"maxDist")
expandStartingCond(startingCond="all",10,"kmeans-")
expandStartingCond(startingCond="all",10,c("maxDist","kmeans-"))

expandStartingCond(startingCond="nearlyAll",10,"")
expandStartingCond(startingCond="nearlyAll",10,"maxDist")
expandStartingCond(startingCond="nearlyAll",10,"kmeans-")
expandStartingCond(startingCond="nearlyAll",10,c("maxDist","kmeans-"))

expandStartingCond(startingCond="kmeans-",10,"")
expandStartingCond(startingCond=c("kmeans-","randomK"),10,"randomK")

cleanProg(cutScreen)
cleanProg(fastOrSlow,,,1) #DISTANCE_METHODS
cleanProg(.clusterLongData.kml)

cld3 <- as.clusterLongData(dn3,time=1:6,timeDataFrame=list(cred=3:8,creq=9:14,croq=c(24:28,NA)))
kml(cld3,toPlot="both")
choice(cld3)
cleanProg(choiceChangeParam,,,1) # choiceChangeParam


my <- gald(
           nbEachClusters=c(20,20,20,20,20),
           functionClusters=list(
           function(t){return(c(0,0))},
           function(t){return(c(10,10))},
           function(t){return(c(0,10))},
           function(t){return(c(10-t,10-t))},
           function(t){return(c(10,10-t))}
           )
           )
kml(my)

my2 <- gald(
           nbEachClusters=c(150,150,150,150,150),
           functionClusters=list(
           function(t){return(c(0,0))},
           function(t){return(c(10,10))},
           function(t){return(c(0,10))},
           function(t){return(c(10-t,10-t))},
           function(t){return(c(10,10-t))}
           )
           )
kml(my2)


mi <- gald(
           nbEachClusters=c(20,20,20,20,20),
           functionClusters=list(
           function(t){return(c(0,0))},
           function(t){return(c(10,10))},
           function(t){return(c(0,10))},
           function(t){return(c(t,t))},
           function(t){return(c(10,10-t))}
           )
           )
kml(mi)

mi2 <- gald(
           nbEachClusters=c(150,30,40,25,15),
           functionClusters=list(
           function(t){return(c(0,0))},
           function(t){return(c(10,10))},
           function(t){return(c(0,10))},
           function(t){return(c(t,t))},
           function(t){return(c(10,10-t))}
           )
           )
kml(mi2)


ma <- gald(time=5:15,
           nbEachClusters=c(30,30,30),
           functionClusters=list(
           function(t){return(c(0,dnorm(10-t)*30))},
           function(t){return(c(dnorm(10-t)*30,0))},
           function(t){return(c(dnorm(10-t)*30,dnorm(10-t)*30))}
           )
           )
pa <- partition(rep(1:3,each=30))
plot(ma,pa)
kml(ma)




