source("../R/clusterLongData3d.r")

cat("
  ####################################################################
 ######################## Test ClusterLongData ######################
####################################################################\n")


data <- data.frame(id=1:5,A11=c(12,13,12,13,12),A13=c(13,13,15,13,15),A14=c(15,12,15,12,15),
                   P11=c(45,46,46,46,46),P12=c(47,48,49,48,49),P14=c(49,49,48,49,49))

detectGlobal(ClusterLongData3d_validity)

### Constructeurs
detectGlobal(clusterLongData3d,1) # all

new("ClusterLongData")

clusterLongData3d(
    traj=tr1,
    idAll=c(100,102,104),
    time=c(1,2,4,8,16),
    varNames=c("P","A"),
    maxNA=3
    )

new("ClusterLongData3d",
    traj=tr2,
    idFewNA=c("i1","i2","i3","i4"),
    idAll=c("i1","i2","i3","i4"),
    time=c(1,2,4),
    varNames=c("P","A","E","R"),
    maxNA=2
    )

new("ClusterLongData3d",
    traj=tr3n,
    idAll=c("i1","i2","i3","i4","i5","i6","i7","i8"),
    idFewNA=c("i1","i2","i3","i4","i5","i6","i7","i8"),
    time=c(1,2,4),
    varNames=c("P","A"),
    maxNA=2
    )


new("ClusterLongData3d",
    traj=tr4n,
    idAll=as.character(1:200),
    idFewNA=as.character(1:200),
    time=c(1,2,4,5,8,10),
    varNames=c("P","A"),
    maxNA=c(1,1)
    )



clusterLongData3d()
#longData(traj=array(c(1,2,3,1,4,6,1,8,10),dim=c(3,3)))

clusterLongData3d(traj=tr1,idAll=as.character(c(101,102,104)),time=c(1,2,4,8,16),varNames=c("P","A"),maxNA=3)
cld3d(traj=tr2,idAll=as.character(c(1,2,3,4)),time=c(1,2,4),varNames=c("P","A","E","R"),maxNA=2)
cld3d(traj=tr3n,idAll=as.character(c(1,2,3,4,5,6,7,8)),time=c(1,2,4),varNames=c("P","A"),maxNA=2)
cld3d(traj=tr4n,idAll=101:300,time=c(1,2,4,5,8,9),varNames=c("Az","er"),maxNA=2)

### Vérification de l'exclusion des manquantes
clusterLongData3d(traj=tr4n,idAll=101:300,time=1:6,varNames=c("P","A"),maxNA=1)
clusterLongData3d(traj=tr4n,idAll=1:200,time=1:6,varNames=c("P","A"),maxNA=2)
clusterLongData3d(traj=tr4n,idAll=1:200,time=1:6,varNames=c("P","A"),maxNA=c(1,1))
clusterLongData3d(traj=tr4n,idAll=1:200,time=1:6,varNames=c("P","A"),maxNA=c(2,2))
clusterLongData3d(traj=tr4n,idAll=1:200,time=1:6,varNames=c("P","A"),maxNA=c(2,1))




### Base de données
detectGlobal(ClusterLongData3d_show)

CLD0 <- clusterLongData3d()
CLD1 <- clusterLongData3d(traj=tr1,idAll=c(101,102,104),time=c(1,2,4,8,16),varNames=c("Pa","Av"),maxNA=3)
CLD1n <- clusterLongData3d(traj=tr1n,idAll=c(101,102,104),time=c(1,2,4,8,16),varNames=c("Pa","Av"),maxNA=3)

CLD2 <- clusterLongData3d(data,timeInData=list(A=c(2,4),P=c(5,7)),time=2:3,varNames=c("Av","Pe"))
CLD2 <- clusterLongData3d(data,timeInData=list(A=c(2,4),P=c(5,7)),time=2:3)
CLD2 <- clusterLongData3d(data,timeInData=list(c(2,4),c(5,7)),time=2:3)
CLD2 <- clusterLongData3d(data,timeInData=list(A=c(2,4),P=c(5,7)),time=2:3,varNames=c("Av","Pe"))
CLD2n <- clusterLongData3d(data,timeInData=list(A=c(2,NA,4),P=c(NA,5,7)),time=2:4)
CLD2 <- clusterLongData3d(data,timeInData=list(V21=c(2,3,4),V4=c(5,6,7)),time=c(11,13,14))

CLD3 <- clusterLongData3d(tr3n,maxNA=2)
CLD3n <- clusterLongData3d(tr3n,maxNA=2)

CLD4 <- clusterLongData3d(tr4)
CLD4n <- clusterLongData3d(tr4n,maxNA=3)

CLD5 <- clusterLongData3d(tra5)
CLD5n <- clusterLongData3d(tra5n,maxNA=4)

CLD6n <- clusterLongData3d(dn6,time=1:6,timeInData=list(cred=3:8,creq=9:14,croq=c(24:28,NA)),maxNA=5)

CLD7 <- clusterLongData3d(tra7)
CLD7n <- clusterLongData3d(tra7n)



cat("
  ###################################################################
 ############################ Get & Set ############################
###################################################################\n")


detectGlobal(ClusterLongData3d_get,2) # CLUSTER_NAMES CRITERION_NAMES
detectGlobal(ClusterLongData3d_set,2) # CLUSTER_NAMES CRITERION_NAMES

############# Set #############
CLD3['add'] <- p3a
CLD3['add'] <- p3b
CLD3['add'] <- p3c
CLD3['add'] <- p3d
CLD3['add'] <- p3e
CLD3['add'] <- p3f

CLD7['add'] <- p4a
CLD7['add'] <- p4b
CLD7['add'] <- p4b
CLD7['add'] <- p4c
CLD7['add'] <- p4d
CLD7['add'] <- p4e
CLD7['add'] <- p4f



############# Get #############
CLD3["varNames"]
CLD3n["idAll"]
CLD3n["idFewNA"]
CLD6n["idAll"]
CLD2n["idFewNA"]
tryBug(CLD3[2])

getClusters(CLD7,3,1,FALSE)
getClusters(CLD7,3,1,TRUE)
getClusters(CLD7,3,2)




cat("
  #####################################################################
 ############################### plot3d ##############################
#####################################################################")

detectGlobal(ClusterLongData3d_num_plot3d)
detectGlobal(ClusterLongData3d_missingY_plot3d)
plot3d(CLD7)
plot3d(CLD7,2)
plot3d(CLD7,c(3,1))
plot3d(CLD7,c(2,1))
plot3d(CLD7,c(3,1))
plot3d(CLD7,c(3,2))

detectGlobal(ClusterLongData3d_plotTrajMeans,1)
plotTraj3d(CLD7,3)
plotTraj3d(CLD7,c(3,1))
plotTraj3d(CLD7,c(4,1))
plotTraj3d(CLD7,c(3,1))
plotTraj3d(CLD7,c(3,2))


detectGlobal(ClusterLongData3d_num_plotMeans3d)
plotMeans3d(CLD7,3)
plotMeans3d(CLD7,c(3,1))
plotMeans3d(CLD7,c(4,1))
plotMeans3d(CLD7,c(3,1))
plotMeans3d(CLD7,c(3,2))




cat("
  #####################################################################
 ############################## plot3dPdf ############################
#####################################################################")

detectGlobal(ClusterLongData3d_num_plot3dPdf)
plot3dPdf(CLD7,2)
plot3dPdf(CLD7,c(3,1))
plot3dPdf(CLD7,c(2,1))
plot3dPdf(CLD7,c(3,1))
plot3dPdf(CLD7,c(3,2))




cat("
  #####################################################################
 ################################ plot ###############################
#####################################################################")

detectGlobal(ClusterLongData3d_num_plot3d)
detectGlobal(ClusterLongData3d_missingY_plot3d)
ClusterLongData3d_plotTrajMeans(CLD7,2)
close.screen(1:3)
ClusterLongData3d_plotTrajMeans(CLD7,2,ajustLegend=-0.2)
close.screen(1:3)
ClusterLongData3d_plotTrajMeans(CLD7,c(3,1))
close.screen(1:3)
ClusterLongData3d_plotTrajMeans(CLD7,c(2,1))
close.screen(1:3)
ClusterLongData3d_plotTrajMeans(CLD7,c(3,1),parMean=parMEAN(type="n"))
close.screen(1:3)
ClusterLongData3d_plotTrajMeans(CLD7,c(3,2))
close.screen(1:3)

detectGlobal(ClusterLongData3d_plotTrajMeans,1)
ClusterLongData3d_plot(CLD7,3)
ClusterLongData3d_plot(CLD3,c(3,1))
ClusterLongData3d_plot(CLD7,c(4,1))
ClusterLongData3d_plot(CLD7,c(3,1))
ClusterLongData3d_plot(CLD3,c(3,2))



detectGlobal(ClusterLongData3d_plotTrajMeans,1)
plot(CLD3)
plot(CLD7,c(4,1))
plot(CLD7,c(3,1),ajustLegend=-0.2)
plot(CLD3,c(3,2))
plot(CLD7,c(4,1),toPlot="traj")
plot(CLD7,c(3,1),toPlot="criterion")
plot(CLD3,c(3,2),toPlot="both")


detectGlobal(ClusterLongData3d_plotTraj,1)
plotTraj(CLD7,c(4,1),ajustLegend=-0.2)
plotTraj(CLD7,c(3,1))
plotTraj(CLD3,c(3,2))

detectGlobal(ClusterLongData3d_plotTraj,1)
plotMeans(CLD7,c(4,1),ajustLegend=-0.2)
plotMeans(CLD7,c(3,1))
plotMeans(CLD3,c(3,2))


cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ Fin Test ClusterLongData +++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
