source("../R/parKml.r")

cat("\n####################################################################
########################### Test ParKml ############################
####################################################################\n")

cleanProg(.ParKml.validity)
cleanProg(parKml,,,2) # DISTANCE_METHODS meanNA
cleanProg(.ParKml.show)
new("ParKml")
parKml()
parKml(saveFreq=300)
pk <- parKml(distanceName="maximum")
pk2 <- parKml(distanceName="manhattan",saveFreq=30,maxIt=50,imputationMethod="LI-LOCBF")
pk3 <- parKml(distanceName="minkowski",power=3,saveFreq=30,maxIt=50,imputationMethod="LI-LOCBF",centerMethod=medianNA,startingCond="randomAll")
pk['distance']


x <- matrix(c(1,2,3,2,3,4),2,byrow=TRUE)
y <- matrix(c(0,0,0,3,3,3),2,byrow=TRUE)
pk['distance'](x,y)
pk2['distance'](x,y)
pk3['distance'](x,y)




cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++ Fin Test ParKml +++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
