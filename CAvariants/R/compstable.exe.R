compstable.exe <-
function(Z){
tZZ <- t(Z) %*% Z
ZtZ <- Z %*% t(Z)
factor <- sum(diag(tZZ))
compsR <- matrix(0, nrow = 5, ncol = 2)
compsC <- matrix(0, nrow = 5, ncol = 2)
###############################
#   row  Category 1 #
###############################

compsR[1, 1] <- ZtZ[1, 1]#Tube Location Component for Category 2#
compsR[1, 2] <- 1 - pchisq(compsR[1, 1], ncol(Z))#P-value of the Location Comp for Category 2#
compsR[2, 1] <- ZtZ[2, 2]#Dispersion Component for Category 2#
compsR[2, 2] <- 1 - pchisq(compsR[2, 1], ncol(Z))#P-value of Dispersion Comp of Category 2#
compsR[3, 1] <- ZtZ[3, 3]#Cubic Component for Category 1#
compsR[3, 2] <- 1 - pchisq(compsR[3, 1], ncol(Z))#P-value for the Cubic Comp of Category 2#
compsR[4, 1] <- factor - (compsR[1, 1] + compsR[2, 1]+compsR[3,1])#Error of Row Components for Category 1#
if(ncol(Z) > 3) {
compsR[4, 2] <- 1 - pchisq(compsR[4, 1], (nrow(Z)-3) * ncol(Z) )#P-value for the Errors of Category 1#
}
else {
compsR[4, 2] <- 0
}
compsR[5, 1] <- factor
compsR[5, 2] <- 1 - pchisq(compsR[5, 1], nrow(Z) * ncol(Z))

###############################
#  column     Category 2 #
###############################

compsC[1, 1] <- tZZ[1, 1]# Location Component for Category 1#
compsC[1, 2] <- 1 - pchisq(compsC[2, 1], nrow(Z))#P-Value for the Location Comp of Category 1#
compsC[2, 1] <- tZZ[2, 2]#Dispersion Component for Category 1#
compsC[2, 2] <- 1 - pchisq(compsC[2, 1], nrow(Z))#P-value for the Dispersion Comp of Category 2#
compsC[3, 1] <- tZZ[3, 3]#Cubic Component for Category 1#
compsC[3, 2] <- 1 - pchisq(compsC[3, 1], nrow(Z))#P-value for the Cubic Comp of Category 2#
compsC[4, 1] <- factor - (compsC[1, 1] + compsC[2, 1]+compsC[3,1])#Error of Row Components for Category 1#
if(ncol(Z) > 3) {
compsC[4, 2] <- 1 - pchisq(compsC[4, 1], nrow(Z) * (ncol(Z) - 3))#P-value for the Errors of Category 1#
}
else {
compsC[4, 2] <- 0
}
compsC[5, 1] <- factor
compsC[5, 2] <- 1 - pchisq(compsC[5, 1], nrow(Z) * ncol(Z))

list(compsR=compsR,compsC=compsC)
#return(comps)

}
