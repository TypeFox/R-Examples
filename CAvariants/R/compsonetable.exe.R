compsonetable.exe <- 
function(Z){
tZZ <- t(Z) %*% (Z)
factor <- sum(diag(tZZ))
comps <- matrix(0, nrow=5, ncol=2)
###############################
#     Column Category         #
###############################
comps[1, 1] <- tZZ[1, 1]#Location Component for Category 2#
comps[1, 2] <- 1 - pchisq(comps[1, 1], ncol(Z))#P-value of the Location Comp for Category 2#
comps[2, 1] <- tZZ[2, 2]#Dispersion Component for Category 2#
comps[2, 2] <- 1 - pchisq(comps[2, 1], ncol(Z))#P-value of Dispersion Comp of Category 2#
comps[3, 1] <-tZZ[3,3]#Cubic Components for Category 2#
comps[3, 2] <- 1 - pchisq(comps[3, 1], ncol(Z))#P-value of Cubic Comp of Category 2#
comps[4, 1] <- factor - (comps[1, 1] + comps[2, 1]+comps[3,1])#Error of Components for Category 2#
if(nrow(Z) > 3) {
comps[4, 2] <- 1 - pchisq(comps[4, 1], (nrow(Z) - 3) * ncol(Z))#P-value of Errors for Category 2#
}
else {
comps[4, 2] <- 0
}
comps[5, 1] <- factor
comps[5, 2] <- 1 - pchisq(comps[4, 1], nrow(Z) * ncol(Z))
return(comps)
}
