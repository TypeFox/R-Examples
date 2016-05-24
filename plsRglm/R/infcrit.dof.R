infcrit.dof <- function(modplsR,naive=FALSE){
if(!(is.null(modplsR$weights) | identical(modplsR$weights,rep(1L,modplsR$nr)) | identical(modplsR$weights,rep(1,modplsR$nr)))){naive=TRUE}
if(modplsR$na.miss.X|modplsR$na.miss.Y){naive=TRUE}
if(!naive){
tempmodplsR_dof <- plsR.dof(modplsR,naive=FALSE)
tempAIC.dof <- aic.dof(modplsR$RSS,modplsR$nr,tempmodplsR_dof$DoF,tempmodplsR_dof$sigmahat)
tempBIC.dof <- bic.dof(modplsR$RSS,modplsR$nr,tempmodplsR_dof$DoF,tempmodplsR_dof$sigmahat)
tempGMDL.dof <- gmdl.dof(tempmodplsR_dof$sigmahat,modplsR$nr,tempmodplsR_dof$DoF,tempmodplsR_dof$yhat)
tempmodplsR_naive <- plsR.dof(modplsR,naive=TRUE)
tempAIC.naive <- aic.dof(modplsR$RSS,modplsR$nr,tempmodplsR_naive$DoF,tempmodplsR_naive$sigmahat)
tempBIC.naive <- bic.dof(modplsR$RSS,modplsR$nr,tempmodplsR_naive$DoF,tempmodplsR_naive$sigmahat)
tempGMDL.naive <- gmdl.dof(tempmodplsR_naive$sigmahat,modplsR$nr,tempmodplsR_naive$DoF,tempmodplsR_naive$yhat)
InfCrit.dof <- t(rbind(tempmodplsR_dof$DoF,tempmodplsR_dof$sigmahat,tempAIC.dof,tempBIC.dof,tempGMDL.dof,tempmodplsR_naive$DoF,tempmodplsR_naive$sigmahat,tempAIC.naive,tempBIC.naive,tempGMDL.naive))
dimnames(InfCrit.dof) <- list(paste("Nb_Comp_",0:modplsR$computed_nt,sep=""), c("DoF.dof","sigmahat.dof","AIC.dof", "BIC.dof", "GMDL.dof","DoF.naive","sigmahat.naive","AIC.naive", "BIC.naive", "GMDL.naive"))
} else {InfCrit.dof <- NULL}
return(InfCrit.dof)
}
