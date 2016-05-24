kfolds2CVinfos_lm <- function(pls_kfolds,MClassed=FALSE,verbose=TRUE) {
if(!(match("dataY",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("dataY", "dataX", "nt", "limQ2set", "modele", "family", "scaleX", "scaleY", "weights", "method", "sparse", "naive", "verbose"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
mf$verbose<-verbose
if(is.null(mf$modele)){mf$modele<-"pls"}
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
(mf[[1]] <- as.name("PLS_lm"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
press_kfolds <- kfolds2Press(pls_kfolds)
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}
if (mf$modele == "pls") {

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "MissClassed", "CV_MissClassed", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}
attr(CVinfos[[nnkk]],"computed_nt") <- computed_nt
    }
}
return(CVinfos)
}



if(!(match("formula",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("formula", "data", "nt", "limQ2set", "modele", "scaleX", "scaleY","weights","subset","contrasts", "method", "sparse", "naive", "verbose"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
mf$verbose<-verbose
if(is.null(mf$modele)){mf$modele<-"pls"}
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
(mf[[1]] <- as.name("PLS_lm_formula"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
press_kfolds <- kfolds2Press(pls_kfolds)
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
if(is.numeric(pls_kfolds$call["limQ2set"])){limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)} else {limQ2=rep(as.numeric(as.character(0.0975)),computed_nt)}
if (mf$modele == "pls") {

    for (nnkk in 1:length(pls_kfolds[[1]])) {
      if(verbose){if(nnkk%%10==1){cat("\n");cat(paste("NK:", nnkk))} else {cat(paste(", ", nnkk))}}
      Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2cum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2_2[1:computed_nt]), c(NA,press_kfolds[[nnkk]][1:computed_nt]), tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt]), tempres$AIC.std[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "MissClassed", "CV_MissClassed", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "AIC.std"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}
attr(CVinfos[[nnkk]],"computed_nt") <- computed_nt
    }
}
if(verbose){cat("\n")};
return(CVinfos)
}
}
