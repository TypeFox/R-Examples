kfolds2CVinfos_beta <- function(pls_kfolds,MClassed=FALSE) {
if(!(match("dataY",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("dataY", "dataX", "nt", "limQ2set", "modele", "family", "scaleX", "scaleY", "weights", "method", "sparse", "naive", "link", "link.phi", "type"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
if (!is.null(mf$family)) {mf$modele <- "pls-glm-family"}
(mf[[1]] <- as.name("PLS_beta"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}

if (as.character(pls_kfolds$call["modele"]) == "pls") {
press_kfolds <- kfolds2Press(pls_kfolds)
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2_2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
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


    }
}

if (as.character(pls_kfolds$call["modele"]) %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
press_kfolds <- kfolds2Press(pls_kfolds)
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
}

if (as.character(pls_kfolds$call["modele"]) == "pls-glm-polr") {
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}



if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
}

if (as.character(pls_kfolds$call["modele"]) %in% c("pls-beta")) {
press_kfolds <- kfolds2Press(pls_kfolds)
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$pseudo.R2[1:computed_nt]), c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "pseudo_R2_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$pseudo.R2[1:computed_nt]), c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "pseudo_R2_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
}

return(CVinfos)
}

if(!(match("formula",names(pls_kfolds$call), 0L)==0L)){
(mf <- pls_kfolds$call)
(m <- match(c("formula", "data", "nt", "limQ2set", "modele", "family", "scaleX", "scaleY", "weights","subset","start","etastart","mustart","offset","control","method","contrasts","method", "sparse", "naive", "link", "link.phi", "type"), names(pls_kfolds$call), 0))
(mf <- mf[c(1, m)])
(mf$typeVC <- "none")
(mf$MClassed <- MClassed)
if (mf$modele %in% c("pls","pls-glm-logistic","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson","pls-glm-polr")){mf$family <- NULL}
(mf[[1]] <- as.name("PLS_beta_formula"))
(tempres <- eval(mf, parent.frame()))
nt <- as.numeric(as.character(pls_kfolds$call["nt"]))
computed_nt <- tempres$computed_nt
if (MClassed==TRUE) {
Mclassed_kfolds <- kfolds2Mclassed(pls_kfolds)
}

if (as.character(pls_kfolds$call["modele"]) == "pls") {
press_kfolds <- kfolds2Press(pls_kfolds)
Q2cum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2_2 <- 1-press_kfolds[[nnkk]][1:min(length(press_kfolds[[nnkk]]),computed_nt)]/tempres$RSS[1:min(length(press_kfolds[[nnkk]]),computed_nt)]
            for (k in 1:min(length(press_kfolds[[nnkk]]),computed_nt)) {Q2cum_2[k] <- prod(press_kfolds[[nnkk]][1:k])/prod(tempres$RSS[1:k])}
            Q2cum_2 <- 1 - Q2cum_2
            if(length(Q2_2)<computed_nt) {Q2_2 <- c(Q2_2,rep(NA,computed_nt-length(Q2_2)))}
            if(length(Q2cum_2)<computed_nt) {Q2_2cum_2 <- c(Q2cum_2,rep(NA,computed_nt-length(Q2cum_2)))}
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


    }
}

if (as.character(pls_kfolds$call["modele"]) %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
press_kfolds <- kfolds2Press(pls_kfolds)
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
}

if (as.character(pls_kfolds$call["modele"]) == "pls-glm-polr") {
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}



if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)]))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
}

if (as.character(pls_kfolds$call["modele"]) %in% c("pls-beta")) {
press_kfolds <- kfolds2Press(pls_kfolds)
preChisq_kfolds <- kfolds2Chisq(pls_kfolds)
Q2Chisqcum_2=rep(NA, nt)
CVinfos <- vector("list",length(pls_kfolds[[1]]))
limQ2 <- rep(as.numeric(as.character(pls_kfolds$call["limQ2set"])),computed_nt)

    for (nnkk in 1:length(pls_kfolds[[1]])) {
            cat(paste("NK:", nnkk, "\n"))
            Q2Chisq_2 <- 1-preChisq_kfolds[[nnkk]][1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]/tempres$ChisqPearson[1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)]

            for (k in 1:min(length(preChisq_kfolds[[nnkk]]),computed_nt)) {Q2Chisqcum_2[k] <- prod(preChisq_kfolds[[nnkk]][1:k])/prod(tempres$ChisqPearson[1:k])}
            Q2Chisqcum_2 <- 1 - Q2Chisqcum_2
            if(length(Q2Chisq_2)<computed_nt) {Q2Chisq_2 <- c(Q2Chisq_2,rep(NA,computed_nt-length(Q2Chisq_2)))}
            if(length(Q2Chisqcum_2)<computed_nt) {Q2Chisqcum_2 <- c(Q2Chisqcum_2,rep(NA,computed_nt-length(Q2Chisqcum_2)))}
            if(length(press_kfolds[[nnkk]])<computed_nt) {press_kfolds[[nnkk]] <- c(press_kfolds[[nnkk]],rep(NA,computed_nt-length(press_kfolds[[nnkk]])))}


if (MClassed==FALSE) {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$pseudo.R2[1:computed_nt]), c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y","pseudo_R2_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
} else {
            CVinfos[[nnkk]] <- t(rbind(tempres$AIC[1:(computed_nt+1)], tempres$BIC[1:(computed_nt+1)], tempres$MissClassed[1:(computed_nt+1)], c(NA,Mclassed_kfolds[[nnkk]][1:computed_nt]), c(NA,Q2Chisqcum_2[1:computed_nt]), c(NA,limQ2[1:computed_nt]), c(NA,Q2Chisq_2[1:computed_nt]), c(NA,preChisq_kfolds[[nnkk]][1:computed_nt]), tempres$ChisqPearson[1:(computed_nt+1)], tempres$RSS[1:(computed_nt+1)], c(NA,tempres$pseudo.R2[1:computed_nt]), c(NA,tempres$R2[1:computed_nt])))
            dimnames(CVinfos[[nnkk]]) <- list(paste("Nb_Comp_",0:computed_nt,sep=""), c("AIC", "BIC",  "MissClassed", "CV_MissClassed", "Q2Chisqcum_Y", "limQ2", "Q2Chisq_Y", "PREChi2_Pearson_Y", "Chi2_Pearson_Y", "RSS_Y", "pseudo_R2_Y", "R2_Y"))
            CVinfos[[nnkk]] <- cbind(CVinfos[[nnkk]],tempres$ic.dof)
}

    }
}

return(CVinfos)
}
}
