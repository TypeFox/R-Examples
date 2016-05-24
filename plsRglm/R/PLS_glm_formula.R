PLS_glm_formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,subset,start=NULL,etastart,mustart,offset,method,control= list(),contrasts=NULL,sparse=FALSE,sparseStop=FALSE,naive=FALSE,verbose=TRUE) {  


##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

if(verbose){cat("____************************************************____\n")}
if(is.null(modele)){naive=FALSE} else {if(modele=="pls"){naive=FALSE} else {if(!missing(naive)){if(verbose){cat(paste("Only naive DoF can be used with PLS GLM\n",sep=""))}}; naive=TRUE}}
if(sparse){sparseStop=TRUE}
if(sparseStop){pvals.expli=TRUE}

if (missing(data)) {data <- environment(formula)}
mf <- mf2 <- match.call(expand.dots = FALSE)
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
m2 <- match(c("formula","weights","subset","start","etastart","mustart","offset","control","method","contrasts"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$na.action <- na.exclude    
mf2$model <- FALSE
mf2[[1L]] <- as.name("glm")
if(match("method",names(mf2), 0L)==0L){mf2$method<-"glm.fit";method<-"glm.fit"}
}
if (modele %in% c("pls-glm-polr")) {
m2 <- match(c("formula","weights","subset","start","method","contrasts"), names(mf2), 0L)
mf2 <- mf2[c(1L, m2)]
mf2$na.action <- na.exclude    
mf2$model <- FALSE
mf2$Hess <- FALSE
mf2[[1L]] <- as.name("polr")
if(match("method",names(mf2), 0L)==0L){mf2$method<-"logistic"} else {if(!(mf2$method %in% c("logistic", "probit", "cloglog", "cauchit"))) {mf2$method<-"logistic"}}
}
m <- match(c("formula", "data", "subset", "weights", "etastart", "mustart", "offset"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf$na.action <- na.pass    
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame(n=sys.nframe()))

if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if (identical(method, "model.frame")) return(mf)
if (!is.character(method) && !is.function(method)) stop("invalid 'method' argument")
if (identical(method, "glm.fit")) control <- do.call("glm.control", control)
}
mt <- attr(mf, "terms")
attr(mt,"intercept")<-0L
dataY <- model.response(mf, "any")
if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}
if(!NoWeights){naive=TRUE; if(verbose){cat(paste("Only naive DoF can be used with weighted PLS or PLS glms\n",sep=""))}} else {NoWeights=TRUE}
if (length(dim(dataY)) == 1L) {
    nm <- rownames(dataY)
    dim(dataY) <- NULL
    if (!is.null(nm)) names(dataY) <- nm
    }
dataX <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    else matrix(, NROW(dataY), 0L)
weights <- as.vector(model.weights(mf))
if (!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
if (!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
offset <- as.vector(model.offset(mf))
if (!is.null(offset)) {
    if (length(offset) != NROW(dataY)) stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(dataY)), domain = NA)
    }

if(any(apply(is.na(dataX),MARGIN=2,"all"))){return(vector("list",0)); cat("One of the columns of dataX is completely filled with missing data"); stop()}
if(any(apply(is.na(dataX),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of dataX is completely filled with missing data"); stop()}
if(identical(dataPredictY,dataX)){PredYisdataX <- TRUE} else {PredYisdataX <- FALSE}
if(!PredYisdataX){
#if(any(apply(is.na(dataPredictY),MARGIN=2,"all"))){return(vector("list",0)); cat("One of the columns of dataPredictY is completely filled with missing data"); stop()}
if(any(apply(is.na(dataPredictY),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of dataPredictY is completely filled with missing data"); stop()}
}

if(any(is.na(dataX))) {na.miss.X <- TRUE} else na.miss.X <- FALSE
if(any(is.na(dataY))) {na.miss.Y <- TRUE} else na.miss.Y <- FALSE
if(any(is.na(dataPredictY))) {na.miss.PredictY <- TRUE} else {na.miss.PredictY <- FALSE}
if(na.miss.X|na.miss.Y){naive=TRUE; if(verbose){cat(paste("Only naive DoF can be used with missing data\n",sep=""))}; if(!NoWeights){if(verbose){cat(paste("Weights cannot be used with missing data\n",sep=""))}};  if(sparse){if(verbose){cat(paste("sparse option cannot be used with missing data\n",sep=""))}; sparse=FALSE}}

if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}
if (is.null(modele) & !is.null(family)) {modele<-"pls-glm-family"}
if (!(modele %in% c("pls","pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson","pls-glm-polr"))) {print(modele);stop("'modele' not recognized")}
if (!(modele %in% "pls-glm-family") & !is.null(family)) {stop("Set 'modele=pls-glm-family' to use the family option")}
if (modele=="pls") {family<-NULL}
if (modele=="pls-glm-Gamma") {family<-Gamma(link = "inverse")}
if (modele=="pls-glm-gaussian") {family<-gaussian(link = "identity")}
if (modele=="pls-glm-inverse.gaussian") {family<-inverse.gaussian(link = "1/mu^2")}
if (modele=="pls-glm-logistic") {family<-binomial(link = "logit")}
if (modele=="pls-glm-poisson") {family<-poisson(link = "log")}
if (modele=="pls-glm-polr") {family<-NULL}
if (!is.null(family)) {
    if (is.character(family)) {family <- get(family, mode = "function", envir = parent.frame(n=sys.nframe()))}
    if (is.function(family)) {family <- family()}
    if (is.language(family)) {family <- eval(family)}
}
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {mf2$family <- family;if(verbose){print(family)}}
if (modele %in% c("pls-glm-polr")) {if(verbose){cat("\nModel:", modele, "\n");cat("Method:", mf2$method, "\n\n")}}
if (modele=="pls") {if(verbose){cat("\nModel:", modele, "\n\n")}}

scaleY <- NULL
if (is.null(scaleY)) {
if (!(modele %in% c("pls"))) {scaleY <- FALSE} else {scaleY <- TRUE}
}
if (scaleY) {if(NoWeights){RepY <- scale(dataY)} else {meanY <- weighted.mean(dataY,weights); stdevY <- sqrt((length(dataY)-1)/length(dataY)*weighted.mean((dataY-meanY)^2,weights)); RepY <- (dataY-meanY)/stdevY; attr(RepY,"scaled:center") <- meanY ; attr(RepY,"scaled:scale") <- stdevY}}
else {
    RepY <- dataY
    attr(RepY,"scaled:center") <- 0
    attr(RepY,"scaled:scale") <- 1
}
if (scaleX) {if(NoWeights){ExpliX <- scale(dataX)} else {meanX <- apply(dataX,2,weighted.mean,weights); stdevX <- sqrt((length(dataY)-1)/length(dataY)*apply((sweep(dataX,2,meanX))^2,2,weighted.mean,weights)); ExpliX <- sweep(sweep(dataX, 2, meanX), 2 ,stdevX, "/"); attr(ExpliX,"scaled:center") <- meanX ; attr(ExpliX,"scaled:scale") <- stdevX}
    if(PredYisdataX){PredictY <- ExpliX} else {PredictY <- sweep(sweep(dataPredictY, 2, attr(ExpliX,"scaled:center")), 2 ,attr(ExpliX,"scaled:scale"), "/")}
}
else {
    ExpliX <- dataX
    attr(ExpliX,"scaled:center") <- rep(0,ncol(dataX))
    attr(ExpliX,"scaled:scale") <- rep(1,ncol(dataX))
    PredictY <- (dataPredictY)
}
if(is.null(colnames(ExpliX))){colnames(ExpliX)<-paste("X",1:ncol(ExpliX),sep=".")}
if(is.null(rownames(ExpliX))){rownames(ExpliX)<-1:nrow(ExpliX)}

XXNA <- !(is.na(ExpliX))
YNA <- !(is.na(RepY))
if(PredYisdataX){PredictYNA <- XXNA} else {PredictYNA <- !is.na(PredictY)}

ExpliXwotNA <- as.matrix(ExpliX)
ExpliXwotNA[!XXNA] <- 0

XXwotNA <- as.matrix(ExpliX)
XXwotNA[!XXNA] <- 0

dataXwotNA <- as.matrix(dataX)
dataXwotNA[!XXNA] <- 0

YwotNA <- as.matrix(RepY)
YwotNA[!YNA] <- 0

dataYwotNA <- as.matrix(dataY)
dataYwotNA[!YNA] <- 0

if(PredYisdataX){PredictYwotNA <- XXwotNA} else {
PredictYwotNA <- as.matrix(PredictY)
PredictYwotNA [is.na(PredictY)] <- 0
}

if (modele == "pls-glm-polr") {
dataY <- as.factor(dataY)
YwotNA <- as.factor(YwotNA)}

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,RSS=rep(NA,nt),RSSresidY=rep(NA,nt),R2=rep(NA,nt),R2residY=rep(NA,nt),press.ind=NULL,press.tot=NULL,Q2cum=rep(NA, nt),family=family,ttPredictY = NULL,typeVC=typeVC,dataX=dataX,dataY=dataY) 
if(NoWeights){res$weights<-rep(1L,res$nr)} else {res$weights<-weights}
res$temppred <- NULL

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (scaleY) {res$YChapeau=rep(attr(RepY,"scaled:center"),nrow(ExpliX))
res$residYChapeau=rep(0,nrow(ExpliX))}
else
{res$YChapeau=rep(mean(RepY),nrow(ExpliX))
res$residYChapeau=rep(mean(RepY),nrow(ExpliX))}
}





################################################
################################################
##                                            ##
##  Beginning of the loop for the components  ##
##                                            ##
################################################
################################################

res$computed_nt <- 0
break_nt <- FALSE
break_nt_sparse <- FALSE
break_nt_sparse1 <- FALSE
break_nt_vc <- FALSE

for (kk in 1:nt) {
XXwotNA <- as.matrix(res$residXX)
XXwotNA[!XXNA] <- 0
YwotNA <- as.matrix(res$residY)
YwotNA[!YNA] <- 0
tempww <- rep(0,res$nc)


temptest <- sqrt(colSums(res$residXX^2, na.rm=TRUE))
if(any(temptest<tol_Xi)) {
break_nt <- TRUE
if (is.null(names(which(temptest<tol_Xi)))) {
  if(verbose){cat(paste("Warning : ",paste(names(which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))}
} else {
  if(verbose){cat(paste("Warning : ",paste((which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))}
}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
rm(temptest)
break
}

res$computed_nt <- kk

##############################################
#                                            #
#     Weight computation for each model      #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if(NoWeights){
tempww <- t(XXwotNA)%*%YwotNA/(t(XXNA)%*%YwotNA^2)
}
if(!NoWeights){
tempww <- t(XXwotNA*weights)%*%YwotNA/(t(XXNA*weights)%*%YwotNA^2)
}
if (pvals.expli) {
tempvalpvalstep <- 2 * pnorm(-abs(tempww)) 
temppvalstep <- (tempvalpvalstep < alpha.pvals.expli)
if(sparse){if(sum(temppvalstep)>1L){tempww[!temppvalstep] <- 0}}
if(sparseStop){if(sum(temppvalstep)==0L){break_nt_sparse <- TRUE}}
res$valpvalstep <- cbind(res$valpvalstep,tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,temppvalstep)
}
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if (!pvals.expli) {
XXwotNA[!XXNA] <- NA
for (jj in 1:(res$nc)) {
    mf2[[2]]<-YwotNA~cbind(res$tt,XXwotNA[,jj])
    tempww[jj] <- coef(eval(mf2, parent.frame(n=sys.nframe())))[kk+1]
}
XXwotNA[!XXNA] <- 0
rm(jj)}
else {
XXwotNA[!XXNA] <- NA
tempvalpvalstep <- rep(0,res$nc)
temppvalstep <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
    mf2[[2]]<-YwotNA~cbind(res$tt,XXwotNA[,jj])
    tmww <- summary(eval(mf2, parent.frame(n=sys.nframe())))$coefficients[kk+1,]
    tempww[jj] <- tmww[1]
    tempvalpvalstep[jj] <- tmww[4]
    temppvalstep[jj] <- (tmww[4] < alpha.pvals.expli)
}
if(sparse){if(sum(temppvalstep)>1L){tempww[!temppvalstep] <- 0}}
if(sparseStop){if(sum(temppvalstep)==0L){break_nt_sparse <- TRUE}}
XXwotNA[!XXNA] <- 0
rm(jj)
res$valpvalstep <- cbind(res$valpvalstep,tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,temppvalstep)
}
}

##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
YwotNA <- as.factor(YwotNA)
if (!pvals.expli) {
XXwotNA[!XXNA] <- NA
library(MASS)
tts <- res$tt
for (jj in 1:(res$nc)) {
    mf2[[2]]<-YwotNA~cbind(tts,XXwotNA[,jj])
    tempww[jj] <- -1*eval(mf2, parent.frame(n=sys.nframe()))$coef[kk]
}
XXwotNA[!XXNA] <- 0
rm(jj,tts)}
else {
XXwotNA[!XXNA] <- NA
library(MASS)
tts <- res$tt
tempvalpvalstep <- rep(0,res$nc)
temppvalstep <- rep(0,res$nc)
mf2$Hess <- TRUE
for (jj in 1:(res$nc)) {
    mf2[[2]]<-YwotNA~cbind(tts,XXwotNA[,jj])
    tmww <- -1*summary(eval(mf2, parent.frame(n=sys.nframe())))$coefficients[kk,]
    tempww[jj] <- tmww[1]
    tempvalpvalstep[jj] <- 2 * pnorm(-abs(tmww[3])) 
    temppvalstep[jj] <- (tempvalpvalstep[jj] < alpha.pvals.expli)
}
if(sparse){if(sum(temppvalstep)>1L){tempww[!temppvalstep] <- 0}}
if(sparseStop){if(sum(temppvalstep)==0L){break_nt_sparse <- TRUE}}
XXwotNA[!XXNA] <- 0
rm(jj,tts)
mf2$Hess <- FALSE
res$valpvalstep <- cbind(res$valpvalstep,tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,temppvalstep)
}
}




##############################################
#                                            #
# Computation of the components (model free) #
#                                            #
##############################################
if((break_nt_sparse)&(kk==1L)){
  if(verbose){cat(paste("No significant predictors (<",alpha.pvals.expli,") found\n",sep=""))}
  if(verbose){cat(paste("Warning only one standard component (without sparse option) was thus extracted\n",sep=""))}
break_nt_sparse1 <- TRUE
}
if((break_nt_sparse)&!(kk==1L)){
res$computed_nt <- kk-1
if(!(break_nt_sparse1)){
  if(verbose){cat(paste("No more significant predictors (<",alpha.pvals.expli,") found\n",sep=""))}
  if(verbose){cat(paste("Warning only ",res$computed_nt," components were thus extracted\n",sep=""))}
}
break}

tempwwnorm <- tempww/sqrt(drop(crossprod(tempww)))

temptt <- XXwotNA%*%tempwwnorm/(XXNA%*%(tempwwnorm^2))

temppp <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
     temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod(XXNA[,jj],temptt^2))
}
res$residXX <- XXwotNA-temptt%*%temppp

if (na.miss.X & !na.miss.Y) {
  if(sparse==FALSE){
    for (ii in 1:res$nr) {
if(rcond(t(cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])<tol_Xi) {
break_nt <- TRUE; res$computed_nt <- kk-1
if(verbose){cat(paste("Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE] < 10^{-12}\n",sep=""))}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
break
}
}
rm(ii)
if(break_nt==TRUE) {res$computed_nt <- kk-1;break}
}
}

if(!PredYisdataX){
  if(sparse==FALSE){
if (na.miss.PredictY & !na.miss.Y) {
for (ii in 1:nrow(PredictYwotNA)) {
if(rcond(t(cbind(res$pp,temppp)[PredictYNA[ii,],,drop=FALSE])%*%cbind(res$pp,temppp)[PredictYNA[ii,],,drop=FALSE])<tol_Xi) {
break_nt <- TRUE; res$computed_nt <- kk-1
if(verbose){cat(paste("Warning : reciprocal condition number of t(cbind(res$pp,temppp)[PredictYNA[",ii,",,drop=FALSE],])%*%cbind(res$pp,temppp)[PredictYNA[",ii,",,drop=FALSE],] < 10^{-12}\n",sep=""))}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
break
}
}
rm(ii)
if(break_nt==TRUE) {res$computed_nt <- kk-1;break}
}
}
}


res$ww <- cbind(res$ww,tempww)
res$wwnorm <- cbind(res$wwnorm,tempwwnorm)
res$tt <- cbind(res$tt,temptt)       
res$pp <- cbind(res$pp,temppp)   




##############################################
#                                            #
#      Computation of the coefficients       #
#      of the model with kk components       #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
if (kk==1) {
tempCoeffC <- solve(t(res$tt[YNA])%*%res$tt[YNA])%*%t(res$tt[YNA])%*%YwotNA[YNA]
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- 0
} else {
if (!(na.miss.X | na.miss.Y)) {
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
else
{
tempCoeffC <- c(rep(0,kk-1),solve(t(res$tt[YNA,kk])%*%res$tt[YNA,kk])%*%t(res$tt[YNA,kk])%*%YwotNA[YNA])  
tempCoeffConstante <- 0
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- diag(res$CoeffCFull)
res$CoeffConstante <- tempCoeffConstante
res$Std.Coeffs <- rbind(tempCoeffConstante,res$wwetoile%*%res$CoeffC)
rownames(res$Std.Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
if (kk==1) {
mf2[[2]]<-YwotNA~1
mf2$model<-TRUE
tempconstglm <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model<-FALSE
res$AIC <- AIC(tempconstglm)
res$BIC <- AIC(tempconstglm, k = log(res$nr))
res$Coeffsmodel_vals <- rbind(summary(tempconstglm)$coefficients,matrix(rep(NA,4*nt),ncol=4))
res$ChisqPearson <- crossprod(residuals.glm(tempconstglm,type="pearson"))  
#if ((modele %in% c("pls-glm-logistic"))|(family$family=="binomial")) {
res$MissClassed <- sum(unclass(res$RepY)!=ifelse(predict(tempconstglm,type="response") < 0.5, 0,1))
#}
rm(tempconstglm)
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
mf2$data <- tttrain  
mf2$model <- TRUE
mf2[[2]]<-YwotNA~.
tempregglm <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
#if ((modele %in% c("pls-glm-logistic"))|(family$family=="binomial")) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(res$RepY)!=ifelse(predict(tempregglm,type="response") < 0.5, 0,1)))
#}
tempCoeffC <- as.vector(coef(tempregglm))
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- tempCoeffConstante
tempCoeffC <- tempCoeffC[-1]
} else {
if (!(na.miss.X | na.miss.Y)) {
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
mf2$data <- tttrain  
mf2$model <- TRUE
mf2[[2]]<-YwotNA~.
tempregglm <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
#if ((modele %in% c("pls-glm-logistic"))|(family$family=="binomial")) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(res$RepY)!=ifelse(predict(tempregglm,type="response") < 0.5, 0,1)))
#}
tempCoeffC <- as.vector(coef(tempregglm))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
tempCoeffC <- tempCoeffC[-1]
}
else
{
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
mf2$data <- tttrain  
mf2$model <- TRUE
mf2[[2]]<-YwotNA~.
tempregglm <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregglm))
res$BIC <- cbind(res$BIC,AIC(tempregglm, k = log(res$nr)))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
res$ChisqPearson <- c(res$ChisqPearson,crossprod(residuals.glm(tempregglm,type="pearson")))
#if ((modele %in% c("pls-glm-logistic"))|(family$family=="binomial")) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(res$RepY)!=ifelse(predict(tempregglm,type="response") < 0.5, 0,1)))
#}
tempCoeffC <- as.vector(coef(tempregglm))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
tempCoeffC <- tempCoeffC[-1]
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- rbind(tempCoeffConstante,res$wwetoile%*%res$CoeffC)
rownames(res$Std.Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######           PLS-GLM-POLR           ######
##############################################


if (modele %in% c("pls-glm-polr")) {
            Varyy <- function(piVaryy) {
            diag(piVaryy[-length(piVaryy)])-piVaryy[-length(piVaryy)]%*%t(piVaryy[-length(piVaryy)])
            }
            Chisqcomp <- function(yichisq,pichisq) {#change 2010/12/29 solve -> ginv
            t(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])%*%MASS::ginv(Varyy(pichisq))%*%(yichisq[-length(yichisq)]-pichisq[-length(pichisq)])
            }
            Chiscompmatrix <- function(rowspi,rowsyi) {
            sum(mapply(Chisqcomp,rowsyi,rowspi))
            }
if (kk==1) {
mf2$model <- TRUE
mf2$Hess <- TRUE
mf2[[2]]<-YwotNA~1
tempconstpolr <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
mf2$Hess <- FALSE
res$AIC <- AIC(tempconstpolr)
res$BIC <- AIC(tempconstpolr, k = log(res$nr))
res$MissClassed <- sum(!(unclass(predict(tempconstpolr,type="class"))==unclass(res$RepY)))
res$Coeffsmodel_vals <- rbind(summary(tempconstpolr)$coefficients,matrix(rep(NA,3*nt),ncol=3))
tempmodord <- predict(tempconstpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempconstpolr,type="probs")))),as.list(as.data.frame(t(tempmat)))))
rm(tempconstpolr)
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
mf2$model <- TRUE
mf2$Hess <- TRUE
mf2[[2]]<-YwotNA~.
mf2$data <- tttrain
tempregpolr <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
mf2$Hess <- FALSE
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- matrix(c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)),ncol=1)
res$CoeffConstante <- tempCoeffConstante
} else {
if (!(na.miss.X | na.miss.Y)) {
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
mf2$model <- TRUE
mf2$Hess <- TRUE
mf2[[2]]<-YwotNA~.
mf2$data <- tttrain
tempregpolr <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
mf2$Hess <- FALSE
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
else
{
rm(tempconstpolr)
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
mf2$model <- TRUE
mf2$Hess <- TRUE
mf2[[2]]<-YwotNA~.
mf2$data <- tttrain
tempregpolr <- eval(mf2, parent.frame(n=sys.nframe()))
mf2$model <- FALSE
mf2$Hess <- FALSE
rm(tttrain)
res$AIC <- cbind(res$AIC,AIC(tempregpolr))
res$BIC <- cbind(res$BIC,AIC(tempregpolr, k = log(res$nr)))
res$MissClassed <- cbind(res$MissClassed,sum(!(unclass(predict(tempregpolr,type="class"))==unclass(res$RepY))))
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempmodord <- predict(tempregpolr,type="class")
tempfff <- ~tempmodord-1
tempm <- model.frame(tempfff, tempmodord)
tempmat <- model.matrix(tempfff, model.frame(tempfff, tempmodord))
res$ChisqPearson <- c(res$ChisqPearson,sum(Chiscompmatrix(as.list(as.data.frame(t(predict(tempregpolr,type="probs")))),as.list(as.data.frame(t(tempmat))))))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- as.matrix(rbind(as.matrix(tempCoeffConstante),res$wwetoile%*%res$CoeffC))
rownames(res$Std.Coeffs) <- c(names(tempregpolr$zeta),colnames(ExpliX))
}




##############################################
#                                            #
#       Prediction of the components         #
#     as if missing values (model free)      #
#       For cross-validating the GLM         #
#                                            #
##############################################





if (!(na.miss.X | na.miss.Y)) {

##############################################
#                                            #
#             Cross validation               #
#           without missing value            #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}
##############################################


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$residYChapeau <- tempregglm$linear.predictors
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values          
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}


##############################################
######              PLS-GLM-POLR         ######
##############################################
if (modele %in% c("pls-glm-polr")) {
tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*tempCoeffConstante
res$Coeffs <- rbind(as.matrix(tempConstante),tempCoeffs)
rownames(res$Coeffs) <- rownames(res$Std.Coeffs)
}                                                                                                                             
##############################################
}

else {
if (na.miss.X & !na.miss.Y) {


##############################################
#                                            #
#             Cross validation               #
#           with missing value(s)            #
#                                            #
##############################################


if (kk==1) {
  if(verbose){cat("____There are some NAs in X but not in Y____\n")}
}

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}
##############################################



##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$residYChapeau <- tempregglm$linear.predictors
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY))
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY-mean(RepY),weights*(RepY-mean(RepY)))
}
}
if(NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau))
}
if(!NoWeights){
res$RSSresidY <- cbind(res$RSSresidY,crossprod(res$residY-res$residYChapeau,weights*(res$residY-res$residYChapeau)))
}

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values                      
res$Yresidus <- dataY-res$YChapeau
if (kk==1) {
if(NoWeights){
res$RSS <- crossprod(dataY-mean(dataY))
}
if(!NoWeights){
res$RSS <- crossprod(dataY-mean(dataY),weights*(dataY-mean(dataY)))
}
}
if(NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus))
}
if(!NoWeights){
res$RSS <- cbind(res$RSS,crossprod(res$Yresidus,weights*res$Yresidus))
}
}


##############################################
######              PLS-GLM-POLR         ######
##############################################
if (modele %in% c("pls-glm-polr")) {
tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")* tempCoeffConstante
res$Coeffs <- rbind(as.matrix(tempConstante),tempCoeffs)
rownames(res$Coeffs) <- rownames(res$Std.Coeffs)
}
##############################################
}

else {
if (kk==1) {
  if(verbose){cat("____There are some NAs both in X and Y____\n")}
}
}
}


##############################################
#                                            #
#      Update and end of loop cleaning       #
#        (Especially useful for PLS)         #
#                                            #
##############################################


##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$uscores <- cbind(res$uscores,res$residY/res$CoeffC[kk])
res$residY <- res$residY - res$tt%*%tempCoeffC 
res$residusY <- cbind(res$residusY,res$residY)


if (kk==1) {
res$AIC.std <- AIC(lm(res$RepY~1,weights=res$weights))
res$AIC.std <- cbind(res$AIC.std,AICpls(kk,res$residY,weights=res$weights))
res$AIC <- AIC(lm(dataY~1,weights=res$weights))
res$AIC <- cbind(res$AIC,AICpls(kk,res$Yresidus,weights=res$weights))
if (MClassed) {
res$MissClassed <- sum(unclass(dataY)!=ifelse(predict(lm(dataY~1,weights=res$weights)) < 0.5, 0,1))
res$MissClassed <- cbind(res$MissClassed,sum(unclass(dataY)!=ifelse(res$YChapeau < 0.5, 0,1)))
tempprob <- res$Probs <- predict(lm(dataY~1,weights=res$weights))
tempprob <- ifelse(tempprob<0,0,tempprob)
res$Probs.trc <- ifelse(tempprob>1,1,tempprob)
res$Probs <- cbind(res$Probs,res$YChapeau)
tempprob <- ifelse(res$YChapeau<0,0,res$YChapeau)
tempprob <- ifelse(tempprob>1,1,tempprob)
res$Probs.trc <- cbind(res$Probs.trc,tempprob)
}
} else {
res$AIC.std <- cbind(res$AIC.std,AICpls(kk,res$residY,weights=res$weights))
res$AIC <- cbind(res$AIC,AICpls(kk,res$Yresidus,weights=res$weights))
if (MClassed) {
res$MissClassed <- cbind(res$MissClassed,sum(unclass(dataY)!=ifelse(res$YChapeau < 0.5, 0,1)))
res$Probs <- cbind(res$Probs,res$YChapeau)
tempprob <- ifelse(res$YChapeau<0,0,res$YChapeau)
tempprob <- ifelse(tempprob>1,1,tempprob)
res$Probs.trc <- cbind(res$Probs.trc,tempprob)
}
}


rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}

##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}

##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
rm(tempCoeffs)
rm(tempConstante)
}
if(verbose){cat("____Component____",kk,"____\n")}
}




##############################################
##############################################
##                                          ##
##    End of the loop on the components     ##
##                                          ##
##############################################
##############################################

if(res$computed_nt==0){
cat("No component could be extracted please check the data for NA only lines or columns\n"); stop()
}


if (pvals.expli&!(modele=="pls")) {
res$Coeffsmodel_vals<-res$Coeffsmodel_vals[1:(dim(res$Coeffsmodel_vals)[1]-(nt-res$computed_nt)),]
}



##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.PredictY | na.miss.Y)) {
  if(verbose){cat("____Predicting X without NA neither in X or Y____\n")}
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- NULL
}
else {
if (na.miss.PredictY & !na.miss.Y) {
  if(verbose){cat("____Predicting X with NA in X and not in Y____\n")}
for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%res$pp[PredictYNA[ii,],,drop=FALSE])%*%t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}
colnames(res$ttPredictY) <- NULL
}
else {
  if(verbose){cat("____There are some NAs both in X and Y____\n")}
}
}




##############################################
#                                            #
#          Computing RSS, PRESS,             #
#           Chi2, Q2 and Q2cum               #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele == "pls") {
res$R2residY <- 1-res$RSSresidY[2:(res$computed_nt+1)]/res$RSSresidY[1]
res$R2 <- 1-res$RSS[2:(res$computed_nt+1)]/res$RSS[1]
if (MClassed==FALSE) {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "AIC.std"))
res$ic.dof<-infcrit.dof(res,naive=naive)
res$InfCrit <- cbind(res$InfCrit,res$ic.dof)
} else {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(NA,res$R2), res$MissClassed, c(NA,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "AIC.std"))
res$ic.dof<-infcrit.dof(res,naive=naive)
res$InfCrit <- cbind(res$InfCrit,res$ic.dof)
}
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$R2residY <- 1-res$RSSresidY[2:(res$computed_nt+1)]/res$RSSresidY[1]
res$R2 <- 1-res$RSS[2:(res$computed_nt+1)]/res$RSS[1]
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC, res$ChisqPearson, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}
if ((modele %in% c("pls-glm-logistic"))|(family$family=="binomial")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC, res$MissClassed, res$ChisqPearson, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC", "Missclassed", "Chi2_Pearson_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY"))
}
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele == "pls-glm-polr") {

res$InfCrit <- t(rbind(res$AIC, res$BIC, res$MissClassed, res$ChisqPearson))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC", "Missclassed", "Chi2_Pearson_Y"))
}




##########################################
#                                        #
#          Predicting responses          #
#                                        #
##########################################


##############################################
######               PLS                ######
##############################################
if (modele == "pls") {
res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
rownames(res$YChapeau) <- rownames(ExpliX)

res$Std.ValsPredictY <- res$ttPredictY%*%res$CoeffC
res$ValsPredictY <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$ttPredictY%*%res$CoeffC

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
if (EstimXNA) {
res$XChapeau <- sweep(sweep(res$Std.XChapeau,2,attr(res$ExpliX,"scaled:scale"),FUN="*"),2,attr(res$ExpliX,"scaled:center"),FUN="+")
rownames(res$XChapeau) <- rownames(ExpliX)
colnames(res$XChapeau) <- colnames(ExpliX)

res$XChapeauNA <- sweep(sweep(res$Std.XChapeau,2,attr(res$ExpliX,"scaled:scale"),FUN="*"),2,attr(res$ExpliX,"scaled:center"),FUN="+")*!XXNA
rownames(res$XChapeau) <- rownames(ExpliX)
colnames(res$XChapeau) <- colnames(ExpliX)
}
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {
res$YChapeau <- as.matrix(tempregglm$fitted.values)            
rownames(res$YChapeau) <- rownames(ExpliX)

ttpred <- data.frame(tt=res$ttPredictY)
res$Std.ValsPredictY <- predict(tempregglm,newdata=ttpred)
res$ValsPredictY <- predict(tempregglm,newdata=ttpred,type = "response")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
res$FinalModel <- tempregglm
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
res$YChapeau <- tempregpolr$fitted.values
res$YChapeauCat <- predict(tempregpolr,type="class")
rownames(res$YChapeau) <- rownames(ExpliX)

ttpred <- data.frame(tt=res$ttPredictY)
res$ValsPredictY <- predict(tempregpolr, data=ttpred,type="probs")
res$ValsPredictYCat <- predict(tempregpolr, data=ttpred,type="class")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt,sep="")
res$FinalModel <- tempregpolr
}


colnames(res$ttPredictY) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$pp) <- colnames(ExpliX)
colnames(res$pp) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$ww) <- colnames(ExpliX)
colnames(res$ww) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$wwnorm) <- colnames(ExpliX)
colnames(res$wwnorm) <- paste("Comp_",1:res$computed_nt,sep="")
rownames(res$wwetoile) <- colnames(ExpliX)
colnames(res$wwetoile) <- paste("Coord_Comp_",1:res$computed_nt,sep="")
rownames(res$tt) <- rownames(ExpliX)
colnames(res$tt) <- paste("Comp_",1:res$computed_nt,sep="")
res$XXwotNA <- XXwotNA
if(verbose){cat("****________________________________________________****\n")}
if(verbose){cat("\n")}
return(res)
}
