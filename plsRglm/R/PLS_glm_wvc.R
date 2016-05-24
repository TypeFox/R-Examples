PLS_glm_wvc <- function(dataY,dataX,nt=2,dataPredictY=dataX,modele="pls",family=NULL,scaleX=TRUE,scaleY=NULL,keepcoeffs=FALSE,keepstd.coeffs=FALSE,tol_Xi=10^(-12),weights,method="logistic",verbose=TRUE) {

##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

if(verbose){cat("____************************************************____\n")}
if(any(apply(is.na(dataX),MARGIN=2,"all"))){return(vector("list",0)); cat("One of the columns of dataX is completely filled with missing data"); stop()}
if(any(apply(is.na(dataX),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of dataX is completely filled with missing data"); stop()}
if(identical(dataPredictY,dataX)){PredYisdataX <- TRUE} else {PredYisdataX <- FALSE}
if(!PredYisdataX){
if(any(apply(is.na(dataPredictY),MARGIN=2,"all"))){return(vector("list",0)); cat("One of the columns of dataPredictY is completely filled with missing data"); stop()}
if(any(apply(is.na(dataPredictY),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of dataPredictY is completely filled with missing data"); stop()}
}
if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}
if(any(is.na(dataX))) {na.miss.X <- TRUE} else na.miss.X <- FALSE
if(any(is.na(dataY))) {na.miss.Y <- TRUE} else na.miss.Y <- FALSE
if(any(is.na(dataPredictY))) {na.miss.PredictY <- TRUE} else {na.miss.PredictY <- FALSE}
if(na.miss.X|na.miss.Y){naive=TRUE; if(verbose){cat(paste("Only naive DoF can be used with missing data\n",sep=""))}; if(!NoWeights){if(verbose){cat(paste("Weights cannot be used with missing data\n",sep=""))}}}
if(!NoWeights){naive=TRUE; if(verbose){cat(paste("Only naive DoF can be used with weighted PLS\n",sep=""))}}


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
    if (modele %in% c("pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-logistic","pls-glm-poisson")) {if(verbose){print(family)}}
    if (modele %in% c("pls-glm-polr")) {if(verbose){cat("\nModel:", modele, "\n");cat("Method:", method, "\n\n")}}
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

if (modele %in% "pls-glm-polr") {
dataY <- as.factor(dataY)
YwotNA <- as.factor(YwotNA)}

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,RSS=rep(NA,nt),RSSresidY=rep(NA,nt),R2=rep(NA,nt),R2residY=rep(NA,nt),press.ind=NULL,press.tot=NULL,Q2cum=rep(NA, nt),family=family,ttPredictY = NULL,typeVC="none",listValsPredictY=NULL) 
if(NoWeights){res$weights<-rep(1L,res$nr)} else {res$weights<-weights}
res$temppred <- NULL

##############################################
######                PLS               ######
##############################################
if (modele %in% "pls") {
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
break_nt_vc <- FALSE

for (kk in 1:nt) {

temptest <- sqrt(colSums(res$residXX^2, na.rm=TRUE))
if(any(temptest<tol_Xi)) {
break_nt <- TRUE
if (is.null(names(which(temptest<tol_Xi)))) {
  if(verbose){cat(paste("Warning : ",paste(names(which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))}
} else {
  if(verbose){cat(paste("Warning : ",paste((which(temptest<tol_Xi)),sep="",collapse=" ")," < 10^{-12}\n",sep=""))}
}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
break
}

res$computed_nt <- kk

XXwotNA <- as.matrix(res$residXX)
XXwotNA[!XXNA] <- 0
YwotNA <- as.matrix(res$residY)
YwotNA[!YNA] <- 0
tempww <- rep(0,res$nc)


##############################################
#                                            #
#     Weight computation for each model      #
#                                            #
##############################################

##############################################
######                PLS               ######
##############################################
if (modele %in% "pls") {
if(NoWeights){
tempww <- t(XXwotNA)%*%YwotNA/(t(XXNA)%*%YwotNA^2)
}
if(!NoWeights){
tempww <- t(XXwotNA*weights)%*%YwotNA/(t(XXNA*weights)%*%YwotNA^2)
}
}

##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
XXwotNA[!XXNA] <- NA
for (jj in 1:(res$nc)) {
    tempww[jj] <- coef(glm(YwotNA~cbind(res$tt,XXwotNA[,jj]),family=family))[kk+1]
}
XXwotNA[!XXNA] <- 0
rm(jj)}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
YwotNA <- as.factor(YwotNA)
XXwotNA[!XXNA] <- NA
library(MASS)
tts <- res$tt
for (jj in 1:(res$nc)) {
    tempww[jj] <- -1*MASS::polr(YwotNA~cbind(tts,XXwotNA[,jj]),na.action=na.exclude,method=method)$coef[kk]
}
XXwotNA[!XXNA] <- 0
rm(jj,tts)}




##############################################
#                                            #
# Computation of the components (model free) #
#                                            #
##############################################

tempwwnorm <- tempww/sqrt(drop(crossprod(tempww)))

temptt <- XXwotNA%*%tempwwnorm/(XXNA%*%(tempwwnorm^2))

temppp <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
     temppp[jj] <- crossprod(temptt,XXwotNA[,jj])/drop(crossprod(XXNA[,jj],temptt^2))
}
res$residXX <- XXwotNA-temptt%*%temppp

if (na.miss.X & !na.miss.Y) {
for (ii in 1:res$nr) {
if(rcond(t(cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])<tol_Xi) {
break_nt <- TRUE; res$computed_nt <- kk-1
if(verbose){cat(paste("Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE] < 10^{-12}\n",sep=""))}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
break
}
}
rm(ii)
if(break_nt) {break}
}

if(!PredYisdataX){
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
if(break_nt) {break}
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
if (modele %in% c("pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
if (kk==1) {
tempconstglm <- glm(YwotNA~1,family=family)
res$Coeffsmodel_vals <- rbind(summary(tempconstglm)$coefficients,matrix(rep(NA,4*nt),ncol=4))
rm(tempconstglm)
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
tempCoeffC <- as.vector(coef(tempregglm))
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- tempCoeffConstante
tempCoeffC <- tempCoeffC[-1]
} else {
if (!(na.miss.X | na.miss.Y)) {
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
tempCoeffC <- as.vector(coef(tempregglm))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
tempCoeffConstante <- tempCoeffC[1]
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
tempCoeffC <- tempCoeffC[-1]
}
else
{
tt<-res$tt
tempregglm <- glm(YwotNA~tt,family=family)
rm(tt)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregglm)$coefficients,matrix(rep(NA,4*(nt-kk)),ncol=4)))
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
if (kk==1) {
tempconstpolr <- MASS::polr(YwotNA~1,na.action=na.exclude,Hess=TRUE,method=method)
res$Coeffsmodel_vals <- rbind(summary(tempconstpolr)$coefficients,matrix(rep(NA,3*nt),ncol=3))
rm(tempconstpolr)
tts <- res$tt
tempregpolr <- MASS::polr(YwotNA~tts,na.action=na.exclude,Hess=TRUE,method=method)
rm(tts)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- matrix(c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)),ncol=1)
res$CoeffConstante <- tempCoeffConstante
} else {
if (!(na.miss.X | na.miss.Y)) {
tts <- res$tt
tempregpolr <- MASS::polr(YwotNA~tts,na.action=na.exclude,Hess=TRUE,method=method)
rm(tts)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
tempCoeffC <- -1*as.vector(tempregpolr$coef)  
tempCoeffConstante <- as.vector(tempregpolr$zeta)
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffConstante,tempCoeffC,rep(NA,nt-kk)))
res$CoeffConstante <- cbind(res$CoeffConstante,tempCoeffConstante)
}
else
{
tts <- res$tt
tempregpolr <- MASS::polr(YwotNA~tts,na.action=na.exclude,Hess=TRUE,method=method)
rm(tts)
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregpolr)$coefficients,matrix(rep(NA,3*(nt-kk)),ncol=3)))
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


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC             
res$Yresidus <- dataY-res$YChapeau
}
##############################################


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
res$residYChapeau <- tempregglm$linear.predictors


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))+attr(res$RepY,"scaled:scale")*res$Std.Coeffs[1]
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values          
res$Yresidus <- dataY-res$YChapeau
}


##############################################
######              PLS-POLR             ######
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


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$tt%*%res$CoeffC            
res$Yresidus <- dataY-res$YChapeau
}
##############################################



##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
res$residYChapeau <- tempregglm$linear.predictors


tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
tempConstante <- attr(res$RepY,"scaled:center")-sum(tempCoeffs*attr(res$ExpliX,"scaled:center"))
res$Coeffs <- rbind(tempConstante,tempCoeffs)

res$YChapeau <- tempregglm$fitted.values                      
res$Yresidus <- dataY-res$YChapeau
}


##############################################
######              PLS-POLR            ######
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
if (modele %in% c("pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
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
}

if(res$computed_nt==0){
  if(verbose){cat("No component could be extracted please check the data for NA only lines or columns\n")}; stop()
}


##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.PredictY | na.miss.Y)) {
if(kk==1){
  if(verbose){cat("____Predicting X without NA neither in X nor in Y____\n")}
}
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- paste("tt",1:kk,sep="")
}
else {
if (na.miss.PredictY & !na.miss.Y) {
if(kk==1){
  if(verbose){cat("____Predicting X with NA in X and not in Y____\n")}
}
res$ttPredictY <- NULL

for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%res$pp[PredictYNA[ii,],,drop=FALSE])%*%t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}

colnames(res$ttPredictY) <- paste("tt",1:kk,sep="")
}
else {
if(kk==1){
  if(verbose){cat("____There are some NAs both in X and Y____\n")}
}
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


##############################################
######              PLS-GLM             ######
##############################################


##############################################
######           PLS-GLM-POLR           ######
##############################################





##########################################
#                                        #
#          Predicting responses          #
#                                        #
##########################################


##############################################
######               PLS                ######
##############################################
if (modele == "pls") {
res$listValsPredictY <- cbind(res$listValsPredictY,attr(res$RepY,"scaled:center")+attr(res$RepY,"scaled:scale")*res$ttPredictY%*%res$CoeffC)
}


##############################################
######              PLS-GLM             ######
##############################################
if (modele %in% c("pls-glm-logistic","pls-glm-family","pls-glm-Gamma","pls-glm-gaussian","pls-glm-inverse.gaussian","pls-glm-poisson")) {
tt <- res$ttPredictY
res$listValsPredictY <- cbind(res$listValsPredictY,predict(object=tempregglm,newdata=data.frame(tt),type = "response"))
}


##############################################
######           PLS-GLM-POLR           ######
##############################################
if (modele %in% c("pls-glm-polr")) {
tts <- res$ttPredictY
if(kk==1){
res$listValsPredictY <- list(predict(tempregpolr,predict(tempregpolr, data.frame(tts=I(tts))),type="probs")) 
} else {
res$listValsPredictY <- c(res$listValsPredictY,list(predict(tempregpolr,predict(tempregpolr, data.frame(tts=I(tts))),type="probs"))) 
}
attr(res$listValsPredictY,"numberlevels") <- nlevels(dataY)
attr(res$listValsPredictY,"modele") <- modele
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


if(verbose){cat("****________________________________________________****\n")}
if(verbose){cat("\n")}
if (!keepcoeffs) {
if (!keepstd.coeffs) {return(list(valsPredict=res$listValsPredictY))} else {return(list(valsPredict=res$listValsPredictY, std.coeffs=res$Std.Coeffs))}}
else {
if (!keepstd.coeffs) {return(list(valsPredict=res$listValsPredictY, coeffs=res$Coeffs))} else {return(list(valsPredict=res$listValsPredictY, coeffs=res$Coeffs, std.coeffs=res$Std.Coeffs))}
}
}
