PLS_lm_formula <- function(formula,data=NULL,nt=2,limQ2set=.0975,dataPredictY=dataX,modele="pls",family=NULL,typeVC="none",EstimXNA=FALSE,scaleX=TRUE,scaleY=NULL,pvals.expli=FALSE,alpha.pvals.expli=.05,MClassed=FALSE,tol_Xi=10^(-12),weights,subset,contrasts=NULL,sparse=FALSE,sparseStop=FALSE,naive=FALSE,verbose=TRUE) {

##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

if(verbose){cat("____************************************************____\n")}
if(sparse){sparseStop=TRUE}
if(sparseStop){pvals.expli=TRUE}

if (missing(data)) {data <- environment(formula)}
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf$na.action <- na.pass    
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
#if (identical(method, "model.frame")) return(mf)
mt <- attr(mf, "terms")
attr(mt,"intercept")<-0L
dataY <- model.response(mf, "any")
if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}
if(!NoWeights){naive=TRUE; if(verbose){cat(paste("Only naive DoF can be used with weighted PLS\n",sep=""))}} else {NoWeights=TRUE}
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
if (!(modele %in% c("pls"))) {break}
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

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,press.ind=NULL,press.tot=NULL,family=family,ttPredictY = NULL,typeVC=typeVC,dataX=dataX,dataY=dataY)
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
if (typeVC == "none") {} else {
if(NoWeights){
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
}
rm(i)
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2) 
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
if (typeVC == "missingdata") {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                if(rcond(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])<tol_Xi) {
                break_nt_vc <- TRUE; res$computed_nt_vc <- kk-1
                if(verbose){cat(paste("Warning : reciprocal condition number of t(temppp.cv[XXNA[",i,",],,drop=FALSE])%*%temppp.cv[XXNA[",i,",],,drop=FALSE]) < 10^{-12}\n",sep=""))}
                if(verbose){cat(paste("Please choose a smaller number of components to cross validate\n",sep=""))}
                break
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])%*%t(temppp.cv[XXNA[i,],,drop=FALSE])%*%(ExpliXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
rm(i)
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
}
else {
if (typeVC %in% c("standard","adaptative")) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ]*weights[-i])%*%YwotNA[-i, ]/(t(XXNA[-i, ]*weights[-i])%*%YwotNA[-i, ]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
}
rm(i)
res$press.ind <- cbind(res$press.ind,weights*(YwotNA-temppred)^2) 
res$press.ind2 <- cbind(res$press.ind2,weights*(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
if (typeVC == "missingdata") {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ]*weights[-i])%*%YwotNA[-i, ]/(t(XXNA[-i, ]*weights[-i])%*%YwotNA[-i, ]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                if(rcond(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])<tol_Xi) {
                break_nt_vc <- TRUE; res$computed_nt_vc <- kk-1
                if(verbose){cat(paste("Warning : reciprocal condition number of t(temppp.cv[XXNA[",i,",],,drop=FALSE])%*%temppp.cv[XXNA[",i,",],,drop=FALSE]) < 10^{-12}\n",sep=""))}
                if(verbose){cat(paste("Please choose a smaller number of components to cross validate\n",sep=""))}
                break
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])%*%t(temppp.cv[XXNA[i,],,drop=FALSE])%*%(ExpliXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
rm(i)
res$press.ind <- cbind(res$press.ind,weights*(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,weights*(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
}
if (!(typeVC %in% c("standard","adaptative","missingdata"))) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____unknown____\n"))}
}
}
}
res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY)
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY,weights*RepY)
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
if (typeVC == "none") {} else {
if(NoWeights){
if (typeVC %in% c("standard","missingdata")) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                if(rcond(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])<tol_Xi) {
                break_nt_vc <- TRUE; res$computed_nt_vc <- kk-1
                if(verbose){cat(paste("Warning : reciprocal condition number of t(temppp.cv[XXNA[",i,",],,drop=FALSE])%*%temppp.cv[XXNA[",i,",],,drop=FALSE]) < 10^{-12}\n",sep=""))}
                if(verbose){cat(paste("Please choose a smaller number of components to cross validate\n",sep=""))}
                break
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])%*%t(temppp.cv[XXNA[i,],,drop=FALSE])%*%(XXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
if (typeVC == "adaptative") {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
if (all(XXNA[i,])) {
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
}
else {
                tempww.cv <- t(XXwotNA[-i, ])%*%YwotNA[-i]/(t(XXNA[-i, ])%*%YwotNA[-i]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                if(rcond(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])<tol_Xi) {
                break_nt_vc <- TRUE; res$computed_nt_vc <- kk-1
                if(verbose){cat(paste("Warning : reciprocal condition number of t(temppp.cv[XXNA[",i,",],,drop=FALSE])%*%temppp.cv[XXNA[",i,",],,drop=FALSE]) < 10^{-12}\n",sep=""))}
                if(verbose){cat(paste("Please choose a smaller number of components to cross validate\n",sep=""))}
                break
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])%*%t(temppp.cv[XXNA[i,],,drop=FALSE])%*%(XXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
                }
}
res$press.ind <- cbind(res$press.ind,(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
if (typeVC %in% c("standard","missingdata","adaptative")) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____unknown____\n"))}
}
}
}
else {
if (typeVC %in% c("standard","missingdata")) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
                tempww.cv <- t(XXwotNA[-i, ]*weights[-i])%*%YwotNA[-i, ]/(t(XXNA[-i, ]*weights[-i])%*%YwotNA[-i, ]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                if(rcond(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])<tol_Xi) {
                break_nt_vc <- TRUE; res$computed_nt_vc <- kk-1
                if(verbose){cat(paste("Warning : reciprocal condition number of t(temppp.cv[XXNA[",i,",],,drop=FALSE])%*%temppp.cv[XXNA[",i,",],,drop=FALSE]) < 10^{-12}\n",sep=""))}
                if(verbose){cat(paste("Please choose a smaller number of components to cross validate\n",sep=""))}
                break
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])%*%t(temppp.cv[XXNA[i,],,drop=FALSE])%*%(XXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
}
res$press.ind <- cbind(res$press.ind,weights*(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,weights*(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
if (typeVC == "adaptative") {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____\n"))}
}
temppp.cv <- res$pp  
temppred <- rep(0, res$nr)
for (i in 1:(res$nr)) {     
if (all(XXNA[i,])) {
                tempww.cv <- t(XXwotNA[-i, ]*weights[-i])%*%YwotNA[-i, ]/(t(XXNA[-i, ]*weights[-i])%*%YwotNA[-i, ]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])
                tempc.cv <- as.vector(tempc.cv)
                temppred[i] <- tempc.cv * XXwotNA[i, ] %*% tempwwnorm.cv
}
else {
                tempww.cv <- t(XXwotNA[-i, ]*weights[-i])%*%YwotNA[-i, ]/(t(XXNA[-i, ]*weights[-i])%*%YwotNA[-i, ]^2)
                tempwwnorm.cv <- tempww.cv/sqrt(drop(crossprod(tempww.cv)))
                temptt.cv <- XXwotNA[-i, ]%*%tempwwnorm.cv/(XXNA[-i, ]%*%(tempwwnorm.cv^2))
                tempc.cv <- solve(t(temptt.cv)%*%temptt.cv)%*%t(temptt.cv)%*%(YwotNA[-i])    
                for (jj in 1:(res$nc)) {
                     temppp.cv[jj,kk] <- crossprod(temptt.cv,(XXwotNA[-i,])[,jj])/drop(crossprod((XXNA[-i, ])[,jj],temptt.cv^2))
                }
                if(rcond(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])<tol_Xi) {
                break_nt_vc <- TRUE; res$computed_nt_vc <- kk-1
                if(verbose){cat(paste("Warning : reciprocal condition number of t(temppp.cv[XXNA[",i,",],,drop=FALSE])%*%temppp.cv[XXNA[",i,",],,drop=FALSE]) < 10^{-12}\n",sep=""))}
                if(verbose){cat(paste("Please choose a smaller number of components to cross validate\n",sep=""))}
                break
                }
                ttPredictY.cv <- (solve(t(temppp.cv[XXNA[i,],,drop=FALSE])%*%temppp.cv[XXNA[i,],,drop=FALSE])%*%t(temppp.cv[XXNA[i,],,drop=FALSE])%*%(XXwotNA[i,])[XXNA[i,]])[kk]
                temppred[i] <- tempc.cv*ttPredictY.cv   
                }
}
res$press.ind <- cbind(res$press.ind,weights*(YwotNA-temppred)^2)
res$press.ind2 <- cbind(res$press.ind2,weights*(dataYwotNA-res$YChapeau-attr(res$RepY,"scaled:scale")*temppred)^2)
}
if (!(typeVC %in% c("standard","missingdata","adaptative"))) {
if (kk==1) {
  if(verbose){cat(paste("____TypeVC____",typeVC,"____unknown____\n"))}
}
}
}
}

res$residYChapeau <- res$tt%*%tempCoeffC
if (kk==1) {
if(NoWeights){
res$RSSresidY <- crossprod(RepY)
}
if(!NoWeights){
res$RSSresidY <- crossprod(RepY,weights*RepY)
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

if(verbose){cat("____Component____",kk,"____\n")}
if(break_nt_vc==TRUE) {break}
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


##############################################
#                                            #
#           Predicting components            #
#                                            #
##############################################

if (!(na.miss.PredictY | na.miss.Y)) {
  if(verbose){cat("____Predicting X without NA neither in X nor in Y____\n")}
res$ttPredictY <- PredictYwotNA%*%res$wwetoile 
colnames(res$ttPredictY) <- paste("Comp_",1:res$computed_nt,sep="")
}
else {
if (na.miss.PredictY & !na.miss.Y) {
  if(verbose){cat("____Predicting X with NA in X and not in Y____\n")}
for (ii in 1:nrow(PredictYwotNA)) {  
      res$ttPredictY <- rbind(res$ttPredictY,t(solve(t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%res$pp[PredictYNA[ii,],,drop=FALSE])%*%t(res$pp[PredictYNA[ii,],,drop=FALSE])%*%(PredictYwotNA[ii,])[PredictYNA[ii,]]))
}
colnames(res$ttPredictY) <- paste("Comp_",1:res$computed_nt,sep="")
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

if (typeVC %in% c("standard","missingdata","adaptative")) {
if(NoWeights){
res$press.tot <- colSums(res$press.ind)
res$press.tot2 <- colSums(res$press.ind2)
}
if(!NoWeights){
res$press.tot <- colSums(res$press.ind*weights)
res$press.tot2 <- colSums(res$press.ind2*weights)
}
res$Q2 <- 1-res$press.tot/res$RSSresidY[-(res$computed_nt+1)]
res$limQ2 <- rep(limQ2set,res$computed_nt)
res$Q2_2 <- 1-res$press.tot2/res$RSS[-(res$computed_nt+1)]


for (k in 1:res$computed_nt) {res$Q2cum[k] <- prod(res$press.tot[1:k])/prod(res$RSSresidY[1:k])}
res$Q2cum <- 1 - res$Q2cum
for (k in 1:res$computed_nt) {res$Q2cum_2[k] <- prod(res$press.tot2[1:k])/prod(res$RSS[1:k])}
res$Q2cum_2 <- 1 - res$Q2cum_2
if (MClassed==FALSE) {
res$InfCrit <- t(rbind(res$AIC,c(NA,res$Q2cum_2), c(NA,res$limQ2), c(NA,res$Q2_2[1:res$computed_nt]), c(NA,res$press.tot2[1:res$computed_nt]), res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY, c(NA,res$press.tot), c(NA,res$Q2), c(NA,res$limQ2), c(NA,res$Q2cum), res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "PRESS_residY", "Q2_residY", "LimQ2", "Q2cum_residY", "AIC.std"))
} else {
res$InfCrit <- t(rbind(res$AIC,c(NA,res$Q2cum_2), c(NA,res$limQ2), c(NA,res$Q2_2[1:res$computed_nt]), c(NA,res$press.tot2[1:res$computed_nt]), res$RSS, c(NA,res$R2), res$MissClassed, c(NA,res$R2residY), res$RSSresidY, c(NA,res$press.tot), c(NA,res$Q2), c(NA,res$limQ2), c(NA,res$Q2cum), res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "Q2cum_Y", "LimQ2_Y", "Q2_Y", "PRESS_Y", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "PRESS_residY", "Q2_residY", "LimQ2", "Q2cum_residY", "AIC.std"))
}
res$ic.dof<-infcrit.dof(res,naive=naive)
res$InfCrit <- cbind(res$InfCrit,res$ic.dof)
} else {
if (MClassed==FALSE) {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(NA,res$R2), c(NA,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "RSS_Y", "R2_Y", "R2_residY", "RSS_residY", "AIC.std"))
} else {
res$InfCrit <- t(rbind(res$AIC, res$RSS, c(NA,res$R2), res$MissClassed, c(NA,res$R2residY), res$RSSresidY, res$AIC.std))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "RSS_Y", "R2_Y", "MissClassed", "R2_residY", "RSS_residY", "AIC.std"))
}
res$ic.dof<-infcrit.dof(res,naive=naive)
res$InfCrit <- cbind(res$InfCrit,res$ic.dof)
}
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
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt,sep="")
rownames(res$Coeffs) <- c("Intercept",colnames(ExpliX))
}

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
