plsRcoxmodel.default <- function(Xplan,time,time2,event,type,origin,typeres="deviance", collapse, weighted, scaleX=TRUE, scaleY=TRUE, nt=min(2,ncol(Xplan)),limQ2set=.0975, dataPredictY=Xplan, pvals.expli=FALSE,alpha.pvals.expli=.05,tol_Xi=10^(-12),weights,control, sparse=FALSE,sparseStop=TRUE,allres=TRUE, verbose=TRUE,...) {

##################################################
#                                                #
#    Initialization and formatting the inputs    #
#                                                #
##################################################

if(verbose){cat("____************************************************____\n")}
dataX<-Xplan
dataY<-time
modele <- "pls-cox"
if(any(apply(is.na(dataX),MARGIN=2,"all"))){return(vector("list",0)); cat("One of the columns of Xplan is completely filled with missing data\n"); stop()}
if(any(apply(is.na(dataX),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of Xplan is completely filled with missing data\n"); stop()}
if(identical(dataPredictY,dataX)){PredYisdataX <- TRUE} else {PredYisdataX <- FALSE}
if(!PredYisdataX){
if(any(apply(is.na(dataPredictY),MARGIN=2,"all"))){return(vector("list",0)); cat("One of the columns of dataPredictY is completely filled with missing data\n"); stop()}
if(any(apply(is.na(dataPredictY),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of dataPredictY is completely filled with missing data\n"); stop()}
}
if(missing(weights)){NoWeights=TRUE} else {if(all(weights==rep(1,length(dataY)))){NoWeights=TRUE} else {NoWeights=FALSE}}
if(any(is.na(dataX))) {na.miss.X <- TRUE} else na.miss.X <- FALSE
if(any(is.na(dataY))) {na.miss.Y <- TRUE} else na.miss.Y <- FALSE
if(any(is.na(dataPredictY))) {na.miss.PredictY <- TRUE} else {na.miss.PredictY <- FALSE}
if(sparse){pvals.expli=TRUE}
if (!is.data.frame(dataX)) {dataX <- data.frame(dataX)}

if ((scaleY & missing(time2))) {if(NoWeights){RepY <- scale(dataY)} else {meanY <- weighted.mean(dataY,weights); stdevY <- sqrt((length(dataY)-1)/length(dataY)*weighted.mean((dataY-meanY)^2,weights)); RepY <- (dataY-meanY)/stdevY; attr(RepY,"scaled:center") <- meanY ; attr(RepY,"scaled:scale") <- stdevY}}
else {
    RepY <- time
    attr(RepY,"scaled:center") <- 0
    attr(RepY,"scaled:scale") <- 1
}
if (scaleX) {if(NoWeights){ExpliX <- scale(dataX)} else {meanX <- apply(dataX,2,weighted.mean,weights); stdevX <- sqrt((length(time)-1)/length(time)*apply((sweep(dataX,2,meanX))^2,2,weighted.mean,weights)); ExpliX <- sweep(sweep(dataX, 2, meanX), 2 ,stdevX, "/"); attr(ExpliX,"scaled:center") <- meanX ; attr(ExpliX,"scaled:scale") <- stdevX}
    if(PredYisdataX){PredictY <- ExpliX} else {PredictY <- sweep(sweep(dataPredictY, 2, attr(ExpliX,"scaled:center")), 2 ,attr(ExpliX,"scaled:scale"), "/")}
}
else {
    ExpliX <- dataX
    attr(ExpliX,"scaled:center") <- rep(0,ncol(dataX))
    attr(ExpliX,"scaled:scale") <- rep(1,ncol(dataX))
    PredictY <- (dataPredictY)
}

try(attachNamespace("survival"),silent=TRUE)
#on.exit(try(unloadNamespace("survival"),silent=TRUE))
mf <- match.call(expand.dots = FALSE)
m <- match(c("time", "time2", "event", "type", "origin"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf[[1L]] <- as.name("Surv")
mf$time <- RepY
YCsurv <- eval(mf, parent.frame())
attr(YCsurv,"scaled:center") <- attr(RepY,"scaled:center")
attr(YCsurv,"scaled:scale") <- attr(RepY,"scaled:scale")
RepY <- YCsurv

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

#YwotNA <- as.matrix(RepY)
#YwotNA[!YNA] <- 0
YwotNA <- YCsurv

#dataYwotNA <- as.matrix(dataY)
#dataYwotNA[!YNA] <- 0
dataYwotNA <- YCsurv

if(PredYisdataX){PredictYwotNA <- XXwotNA} else {
PredictYwotNA <- as.matrix(PredictY)
PredictYwotNA [is.na(PredictY)] <- 0
}

res <- list(nr=nrow(ExpliX),nc=ncol(ExpliX),nt=nt,ww=NULL,wwnorm=NULL,wwetoile=NULL,tt=NULL,pp=NULL,CoeffC=NULL,uscores=NULL,YChapeau=NULL,residYChapeau=NULL,RepY=RepY,na.miss.Y=na.miss.Y,YNA=YNA,residY=RepY,ExpliX=ExpliX,na.miss.X=na.miss.X,XXNA=XXNA,residXX=ExpliX,PredictY=PredictYwotNA,ttPredictY = NULL,dataX=dataX,dataY=dataY) 
if(exists("XplanFormula")){res$XplanFormula=XplanFormula} else {XplanFormula=NULL}
if(NoWeights){res$weights<-rep(1L,res$nr)} else {res$weights<-weights}
res$temppred <- NULL


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
#YwotNA <- as.matrix(res$residY)
YwotNA <- res$residY
#YwotNA[!YNA] <- 0
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
######              PLS-COX             ######
##############################################
if (modele %in% c("pls-cox")) {

mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b[[1L]] <- as.name("coxph")
if (!pvals.expli) {
XXwotNA[!XXNA] <- NA
for (jj in 1:(res$nc)) {
    mf2b$data <- data.frame(cbind(res$tt,XXwotNA[,jj,drop=FALSE]))
    mf2b$formula <- as.formula(YwotNA~.)
    tempww[jj] <- coef(eval(mf2b, parent.frame()))[kk]
}
XXwotNA[!XXNA] <- 0
rm(jj)}
else {
XXwotNA[!XXNA] <- NA
tempvalpvalstep <- rep(0,res$nc)
temppvalstep <- rep(0,res$nc)
for (jj in 1:(res$nc)) {
    mf2b$data <- data.frame(cbind(res$tt,XXwotNA[,jj,drop=FALSE]))
    mf2b$formula <- as.formula(YwotNA~.)
    tmww <- summary(eval(mf2b, parent.frame()))$coefficients[kk,]
    tempww[jj] <- tmww[1]
    tempvalpvalstep[jj] <- tmww[5]
    temppvalstep[jj] <- (tmww[5] < alpha.pvals.expli)
}
if(sparse&sparseStop){
  if(sum(temppvalstep)==0L){
    break_nt_sparse <- TRUE}
  else 
  {tempww[!temppvalstep] <- 0}}
XXwotNA[!XXNA] <- 0
rm(jj)
res$valpvalstep <- cbind(res$valpvalstep,tempvalpvalstep)
res$pvalstep <- cbind(res$pvalstep,temppvalstep)
}
}

#assign("YwotNA", YwotNA, envir=parent.frame(n=sys.nframe()))
#assign("tts", tts, envir=parent.frame(n=sys.nframe()))
#assign("XXwotNA", XXwotNA, envir=parent.frame(n=sys.nframe()))
# assign("jj", jj, envir=parent.frame(n=sys.nframe()))



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
for (ii in 1:res$nr) {
if(rcond(t(cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[ii,],,drop=FALSE])<tol_Xi) {
break_nt <- TRUE; res$computed_nt <- kk-1
if(verbose){cat(paste("Warning : reciprocal condition number of t(cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE])%*%cbind(res$pp,temppp)[XXNA[",ii,",],,drop=FALSE] < 10^{-12}\n",sep=""))}
if(verbose){cat(paste("Warning only ",res$computed_nt," components could thus be extracted\n",sep=""))}
break
}
}
rm(ii)
if(break_nt==TRUE) {break}
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
if(break_nt==TRUE) {break}
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
######              PLS-COX            ######
##############################################
if (modele %in% c("pls-cox")) {
if (kk==1) {
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YCsurv~1)
mf2b[[1L]] <- as.name("coxph")
tempconstcox <- eval(mf2b, parent.frame())
res$AIC <- suppressWarnings(extractAIC(tempconstcox)[2])
res$BIC <- suppressWarnings(extractAIC(tempconstcox, k = log(res$nr))[2])
res$Coeffsmodel_vals <- rbind(rep(NA,5),matrix(rep(NA,5*nt),ncol=5))
rm(tempconstcox)
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YwotNA~.)
#tt<-data.frame(res$tt); colnames(tt) <- paste("Comp_",1:length(tt),sep="")
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
#mf2b$data <- tt
mf2b$data <- tttrain  
mf2b$model <- TRUE
mf2b[[1L]] <- as.name("coxph")
tempregcox <- eval(mf2b, parent.frame())
rm(tttrain)
res$AIC <- cbind(res$AIC,extractAIC(tempregcox)[2])
res$BIC <- cbind(res$BIC,extractAIC(tempregcox, k = log(res$nr))[2])
res$Coeffsmodel_vals <- cbind(rbind(summary(tempregcox)$coefficients,matrix(rep(NA,5*(nt-kk)),ncol=5)))
tempCoeffC <- as.vector(coef(tempregcox))
res$CoeffCFull <- matrix(c(tempCoeffC,rep(NA,nt-kk)),ncol=1)
} else {
if (!(na.miss.X | na.miss.Y)) {
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YwotNA~.)
#tt<-data.frame(res$tt); colnames(tt) <- paste("Comp_",1:length(tt),sep="")
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
#mf2b$data <- tt
mf2b$data <- tttrain  
mf2b$model <- TRUE
mf2b[[1L]] <- as.name("coxph")
tempregcox <- eval(mf2b, parent.frame())
rm(tttrain)
res$AIC <- cbind(res$AIC,extractAIC(tempregcox)[2])
res$BIC <- cbind(res$BIC,extractAIC(tempregcox, k = log(res$nr))[2])
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregcox)$coefficients,matrix(rep(NA,5*(nt-kk)),ncol=5)))
tempCoeffC <- as.vector(coef(tempregcox))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
else
{
mf2b <- match.call(expand.dots = TRUE)
m2b <- match(c(head(names(as.list(args(coxph))),-2),head(names(as.list(args((coxph.control)))),-1)), names(mf2b), 0L)
mf2b <- mf2b[c(1L, m2b)]
mf2b$formula <- as.formula(YwotNA~.)
#tt<-data.frame(res$tt); colnames(tt) <- paste("Comp_",1:length(tt),sep="")
tttrain<-data.frame(YwotNA=YwotNA,tt=res$tt)
#mf2b$data <- tt
mf2b$data <- tttrain  
mf2b$model <- TRUE
mf2b[[1L]] <- as.name("coxph")
tempregcox <- eval(mf2b, parent.frame())
rm(tttrain)
res$AIC <- cbind(res$AIC,extractAIC(tempregcox)[2])
res$BIC <- cbind(res$BIC,extractAIC(tempregcox, k = log(res$nr))[2])
res$Coeffsmodel_vals <- cbind(res$Coeffsmodel_vals,rbind(summary(tempregcox)$coefficients,matrix(rep(NA,5*(nt-kk)),ncol=5)))
tempCoeffC <- as.vector(coef(tempregcox))  
res$CoeffCFull <- cbind(res$CoeffCFull,c(tempCoeffC,rep(NA,nt-kk)))
}
}

res$wwetoile <- (res$wwnorm)%*%solve(t(res$pp)%*%res$wwnorm)
res$CoeffC <- tempCoeffC
res$Std.Coeffs <- res$wwetoile%*%res$CoeffC
rownames(res$Std.Coeffs) <- colnames(ExpliX)
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
######              PLS-COX            ######
##############################################
if (modele %in% c("pls-cox")) {
res$residYChapeau <- tempregcox$linear.predictors

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
res$Coeffs <- tempCoeffs

res$YChapeau <- as.matrix(predict(tempregcox, type='expected'))            
res$Yresidus <- residuals(tempregcox,type="martingale")
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
######              PLS-COX            ######
##############################################
if (modele %in% c("pls-cox")) {
res$residYChapeau <- tempregcox$linear.predictors

tempCoeffs <- res$wwetoile%*%res$CoeffC*attr(res$RepY,"scaled:scale")/attr(res$ExpliX,"scaled:scale")
res$Coeffs <- tempCoeffs

res$YChapeau <- as.matrix(predict(tempregcox, type='expected'))            
res$Yresidus <- residuals(tempregcox,type="martingale")
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
######              PLS-COX            ######
##############################################
if (modele %in% c("pls-cox")) {
res$residY <- res$residY 
res$residusY <- cbind(res$residusY,res$residY)

rm(tempww)
rm(tempwwnorm)
rm(temptt)
rm(temppp)
rm(tempCoeffC)
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
  if(verbose){cat("No component could be extracted please check the data for NA only lines or columns\n")}; stop()
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
  if(verbose){cat("____Predicting X without NA neither in X nor in Y____\n")}
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
######              PLS-COX            ######
##############################################
if (modele %in% c("pls-cox")) {
res$InfCrit <- t(rbind(res$AIC, res$BIC))
dimnames(res$InfCrit) <- list(paste("Nb_Comp_",0:res$computed_nt,sep=""), c("AIC", "BIC"))
}






##########################################
#                                        #
#          Predicting responses          #
#                                        #
##########################################



##############################################
######              PLS-COX            ######
##############################################
if (modele %in% c("pls-cox")) {
res$YChapeau <- as.matrix(predict(tempregcox, type='expected'))            
rownames(res$YChapeau) <- rownames(ExpliX)

ttpred <- data.frame(tt=res$ttPredictY)
res$Std.ValsPredictY <- predict(tempregcox,newdata=ttpred, type = "lp")
res$ValsPredictY <- predict(tempregcox,newdata=ttpred,type = "risk")

res$Std.XChapeau <- res$tt%*%t(res$pp)
rownames(res$Std.XChapeau) <- rownames(ExpliX)
names(res$CoeffC) <- paste("Coeff_Comp_Reg",1:res$computed_nt)
rownames(res$Coeffs) <- colnames(ExpliX)
res$FinalModel <- tempregcox
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
class(res) <- "plsRcoxmodel"
return(res)
}
