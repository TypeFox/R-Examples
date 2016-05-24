
AUC <- function(sens, spec)
{
.tab.res<-data.frame(se=sens, sp=spec)
.tab.res<-.tab.res[!is.na(.tab.res$sp + .tab.res$se),]
.tab.res$sp1 <- 1-.tab.res$sp
.tab.res<-.tab.res[order(.tab.res$sp1, .tab.res$se),]

return( sum((.tab.res$sp1[2:length(.tab.res$sp1)] -
  .tab.res$sp1[1:(length(.tab.res$sp1)-1)]) *
  (.tab.res$se[2:length(.tab.res$se)])) )
}


ROC <- function(status, marker, cut.values)
{

resultats<-data.frame(
 x = 1-sapply(cut.values, 
  function(x) {sum(1*(marker[status==0]<=x))/sum(status==0)} ),
 y = sapply(cut.values,
  function(x) {sum(1*(marker[status==1]>x))/sum(status==1)} ),
 z = cut.values )

resultats1<-resultats[order(resultats$x, resultats$y),]
resultats1[dim(resultats1)[1]+1,1]<-1
resultats1[dim(resultats1)[1],2]<-1

return( list(
cut.values = resultats$z,
TP = resultats$y,
FP = resultats$x,
AUC = sum((resultats1$x[2:length(resultats1$x)] - resultats1$x[1:(length(resultats1$x)-1)]) * 0.5 * (resultats1$y[2:length(resultats1$x)] + resultats1$y[1:(length(resultats1$x)-1)]))  ) )
}


boot.ROCt <- function(times, failures, features, N.boot, precision, prop, pro.time, fold.cv=5, lambda1=NULL)
{

if(!is.null(lambda1)) {

.E <- features

.G <- dim(features)[2]

.N <- dim(features)[1]

.Data <- data.frame(ident = 1:length(times), t=times, f=failures)

cut.off  <- qnorm(precision, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

pen0 <- penalized(Surv(.Data$t, .Data$f), penalized = .E,
 lambda1=lambda1, trace = FALSE)

Coef0  <- coefficients(pen0, "all")

Score0 <- as.vector(.E %*% coefficients(pen0, "all"))

if (sd(Score0)==0) {Score0 <- (Score0 - mean(Score0)) / sd(Score0)}

auc.train <- rep(-99, length(pro.time))
auc.valid <- rep(-99, length(pro.time))
auc.632 <- rep(-99, length(pro.time))
auc.632p <- rep(-99, length(pro.time))

B.train <- matrix(-99, nrow=N.boot, ncol=.G)

sd.train <- rep(-99, N.boot)
mean.train <- rep(-99, N.boot)

fn.train <- matrix(-99, nrow=N.boot, ncol=length(precision))
fp.train <- matrix(-99, nrow=N.boot, ncol=length(precision))

fn.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))
fp.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))

i <- 1

while (i <= N.boot) {
 
num <- sample(.Data$ident, size = .N, replace = TRUE)

# Calcul sur l'échantillon de Bootstrap (train)

Eb <- .E[num,] ; Data.b <- .Data[num,]

pen1 <- try(penalized(Surv(Data.b$t, Data.b$f), penalized = Eb,
 lambda1=lambda1, trace = FALSE), silent = FALSE)

if(inherits(pen1, "try-error")) { i <- i; cat(" test1 ") } else {

B.train[i,]  <- coefficients(pen1, "all")

if(sum(B.train[i,])==0)
{
fn.train[i,]<-precision
fp.train[i,]<-1-precision
fn.valid[i,]<-precision
fp.valid[i,]<-1-precision
i <- i + 1
} else {

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])

sd.train[i] <- sd(Data.b$S.l1)

mean.train[i] <- mean(Data.b$S.l1)

Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(survivalROC(Stime = Data.b$t, status = Data.b$f, marker = Data.b$S.l1, predict.time = pro.time, cut.values = cut.off, method = "NNE", lambda = prop), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test2 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i } else {
 
fn.train[i,]<-1-ROCt.l1$TP[-1]
fp.train[i,]<-ROCt.l1$FP[-1]

# Calcul sur l'échantillon de validation (valid)

Eb <- .E[-unique(sort(num)),] ; Data.b <- .Data[-unique(sort(num)),]

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])
Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(survivalROC(Stime = Data.b$t, status = Data.b$f, marker = Data.b$S.l1, predict.time = pro.time, cut.values = cut.off, method = "NNE", lambda = prop), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test3 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i; cat(" test4 ") } else {
 
fn.valid[i,]<-1-ROCt.l1$TP[-1]
fp.valid[i,]<-ROCt.l1$FP[-1]

cat("   sample ", i, "\n", sep="")

i <- i + 1

} } } } } } }

# Calcul des estimateurs "0.632"

mean.fn.train <- apply(fn.train, FUN="mean", MARGIN=2)
mean.fp.train <- apply(fp.train, FUN="mean", MARGIN=2)

mean.fn.valid <- apply(fn.valid, FUN="mean", MARGIN=2)
mean.fp.valid <- apply(fp.valid, FUN="mean", MARGIN=2)

fn.632 <- 0.368 * mean.fn.train + 0.632 * mean.fn.valid
fp.632 <- 0.368 * mean.fp.train + 0.632 * mean.fp.valid

r.fn <- (mean.fn.valid - mean.fn.train)/(precision - mean.fn.train)
r.fp <- (mean.fp.valid - mean.fp.train)/((1-precision) - mean.fp.train)

r.fn <- pmax(pmin(r.fn, 1), 0)
r.fp <- pmax(pmin(r.fp, 1), 0)

fn.632p<-(1-(0.632/(1-0.368*r.fn)))*mean.fn.train + (0.632/(1-0.368*r.fn)) * mean.fn.valid
fp.632p<-(1-(0.632/(1-0.368*r.fp)))*mean.fp.train + (0.632/(1-0.368*r.fp)) * mean.fp.valid

auc.train <- AUC(sens = 1-mean.fn.train, spec = 1-mean.fp.train)
auc.valid <- AUC(sens = 1-mean.fn.valid, spec = 1-mean.fp.valid)
auc.632   <- AUC(sens = 1-fn.632, spec = 1-fp.632)
auc.632p  <- AUC(sens = 1-fn.632p, spec = 1-fp.632p)

return( list(
Coef=Coef0,
Signature=Score0,
Lambda=lambda1,
Model=pen0,
AUC=data.frame(Apparent=auc.train,
  CV=auc.valid, s632=auc.632, p632=auc.632p),
cut.values=cut.off,
ROC.Apparent = data.frame(FNR = mean.fn.train, FPR = mean.fp.train),
ROC.CV = data.frame(FNR = mean.fn.valid, FPR = mean.fp.valid),
ROC.632 = data.frame(FNR = fn.632, FPR = fp.632),
ROC.632p = data.frame(FNR = fn.632p, FPR = fp.632p) ) )

} else {


.E <- features

.G <- dim(features)[2]

.N <- dim(features)[1]

.Data <- data.frame(ident = 1:length(times), t=times, f=failures)

cut.off  <- qnorm(precision, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

opt0 <- optL1(Surv(.Data$t, .Data$f), penalized = .E,
 fold=fold.cv, trace = FALSE)

pen0 <- penalized(Surv(.Data$t, .Data$f), penalized = .E,
 lambda1=opt0$lambda, trace = FALSE)

Coef0  <- coefficients(pen0, "all")

Score0 <- as.vector(.E %*% coefficients(pen0, "all"))

if (sd(Score0)==0) {Score0 <- (Score0 - mean(Score0)) / sd(Score0)}

auc.train <- rep(-99, length(pro.time))
auc.valid <- rep(-99, length(pro.time))
auc.632 <- rep(-99, length(pro.time))
auc.632p <- rep(-99, length(pro.time))

B.train <- matrix(-99, nrow=N.boot, ncol=.G)

sd.train <- rep(-99, N.boot)
mean.train <- rep(-99, N.boot)

fn.train <- matrix(-99, nrow=N.boot, ncol=length(precision))
fp.train <- matrix(-99, nrow=N.boot, ncol=length(precision))

fn.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))
fp.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))

i <- 1

while (i <= N.boot) {
 
num <- sample(.Data$ident, size = .N, replace = TRUE)

# Calcul sur l'échantillon de Bootstrap (train)

Eb <- .E[num,] ; Data.b <- .Data[num,]

opt1 <- try(optL1(Surv(Data.b$t, Data.b$f), penalized = Eb,
 fold=fold.cv, trace = FALSE), silent = FALSE)

if(inherits(opt1, "try-error")) { i <- i; cat(" test0 ") } else {

B.train[i,]  <- coefficients(opt1$fullfit,"all")

if(sum(B.train[i,])==0)
{
fn.train[i,]<-precision
fp.train[i,]<-1-precision
fn.valid[i,]<-precision
fp.valid[i,]<-1-precision
i <- i + 1
} else {

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])

sd.train[i] <- sd(Data.b$S.l1)

mean.train[i] <- mean(Data.b$S.l1)

Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(survivalROC(Stime = Data.b$t, status = Data.b$f, marker = Data.b$S.l1, predict.time = pro.time, cut.values = cut.off, method = "NNE", lambda = prop), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test2 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i } else {
 
fn.train[i,]<-1-ROCt.l1$TP[-1]
fp.train[i,]<-ROCt.l1$FP[-1]

# Calcul sur l'échantillon de validation (valid)

Eb <- .E[-unique(sort(num)),] ; Data.b <- .Data[-unique(sort(num)),]

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])
Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(survivalROC(Stime = Data.b$t, status = Data.b$f, marker = Data.b$S.l1, predict.time = pro.time, cut.values = cut.off, method = "NNE", lambda = prop), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test3 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i; cat(" test4 ") } else {
 
fn.valid[i,]<-1-ROCt.l1$TP[-1]
fp.valid[i,]<-ROCt.l1$FP[-1]

cat("   sample ", i, "\n", sep="")

i <- i + 1

} } } } } } } 

# Calcul des estimateurs "0.632"

mean.fn.train <- apply(fn.train, FUN="mean", MARGIN=2)
mean.fp.train <- apply(fp.train, FUN="mean", MARGIN=2)

mean.fn.valid <- apply(fn.valid, FUN="mean", MARGIN=2)
mean.fp.valid <- apply(fp.valid, FUN="mean", MARGIN=2)

fn.632 <- 0.368 * mean.fn.train + 0.632 * mean.fn.valid
fp.632 <- 0.368 * mean.fp.train + 0.632 * mean.fp.valid

r.fn <- (mean.fn.valid - mean.fn.train)/(precision - mean.fn.train)
r.fp <- (mean.fp.valid - mean.fp.train)/((1-precision) - mean.fp.train)

r.fn <- pmax(pmin(r.fn, 1), 0)
r.fp <- pmax(pmin(r.fp, 1), 0)

fn.632p<-(1-(0.632/(1-0.368*r.fn)))*mean.fn.train + (0.632/(1-0.368*r.fn)) * mean.fn.valid
fp.632p<-(1-(0.632/(1-0.368*r.fp)))*mean.fp.train + (0.632/(1-0.368*r.fp)) * mean.fp.valid

auc.train <- AUC(sens = 1-mean.fn.train, spec = 1-mean.fp.train)
auc.valid <- AUC(sens = 1-mean.fn.valid, spec = 1-mean.fp.valid)
auc.632   <- AUC(sens = 1-fn.632, spec = 1-fp.632)
auc.632p  <- AUC(sens = 1-fn.632p, spec = 1-fp.632p)

return( list(
Coef=Coef0,
Signature=Score0,
Lambda=opt0$lambda,
Model=pen0,
AUC=data.frame(Apparent=auc.train,
  CV=auc.valid, s632=auc.632, p632=auc.632p),
cut.values=cut.off,
ROC.Apparent = data.frame(FNR = mean.fn.train, FPR = mean.fp.train),
ROC.CV = data.frame(FNR = mean.fn.valid, FPR = mean.fp.valid),
ROC.632 = data.frame(FNR = fn.632, FPR = fp.632),
ROC.632p = data.frame(FNR = fn.632p, FPR = fp.632p) ) )
}

}










boot.ROC <- function(status, features, N.boot, precision, fold.cv=5, lambda1=NULL)
{

if(!is.null(lambda1)) {

.E <- features

.G <- dim(features)[2]

.N <- dim(features)[1]

.Data <- data.frame(ident = 1:length(status), f=status)

cut.off  <- qnorm(precision, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

pen0 <- try(penalized(.Data$f, penalized = .E,
 lambda1=lambda1, model="logistic", trace = FALSE), silent = FALSE)
 
Coef0 <- coefficients(pen0, "all")[-1]

Score0 <- as.vector(.E %*% coefficients(pen0, "all")[-1])

if (sd(Score0)==0) {Score0 <- (Score0 - mean(Score0)) / sd(Score0)}

B.train <- matrix(-99, nrow=N.boot, ncol=.G)

sd.train <- rep(-99, N.boot)
mean.train <- rep(-99, N.boot)

fn.train <- matrix(-99, nrow=N.boot, ncol=length(precision))
fp.train <- matrix(-99, nrow=N.boot, ncol=length(precision))

fn.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))
fp.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))

i <- 1

while (i <= N.boot) {
 
num <- sample(.Data$ident, size = .N, replace = TRUE)

# Calcul sur l'échantillon de Bootstrap (train)

Eb <- .E[num,] ; Data.b <- .Data[num,]

pen1 <- try(penalized(Data.b$f, penalized = Eb,
 lambda1=lambda1, model="logistic", trace = FALSE), silent = FALSE)
 
if(inherits(pen1, "try-error")) { i <- i; cat(" test1 ") } else {

B.train[i,]  <- coefficients(pen1, "all")[-1]

if(sum(B.train[i,])==0)
{
fn.train[i,]<-precision
fp.train[i,]<-1-precision
fn.valid[i,]<-precision
fp.valid[i,]<-1-precision
i <- i + 1
} else {

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])

sd.train[i] <- sd(Data.b$S.l1)

mean.train[i] <- mean(Data.b$S.l1)

Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(ROC(status = Data.b$f, marker = Data.b$S.l1, cut.values = cut.off), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test2 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i } else {
 
fn.train[i,]<-1-ROCt.l1$TP
fp.train[i,]<-ROCt.l1$FP

# Calcul sur l'échantillon de validation (valid)

Eb <- .E[-unique(sort(num)),] ; Data.b <- .Data[-unique(sort(num)),]

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])
Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(ROC(status = Data.b$f, marker = Data.b$S.l1, cut.values = cut.off), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test3 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i; cat(" test4 ") } else {
 
fn.valid[i,]<-1-ROCt.l1$TP
fp.valid[i,]<-ROCt.l1$FP

cat("   sample ", i, "\n", sep="")

i <- i + 1

} } } } } } }

# Calcul des estimateurs "0.632"

mean.fn.train <- apply(fn.train, FUN="mean", MARGIN=2)
mean.fp.train <- apply(fp.train, FUN="mean", MARGIN=2)

mean.fn.valid <- apply(fn.valid, FUN="mean", MARGIN=2)
mean.fp.valid <- apply(fp.valid, FUN="mean", MARGIN=2)

fn.632 <- 0.368 * mean.fn.train + 0.632 * mean.fn.valid
fp.632 <- 0.368 * mean.fp.train + 0.632 * mean.fp.valid

r.fn <- (mean.fn.valid - mean.fn.train)/(precision - mean.fn.train)
r.fp <- (mean.fp.valid - mean.fp.train)/((1-precision) - mean.fp.train)

r.fn <- pmax(pmin(r.fn, 1), 0)
r.fp <- pmax(pmin(r.fp, 1), 0)

fn.632p<-(1-(0.632/(1-0.368*r.fn)))*mean.fn.train + (0.632/(1-0.368*r.fn)) * mean.fn.valid
fp.632p<-(1-(0.632/(1-0.368*r.fp)))*mean.fp.train + (0.632/(1-0.368*r.fp)) * mean.fp.valid

auc.train <- AUC(sens = 1-mean.fn.train, spec = 1-mean.fp.train)
auc.valid <- AUC(sens = 1-mean.fn.valid, spec = 1-mean.fp.valid)
auc.632  <- AUC(sens = 1-fn.632, spec = 1-fp.632)
auc.632p  <- AUC(sens = 1-fn.632p, spec = 1-fp.632p)

return( list(
Coef=Coef0,
Signature=Score0,
Lambda=lambda1,
Model=pen0,
AUC=data.frame(Apparent=auc.train,
  CV=auc.valid, s632=auc.632, p632=auc.632p),
cut.values=cut.off,
ROC.Apparent = data.frame(FNR = mean.fn.train, FPR = mean.fp.train),
ROC.CV = data.frame(FNR = mean.fn.valid, FPR = mean.fp.valid),
ROC.632 = data.frame(FNR = fn.632, FPR = fp.632),
ROC.632p = data.frame(FNR = fn.632p, FPR = fp.632p) ) )

} else {

.E <- features

.G <- dim(features)[2]

.N <- dim(features)[1]

.Data <- data.frame(ident = 1:length(status), f=status)

cut.off  <- qnorm(precision, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

opt0 <- try(optL1(.Data$f, penalized = .E,
 fold=fold.cv, model="logistic", trace = FALSE), silent = FALSE)

pen0 <- try(penalized(.Data$f, penalized = .E,
 lambda1=opt0$lambda, model="logistic", trace = FALSE), silent = FALSE)
 
Coef0 <- coefficients(pen0, "all")[-1]

Score0 <- as.vector(.E %*% coefficients(pen0, "all")[-1])

if (sd(Score0)==0) {Score0 <- (Score0 - mean(Score0)) / sd(Score0)}

B.train <- matrix(-99, nrow=N.boot, ncol=.G)

sd.train <- rep(-99, N.boot)
mean.train <- rep(-99, N.boot)

fn.train <- matrix(-99, nrow=N.boot, ncol=length(precision))
fp.train <- matrix(-99, nrow=N.boot, ncol=length(precision))

fn.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))
fp.valid<-matrix(-99, nrow=N.boot, ncol=length(precision))

i <- 1

while (i <= N.boot) {
 
num <- sample(.Data$ident, size = .N, replace = TRUE)

# Calcul sur l'échantillon de Bootstrap (train)

Eb <- .E[num,] ; Data.b <- .Data[num,]

opt1 <- try(optL1(Data.b$f, penalized = Eb,
 fold=fold.cv, model="logistic", trace = FALSE), silent = FALSE)

if(inherits(opt1, "try-error")) { i <- i; cat(" test0 ") } else {

B.train[i,]  <- coefficients(opt0$fullfit,"all")[-1]

if(sum(B.train[i,])==0)
{
fn.train[i,]<-precision
fp.train[i,]<-1-precision
fn.valid[i,]<-precision
fp.valid[i,]<-1-precision
i <- i + 1
} else {

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])

sd.train[i] <- sd(Data.b$S.l1)

mean.train[i] <- mean(Data.b$S.l1)

Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(ROC(status = Data.b$f, marker = Data.b$S.l1, cut.values = cut.off), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test2 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i } else {
 
fn.train[i,]<-1-ROCt.l1$TP
fp.train[i,]<-ROCt.l1$FP

# Calcul sur l'échantillon de validation (valid)

Eb <- .E[-unique(sort(num)),] ; Data.b <- .Data[-unique(sort(num)),]

Data.b$S.l1 <- as.vector(Eb %*% B.train[i,])
Data.b$S.l1 <- (Data.b$S.l1 - mean.train[i]) / sd.train[i]

ROCt.l1 <- try(ROC(status = Data.b$f, marker = Data.b$S.l1, cut.values = cut.off), silent = TRUE)

if(inherits(ROCt.l1, "try-error")) { i <- i; cat(" test3 ") } else {

if(length(ROCt.l1$TP)==0 | is.na(sum(ROCt.l1$TP[-1])) | is.na(sum(ROCt.l1$FP[-1])))
 { i <- i; cat(" test4 ") } else {
 
fn.valid[i,]<-1-ROCt.l1$TP
fp.valid[i,]<-ROCt.l1$FP

cat("   sample ", i, "\n", sep="")

i <- i + 1

} } } } } } }

# Calcul des estimateurs "0.632"

mean.fn.train <- apply(fn.train, FUN="mean", MARGIN=2)
mean.fp.train <- apply(fp.train, FUN="mean", MARGIN=2)

mean.fn.valid <- apply(fn.valid, FUN="mean", MARGIN=2)
mean.fp.valid <- apply(fp.valid, FUN="mean", MARGIN=2)

fn.632 <- 0.368 * mean.fn.train + 0.632 * mean.fn.valid
fp.632 <- 0.368 * mean.fp.train + 0.632 * mean.fp.valid

r.fn <- (mean.fn.valid - mean.fn.train)/(precision - mean.fn.train)
r.fp <- (mean.fp.valid - mean.fp.train)/((1-precision) - mean.fp.train)

r.fn <- pmax(pmin(r.fn, 1), 0)
r.fp <- pmax(pmin(r.fp, 1), 0)

fn.632p<-(1-(0.632/(1-0.368*r.fn)))*mean.fn.train + (0.632/(1-0.368*r.fn)) * mean.fn.valid
fp.632p<-(1-(0.632/(1-0.368*r.fp)))*mean.fp.train + (0.632/(1-0.368*r.fp)) * mean.fp.valid

auc.train <- AUC(sens = 1-mean.fn.train, spec = 1-mean.fp.train)
auc.valid <- AUC(sens = 1-mean.fn.valid, spec = 1-mean.fp.valid)
auc.632  <- AUC(sens = 1-fn.632, spec = 1-fp.632)
auc.632p  <- AUC(sens = 1-fn.632p, spec = 1-fp.632p)

return( list(
Coef=Coef0,
Signature=Score0,
Lambda=opt0$lambda,
Model=pen0,
AUC=data.frame(Apparent=auc.train,
  CV=auc.valid, s632=auc.632, p632=auc.632p),
cut.values=cut.off,
ROC.Apparent = data.frame(FNR = mean.fn.train, FPR = mean.fp.train),
ROC.CV = data.frame(FNR = mean.fn.valid, FPR = mean.fp.valid),
ROC.632 = data.frame(FNR = fn.632, FPR = fp.632),
ROC.632p = data.frame(FNR = fn.632p, FPR = fp.632p) ) )

}

}
