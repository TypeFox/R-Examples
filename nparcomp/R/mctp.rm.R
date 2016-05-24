mctp.rm <-
function (formula, data, type = c("UserDefined", "Tukey", "Dunnett", "Sequen",
    "Williams", "Changepoint", "AVE", "McDermott", "Marcus","UmbrellaWilliams"),
    control = NULL, conf.level = 0.95, alternative = c("two.sided",
        "lower", "greater"), rounds = 3, correlation = FALSE,
    asy.method = c("fisher", "normal", "mult.t"), plot.simci = FALSE,
    info = TRUE, contrast.matrix = NULL) {

    input.list <- list(formula = formula, data = data, type = type[1],
                   conf.level=conf.level, alternative=alternative,
                   asy.method=asy.method, plot.simci=plot.simci,
                   control=control, info=info, rounds=rounds, 
                   contrast.matrix=contrast.matrix, correlation=correlation)
    conflevel<-conf.level

#-------------------------Quality Checks---------------------------------------#
     if (conflevel >= 1 || conflevel <= 0) {
        stop("The confidence level must be between 0 and 1!")
        if (is.null(alternative)) {
            stop("Please declare the alternative! (two.sided, lower, greater)")
        }
    }
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    asy.method <- match.arg(asy.method)
    if (length(formula) != 3) {
        stop("You can only analyse one-way layouts!")
    }
    dat <- model.frame(formula, data)
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    response <- dat[, 1]

    factorx <- as.factor(dat[, 2])

    fl <- levels(factorx)

    a <- nlevels(factorx)
    if (a <= 2) {
        stop("You want to perform a two-sample test. Please use the function npar.t.test.paired")
    }
    samples <- split(response, factorx)

    n <- sapply(samples, length)
    if (any(n == 1)) {
        warn <- paste("The factor level", fl[n == 1], "has got only one observation!")
        stop(warn)
    }

    a <- length(n)

##############contrast matrix####################
if (type=="UserDefined"){
if(is.null(contrast.matrix)){stop("Please eanter a contrast matrix!")}
ch<-contrast.matrix
rownames(ch)<-paste("C",1:nrow(ch))

colnames(ch)<-fl}

if (type !="UserDefined"){
if (is.null(control)){icon<-1}
if (!is.null(control)){icon<-which(fl==control)}
ch<-contrMat(n=n,type,base=icon) }

nc<-nrow(ch)
connames<-rownames(ch)
Con<-matrix(ch, ncol=a)
rownames(Con)<-connames
colnames(Con)<-colnames(ch)
##################################################

#########estimators#################################
tmp1<-sort(rep(1:a,a))
tmp2<-rep(1:a,a)
pairRanks<-lapply(1:(a^2),function(arg)
rank(c(samples[[tmp1[arg]]],samples[[tmp2[arg]]])))
intRanks<-lapply(1:a, function(arg)
rank(samples[[arg]]))
n<-n[[1]]
plis<-lapply(1:(a^2), function(arg)
1/n*(mean(pairRanks[[arg]][(n+1):(2*n)])-(n+1)/2) )
vec.plis<-as.numeric(as.character(plis))
pd<-c()
for (i in 1:a){ 
pd[i]<-1/a*sum((tmp2==i)*vec.plis)
  }
pd1 <- (pd == 1)
pd0 <- (pd == 0)
pd[pd1] <- 0.999
pd[pd0] <- 0.001
Zlong<-c()
for (i in 1:(a^2)){
Zlong<-cbind(Zlong,1/n*(pairRanks[[i]][1:n]-intRanks[[tmp1[i]]]-pairRanks[[i]][(n+1):(2*n)]+intRanks[[tmp2[i]]]))
}
Zquer<-1/n*colSums(Zlong)
Yd<-matrix(rep(0,a^4),nrow=a^2,ncol=a^2)

for (k in 1:n){
Yd<-Yd+1/(n-1)*(t(Zlong[k,]-t(Zquer))%*%(Zlong[k,]-t(Zquer)))
}
W<- matrix(rep(1/a*diag(a),a),nrow=a)
Vd<-W%*%Yd%*%t(W)

    logit.pd<-log(c(pd/(1-pd)))
    logit.pd.dev<-diag(1/c((pd*(1-pd))))
Lower.logit1 <-logit.pd-qnorm(conf.level)/sqrt(n)*sqrt(c(diag(logit.pd.dev%*%Vd%*%t(logit.pd.dev))))
Upper.logit1 <-logit.pd+qnorm(conf.level)/sqrt(n)*sqrt(c(diag(logit.pd.dev%*%Vd%*%t(logit.pd.dev))))
Lower.logit <- exp(Lower.logit1)/(1+exp(Lower.logit1))
Upper.logit <- exp(Upper.logit1)/(1+exp(Upper.logit1))

    corr.mat <- function(m, nc) {
        rho <- matrix(c(0), ncol = nc, nrow = nc)
        for (i in 1:nc) {
            for (j in 1:nc) {
                rho[i, j] <- m[i, j]/sqrt(m[i, i] * m[j, j])
            }
        }
        return(rho)
      }
dfT<-n-1
Cpd<-Con%*%pd
CV<-Con%*%Vd%*%t(Con)
rhobf<- corr.mat(CV,nc)
########################################################
p.adj<-c()

#------------------------Compute adjusted p-Values and SCI---------------------#
switch(
asy.method,
#----------------------Multi T-DISTRIBUTION------------------------------------#
mult.t={
T<-sqrt(n)*(Con%*%(pd-rep(1,a)*1/2))/sqrt(c(diag(CV)))
AsyMethod <- paste("Multi - T with", round(dfT,rounds), "DF")
switch(alternative,

#----------------------Two-sided alternative-----------------------------------#
two.sided={
text.Output <- paste("True differences of relative effects are less or equal than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvt(lower=-abs(T[pp]), upper= abs(T[pp]), delta=rep(0,nc), df=dfT, corr=rhobf)[1]}
crit<- qmvt(conflevel, corr = rhobf, tail = "both", df=dfT)$quantile
Lower <- Cpd - crit/sqrt(n)*sqrt(c(diag(CV)))
Upper <- Cpd + crit/sqrt(n)*sqrt(c(diag(CV)))
},
#--------------------Alternative= LOWER----------------------------------------#
less={
text.Output <- paste("True differences of relative effects are less than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvt(lower = T[pp] , upper = Inf, df=dfT,
                delta = rep(0, nc), corr = rhobf)}
crit<- qmvt(conflevel, df=dfT, corr = rhobf, tail = "lower")$quantile
Lower <- rep(-1,nc)
Upper <- Cpd + crit/sqrt(n)*sqrt(c(diag(CV)))
},
#--------------------Alternative= GREATER--------------------------------------#
greater={
text.Output <- paste("True differences of relative effects are greater than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvt(lower =-Inf , upper =T[pp],  df=dfT,
                delta = rep(0, nc), corr = rhobf)[1]}
crit<- qmvt(conflevel, corr = rhobf, df=dfT, tail = "lower")$quantile
Lower <-Cpd - crit/sqrt(n)*sqrt(c(diag(CV)))
Upper <- rep(1,nc)}
)
},
#-------------------------------Multi NORMAL-----------------------------------#
normal={
AsyMethod <- "Normal - Approximation"
  T<-sqrt(n)*(Con%*%(pd-rep(1,a)*1/2))/sqrt(c(diag(CV)))
switch(alternative,
#----------------------Two-sided alternative-----------------------------------#
two.sided={
text.Output <- paste("True differences of relative effects are less or equal than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvnorm(lower=-abs(T[pp]), upper= abs(T[pp]), mean=rep(0,nc),corr=rhobf)[1]}
crit<- qmvnorm(conflevel, corr = rhobf, tail = "both")$quantile
Lower <- Cpd - crit/sqrt(n)*sqrt(c(diag(CV)))
Upper <- Cpd + crit/sqrt(n)*sqrt(c(diag(CV)))
},
#--------------------Alternative= LOWER----------------------------------------#
less={
text.Output <- paste("True differences of relative effects are less than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvnorm(lower = T[pp] , upper = Inf,
                mean = rep(0, nc), corr = rhobf)}
crit<- qmvnorm(conflevel, corr = rhobf, tail = "lower")$quantile
Lower <- rep(-1,nc)
Upper <- Cpd + crit/sqrt(n)*sqrt(c(diag(CV)))
},
#--------------------Alternative= GREATER--------------------------------------#
greater={
text.Output <- paste("True differences of relative effects are greater than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvnorm(lower =-Inf , upper =T[pp],
                mean = rep(0, nc), corr = rhobf)}
crit<- qmvnorm(conflevel, corr = rhobf, tail = "lower")$quantile
Lower <-Cpd - crit/sqrt(n)*sqrt(c(diag(CV)))
Upper <- rep(1,nc)}
)
},
#-------------------------------FISHER-TRANS-----------------------------------#
fisher={
AsyMethod <- paste("Fisher with", round(dfT,rounds), "DF")
 Cfisher<-1/2*log((1+Cpd)/(1-Cpd))
Vfisherdev<-diag(c(1/(1-Cpd^2)))
Vfisher<-Vfisherdev%*%CV%*%t(Vfisherdev)
T<-sqrt(n)*Cfisher/sqrt(c(diag(Vfisher)))
switch( alternative,
#----------------------Two-sided alternative-----------------------------------#
two.sided={
text.Output <- paste("True differences of relative effects are less or equal than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvt(lower=-abs(T[pp]), upper= abs(T[pp]), delta=rep(0,nc),corr=rhobf, df=dfT)[1]}
crit<- qmvt(conflevel, corr = rhobf, tail = "both",df=dfT)$quantile
Lower1 <- Cfisher - crit/sqrt(n) * sqrt(c(diag(Vfisher)))
Upper1 <- Cfisher + crit/sqrt(n) * sqrt(c(diag(Vfisher)))
Lower <- (exp(2*Lower1)-1)/(exp(2*Lower1)+1)
Upper <-  (exp(2*Upper1)-1)/(exp(2*Upper1)+1)
},
#--------------------Alternative= LOWER----------------------------------------#
less={
text.Output <- paste("True differences of relative effects are less than 0")
for (pp in 1:nc) {
p.adj[pp]<-pmvt(lower = T[pp] , upper = Inf,delta = rep(0, nc), df=dfT, corr = rhobf)}
crit<- qmvt(conflevel, corr = rhobf, tail = "lower",df=dfT)$quantile
Lower <- rep(-1,nc)
Upper1 <- Cfisher + crit/sqrt(n) * sqrt(c(diag(Vfisher)))
Upper <-  (exp(2*Upper1)-1)/(exp(2*Upper1)+1)
},
#--------------------Alternative= GREATER--------------------------------------#
greater={
text.Output <- paste("True differences of relative effects are greater than 0")
for (pp in 1:nc) {
p.adj[pp]<-1-pmvt(lower =-Inf , upper =T[pp],delta = rep(0, nc), corr = rhobf,df=dfT)}
crit<- qmvnorm(conflevel, corr = rhobf, tail = "lower")$quantile
Lower1 <- Cfisher - crit/sqrt(n) * sqrt(c(diag(Vfisher)))
Lower <-  (exp(2*Lower1)-1)/(exp(2*Lower1)+1)
Upper <- rep(1,nc)}
)
}
)

data.info <- data.frame(Sample=fl, Size=n, Effect = pd,Lower=Lower.logit,Upper=Upper.logit)
Analysis.of.Relative.Effects <- data.frame(Estimator=round(Cpd,rounds), Lower=round(Lower,rounds), Upper=round(Upper,rounds),
 Statistic = round(T,rounds), p.Value=p.adj)
Overall<-data.frame(Quantile=crit, p.Value=min(p.adj))
result<-list(Data.Info=data.info, Contrast=Con, Analysis=Analysis.of.Relative.Effects, Overall=Overall)

if (plot.simci == TRUE) {
text.Ci<-paste(conflevel*100, "%", "Simultaneous Confidence Intervals")
 Lowerp<-"|"
       plot(Cpd,1:nc,xlim=c(-1,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(Lower,1:nc, pch=Lowerp,font=2,cex=2)
              points(Upper,1:nc, pch=Lowerp,font=2,cex=2)
              abline(v=0, lty=3,lwd=2)
              for (ss in 1:nc){
              polygon(x=c(Lower[ss],Upper[ss]),y=c(ss,ss),lwd=2)}
              axis(1, at = seq(-1, 1, 0.1))
              axis(2,at=1:nc,labels=connames)
                box()
 title(main=c(text.Ci, paste("Type of Contrast:",type), paste("Method:", AsyMethod )))


 }
  
 if (info == TRUE) {
        cat("\n", "#----------------Nonparametric Multiple Comparisons for relative effects---------------#", "\n","\n",
        "-", "Alternative Hypothesis: ", text.Output,"\n",
        "-", "Estimation Method: Global Pseudo ranks","\n",
        "-", "Type of Contrast", ":", type, "\n", "-", "Confidence Level:",
            conflevel*100,"%", "\n", "-", "Method", "=", AsyMethod, "\n","\n",
                  "#--------------------------------------------------------------------------------------#","\n",

            "\n")
    }
if (correlation == TRUE){
result$Covariance <- CV
result$Correlation <- rhobf
}
result$input<-input.list
result$text.Output<-text.Output
result$connames<-connames
result$AsyMethod<-AsyMethod
class(result)<-"mctp.rm"
return(result)
}
