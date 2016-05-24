npar.t.test.paired <-
function(formula, data, conf.level = 0.95, alternative = c("two.sided",
    "less", "greater"), nperm=10000, rounds = 3, info = TRUE, plot.simci=TRUE) {

input.list <- list(formula = formula, data = data, plot.simci=plot.simci,
                   conf.level=conf.level, alternative=alternative,
                   info=info, rounds=rounds, nperm=nperm)
  alpha<-1-conf.level

#------------------------------Basic Checks------------------------------------#
    if (alpha >= 1 || alpha <= 0) {
        stop("The confidence level must be between 0 and 1!")
        if (is.null(alternative)) {
            stop("Please declare the alternative! (two.sided, less, greater)")
        }
    }
    alternative <- match.arg(alternative)
    if (length(formula) != 3) {
        stop("You can only analyse one-way layouts!")
    }

#-----------------------Arrange the data---------------------------------------#
    dat <- model.frame(formula, droplevels(data))
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    response <- dat[, 1]
    #print(response)
    factorx <- as.factor(dat[, 2])
    fl <- levels(factorx)
    a <- nlevels(factorx)
if (a > 2) {stop("You want to perform a contrast test (the factor variable has more than two levels). Please use the function mctp.rm!")}
    samples <- split(response, factorx)
    n <- sapply(samples, length)
    n1<-n[1]
    n2<-n[2]
    if (any(n == 1)) {
        warn <- paste("The factor level", fl[n == 1], "has got only one observation!")
        stop(warn)
    }
    if (n1!=n2){
        warn <-"The factor levels have unequal sample sizes. For independent samples use npar.t.test!"
        stop(warn)
    }
    N <- sum(n)
    n<-n1
     #   print(n)
    n1<-n1+1
   # print(N)
  cmpid <- paste("p(", fl[1], ",", fl[2], ")", sep = "")

#-------------Brunner-Munzel test-----#
#tcrit<-qt(0.975,n-1)
x<-c(samples[[1]],samples[[2]])
#print(x)
xinverse<-c(samples[[2]],samples[[1]])
x1<-samples[[1]]
x2<-samples[[2]]
rx<-rank(x)
rxinverse<-rank(xinverse)
rx1<-rx[1:n]
rx2<-rx[n1:N]
rix1<-rank(x1)
rix2<-rank(x2)
BM1<-1/n*(rx1-rix1)
BM2<-1/n*(rx2-rix2)
BM3<-BM1-BM2
BM4<-1/(2*n)*(rx1 - rx2)
pd<-mean(BM2)
#print(pd)
m<-mean(BM3)
v<-(sum(BM3^2)-n*m^2)/(n-1)
v0<-(v==0)
v[v0]<-1/n
T<-sqrt(n)*(pd-1/2)/sqrt(v)

#-------studentized permutation test---#
if(n<=13){
  nperm=2^n
  p<-0
  for (i in 1:n){
    a<-rep(c(rep(c(i,i+n),nperm/(2^i)),rep(c(i+n,i),nperm/(2^i))),2^(i-1))
    p<-rbind(p,a)
  }
  p<-p[2:(n+1),]
  P<-matrix(p,ncol=nperm)

  xperm<-matrix(x[P],nrow=N,ncol=nperm)
  rxperm<-matrix(rx[P],nrow=N,ncol=nperm)
}
else{
  P<-matrix(nrow=n,ncol=10000)
  permu<-function(x){
    n<-length(x)
    result<-sample(c(0,1),size=n,replace=TRUE)
    return(result)
    }
  P1<-apply(P,2,permu)
  P2<-rbind(P1,P1)
  xperm<-x*P2+xinverse*(1-P2)
  rxperm<-rx*P2+rxinverse*(1-P2)
}
xperm1<-xperm[1:n,]
xperm2<-xperm[n1:N,]
rperm1<-rxperm[1:n,]
rperm2<-rxperm[n1:N,]
riperm1<-apply(xperm1,2,rank)
riperm2<-apply(xperm2,2,rank)
BMperm2<-1/n*(rperm2-riperm2)
BMperm3<-1/n*(rperm1-riperm1)-BMperm2
pdperm<-colMeans(BMperm2)
mperm3<-colMeans(BMperm3)
vperm3<-(colSums(BMperm3^2)-n*mperm3^2)/(n-1)
vperm30<-(vperm3==0)
vperm3[vperm30]<-1/n
#print(vperm3)
Tperm<-sqrt(n)*(pdperm-1/2)/sqrt(vperm3)
p1perm<-mean(Tperm<=T)
pq1<-sort(Tperm)[(floor((1-alpha/2)*nperm)+1)]
pq2<-sort(Tperm)[(floor((1-alpha)*nperm)+1)]

switch(alternative,
#----------------------Two-sided alternative-----------------------------------#
two.sided={
text.H0<-paste("p=1/2")
text.Output <- paste("True relative effect is less or greater than 1/2")
BM<-(2*min(pt(T,n-1),1-pt(T,n-1)))
PERM<-min(2*p1perm,2*(1-p1perm))
BM.lower<-pd-qt(1-alpha/2,n-1)*sqrt(v/n)
BM.upper<-pd+qt(1-alpha/2,n-1)*sqrt(v/n)
PERM.lower<-pd-pq1*sqrt(v/n)
PERM.upper<-pd+pq1*sqrt(v/n)
},
#--------------------Alternative= LOWER----------------------------------------#
less={
text.H0<-paste("p>=1/2")
text.Output <- paste("True relative effect is less than 1/2")
BM<-pt(T,n-1)
PERM<-p1perm
BM.lower<-0
BM.upper<-pd+qt(1-alpha,n-1)*sqrt(v/n)
PERM.lower<-0
PERM.upper<-pd+pq2*sqrt(v/n)
},
#--------------------Alternative= GREATER--------------------------------------#
greater={
text.H0<-paste("p<=1/2")
text.Output <- paste("True relative effect is greater than 1/2")
BM<-1-pt(T,n-1)
PERM<- 1-p1perm
BM.lower<-pd-qt(1-alpha,n-1)*sqrt(v/n)
BM.upper<-1
PERM.lower<-pd-pq2*sqrt(v/n)
PERM.upper<-1
}
)
#------------------------------------------------------------------------------#

if (plot.simci == TRUE) {
par(mfrow=c(2,1),oma=c(0,0,0,0))
plotz<-1
text.Ci<-paste((1-alpha)*100, "%", "Confidence Interval for p")
 Lowerp<-"|"
       plot(rep(pd,plotz),1:plotz,xlim=c(0,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(BM.lower,1:plotz, pch=Lowerp,font=2,cex=2)
              points(BM.upper,1:plotz, pch=Lowerp,font=2,cex=2)
              abline(v=0.5, lty=3,lwd=2)
              for (ss in 1:plotz){
              polygon(x=c(BM.lower[ss],BM.upper[ss]),y=c(ss,ss),lwd=2)}
              axis(1, at = seq(0, 1, 0.1))
              axis(2,at=1:plotz,labels="BM",font=2)
                box()
 title(main=c(text.Ci, " Method: Brunner-Munzel (BM), Permutation (PERM)" ))
       plot(rep(pd,plotz),1:plotz,xlim=c(0,1), pch=15,axes=FALSE,xlab="",ylab="")
       points(PERM.lower,1:plotz, pch=Lowerp,font=2,cex=2)
              points(BM.upper,1:plotz, pch=Lowerp,font=2,cex=2)
              abline(v=0.5, lty=3,lwd=2)
              for (ss in 1:plotz){
              polygon(x=c(PERM.lower[ss],PERM.upper[ss]),y=c(ss,ss),lwd=2)}
              axis(1, at = seq(0, 1, 0.1))
              axis(2,at=1:plotz,labels="PERM",font=2)
                box()
}



 if (info == TRUE) {
 cat("\n", "#----------------Nonparametric Paired t Test-------------------------------------------#", "\n","\n",
     "-", "Sample Size: ", n,"\n",
     "-", "Factor Levels: ", fl,"\n",
     "-", "H0: ", text.H0,"\n",
     "-", "Alternative Hypothesis: ", text.Output,"\n",
     "-", "Confidence Level:", conf.level*100,"%", "\n", 
     "-", "Method:", "Brunner-Munzel (BM), Permutation (PERM)", 
     
     "\n","\n",
     "#--------------------------------------------------------------------------------------#","\n",

            "\n")
    }


result.matrix<-round(matrix(c(BM.lower,pd,BM.upper,T,BM,PERM.lower,pd,PERM.upper,T,PERM),nrow=2,ncol=5,byrow=TRUE),rounds)
colnames(result.matrix)<-c("Lower","p.hat","Upper","T","p.value")
rownames(result.matrix)<-c("BM","PERM")

    Method<-"Brunner-Munzel Test (BM), Studentized Permutation Test (PERM)"
    methodvec<-c("BM","PERM")
    data.info <- data.frame(Sample=fl, Size=c(n,n))
    result<-list(Info=data.info, Analysis=result.matrix) 

result$input<-input.list
result$text.Output<-text.Output
result$methodvec<-methodvec
result$Method<-Method
class(result)<-"nparttestpaired"
return(result)
}

