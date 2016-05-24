Tukey1df<-function(data) {
y<-data[,1]
Afactor<-data[,2]
Bfactor<-data[,3]
tst1<-is.factor(Afactor)
tst2<-is.factor(Bfactor)
tst3<-is.numeric(y)
if (tst1&tst2&tst3) {
a<-nlevels(Afactor)
b<-nlevels(Bfactor)
  }  else {stop("The first column of the data frame is the numeric response, the 2nd and 3rd columns should be coded as factors") }
tst4<-max(a,b)>2
tst5<-length(y)==a*b
if (tst4&tst5) {
ybb<-with(data, tapply(y, Bfactor, mean))
yba<-with(data, tapply(y, Afactor, mean))
sbb<-with(data, tapply(y, Bfactor, sum))
sba<-with(data, tapply(y, Afactor, sum))
ybardd<-mean(y)
CT<-(sum(y)^2)/(a*b)
ssA<-sum(sba^2/b)-CT
ssB<-sum(sbb^2/a)-CT
ssE<-sum(y^2)-CT-ssA-ssB
ybdj<-rep(ybb,a)
prody<-y*ybdj
sumprod<-tapply(prody,Afactor,sum)
leftsum<-sum(sumprod*yba)
ssAB<-(a*b*(leftsum-(ssA+ssB+a*b*ybardd^2)*ybardd)^2/(ssA*ssB))
ssR<-ssE-ssAB
F<-ssAB/(ssR/((a-1)*(b-1)-1))
Pval<-1-pf(F,1,((a-1)*(b-1)-1))
cat("Source           df     SS        MS        F     Pr>F","\n")
cat("A            ",paste(format(a-1, width=6)," ", format(round(ssA,4),justify="right"),"  ",format(round(ssA/(a-1),4), justify="right"),"\n"),sep="")
cat("B            ",paste(format(b-1, width=6)," ", format(round(ssB,4),justify="right"),"  ",format(round(ssB/(b-1),4), justify="right"),"\n"),sep="")
cat("Error        ",paste(format((b-1)*(a-1), width=6)," ", format(round(ssE,4),justify="right"),"  ",format(round(ssE/(a-1)*(b-1),4), justify="right"),"\n"),sep="")
cat("NonAdditivity",paste(format(1, width=6)," ", format(round(ssAB,4),justify="right"),"  ",format(round(ssAB,4),justify="right"),"  ",format(round(F,2),justify="right"),"  ",format(round(Pval,4),justify="right"),"\n"),sep="")
cat("Residual     ",paste(format((b-1)*(a-1)-1, width=6)," ", format(round(ssR,4),justify="right"),"  ",format(round(ssR/((a-1)*(b-1)-1),4), justify="right"),"\n"),sep="")
   } else {stop("This function only works for unreplicated 2-factor factorials with >2 levels for one of the factors")}
}
