ROCscore<-function(xobs, yobs, type="output")
{
# verification
 n1<-nrow(xobs)    # number of observations
 p.input<-ncol(xobs)     # number of inputs

# verification of some conditions
if(n1!=nrow(yobs)) stop("xobs and yobs have not the same number of observations")
match.arg(type,c("output","input","hyper"))

# initialisation
alpha=seq(0.01,1,by=0.01)
m=c(1,5,8,11,14,17,20,23,26,29,32,35,40,45,50,55,60,65,70,75,80,85,90,95,seq(100,2500,length.out=76))

 # Representation of the alpha quantile efficiency frontier 
  # 1- Computation of the alpha-quantile score
 perf1=NULL
 perf2=NULL 
  for(k in 1:length(alpha))
  {
  res1<-alphascore(xobs,yobs, alpha=alpha[k])[,type]  # score computed on (xtab,ytab)
  res2<-ordermscore(xobs,yobs, m=m[k])[,type]   # score computed on (xtab,ytab)
  perf1=c(perf1,length(which(res1>1))/n1)
  perf2=c(perf2,length(which(res2>1))/n1)
  }
  plot(alpha,perf1,type='l',xlab="alpha",ylab="Percentage of super-efficiency firms",ylim=c(0,1),col='royalblue')
  par(new=TRUE)
  plot(m, perf2,type='l',lty=2, ann = FALSE, yaxt = "n", xaxt="n",yaxt="n",ylim=c(0,1),col='red')
  axis(3)
  mtext(side=3,"m")
  legend("topright",legend=c("f(alpha)","f(m)"),lty=1:2,col=c("royalblue","red"))

 res<-data.frame(alpha,perf1,m,perf2)
 names(res)=c("alpha","f(alpha)","m","f(m)")
 return(res)
}  