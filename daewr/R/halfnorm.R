halfnorm<-function(effects,labs,alpha=.05,refline="TRUE") {
crit<-LenthPlot(effects,alpha=alpha,plt=FALSE)["ME"]
effp<-abs(effects)
effp<-sort(effp)
names<-names(effp)
names<-gsub(':','',names)
names<-gsub('1','',names)
le<-length(effects)
r<- c(1:le)
zscore<-c(rep(0,le))
 for (i in 1:le) {
     zscore[i]<-qnorm( ( ( r[i]-.5)/le+1)/2)
     logc<-(abs(effp[i])<=crit)
     if (logc) {names[i]<-" "
               }
                  }
plot(zscore,effp,xlab="Half Normal scores", ylab="abs(effects)")
text(zscore,effp,names,pos=1)
cat("zscore=", zscore)
cat("effp=",effp) 
# calculate pse statistic
ahe<-abs(effects)
s0<-1.5*median(ahe)
selhe<-ahe<(2.5*s0)
pse=1.5*median(ahe[selhe])
if (refline) {
# add reference line to plot
abline(0,pse)
             }
                   }
