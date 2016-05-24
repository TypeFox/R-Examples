fullnormal<-function(effects,labs,alpha=.05,refline="TRUE") {
crit<-LenthPlot(effects,alpha=alpha,plt=FALSE)["ME"]
names<-names(effects)
names<-gsub(':','',names)
names<-gsub('1','',names)
le<-length(effects)
 for (i in 1:le) {
     logc<-(abs(effects[i])<=crit)
     if (logc) {names[i]<-" "}
                  }
qqnorm(effects, ylab="Estimated Effects", xlab="Normal Scores")
x<-qqnorm(effects,plot=FALSE)
zscr<-(x$x)
# Splits effects into positive and negative for labeling
effp<-effects[zscr>0]
zp<-zscr[zscr>0]
namep<-names[zscr>0]
effn<-effects[zscr<0]
zn<-zscr[zscr<0]
namen<-names[zscr<0]
text(zp,effp,namep,pos=1) 
text(zn,effn,namen,pos=3)  
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
