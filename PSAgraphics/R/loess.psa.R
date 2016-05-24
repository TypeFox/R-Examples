loess.psa<-function(response, treatment = NULL, 
    propensity = NULL, family = "gaussian", span = .7, degree=1, minsize=5,
    xlim=c(0,1),  colors = c('dark blue','dark green','blue','dark green'), 
    legend.xy = "topleft", legend = NULL, int = 10, 
    lines = TRUE, strata.lines = TRUE, rg = TRUE, 
    xlab="Estimated Propensity Scores", ylab="Response", pch = c(16,1), ...){

#Creating a data frame to allow subsetting. 

if(is.vector(response)){
    rtp<-data.frame(response,treatment,propensity)}else{
    rtp<-as.data.frame(response)
    dimnames(rtp)[2]<-list(c("response","treatment","propensity"))}

if(is.null(treatment)){treatment<-response[,2]}
sut<-sort(unique(treatment))

response<-rtp[,1]    
treatment<-rtp[,2]
propensity<-rtp[,3]

#Getting the loess estimates of response
loess.0<-loess(rtp$r~rtp$p, subset = treatment == sut[1], family = family, span = span, degree=degree)
loess.1<-loess(rtp$r~rtp$p, subset = treatment == sut[2], family = family, span = span, degree=degree)

#Plotting points (propensity,response) by treatment levels 0 and 1
plot(rtp$p,rtp$r,type="n",xlim=range(rtp$p),xlab=xlab,ylab=ylab)
points(rtp$p[treatment==sut[1]],rtp$r[treatment==sut[1]],col=colors[1],pch = pch[1], ... )
points(rtp$p[treatment==sut[2]],rtp$r[treatment==sut[2]],col=colors[2],pch = pch[2], ...)
li<-length(int)
if(strata.lines){int2 <- int 
	if(li == 1){int2 <- quantile(propensity, seq(0, 1, 1/int))
		        li <- length(int2)
		       }
	for(i in 1:li){abline(v = int2[i], lwd = .5, lty =  3, col = "dark grey")}
		       }


n0<-length(rtp$p[treatment==sut[1]])
n1<-length(rtp$p[treatment==sut[2]])
rtp.p.0<-rtp$p[treatment==sut[1]]
rtp.p.1<-rtp$p[treatment==sut[2]]

if(lines){ord.0<-order(rtp.p.0)
          ord.1<-order(rtp.p.1)
          lines(rtp.p.0[ord.0],loess.0$f[ord.0],cex=.6,col=colors[3], lwd = 1.5, lty = 1)
          lines(rtp.p.1[ord.1],loess.1$f[ord.1],col=colors[4], lwd = 1.5, lty = 2)}
     else{points(rtp$p[treatment==sut[1]],loess.0$f,pch=3,cex=.6,col=colors[3])
          points(rtp$p[treatment==sut[2]],loess.1$f,pch=4,cex=.6,col=colors[4])}

if(is.null(legend)){legend<-sut}       
legend(x = legend.xy, y = NULL, legend = legend, fill = colors[1:2], bty = "n")

#Rug plots of propensity scores
if(rg){
    rug(rtp$p[treatment==sut[1]], ticksize=.02, side=1, col=colors[1])
    rug(rtp$r[treatment==sut[1]], ticksize=.02, side=2, col=colors[1])
    rug(rtp$p[treatment==sut[2]], ticksize=.02, side=3, col=colors[2]) 
    rug(rtp$r[treatment==sut[2]], ticksize=.02, side=4, col=colors[2])
    box()                        
      }                              

##Getting Loess based estimate of effect size.

#Divide up (0,1) into int subintervals
ed <- length(int)
if(ed==1){prop.labels<-as.numeric(cut(propensity, quantile(propensity, seq(0, 1, 1/int)),
      include.lowest = TRUE, labels = FALSE)); nint<-int; dp<-0}
if(ed>1){
  if(int[1]==0 & int[ed]==1){subints<-int; nint<-ed-1; dp<-0}
  if(int[1] >0 & int[ed]==1){subints<-c(0,int); nint<-ed; dp<-1}
  if(int[1]==0 & int[ed] <1){subints<-c(int,1); nint<-ed; dp<-2}
  if(int[1] >0 & int[ed] <1){subints<-c(0,int,1); nint<-ed+1; dp<-3}
  prop.labels<-as.numeric(cut(propensity,breaks=subints,include.lowest = TRUE, labels = FALSE))
        }     

#Checking for empty T/C strata
table.plt<-table(prop.labels,treatment)
ttable.plt<-table.plt[,1]*table.plt[,2]
flag.0<-1
for(i in 1:nint){flag.0<-flag.0*ttable.plt[i]}
if(flag.0 == 0) print("Warning: Some strata-treatment levels have no cases.  Redefine 'int'.")

#Labels propensities as belonging to one of these subintervals
o<-order(treatment)
prop.labels<-prop.labels[o]
resp.o<-response[o]
ncontrol<-table(treatment)[1]
ntreat<-table(treatment)[2]

#Means for each of the subintervals; differences
means.0<-tapply(loess.0$fitted,prop.labels[1:ncontrol],mean)
var.0<-tapply(resp.o[1:ncontrol],prop.labels[1:ncontrol],var)
ncp1<-ncontrol+1
ncpnt<-ncontrol+ntreat
means.1<-tapply(loess.1$fitted,prop.labels[ncp1:ncpnt],mean)
var.1<-tapply(resp.o[ncp1:ncpnt],prop.labels[ncp1:ncpnt],var)
diff.means<-means.1-means.0

#Creating the weights
#Counts in each subinterval, treatment and control
counts.0<-table(prop.labels[1:ncontrol])
counts.1<-table(prop.labels[(ncontrol+1):(ncontrol+ntreat)])
#Indicator variable for nonempty subintervals by treatments
indicator.0<-NULL
indicator.1<-NULL
for(i in 1:nint)if(as.logical(counts.0[i]< minsize)){indicator.0<-c(indicator.0,0)}else{indicator.0<-c(indicator.0,1)}
for(i in 1:nint)if(as.logical(counts.1[i]< minsize)){indicator.1<-c(indicator.1,0)}else{indicator.1<-c(indicator.1,1)}
#Now add counts and multiply by each indicator to get the weights of each subinterval
wts<-(counts.0+counts.1)*indicator.0*indicator.1

#Working on the sd estimator
#First, get rid of potential NAs from var applied to strata with 1 case
mash<-cbind(var.0,var.1,counts.0,counts.1,indicator.0,indicator.1)
mash<-na.omit(mash)
sum.wtd.var<-(mash[,1]/mash[,3] + mash[,2]/mash[,4])*mash[,5]*mash[,6]
nused.int<-sum(indicator.0*indicator.1)

#And finally the calculator of the direct effect estimator

out.table<-cbind(counts.0,counts.1,means.0,means.1,diff.means)
colnames(out.table)<-c(paste("counts.",sut[1],sep=""), paste("counts.",sut[2],sep=""), 
paste("means.",sut[1],sep=""), paste("means.",sut[2],sep=""), "diff.means")

if(dp==0){dee <- sum((wts*diff.means),na.rm=TRUE)/sum(wts)
          sd.wt <- ((sum(sum.wtd.var))^.5)/nused.int}
if(dp==1){dee <- sum((wts[-1]*diff.means[-1]),na.rm=TRUE)/sum(wts[-1])
          sd.wt <- ((sum(sum.wtd.var[-1]))^.5)/nused.int
          out.table<-out.table[-1,] 
          rownames(out.table)<-1:(nint-1)}
if(dp==2){dee <- sum((wts[-nint]*diff.means[-nint]),na.rm=TRUE)/sum(wts[-nint])
          sd.wt <- ((sum(sum.wtd.var[-nint]))^.5)/nused.int
          out.table <- out.table[-nint,]
          rownames(out.table)<-1:(nint-1)}
if(dp==3){dee <- sum((wts[-c(1,nint)]*diff.means[-c(1,nint)]),na.rm=TRUE)/sum(wts[-c(1,nint)])
          sd.wt <- ((sum(sum.wtd.var[-c(1,nint)]))^.5)/nused.int
          out.table<-out.table[-c(1,nint),]
          rownames(out.table)<-1:(nint-2)}
               
#Flag: Testing for large changes in dee for various strata sizes; only done with default int = 15.
#Much of the code that follows is similar to above.
int.test<-c(5,10,15)
effect.test<-c(1,1,dee)
res<-cbind(int.test,effect.test)
it<-1
if(ed == 1) it<-int[1]
if(it == 15){
indx<-0
for(i in c(5,10)){
indx<-indx+1
prop.labels<-as.numeric(cut(propensity, quantile(propensity, seq(0, 1, 1/i)),
      include.lowest = TRUE, labels = FALSE))
      
nint<-i

#Flag: Labels propensities as belonging to one of these subintervals
o<-order(treatment)
prop.labels<-prop.labels[o]
resp.o<-response[o]
ncontrol<-table(treatment)[1]
ntreat<-table(treatment)[2]

#Flag: Means for each of the subintervals; differences
means.0<-tapply(loess.0$fitted,prop.labels[1:ncontrol],mean)
ncp1<-ncontrol+1
ncpnt<-ncontrol+ntreat
means.1<-tapply(loess.1$fitted,prop.labels[ncp1:ncpnt],mean)
difff.means<-means.1-means.0

#Flag: Creating the weights
#Counts in each subinterval, treatment and control
counts.0<-table(prop.labels[1:ncontrol])
counts.1<-table(prop.labels[(ncontrol+1):(ncontrol+ntreat)])
#Indicator variable for nonempty subintervals by treatments
indicator.0<-NULL
indicator.1<-NULL
for(i in 1:nint)if(as.logical(counts.0[i]< minsize)){indicator.0<-c(indicator.0,0)}else{indicator.0<-c(indicator.0,1)}
for(i in 1:nint)if(as.logical(counts.1[i]< minsize)){indicator.1<-c(indicator.1,0)}else{indicator.1<-c(indicator.1,1)}
#Now add counts and multiply by each indicator to get the weights of each subinterval
wtts<-(counts.0+counts.1)*indicator.0*indicator.1
res[indx,2] <- sum((wtts*difff.means),na.rm=TRUE)/sum(wtts)
}
#Flag!
flag<-FALSE
err<-.05*(sd(response))
if(abs(res[1,2]-res[2,2]) > err) flag<-TRUE
if(abs(res[1,2]-res[3,2]) > err) flag<-TRUE
if(abs(res[2,2]-res[3,2]) > err) flag<-TRUE
if(flag) print("Warning: Effect size estimate unstable with changes in number of strata")
if(flag) print(res)
}

#Rough 95 CI
CI95<-c(dee-2*sd.wt,dee+2*sd.wt)

out<-list(dee, sd.wt, CI95, out.table)
names(out)<-c("ATE", "se.wtd", "CI95", "summary.strata")
return(out)
                              }

