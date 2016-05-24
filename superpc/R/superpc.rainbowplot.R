superpc.rainbowplot=function(data, pred, sample.labels,  competing.predictors, call.win.metafile=FALSE){



extrapolate.surv<- function(y,ic,ncat=100){
cut<-seq(min(y),max(y),length=ncat)
b<-surv.to.class2(y,ic,cutoffs=cut)
yhat<-b$prob[,-1]%*%cut
return(yhat)
}





yhat= extrapolate.surv(data$y,data$censoring.status)
yhat[data$censoring.status==1]=data$y[data$censoring.status==1]


ncomp.predictors=length(competing.predictors)

 if (call.win.metafile) {
        win.metafile(width=10,height=2+ncomp.predictors)
    }


npanels=2*(ncomp.predictors+4)
layout(matrix(1:npanels,ncomp.predictors+4,2, byrow = TRUE),widths=c(.8,.2),heights=c(.2,rep(.1, (ncomp.predictors+2)), .05))
#layout.show(14)
par(mar=c(2,0,1,0))
par(cex=.8)
o=order(pred)
n=length(pred)
plot(0,0, xlim=c(0,1),ylim=c(0,1),type="n",axes=F)
for(i in 1:n){
  text(i/n,.5,labels=sample.labels[i],srt=90,cex=.5)
}

cols=rep(c("green","blue"),ncomp.predictors)
plot(0,0, xlim=c(0,1),ylim=c(0,1),type="n",axes=F)
my.barplot(yhat[o], label="survival", type="continuous", col="gray")
my.barplot(pred[o],label="supervised PC", type="continuous", col="red")
for(ii in 1:length(competing.predictors)){
if(!is.factor(competing.predictors[[ii]])){type="continuous"}
if(is.factor(competing.predictors[[ii]])){type="discrete"}

my.barplot(competing.predictors[[ii]][o], label=names(competing.predictors)[[ii]], type=type, col=cols[ii])
}

if (call.win.metafile) {
        dev.off()
    }

return()
}

my.barplot=function(x,label, type=c("continuous", "discrete"), col=c("red","green","blue", "gray")){
n=length(x)
 reds <- rgb(red=(0:n)/n, green=0,blue=0, names=paste("red",0:n,sep="."))
 greens <- rgb(green=(0:n)/n, red=0,blue=0, names=paste("green",0:n,sep="."))
 blues <- rgb(blue=(0:n)/n, green=0,red=0, names=paste("blue",0:n,sep="."))

if(type=="continuous"){
if(col=="gray"){
 palette(gray(seq(0,.9,length=n)))
 cols=rank(x)
 }

if(col=="red"){cols=reds[rank(x)]}
if(col=="green"){cols=greens[rank(x)]}
if(col=="blue"){cols=blues[rank(x)]}
 nc=4 
 temp=quantile(x,c(0.25, .5, .75, 1))
 values.legend=paste("<",round(temp,2),sep="")
 cols.legend=sort(cols[trunc(n*c(0.25,.5,.75,1))])

 if(length(unique(x))<6){
# we guess that x is an ordered discrete variable
   nc= length(unique(x))
   values.legend=sort(unique(x))
  }
  
}

if(type=="discrete"){
 palette("default") 
dd=sort(names(table(x)))
nc=length(dd)
cols=match(x,dd)
values.legend=dd
cols.legend=1:nc
}

par(mar=c(0,0,0,0))
 
plot(0,0,xlim=c(0,1),ylim=c(0,1/n),type="n",axes=F,xlab="",ylab="")
xval=c(0,1/n,1/n,0)
yval=c(0,0,1/n,1/n)
for(i in 1:n){
  polygon(xval,yval,col=cols[i])
  xval=xval+1/n
}
plot(0,0,xlim=c(0,1.25),ylim=c(0,1),type="n",axes=F,xlab="",ylab="")


h=.25

text(.5,.7,label=label, cex=.8)
xval=c(0,1/nc,1/nc,0)
yval=c(0,0,h,h)+h
for(i in 1:nc){
   polygon(xval,yval,col=cols.legend[i])
   text((xval[1]+xval[2])/2,yval[1]-.2, labels=values.legend[i], cex=.5)
  xval=xval+1/nc
}

return()
}




surv.to.class2 <- function (y, icens, cutoffs=NULL, n.class=NULL,  class.names=NULL, newy=y, newic=icens) 

# this is the function "pamr.surv.to.class2" from the pamr libarary
# the auxiliary functions below are also from pamr

# Splits patients into classes based on their survival times
# The user can either specify the number of classes or the survival
# time cutoffs.
#
# y - vector of survival times
# icens - censoring indicator
# cutoffs - survival time cutoffs
# n.class - number of classes to create
# class.names - optional vector of names for each class
{
        require(survival)
        if (is.null(cutoffs) & is.null(n.class)) {
                stop("Must specify either cutoffs or n.class")
        }
        if (!is.null(cutoffs) & !is.null(n.class)) {
                stop("Can't have both cutoffs and n.class specified")
        }
        data.sfit <- survfit(Surv(y,icens)~1)
        if (!is.null(cutoffs)) {
                if (is.null(class.names)) {
                        class.names <- 1:(length(cutoffs)+1)
                }
                cur.mat <- gen.y.mat2(list(y=y, icens=icens), cutoffs, class.names,                                              newdata=list(y=newy, icens=newic))
        }
        else {
                if (n.class==1) {
                        stop("Must have at least two classes")
                }
                if (is.null(class.names)) {
                        class.names <- 1:n.class
                }
                cur.quantiles <- seq(from=0, to=1, length=n.class+1)
                cur.quantiles <- cur.quantiles[2:n.class]
                cutoffs <- quantile(y[icens==1], cur.quantiles)
                cur.mat <- gen.y.mat2(list(y=y, icens=icens), cutoffs, class.names,
                                newdata=list(y=newy, icens=newic))
        }
        mle.classes <- apply(cur.mat, 1, get.mle.class)
         return(list(class=as.numeric(mle.classes), prob=cur.mat, cutoffs=cutoffs))
}
gen.y.mat2 <- function(surv.data, cutoffs, class.names=NULL, newdata=surv.data)
# Calculates the probability that a given patient belongs to a given
# class.  Returns a matrix where entry (i,j) is the probability that
# patient i belongs to class j.  The function calculates the
# probability that a given patient dies between two given cutoffs,
# and uses this information to calculate the probability that
# a patient with a censored survival time died in a given interval.
{
         data.sfit <- survfit(Surv(surv.data$y,surv.data$icens)~1)
         surv.ndx <- find.surv.ndx(cutoffs, data.sfit$time)
         surv.probs <- c(0, 1-data.sfit$surv[surv.ndx],1)
         surv.probs <- c(rep(0, sum((surv.ndx==0))), surv.probs)
         cutoffs <- c((min(surv.data$y)-1), cutoffs, (max(surv.data$y)+1))
         y.size <- length(cutoffs)
         y.mat <- matrix(0,nrow=length(newdata$y), ncol=(y.size-1))
         for (i in 2:y.size) {
                 cur.int.prob <- surv.probs[i] - surv.probs[i-1]
                 y.mat[((newdata$y<=cutoffs[i])&(newdata$y>cutoffs[i-1])&
                         (newdata$icens==1)),i-1] <- 1
                 which.x <- ((newdata$icens==0)&(newdata$y<=cutoffs[i-1]))
                 if (sum(which.x)>0) {
                         which.x.vals <- newdata$y[which.x]
                         surv.ndx <- find.surv.ndx(which.x.vals,
                                 data.sfit$time)
                         y.mat[which.x,i-1][surv.ndx==0] <- cur.int.prob
                         y.mat[which.x,i-1][surv.ndx!=0] <- cur.int.prob /
                                 data.sfit$surv[surv.ndx]
                 }
                 which.x <- ((newdata$icens==0)&(newdata$y>cutoffs[i-1])&
                         (newdata$y<=cutoffs[i]))
                 if (sum(which.x>0)) {
                         which.x.vals <- newdata$y[which.x]
                         surv.ndx <- find.surv.ndx(which.x.vals,
                                 data.sfit$time)
                         y.mat[which.x,i-1][surv.ndx==0] <- surv.probs[i]
                         y.mat[which.x,i-1][surv.ndx!=0] <- 1 -
                                 (1 - surv.probs[i]) / data.sfit$surv[surv.ndx]
                 }
         }
         if (!is.null(class.names)) {
                 y.mat <- as.data.frame(y.mat)
                 names(y.mat) <- class.names
                 y.mat <- as.matrix(y.mat)
         }
         y.mat
}

get.surv.q <- function(surv.obj, quantile) 
{
    ndx <- sum(surv.obj$surv > quantile)
    if (ndx==0)
        return(0)
    else
        return(surv.obj$time[ndx])
}
find.surv.ndx <- function(newtimes, oldtimes) 
{
	first <- apply(as.matrix(newtimes), 1, function(e1,e2) (e1>=e2), e2=oldtimes)
	as.vector(apply(first, 2, sum))
}
get.mle.class <- function(y.row) 
{
	i <- 1+sum((max(y.row)>cummax(y.row)))
	if (!is.null(names(y.row)[i])) {
		return(names(y.row)[i])
	}
	else return(i)
}

