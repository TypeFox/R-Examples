`lorenz.curve` <-
function(data, weight=NULL, mul=FALSE, plot.it=TRUE, 
   main=NULL, xlab=NULL, ylab=NULL, xlim=c(0,1), ylim=c(0,1), ... )
{

### Function for the lorenz curve. Lorenz Single is the default and
### Lorenz Multiple is an option. The input data should be 
### a data frame with 2 columns. The first column will be treated
### as a data vector, and the second column to be treated as weight vector.
### Alternatively, data and weight can be entered as separate one-column vectors.


### Check length of data and weight vectors ###

if(any(is.na(data))){stop("NAs in data. Please try again.")}

if( is.vector(data) & !( is.list(data) & length(data)>1 ) | (is.data.frame(data)&length(data)==1) )
   {
   if(is.null(weight)){stop("A single-column weight vector is required. Please try again.")}

   if( !(is.vector(weight)) | ( is.list(weight) & length(weight)>1 ) ){
      stop("The weight input is not a single-column vector. Please try again.")}

   if(any(is.na(weight))){stop("NAs in the weight vector. Please try again.")}
         
   dframe=data.frame(data,weight)   
   names(dframe)[1]=deparse(substitute(data))
   names(dframe)[2]=deparse(substitute(weight))
   }

else{dframe=data.frame(data)[1:2]}

if( any(is.na(dframe[,1])) | is.factor(dframe[,1]) | is.character(dframe[,1]) ){
   stop("The first column contains invalid input. Please try again.")}

if( any(is.na(dframe[,2])) | is.factor(dframe[,2]) | is.character(dframe[,2]) ){
   stop("The second column contains invalid input. Please try again.")}
     

      
### Process the data vector based on weighting###


if(mul)
   {
   vv=NULL
  
   for (k in 1:nrow(dframe))
      {
      if(dframe[k,2]>1){ vv=c(vv,rep(dframe[k,2]/dframe[k,1], dframe[k,2])) }
      else{ vv=c(vv,dframe[k,2]/dframe[k,1]) }
      }
   if(is.null(main)){main="Lorenz Curve Multiple"}
   }

else{vv=dframe[,2]/dframe[,1]; if(is.null(main)){main="Lorenz Curve"}}


nn<-length(vv)
relative.mean.dev<-1/nn*sum(abs(vv-mean(vv)))/mean(vv)

d<-0
for (i in 1:nn)
{
  for (j in 1:nn) d<-d+abs(vv[i]-vv[j])
}

### gini index ###
gini<-d/nn/(nn-1)/2/mean(vv)
# RDR<-(max(vv)-min(vv))/mean(vv)  
case<-sort(vv)

tot<-sum(case)
qs<-case/tot
qa<-c(0)
qscomp<-c(qa,qs)
y<-cumsum(qscomp)
x<-seq(0,1,1/(length(qs)))
par(mfrow=c(1,1))


if(plot.it){

  if(is.null(xlab)){xlab=paste("Cumulative fraction of",names(dframe)[1]) }
  if(is.null(ylab)){ylab=paste("Cumulative fraction of",names(dframe)[2]) }

  plot(x,y, type="l",xaxs="i",yaxs="i",main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...) ##Lorena Curve
  abline(0,1, ...)
  xx<-x
  yy1<-x
  yy2<-y
  segments(xx,yy1,xx,yy2)    ################ add the lines in the curve area
  legend(0.05,0.8,"rel.mean.de=",bty="n", cex=0.6)
  legend(0.18,0.8,round(relative.mean.dev,3),bty="n", cex=0.6)
  legend(0.05,0.75,"gini.index=",bty="n", cex=0.6)
  legend(0.18,0.75,round(gini,3),bty="n", cex=0.6)
  # legend(0.05,0.7,"RDR=",bty="n", cex=0.6)
  # legend(0.18,0.7,round(RDR,3),bty="n", cex=0.6)
  legend(0.05,0.65,"L(1/2)=",bty="n", cex=0.6)
  legend(0.18,0.65,round(median(y),3),bty="n", cex=0.6)
  }

}

