speedprof2 <-
function(veg,timescale,orders,y=1,adjust=FALSE) {
  deford<- is.null(orders)
  if(deford == TRUE) orders<- c(1,2,3,4)
  maxord<- length(orders)
  defy<- is.null(y)
  if(defy == TRUE) y<- 1
# adjusting scores to 100% if requested
  if(adjust == TRUE) {
      adj<- function(x) {100*x/sum(x)}
      veg<- t(apply(veg,1,adj))
  }
#  par(mfrow=c(1,1),omi=c(2,0,0,0))
# distance matrix of releves, full mode
  mde <- vegdist(veg^y,method = "euclidean",diag=TRUE,upper=TRUE)  
  dimde <- dim(veg)[1]                          # matrix dimension
  mde<- as.matrix(mde,dim=c(dimde,dimde))
 o.speedprof<- list(nrel=dimde,dmatrix=mde,timescale=timescale,orders=orders)
 }
