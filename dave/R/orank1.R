orank1 <-
function(veg,use=c("columns","rows"),rlimit=5,y=1,x.axis=NULL,y.axis=NULL) {
# revised 25. 6. 2014 (ordination now with euclidean distance)
# default column names if missing
  defn<- is.null(names(veg))
  if(defn == TRUE) colnames(veg)<- seq(1,length(veg[1,]),1)
  defr<- is.null(rownames(veg))
  if(defn == TRUE) row.names(veg)<- seq(1,length(veg[,1]),1)

#
# erase empty species vectors
#
  f.s<- apply(abs(veg),2,sum)
  vegr<- veg[,f.s > 0]
  if(use == "columns") ranktemp <- vegr
  if(use == "rows") ranktemp <- t(vegr)
  vnames<- dimnames(ranktemp)[[2]]        # the variable names
  vnames<- substr(vnames,1,20)            # cutting to 20 chrs
  r <- cor(ranktemp^y)
  n <- length(r[1,])
  nm <- n-1           # number of ranks computed
  nm<- min(nm,rlimit)
  small<- 0.000001
  hrank <- rep(0,n)   # final ranking
  hvar <- rep(0,n)    # final variance of the same
  cumvar <- rep(0,n)  # cumulative variance of the same
  proc <- rep("TRUE",n)
#
  inertia <- sum(diag(r))
# in case there are no coordinates given, use ordination axes
  xnull<- is.null(x.axis)
  if(xnull == TRUE) {
     if(use == "rows") ranktemp <- vegr
     if(use == "columns") ranktemp <- t(vegr)
     db<- vegdist(ranktemp,method="euclidean")
     o.pco<- pco(db,k=2)
     x.axis<- o.pco$points[,1]
     y.axis<- o.pco$points[,2]
  }
# loop for final ranks
 r1 <- r
 for (k in 1:nm) {
     r <- r1
     ss <- rep(0.0,n)
     dr<-diag(r)
     if(sum(dr) > 0) {  
        sst<-ss[proc==T]
        sr2<-colSums(r[proc==T,proc==T]^2)
        drt<-dr[proc==T]
        drt[drt<small]<-small
	    ss[proc==T]<-(sst+sr2)/drt
        hr<- which.max(ss)
        hrank[k] <- hr
        hvar[k] <- ss[hr]
        cumvar[k:nm]<-  cumvar[k:nm]+hvar[k]
        proc[hr] <- FALSE           
        for (j in 1:n) {           # reducing r in r1
          trh<-r[hr,]*r[hr,j]/max(r[hr,hr],small)
          r1[j,]<-r[j,]-trh           
        }
     } else {                      # do nothing
     }
  }
  orth <- length(hrank[hrank])     # no. of orthogonal components
  as.double(hvar)
  relvar<- hvar/inertia*100
  cumrelvar<- cumvar/inertia*100
  orank<- list(use=use,n.ranks=orth,var.names=vnames[hrank[1:orth]],var.explained=hvar[1:orth],var.percent=hvar[1:orth]/inertia*100,cum.var=cumrelvar[1:orth],x.axis=x.axis,y.axis=y.axis,all.rownam=rownames(veg),all.colnam=colnames(veg))
}
