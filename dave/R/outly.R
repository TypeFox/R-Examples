outly <-
function(veg,thresh=0.2,y=0.5){
  dmm<-ncol(veg)
  dc<- as.matrix(as.dist((1-cor(t(veg^0.5)))/2))
  diag(dc)<- 100
  neighb<- apply(dc,1,min)
  neighindex<- apply(dc,1,which.min)
  relevenames<- rownames(veg)
  threshv<- rep(thresh,nrow(veg))
  neighname<- relevenames[neighindex]
  diag(dc)<- 0
  o.pco<- pco(dc)
# reduce data 
  veg.red<- veg[neighb <= thresh,]
# remove empty species
  rf<- apply(sign(veg.red),2,sum)
  veg.red<- veg.red[,rf > 0]
  o.outlier<- list(threshold=thresh,y=y,rel.names=relevenames,neigh.names=neighname,neigh.dist=neighb,olddim=dim(veg),newdim=dim(veg.red),new.data=veg.red,pco.points=o.pco$points)
  }
