fspa2 <-
function(veg,method,d.rev=0.5,n.groups=3) {
  toolarge<- 1e12
# settig defaults
  def.rev<- is.null(d.rev)
  if(def.rev == TRUE) d.rev<- 0.5
  rev<- d.rev
# default for distance measure
  def.meth<- is.null(method)
  if(def.meth == TRUE) method<- "euclidean"
# default for group number
  def.gr<- is.null(n.groups)
  if(def.gr == TRUE) n.groups<- 3
#
# end default
#
  order<- length(veg[,1])
  mde <- vegdist(veg,method,diag=TRUE,upper=TRUE)
  hclust.r <-hclust(mde,method="ward.D")
  memb.r<- cutree(hclust.r,k = n.groups)
  den<- as.dendrogram(hclust.r)
#  plot(den) if required here
  md<- as.matrix(mde)
# ordination with pco
  out.pco<- pco(mde,k=6)
  oldscores<- out.pco$points
  lambda1<- round(out.pco$eig[1]/sum(out.pco$eig),digits=2)
  lambda2<- round(out.pco$eig[2]/sum(out.pco$eig),digits=2)

# ordering increasing
  omd<- order(md)
# finding treshold tresh from revised portion rev
  nrev<- as.integer(length(as.vector(md))*(1-rev))
  tresh<- md[omd][nrev]
# connecting distant points:
# nline counts the lines
# iindex is index of starting point
# jindex is index of ending point
#
  nline<- 0
  iindex<- rep(0,order*order)
  jindex<- rep(0,order*order)
  for(i in 1:order) for(j in i:order) {
   if(md[i,j] < tresh) {
        nline<- nline+1
        iindex[nline]<- i
        jindex[nline]<- j
        }
  }
  sel<- which(md >= tresh)
# replacing too large distances by toolarge
  md[sel]<- toolarge
# replacing all distances by nearest path
  for(kk in 1:3){
  for(i in 1:order) for(j in 1:order) {
      md[i,j]<- min(md[i,]+md[j,])
  }
  }
#
  out.pco<- pco(md,k=3)
#  par(op)
  lambda1<- round(out.pco$eig[1]/sum(out.pco$eig),digits=2)
  lambda2<- round(out.pco$eig[2]/sum(out.pco$eig),digits=2)
# output list
  outfspa<- list(oldpoints=oldscores,symbols=memb.r,nline=nline,startline=iindex,endline=jindex,newpoints=out.pco$points,dmat.before=mde,dmat.after=md,d.rev=rev,eig=out.pco$eig)
  }
