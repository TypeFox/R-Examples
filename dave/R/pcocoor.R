pcocoor <-
function(veg,method,y=1) {
#
# Computing and plotting PCOA using vegan() and labdsv()
# followed by calculating and plotting correlations with
# species. sveg is a data frame.
#
# sel.spec<- sel.sp
  vdb<- vegdist(veg^y,method)
  outc<- pco(vdb,k=6)
  ssveg<- scale(veg,center=TRUE, scale=TRUE)
  sspoints<- scale(outc$points,center=TRUE, scale=TRUE)
  Eig.pco<- (t(sspoints) %*% ssveg)/(ncol(veg)-1)
# list
  o.pcobiplot<- list(nrel=nrow(veg),nspec=ncol(veg),naxes=ncol(outc$points),rpoints=outc$points,spoints=Eig.pco,allspnames=names(veg))
  }
