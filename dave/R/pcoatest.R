pcoatest <-
function(veg,y=1) {
  sveg<- veg^y
#
# coordinates of all 6 methods
#
  vdm<- vegdist(sveg,method="euclidean")
  oute<- pco(vdm)

  vdm<- vegdist(sveg,method="manhattan")
  outm<- pco(vdm)

# chord distance
  vdm<- dist(decostand(sveg,"norm"),method="euclidean")
  outco<- pco(vdm)

  vdm<- vegdist(sveg,method="canberra")
  outc<- pco(vdm)

  vdm<- vegdist(sveg,method="bray")
  outb<- pco(vdm)

  vdm<- as.dist((1 - cor(t(sveg)))/2)
  outcor<- pco(vdm)

# PCA
  outpca<- prcomp(sveg)
  pcscor<- scores(outpca,display="sites",choices=1:2)
#
# procrustes analysis for all 6 methods
#
# euclid
proce<- procrustes(pcscor,scores(oute))

# manhattan
procm<- procrustes(pcscor,scores(outm))

# chord
procco<- procrustes(pcscor,scores(outco))

# canberra
procca<- procrustes(pcscor,scores(outc))

# bray
procb<- procrustes(pcscor,scores(outb))

# 1-correlation
procc<- procrustes(pcscor,scores(outcor))
#
# output list
#
 pcovar<- list(nrel=nrow(sveg),nspec=ncol(sveg),y=y,euclidpca=proce$Yrot,euclidpco=proce$X,manhpca=procm$Yrot,manhpco=procm$X,cordpca=procco$Yrot,cordpco=procco$X,canpca=procca$Yrot,canpco=procca$X,bpca=procb$Yrot,bpco=procb$X,corpca=procc$Yrot,corpco=procc$X)
}
