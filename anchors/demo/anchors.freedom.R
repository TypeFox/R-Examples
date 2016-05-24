#######################################################################
## 
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2007-02-01
##
## Modified
## - 2008-05-09 : JW
##   anchors 3.0 syntax
##
#######################################################################
cat("Replication of Wand et al 2007 rank analysis of freedom data\n")

data(freedom)

## ordering of vignettes
z1 <- anchors.order( ~ vign1+vign2+vign3+vign4+vign5+vign6, data=freedom)
summary(z1,top=10,digits=2)

z2 <- anchors.order( ~ vign2+vign1+vign3+vign5+vign4+vign6, data=freedom)
summary(z2,top=10,digits=2)


## anchors
a6 <- anchors(self ~ vign2+vign1+vign3+vign5+vign4+vign6, freedom, method="C")
summary(a6)

## parsimonious, with few ties
fo3 <- self ~ vign1+vign3+vign6
a3 <- anchors(fo3 , freedom, method="C")
summary(a3)

a3a <- anchors(fo3, freedom, method="C", subset=country=="Eastasia")
a3b <- anchors(fo3, freedom, method="C", subset=country=="Oceania"  )
summary(a3a)
summary(a3b)

## note can also use the anchors() expanded/list syntax for formulae
## when using method = ("C","order","minentropy")
fo2 <- list(self = self ~ 1,
            vign = cbind(vign2,vign4) ~ 1 )
a2a <- anchors(fo2, freedom, method="C", subset=country=="Eastasia")
a2b <- anchors(fo2, freedom, method="C", subset=country=="Oceania"  )
summary(a2a)
summary(a2b)

tb <- xtabs( ~ self, freedom, subset=country=="Eastasia"); tb/sum(tb)
tb <- xtabs( ~ self, freedom, subset=country=="Oceania");   tb/sum(tb)
