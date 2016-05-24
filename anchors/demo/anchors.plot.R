#######################################################################
## 
## Author:    Jonathan Wand <wand(at)stanford.edu>
## Created:   2008-05-09 
##
#######################################################################
cat("Demo of anchors library plot functions\n")
## 
data(mexchn)
fo <- list(self = xsayself ~ 1,
           vign = cbind(xsay5,xsay3,xsay1)    ~ 1,
           tau  = ~ age + china,
           tau1 = ~ age + china + male,
           cpolr= ~ 1)

## 
o1 <- anchors(fo, data = mexchn, method="order", options=anchors.options(ties="interval"))
summary(o1)
barplot(o1)
o2 <- anchors(fo, data = mexchn, method="order", options=anchors.options(ties="nominal"))
summary(o2)
barplot(o2)

##
r1m <- anchors(fo, data = mexchn, method="B", subset=china==0)
r1c <- anchors(fo, data = mexchn, method="B", subset=china==1)
summary(r1m)
barplot(r1m,r1c)
barplot(r1m,r1c , ties="cpolr")
barplot(r1m,r1c , ties="minentropy")
barplot(r1m,r1c , ties="omit")

r2 <- anchors(fo, data = mexchn, method="C")
summary(r2)
barplot(r2 , ties="uniform")
barplot(r2 , ties="cpolr")

##
r3 <- anchors(fo, data = mexchn, method="C", combn=TRUE)
plot( r3, xy=c("estimated","minimum"))
plot.anchors.combn( r3$combn, xy=c("minimum","interval"))

