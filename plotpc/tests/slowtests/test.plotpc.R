# test.plotpc: test earth with a biggish model

library(plotpc)
library(grid)
if(!interactive())
    postscript(width=6, height=6)
options(warn=1) # print warnings as they occur

example(plotpc)
example(plotld)

data(iris)

# reproduce pairs plot from help page
pc <- princomp(iris[,-5]) # -5 to drop Species
pairs(pc$scores, col=c(2,3,4)[unclass(iris$Species)])

x <- iris[,c(4,3)] # select Petal.Width and Petal.Length
plotpc(x, main="test plotpc",
       xrange=14,
       breaks1=5, gp.hist1=gpar(col="gray", fill=0), height1=1,
       offset1=NULL,
       breaks2=6, gp.hist2=gpar(col=1, fill="gray"),
       offset2=5, height2=-.5, flip2=TRUE,
       breaksx=20, gp.histx=gpar(col="red", fill=0), heightx=2,
       heighty=0)

library(alr3) # get banknote data
data(banknote)
plotld(banknote[,-7], npc=10) # -7 to drop Y
plotld(banknote[,-7], abs.=TRUE)

# Flury and Riedwyl "Multivariate Statistics A Practical Approach" figure 10.2
x <- banknote[101:200,4:5] # select Forged, Top and Bottom
plotpc(x, main="Figure 10.2: The projection\nwith the largest variance\n",
       breaks=12, height1=-2, offset1=5, height2=0)

# Flury and Riedwyl figure 10.3
x <- banknote[101:200,4:5] # select Forged, Top and Bottom
plotpc(x, main="Figure 10.3\nProjections associated with principal\ncomponents in two dimensions",
       heightx=0, heighty=0, breaks=12,
       axis.len1=0, height1=-1,
       axis.len2=0, height2=-4, offset2=7.5)

# Flury and Riedwyl figure 10.4
x <- banknote[101:200,4:5] # select Forged, Top and Bottom
vp <- plotpc(x, main="Figure 10.4: Principal axes", heightx=0, heighty=0, breaks=12,
             axis.lenx=5, axis.leny=5,
             axis.len1=4, gp.axis1=gpar(col=1), height1=0,
             axis.len2=6, gp.axis2=gpar(col=1), height2=0)
pushViewport(vp)
grid.text("U1", x=unit(5.8, "native"), y=unit(13.6, "native"), gp=gpar(cex=.8, font=2))
grid.text("U2", x=unit(12, "native"), y=unit(14, "native"), gp=gpar(cex=.8, font=2))
grid.text("X4", x=unit(16.8, "native"), y=unit(11.2, "native"), gp=gpar(cex=.8, font=2))
grid.text("X5", x=unit(10.6, "native"), y=unit(14.8, "native"), gp=gpar(cex=.8, font=2))
popViewport()

# Flury and Riedwyl figure 10.13
x <- banknote[101:200,4:5] # select Forged, Top and Bottom
vp <- plotpc(x, main="Figure 10.13\nRegression lines and\nfirst principal component",
             heightx=0, heighty=0, height1=0, height2=0,
             gp.points=gpar(cex=.6, col="gray"), pch=1,
             axis.len1=4, gp.axis1=gpar(col=1),
             axis.len2=.2, gp.axis2=gpar(col=1),
             yonx=TRUE, xony=TRUE)
pushViewport(vp)
grid.text("PC1", x=unit(6.5, "native"), y=unit(13.4, "native"), gp=gpar(cex=.8, font=2))
popViewport()

# Flury and Riedwyl figure 7.1 simple version
x <- banknote[,4:5] # select Top and Bottom
plotpc(x, xrange=24, main="Figure 7.1: Various projections\n(simple version)\n",
       height=-2,  # reverse directions of histograms to match Flury and Riedwyl
       breaks=12, heightx=0, heighty=0, height1=0, height2=0, axis.len1=0, axis.len2=0,
       angle3=30, angle4=102, angle5=174, angle6=246, angle7=318)

if(!interactive()) {
    dev.off()         # finish postscript plot
    q(runLast=FALSE)  # needed else R prints the time on exit (R2.5 and higher) which messes up the diffs
}
