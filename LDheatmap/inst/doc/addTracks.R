### R code from vignette source 'addTracks.Rnw'

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("LDheatmap")


###################################################
### code chunk number 2: load
###################################################
library(LDheatmap)
data(GIMAP5.CEU)
load(system.file("extdata/addTracks.RData",package="LDheatmap"))


###################################################
### code chunk number 3: fig1com
###################################################
ll<-LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)


###################################################
### code chunk number 4: fig1
###################################################
ll<-LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)


###################################################
### code chunk number 5: addTracks.Rnw:93-94 (eval = FALSE)
###################################################
## llGenes <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")


###################################################
### code chunk number 6: addGenes (eval = FALSE)
###################################################
## grid.newpage()
## grid.draw(llGenes$LDheatmapGrob)


###################################################
### code chunk number 7: fig2
###################################################
grid.newpage()
grid.draw(llGenes$LDheatmapGrob)


###################################################
### code chunk number 8: addRecomb (eval = FALSE)
###################################################
## llGenesRecomb <- LDheatmap.addRecombRate(llGenes, chr="chr7", genome="hg18")


###################################################
### code chunk number 9: addRecomb (eval = FALSE)
###################################################
## grid.newpage()
## grid.draw(llGenesRecomb$LDheatmapGrob)


###################################################
### code chunk number 10: fig3
###################################################
grid.newpage()
grid.draw(llGenesRecomb$LDheatmapGrob)


###################################################
### code chunk number 11: addTracks.Rnw:144-148
###################################################
set.seed(1)
atests<-runif(nrow(GIMAP5.CEU$snp.support))
names(atests)<-rownames(GIMAP5.CEU$snp.support)
atests["rs6598"]<-1e-5


###################################################
### code chunk number 12: addScatter
###################################################
llGenesRecombScatter<-LDheatmap.addScatterplot(llGenesRecomb,-log10(atests),
ylab="-log10(p-values)")


###################################################
### code chunk number 13: fig4
###################################################
grid.newpage()
grid.draw(llGenesRecombScatter$LDheatmapGrob)


###################################################
### code chunk number 14: addQplot
###################################################
require("ggplot2")
posn<-GIMAP5.CEU$snp.support$Position
manhattan2<-ggplotGrob(
{
qplot(posn,-log10(atests),xlab="", xlim=range(posn),asp=1/10)
last_plot() + theme(axis.text.x=element_blank(),
              axis.title.y = element_text(size = rel(0.75)))
}
)
llQplot<-LDheatmap.addGrob(ll,manhattan2,height=.7)


###################################################
### code chunk number 15: fig5
###################################################
grid.newpage()
grid.draw(llQplot$LDheatmapGrob)


###################################################
### code chunk number 16: addQplot2
###################################################
manhattan2<-ggplotGrob(
{
qplot(posn,-log10(atests),xlab="", xlim=range(posn),asp=1/10)
last_plot() + theme(axis.text.x=element_blank(),
              axis.title.y = element_text(size = rel(0.75)))
}
)
llQplot2<-LDheatmap.addGrob(ll,rectGrob(gp=gpar(col="white")),height=.2)
pushViewport(viewport(x=.48,y=.76,width=.99,height=.2))
grid.draw(manhattan2)
popViewport(1)
dev.off()


###################################################
### code chunk number 17: addGrob
###################################################
llImage<-LDheatmap.addGrob(ll,rasterGrob(GIMAP5ideo))


###################################################
### code chunk number 18: fig6
###################################################
grid.newpage()
grid.draw(llImage$LDheatmapGrob)


###################################################
### code chunk number 19: pl
###################################################
names(ll$LDheatmapGrob$children)


###################################################
### code chunk number 20: addTracks.Rnw:295-296
###################################################
names(llGenesRecombScatter$LDheatmapGrob$children)


###################################################
### code chunk number 21: addTracks.Rnw:306-308
###################################################
names(llQplot$LDheatmapGrob$children)
names(llImage$LDheatmapGrob$children)


