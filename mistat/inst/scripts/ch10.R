###################################################
### Chap10Start
###################################################
library(mistat)
library(qcc)

qcc.options(bg.margin = "white",
            violating.runs = list(pch=15, col=1),
            beyond.limits = list(pch=15, col=1))


###################################################
### MultivariateCC
###################################################
data(ALMPIN)

Base <- ALMPIN[1:30,]                 

MeansBase <- colMeans(Base)           

CovBase <- cov(Base)                  

AlmpinSubset <- ALMPIN[-(1:30),]              

library(qcc)                          

Mqcc <- mqcc(data=AlmpinSubset,               
             type="T2.single",        
             center=MeansBase,        
             cov=CovBase,             
             plot=T)                  

summary(Mqcc)                         


###################################################
### PlotT2ChartAluminiumPins
###################################################

invisible(mqcc(data=AlmpinSubset, 
               type="T2.single", 
               center=MeansBase, 
               cov=CovBase, 
               plot=TRUE, main=""))

rm(Base, CovBase, MeansBase, AlmpinSubset, Mqcc)


###################################################
### PlotScatterMatrixPlace
###################################################
library(ggplot2)
library(grid)

# Thanks to Stephen Turner:
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

data(PLACE)

PLACE$code <- factor(c(rep("lDev", 9*16),
                       rep("mDev", 3*16),
                       rep("hDev", 14*16)))

P1 <- ggplot(PLACE, aes(xDev, yDev, colour=code))

P1 <- P1 + geom_point(alpha=8/10, size=2) + geom_density2d(alpha=8/10, size=0.2) + 
  scale_color_grey(start=0.4, end=0.8) + theme_bw() +theme(legend.position = "none")

P2 <- ggplot(PLACE, aes(xDev, tDev, colour=code))

P2 <- P2 + geom_point(alpha=8/10, size=2) + geom_density2d(alpha=8/10, size=0.2) + 
  scale_color_grey(start=0.4, end=0.8) + theme_bw() + theme(legend.position = "none")

P3 <- ggplot(PLACE, aes(yDev, tDev, colour=code))

P3 <- P3 + geom_point(alpha=8/10, size=2) + geom_density2d(alpha=8/10, size=0.2) + 
  scale_color_grey(start=0.4, end=0.8) + theme_bw() + theme(legend.position = "none")

arrange(P1, P3, P2, nrow=2, ncol=2)
rm(P1, P2, P3)


###################################################
### PlotConditionalHistogramPlace
###################################################
PLACE$code <- cut(PLACE$xDev, breaks=c(-1, 0.001, 1))

P1 <- ggplot(PLACE, aes(x = xDev, fill=code))

P1 <- P1 + geom_histogram(binwidth=0.0005, colour="black") + 
  scale_fill_grey(start=0.5, end=0.8) + 
  theme_bw() + 
  theme(legend.position = "none") + coord_flip()

P2 <- ggplot(PLACE, aes(x = yDev, fill=code))

P2 <- P2 + geom_histogram(binwidth=0.0005, colour="black") + 
  scale_fill_grey(start=0.5, end=0.8) + 
  theme_bw() + 
  theme(legend.position = "none") + coord_flip()

P3 <- ggplot(PLACE, aes(x = tDev, fill=code))

P3 <- P3 + geom_histogram(binwidth=0.02, colour="black") + 
  scale_fill_grey(start=0.5, end=0.8) + 
  theme_bw() + 
  theme(legend.position = "none") + coord_flip()


arrange(P1, P2, P3, nrow=1)

rm(P1, P2, P3, arrange, vp.layout)


###################################################
### PlotPlaceXyDev
###################################################
PLACE$code <- cut(PLACE$xDev, breaks=c(-1, 0.001, 1))

PLACE2 <- reshape(data=PLACE[, !names(PLACE) %in% "tDev"], 
                  varying=list(c("xDev", "yDev")), 
                  direction="long")

PLACE2$Dev <- factor(PLACE2$time)

levels(PLACE2$Dev) <- c("xDev", "yDev")

PLACE2$Deviation <- PLACE2$xDev

ggplot(data=PLACE2, 
       aes(x=crcBrd, y=Deviation, colour=code, shape=Dev)) + 
  geom_point(alpha=8/10, size=4) + 
  scale_color_grey(start=0.2, end=0.8) + 
  theme_bw() + 
  guides(shape = guide_legend("Error in placement along axis"), 
         colour=FALSE) +
  theme(legend.position = "top", 
        legend.box = "horizontal")

rm(PLACE2)


###################################################
### MultivariateCCPlace
###################################################
data(PLACE)

Base <- subset(x=PLACE, subset=crcBrd <= 9, select=c("xDev", "yDev", "tDev"))                

MeansBase <- colMeans(Base)           

CovBase <- cov(Base)    

Next <- subset(x=PLACE, 
               subset=crcBrd > 9, 
               select=c("xDev", "yDev", "tDev"))                             

library(qcc)           

Mqcc <- mqcc(data=PLACE[, c("xDev", "yDev", "tDev")], 
             confidence.level=(1-0.0000152837)^3, # default to (1-0.0027)^Nvariables
             type="T2.single",        
             center=MeansBase,        
             cov=CovBase,             
             plot=FALSE)


#summary(Mqcc)


###################################################
### PlotT2ChartPlace
###################################################

invisible(mqcc(data=PLACE[, c("xDev", "yDev", "tDev")], 
               confidence.level=(1-0.0000152837)^3, # default to (1-0.0027)^Nvariables
               type="T2.single",        
               center=MeansBase,        
               cov=CovBase, 
               plot=TRUE))

rm(Base, CovBase, MeansBase, Next, Mqcc)


###################################################
### MahalanobisT2Diss
###################################################
data(DISS)

mahalanobisT2(DISS[, c("batch", "min15", "min90")], 
              factor.name="batch", 
              compare.to=c(15, 15))


###################################################
### PlotScatterDiss
###################################################
ggplot(data=DISS, 
       aes(x=min15, y=min90, colour=batch)) + 
  geom_point(alpha=8/10, size=4) + 
  scale_color_grey(start=0.2, end=0.8) + 
  theme_bw() +
  theme(legend.position = "top", 
        legend.box = "horizontal")


###################################################
### PlotMahalanobisT2Diss
###################################################
invisible(mahalanobisT2(DISS[, c("batch", "min15", "min90")], 
                        factor.name="batch", 
                        plot=TRUE,
                        compare.to=c(15, 15)))


###################################################
### Chap10End
###################################################
rm(ALMPIN, PLACE, DISS)
detach(package:qcc)
detach(package:mistat)
detach(package:ggplot2)
detach(package:grid)
