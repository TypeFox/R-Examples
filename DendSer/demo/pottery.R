
# Pottery example as in Section 4.1 of Advances in dendrogram seriation for application to visualization, D. Earle and C. Hurley

library(DendSer)

if (!("PairViz" %in% rownames(installed.packages()))) install.packages("PairViz")
if (!("HSAUR2" %in% rownames(installed.packages()))) install.packages("HSAUR2")


library(PairViz)
library(HSAUR2)

data(pottery)
sitename <- c("Gloucester","Llanedryn","Caldicot","IsleThorns","AshleyRails")[pottery[,"kiln"]]
pottery$Site <- factor(sitename)
reg <- c("Gloucester","Wales","Wales","NewForest","NewForest")[pottery[,"kiln"]]
pottery$Region <- factor(reg)


od <- c(2,9,6,1,3,7,4,5,8) # random order
dat <- scale(pottery[,od])

d <- dist(dat) 
h <- hclust(d, method="average")
o <- DendSer(h, d, costBAR)
# o<- dser(d,cost=costBAR)  shortcut for previous 2 steps


# Heatmaps of distance matix
cols  <- sequential_hcl(50)
dev.new(width=3, height=3); par(mar=c(1,1,1,1))
plotAsColor(as.matrix(d), h$order, col=cols, main="",rank=T,useRaster=T)
  
dev.new(width=3, height=3); par(mar=c(1,1,1,1))
plotAsColor(as.matrix(d), o, col=cols, main="",rank=T,useRaster=T) # Figure 9(a)



# Glyph array

clus <- cutree(h, 3)
table(clus,pottery$Site)
table(clus,pottery$Region)

gcols <- c("#89b7e5","#B2F0A2","red3")  


dev.new(width=3, height=3)
par(mar=c(0,0,0,0))
stars(dat[o,], col.stars=gcols[clus[o]],  lwd=.05,labels=NULL, cex=1,radius=F) # Figure 9(b)

# pcp, no ordering of vars

dev.new(width=7, height=2.25)
pcp(dat,  col=gcols[clus], horiz = TRUE, mar=c(3,2,1,2),
        main="", xaxs = "i") # Figure 10(a)
  
# Order vars to highlight cluster separation


panelMerit <- function(dats,clus,...){

  gmean <- aggregate(dats,list(clus),mean)[,-1]
  nv <- ncol(dats)
  dv <- matrix(0,nv,nv)

  for (i in 2:nv) {
	  for (j in 1:(i-1)){
	    cij <- gmean[,c(i,j)]
	    dv[i,j] <- sum(dist(cij))
	    dv[j,i] <- dv[i,j]
	}
	}
  dv
 	}

dats <-   apply(dat, 2, function(x) (x - min(x))/(max(x) - 
            min(x)))
            
merit <- panelMerit(dats,clus)            
ov<- dser(as.dist(max(merit)-merit),cost=costLPL)



dev.new(width=7, height=2.25)
pcp(dat[,ov],  col=gcols[clus], horiz = TRUE, mar=c(3,2,1,2),
        main="", xaxs = "i")  # Figure 10(b)

