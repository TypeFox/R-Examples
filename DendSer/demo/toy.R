
# Toy examples as in Section 3.3 of Advances in dendrogram seriation for application to visualization, D. Earle and C. Hurley
# Constructs toy figs

library(seriation)


if (!("mvtnorm" %in% rownames(installed.packages()))) install.packages("mvtnorm")

library(mvtnorm)




# like Chen seriation in package seriation, except spilt point found as described by Chen
r2e <- function(dis)
{
     dis <- as.matrix(dis)
    while (qr(dis)$rank > 2) dis <- cor(dis)

    e <- eigen(dis)$vectors[,1:2]
    p.group <- which(e[,1] >= 0)
    p.order <- p.group[order(e[p.group, 2], decreasing=TRUE)]

    m.group <- which(e[,1] < 0)
    m.order <- m.group[order(e[m.group, 2])]

    n1 <- length(m.order); n2 <- length(p.order)

    dis1 <- dis[m.order[1], p.order[n2]]
    dis2 <- dis[m.order[n1], p.order[1]]

    a <- min(dis1, dis2)#find best joining point

    if (a==dis1) o <- c(m.order, p.order)
    else o <- c(m.order[n1:1], p.order[n2:1])

    return(o)
}


calcorders <- function(xd){
  xd1 <- dist(xd)
  xd2 <- as.matrix(xd1)
  hxd <- hclust(xd1,"average")
  xo <- hxd$order
  xolo <- get_order(seriate(xd1,method="OLO", control=list(method="average")))
  # xpl <- DendSer(hxd,as.matrix(xd1), cost=costPL, node_op="c0",direction="down",GW=T)
  xgw <- reorder.hclust(hxd,xd1)$order
  xarc <- DendSer(hxd,as.matrix(xd1), cost=costARc, node_op="r01",direction="down",GW=T)
  xbar <- DendSer(hxd,as.matrix(xd1), cost=costBAR, node_op="r01",direction="down",GW=T)
  cs <- cmdscale(xd1,1)

  xlsmds <- DendSer(hxd, ser_weight=cs,cost=costLS)
  xarsa<- get_order(seriate(xd1,method="ARSA"))
  xtsp <- get_order(seriate(xd1,method="TSP",control=list(method="farthest_insertion")))
  xchen<- r2e(xd1)
  xmds <- order(cs)
  xx <- cbind(xgw, xtsp,xolo,xbar,xarc,xlsmds,xarsa,xchen,xmds)

  colnames(xx) <- c("gw","tsp","olo","bar","arc","ls-mds","arsa","r2e","mds")
  xx
}


plotorders1 <- function(xd,ords,xname="",col="grey70",f=2,main=T){
lcol <- "black"

dev.new(width=ncol(ords)*f,height=f)
par(mar=c(1,1,1,1))
par(mar=c(.2,.2,.2,.2))

par(mfrow=c(1,ncol(ords)))


for (i in 1:ncol(ords)){
	o <- ords[,i]
	m<- paste(xname, colnames(ords)[i])
	print(plot(xd,col=col,xaxt="n", yaxt="n",bty="o",ann=F,asp=1,fg = gray(0.7))); lines(xd[o,], col=lcol)
	if(main) title(m)
}
}


plotorders2 <- function(xd,ords,xname="",f=2,main=T,...){
	
xd1 <- dist(xd)
xd2 <- as.matrix(xd1)

dev.new(width=ncol(ords)*f,height=f)
par(mar=c(1,1,1,1))
par(mar=c(.2,.2,.2,.2))

par(mfrow=c(1,ncol(ords)))


for (i in 1:ncol(ords)){
	o <- ords[,i]
	m<- paste(xname, colnames(ords)[i])
	print(plotAsColor(xd2, ann=F, order=o,useRaster=T,rank=T,...))
	if (main) title(m)
}
}




#-------------------------


hcol <-sequential_hcl(50)

#  linear pattern
xd1 <-rmvnorm(100, mean=c(0,0), sigma=matrix(c(1,.95,.95,1),nrow=2))

ords1 <- calcorders(xd1)
m1 <- c(3:6)
m2 <- c("tsp","arsa","r2e")

plotorders1(xd1,ords1[,m1],"biv",f=2,main=T)
plotorders2(xd1,ords1[,m1],"biv",col=hcol,f=2,main=T)

plotorders1(xd1,ords1[,m2],"biv",f=2,main=T)
plotorders2(xd1,ords1[,m2],"biv",col=hcol,f=2,main=T)



#-------------------------
# u curve. pl/olo/bar/arc follow curve, mds does not

xd2 <-matrix(runif(200),ncol=2)
xd2[,1]<- 3*xd2[,1]
xd2[,2]<- 1.4*xd2[,2]+5*sin(xd2[,1])
xd2<-scale(xd2)

ords2 <- calcorders(xd2)


# Uncomment next 4 lines for plots of curve data 
# plotorders1(xd2,ords2[,m1],"curve",f=2,main=T)
# plotorders2(xd2,ords2[,m1],"curve",col=hcol,f=2,main=T)

# plotorders1(xd2,ords2[,m2],"curve",f=2,main=T)
# plotorders2(xd2,ords2[,m2],"curve",col=hcol,f=2,main=T)


#-------------------------



# circular data

xd3 <-rmvnorm(100, mean=c(0,0), sigma=0.3*diag(2))
inc <- seq(0,2*pi, length.out= 100)
xd3[,1] <-  10*cos(inc) - 4*abs(xd3[,1])
xd3[,2] <- 10*sin(inc) - 4*abs(xd3[,2])

ords3 <- calcorders(xd3)

# Uncomment next 4 lines for plots of circle data 

# plotorders1(xd3,ords3[,m1],"circle",main=T)
# plotorders2(xd3,ords3[,m1],"circle",col=hcol,main=T)
# plotorders1(xd3,ords3[,m2],"circle",main=T)
# plotorders2(xd3,ords3[,m2],"circle",col=hcol,main=T)




#-------------------------

# 4 clusters


nc<-4
means <- 4*matrix(c(1,1,1,0,0,1,0,0),byrow=T,nrow=nc) +0.1*matrix(runif(nc*2),byrow=T,nrow=nc)
xd4 <-rmvnorm(100, mean=c(0,0), sigma=.3*diag(2))

mem <- rep(1:4,each=25)
xd4 <- xd4+ means[mem,]



ords4 <- calcorders(xd4)

# Uncomment next 4 lines for plots of cluster data 


# plotorders1(xd4,ords4[,m1],"clusters",main=T)
# plotorders2(xd4,ords4[,m1],"clusters",col=hcol,main=T)

# plotorders1(xd4,ords4[,m2],"clusters",main=T)
# plotorders2(xd4,ords4[,m2],"clusters",col=hcol,main=T)



