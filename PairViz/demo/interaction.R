####### Code for constructing figures of Section 3.2 (Hurley and Oldford, 2008) ###################

#Rat survival data from Box G and Cox D (1964) "An analysis of transformations" J. Roy. #Stat. Soc. Series B. 26 211.


surv <- c(3.1,4.5,4.6,4.3,8.2,11,8.8,7.2,4.3,4.5,6.3,7.6,4.5,7.1,6.6,6.2,
3.6,2.9,4.0,2.3,9.2,6.1,4.9,12.4,4.4,3.5,3.1,4.0,5.6,10,7.1,3.8,2.2,2.1,1.8,2.3,3.0,3.7,3.8,2.9,2.3,2.5,2.4,2.2,3.0,3.6,3.1,3.3)

treat <- rep(c("A","B","C","D"),times=3,each=4)
poison <- rep(c("P1","P2","P3"),each=16)
d <- data.frame(surv,treat,poison )
a <- aov(surv~treat+poison+treat:poison,data=d)

summary(a)


# Produce matrix of cell means- 

m <- tapply(surv, list(poison,treat), mean)

require(PairViz)	 			

####### Function to construct a plot ###################



inter_plot <- function(x,m,cols,xaxt="n",ylab="Mean response", xlab="Treatment",lx,ly,...){
	par(mar=c(3,3,3,3))
	par(cex=.8)
	matplot(x,t(m),type="o",col=cols,xaxt=xaxt,ylab=ylab, 
	pch=1,  lty=1, xlab=xlab,...)
    axis(1,x,colnames(m))
   legend(lx,ly,legend=rownames(m), col=cols,fill=cols,box.lty=0,cex=.8)
}


cols <- c("red","green","blue","purple")

#------------------------------------

dev.new(width=3.5,height=3)
inter_plot(1:ncol(m),m,cols,main="Dataset order: h_0",lx=3.3,ly=9)

# Could construct this one with interaction.plot but difficult for plots with repeated treatments
interaction.plot(d$treat,d$poison,d$surv, col=cols,  lty=c(1,1,1),ylab="Mean response", xlab="Treatment", main="Dataset order: h_0")
# add the points
interaction.plot(treat,poison,surv, col=cols, type="p", pch=22,legend=FALSE,add=TRUE,ann=FALSE, axes=FALSE)


#------------------------------------
o <- c(1,4,3,2,NA,4,2,1,3)
# o <- hpaths(1:ncol(m)); o <- c(o[1,],NA,o[2,])

x <- 1:ncol(m)
x <- c(x,ncol(m)+1,x+ncol(m))
dev.new(width=5.5,height=3)
inter_plot(x,m[,o],cols,main="Hamiltonian decomposition: h_1 and h_2",lx=7.5,ly=8)
#------------------------------------

o <-hpaths(1:4,matrix=FALSE)

dev.new(width=5.5,height=3)
inter_plot(1:8,m[,o],cols,main="Eulerian path",lx=7.5,ly=8.5)
#------------------------------------
mm <- scale(m, apply(m,2,mean), FALSE)

o <- c(1,4,3,2,NA,4,2,1,3)
# o <- hpaths(1:ncol(m)); o <- c(o[1,],NA,o[2,])

x <- 1:ncol(m)
x <- c(x,ncol(m)+1,x+ncol(m))

dev.new(width=5.5,height=3)
inter_plot(x,mm[,o],cols,main="Adjusted interaction plot: h_1 and h_2",lx=7.5,ly=-2,ylab="Adjusted mean  response")

#------------------------------------

mdiff <- m[c(1,1,2),]- m[c(2,3,3),]
rownames(mdiff) <- c("P1-P2","P1-P3", "P2-P3")


dcols <- c("sienna"  ,"magenta","cyan")
o <- c(1,4,3,2,NA,4,2,1,3)
# o <- hpaths(1:ncol(m)); o <- c(o[1,],NA,o[2,])

x <- 1:ncol(m)
x <- c(x,ncol(m)+1,x+ncol(m))


dev.new(width=5.5,height=3)

inter_plot(x,mdiff[,o],dcols,main="Differenced interaction plot: h_1 and h_2",lx=7,ly=5.4,ylab="Mean  response differences")

#------------------------------------

a <- aov(surv ~ treat+poison+treat*poison,data=d)
hsd <- TukeyHSD(a,conf.level = .99)[[3]]
#yusr <- range(pretty(hsd[,2:3]))
yusr <- c(-10,10)
bcols <- rep(c("red","green","blue"), each=4)


dev.new(width=4,height=7)
par(mar=c(2,4,2,4))
par(mfrow=c(3,1))
par(mgp=c(2,1,0))


dd <- split(d, list(treat,poison))
dd <- lapply(dd, function(x) x[,1])
treats <- levels(d[,"treat"])
poisons <- levels(d[,"poison"])
ot <- hpaths(1:4, matrix=FALSE)

for (i in 1:3) {
  o <- ot+ 4*(i-1)
  if (i==1) main <- "Rat survival by treatment within poison" else main <- ""
  mc_plot(dd,a,path=o,ylim=range(d$surv),ylab=poisons[i], 
  col=bcols,   main=main, names= treats[ot],lty=1,cifunction=function(a,level) TukeyHSD(a,conf.level= level)[[3]], ci.yusr=yusr,ci.ex=2,sig.lwd = 1)
}
