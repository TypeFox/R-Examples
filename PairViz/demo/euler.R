require(PairViz)

# Comparison of algorithms for generating tours on K9	
dev.new(width=6,height=2)
par(mar=c(2,2,2,1))
par(mfcol=c(1,2))
par(tcl = -.2, cex.axis=.7,mgp=c(3,.3,0),cex.main=.8,cex=.8)




h8 <- hpaths(8,matrix=FALSE)
h9 <- hpaths(9,matrix=FALSE)



plot(1:length(h8),h8,main= "hpaths(8)",pch=20,ylim=c(1,9),xlim=c(0,38),panel.first=grid(0,NULL,col="grey70"))
lines(1:length(h8),h8,col="grey70")

lines(8:9,h8[8:9],col="grey20")
lines(16:17,h8[16:17],col="grey20")
lines(24:25,h8[24:25],col="grey20")

plot(1:length(h9),h9,main= "hpaths (9)" ,pch=20,ylim=c(1,9),xlim=c(0,38),panel.first=grid(0,NULL,col="grey70"))
lines(1:length(h9),h9,col="grey70")
lines(9:11,h9[9:11],col="grey20")
lines(18:20,h9[18:20],col="grey20")
lines(27:29,h9[27:29],col="grey20")
lines(1:2,h9[1:2],col="grey20")
lines(36:37,h9[36:37],col="grey20")


#----------------------------


dev.new(width=6,height=3.8)
par(mar=c(2,2,2,1))
par(mfrow=c(2,2))
par(tcl = -.2, cex.axis=.7,mgp=c(3,.3,0),cex.main=.8)

n <- 8	

e1 <- eseq(n)
e2 <- eulerian(n)
e3 <- hpaths(n,matrix=FALSE)
x <- 1:length(e1)


plot(x,e1,main="eseq(8)",pch=20,ylim=c(1,9),xlim=c(0,38),panel.first=grid(0,NULL,col="grey70"))
lines(x,e1,col="grey70")

plot(x,e2,main="eulerian(8)" ,pch=20,ylim=c(1,9),xlim=c(0,38),panel.first=grid(0,NULL,col="grey70"))
lines(x,e2,col="grey70")


n <- 9	

e1 <- eseq(n)
e2 <- eulerian(n)
e3 <- hpaths(n,matrix=FALSE)
x <- 1:length(e1)
plot(x,e1,main="eseq(9)",pch=20,ylim=c(1,9),xlim=c(0,38),panel.first=grid(0,NULL,col="grey70"))
lines(x,e1,col="grey70")
lines(x[1:7],e1[1:7],col="grey20")
#lines(x[30:37],e1[30:37],col="grey20")

plot(x,e2,main="eulerian(9)" ,pch=20,ylim=c(1,9),xlim=c(0,38),panel.first=grid(0,NULL,col="grey70"))
lines(x,e2,col="grey70")
lines(x[1:7],e2[1:7],col="grey20")
#lines(x[30:37],e2[30:37],col="grey20")


#----------------------------


# Real data example
d <- eurodist
e1 <- eseq(length(labels(d)))
e2 <- eulerian(d)

require(gclus)
path1 <- order.hclust(-d,method="average")
#require(seriation)
#path1 <- get_order(seriate(d,method="gw",control=(list(method="average"))),1) # same as order.hclust
#path1 <- get_order(seriate(d,method="olo",control=(list(method="average"))),1))) # alternative

e4 <- weighted_hpaths(d,path1= path1,matrix=TRUE)
e3 <- weighted_hpaths(d,path1= path1,matrix=FALSE)



dev.new(width=6,height=6)
par(mar=c(2,4,2,1))


layout(matrix(1:3,nrow=3),heights=c(1,1,1))

par(tcl = -.2, cex.axis=1,mgp=c(2,.3,0),cex.main=1.2,cex=.5)

pw <- path_weights(dist2edge(d),e1)


plot(1:length(pw),pw, main= "Algorithm eseq: Eurodist edge weights",ylab="distance",col="grey70",type="l")
points(1:length(pw),pw,pch=20,col="grey30")

pw <- path_weights(dist2edge(d),e2)
plot(1:length(pw),pw,type="l",col="grey70",main= "Weighted eulerian on Eurodist",ylab="distance")
points(1:length(pw),pw,pch=20,col="grey30")

pw <- path_weights(dist2edge(d),e3)
plot(1:length(pw),pw,type="l",col="grey70",main= "Weighted hamiltonians on Eurodist ",ylab="distance")
points(1:length(pw),pw,pch=20,col="grey30")


i <- 21
avepw <- rep(0,length(pw))
for (e in i*(1:10)){
	s <- e-i+1
   lines(c(e,e),c(0,5000),col="grey70",lty=3)
   avepw[s:e] <- mean(pw[s:e])
}


lines(1:length(pw),avepw,col="grey10")



