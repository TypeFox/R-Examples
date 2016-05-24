y <- c( seq(0, 1, by=0.01) )
n <- c(1,5,10,20)
n <- rep(n,times=length(y))
y <- rep(y,each=4)

density <-  n * (1-y)^{n-1}
groups <- paste("n =",n)
groups <- factor(groups,levels=unique(groups))
myplot <- xyplot(density~y, groups=groups, type="l", 
            main="Pdf of the mininum of a sample from Unif(0,1)",
            key=simpleKey(levels(groups),columns=2,lines=TRUE,points=FALSE),
            xlim=c(0,0.20))
