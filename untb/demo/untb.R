if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))


# Following plots and demos illustrate various functionality from the package.
# Some of the parameters have been altered in the interest of speed.

require(untb)


###################################################
### chunk number 3: SaundersSummary
###################################################
data(saunders)
summary(saunders.tot)


###################################################
### chunk number 4: prestonSaunders
###################################################
preston(saunders.tot,n=10)

###################################################
### chunk number 6: calculate_uncertainty_Saunders
###################################################
plot(extant(extractor(saunders.exposed,1)),uncertainty=TRUE)


###################################################
### chunk number 9: supportTheta
###################################################
S <- no.of.spp(saunders.tot)
J <- no.of.ind(saunders.tot)
theta <- seq(from=25,to=39,len=55)
jj <- theta.likelihood(theta=theta,S=S,J=J,give.log=TRUE)
support <- jj-max(jj)
plot(theta,support,xlab=paste("Biodiversity parameter",expression(theta)),ylab="support")
abline(h= -2)


###################################################
### chunk number 14: e.lowandhigh
###################################################
n <- 10
x <- expected.abundance(J=n, theta=3)
e.low  <- expected.abundance(J=n,theta=4)
e.high <- expected.abundance(J=n,theta=2)
plot(x)
segments(x0=1:n,x1=1:n,y0=e.low,y1=e.high)


###################################################
### chunk number 16: calculate_thirdRank
###################################################
rank3 <- table(replicate(100,rand.neutral(J=20,theta=2)[3]))
plot(rank3,xlab="abundance of third ranked species",ylab="frequency")


###################################################
### chunk number 18: calculate_species_table
###################################################
 {
set.seed(0);
a <- species.table(untb(start=rep(1,60),prob=0.02, gens=4000,keep=TRUE))
}
matplot(a,type="l",lty=1,xlab="time (generation)",ylab="abundance")


###################################################
### chunk number 20: SampleTenThousand
###################################################
set.seed(0)
rn <- rand.neutral(1e5, theta=50)
jj <- isolate(rn,size=1000)
a <- untb(start=jj, prob=0.01, D=1000, gens=1000, meta=rn)
a.logkda <- logkda(a)
op <- optimal.params(a,log.kda=a.logkda)
v.opt <- volkov(no.of.ind(a), op, bins=TRUE)
v.true <- volkov(no.of.ind(a), c(100,0.01), bins=TRUE)


###################################################
### chunk number 21: PlotSampleTenThousand
###################################################
pa <- preston(a,n=12)
pa.names <- sub(" ", "", names(pa))
jj <- plot(pa,ylim=c(0,15),axisnames=FALSE,
ylab="Number of species",xlab="Abundance class")
axis(1, at=jj, labels=FALSE, lty=0)
text(jj, par("usr")[3]-0.65, srt=90, cex=0.8, adj=1, labels=pa.names,xpd=TRUE)

points(jj, v.opt[1:12], type="b",col="red",pch=1)
points(jj, v.true[1:12], type="b",col="blue",pch=4)
par(xpd=2)
legend("topright", c("best estimate","true"), pch=c(1,4), col=c("red","blue"), lty=c(1,1))





lof <- function(gens,p,size){
  untb(start=isolate(rn,size=size), prob=p, gens=gens, D=100, meta=rn)
}

do.stuff <- function(n,p,gens,size){
  show.stuff <- FALSE
  ans <- matrix(0,n,2)
  colnames(ans) <- c("theta","m")
  for(i in 1:n){
    jj <- lof(gens,p,size)
    if(show.stuff){print(jj)}
    l <- logkda.pari(jj)
    ans[i,] <- optimal.params(D=jj, log.kda=l)
  }
  return(ans)
}


x100 <- do.stuff(10,0.01,300,100)
x200 <- do.stuff(10,0.01,300,200)
x300 <- do.stuff(10,0.01,300,300)

plot(rbind(x100,x200,x300),type="n",log="xy",xlab="theta",ylab="m",main="Maximum likelihood estimates of theta and m")
points(x100,col="black",pch=1)
points(x200,col="red",pch=2)
points(x300,col="blue",pch=3)
points(50,0.01,pch=4,lwd=3, cex=2)
legend("topright" , c("100","200","300"), col=c("black","red","blue"),pch=1:3,title="local community size")



###################################################
### chunk number 22: differentThetas
###################################################
set.seed(0)
f <- function(gens,p){
  display.untb(start=sample(as.census(untb(start=1:100,gens=gens,D=100,prob=p))),gens=0,main="",cex=1.7, asp=1)
}

g <- function(u="title", ...){
  par(srt=0)
  par(mai=c(0,0,0,0))
  plot.new()
  text(0.5,0.5,u,...)
}

h <- function(u="title", ...){
  par(mai=c(0,0,0,0))
  par(srt=90)
  plot.new()
  text(0.5,0.5,u, ...)
}

nf <- layout(matrix(
                    c(00,01,00,02,00,03,
                      04,05,00,06,00,07,
                      00,00,00,00,00,00,
                      08,09,00,10,00,11,
                      00,00,00,00,00,00,
                      12,13,00,14,00,15,
                      00,00,00,00,00,00,
                      16,17,00,18,00,19),8,6,byrow=TRUE),
             c(1,4, 1,4, 1,4),
             c(1,4, 1,4, 1,4, 1,4),
             TRUE)

g(expression(t==10))
g(expression(t==50))
g(expression(t==100))

h(expression(theta==0))
f(10,0)
f(50,0)
f(100,0)
h(expression(theta==0.1))
f(10,0.001)
f(50,0.001)
f(100,0.001)
h(expression(theta==1))
f(10,0.01)
f(50,0.01)
f(100,0.01)
h(expression(theta==10))
f(10,0.1)
f(50,0.1)
f(100,0.1)


