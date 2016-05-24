### R code from vignette source 'untbpaper.Rnw'

###################################################
### code chunk number 1: overalloptions
###################################################



###################################################
### code chunk number 2: untbpaper.Rnw:90-91
###################################################
require("untb",quietly=TRUE)


###################################################
### code chunk number 3: time_saver
###################################################
calc_from_scratch <- FALSE


###################################################
### code chunk number 4: generate_rn
###################################################
if(calc_from_scratch){
  set.seed(0)
  rn <- rand.neutral(5e6, theta=50)
} else {
  load("rn.Rdata")
}


###################################################
### code chunk number 5: SaundersSummary
###################################################
data("saunders")
summary(saunders.tot)


###################################################
### code chunk number 6: prestonSaunders
###################################################
preston(saunders.tot,n=9)


###################################################
### code chunk number 7: prestonSaundersTemp
###################################################
jj.preston <- preston(saunders.tot)
jj.summary <- summary(saunders.tot)


###################################################
### code chunk number 8: calculate_uncertainty_Saunders
###################################################
n <- 10
J <- no.of.ind(saunders.tot)
unc <- list()
theta_hat <- optimal.theta(saunders.tot)
for(i in 1:n){
  unc[[i]] <- rand.neutral(J=J, theta=theta_hat)
}


###################################################
### code chunk number 9: plotSaunders
###################################################
plot(saunders.tot,uncertainty=FALSE)
for(i in 1:n){
  points(seq_along(unc[[i]]),unc[[i]],type="l",col="grey")
}


###################################################
### code chunk number 10: optimalThetaSaunders
###################################################
optimal.theta(saunders.tot)


###################################################
### code chunk number 11: supportTheta
###################################################
S <- no.of.spp(saunders.tot)
J <- no.of.ind(saunders.tot)
theta <- seq(from=25,to=39,len=55)
jj <- theta.likelihood(theta=theta,S=S,J=J,give.log=TRUE)
support <- jj-max(jj)
plot(theta,support,xlab=paste("Biodiversity parameter",expression(theta)),ylab="support")
abline(h= -2)


###################################################
### code chunk number 12: calculate_mle
###################################################
f <- function(local_size,gens){
  jj <- isolate(rn,size=local_size)
  a <- untb(start=jj, prob=0.01, D=local_size, gens=gens, meta=rn)
  optimal.params(a)
}

if(calc_from_scratch){
  x100   <- t(replicate(100,f(100   ,999)))
  x1000  <- t(replicate(100,f(1000  ,999)))
  x10000 <- t(replicate(100,f(10000 ,999)))
} else {
  load("dat.Rdata")
  x100   <- dat[,,1]
  x1000  <- dat[,,2]
  x10000 <- dat[,,3]
}
  


###################################################
### code chunk number 13: estimateMandTheta
###################################################
plot(x100,log="xy",xlim=c(1,80),ylim=c(0.001,1),col="black",pch=1,
main=paste("Maximum likelihood estimates of ", expression(m), " and ",expression(theta))
     )
points(x1000,col="red",pch=2) 
points(x10000,col="blue",pch=3) 

points(50,0.01,pch=4,lwd=3,cex=2)

legend( "bottomleft", c("100","1000","10000"),col=c("black","red","blue"), pch=1:3, title="Local community size")


###################################################
### code chunk number 14: e.lowandhigh
###################################################
n <- 20
x <- expected.abundance(J=n, theta=3)
e.low  <- expected.abundance(J=n,theta=4)
e.high <- expected.abundance(J=n,theta=2)


###################################################
### code chunk number 15: expectedAbundance
###################################################
plot(x)
segments(x0=1:n,x1=1:n,y0=e.low,y1=e.high)


###################################################
### code chunk number 16: calculate_thirdRank
###################################################
rank3 <- table(replicate(1000,rand.neutral(J=20,theta=2)[3]))


###################################################
### code chunk number 17: plot_thirdRank
###################################################
plot(rank3,xlab="abundance of third ranked species",ylab="frequency")


###################################################
### code chunk number 18: calculate_species_table
###################################################
if(calc_from_scratch){
  set.seed(0);
  synthetic.spp <- species.table(untb(start=rep(1,60),prob=0.002, gens=40000, keep=TRUE))
} else {
  load("synthetic_spp.Rdata")
}


###################################################
### code chunk number 19: matplot_species_table
###################################################
plot(1:10,xlim=c(0,40000),ylim=c(0,60),type="n",xlab="time (generation)",ylab="abundance")
"showabundance" <- function(x, ...){
  jj <- rle(x)
  x <- cumsum(jj$lengths)
  y <- jj$values
  segments(x0=c(1,x),y0=y,x1=c(x,x[length(x)]),y1=y, ...)  #horizontal ones
  segments(x0=x[-length(x)],y0=y[-1],x1=x[-length(x)],y1=y[-length(y)], ...)  # vertical ones
}
showabundance(synthetic.spp[,1],col="black")
showabundance(synthetic.spp[,2],col="red")
showabundance(synthetic.spp[,3],col="green")
showabundance(synthetic.spp[,4],col="blue")
showabundance(synthetic.spp[,5],col="cyan")
showabundance(synthetic.spp[,6],col="magenta")


if(FALSE){matplot(synthetic.spp,type="l",lty=1,xlab="time (generation)",ylab="abundance")}


###################################################
### code chunk number 20: SampleTenThousand
###################################################
set.seed(0)
jj <- isolate(rn,size=10000)
a <- untb(start=jj, prob=0.01, D=10000, gens=1000, meta=rn)
if(calc_from_scratch){
  a.logkda <- logkda(a)
} else {
  load("a_logkda.Rdata")
}
op <- optimal.params(a,log.kda=a.logkda)
v.opt <- volkov(no.of.ind(a), op, bins=TRUE)
v.true <- volkov(no.of.ind(a), c(100,0.01), bins=TRUE)


###################################################
### code chunk number 21: PlotSampleTenThousand
###################################################
pa <- preston(a,n=12)
pa.names <- sub(" ", "", names(pa))
jj <- plot(pa,ylim=c(0,27),axisnames=FALSE,
ylab="Number of species",xlab="Abundance class")
axis(1, at=jj, labels=FALSE, lty=0)
text(jj, par("usr")[3]-0.65, srt=90, cex=0.8, adj=1, labels=pa.names,xpd=TRUE)

points(jj, v.opt[1:12], type="b",col="red",pch=1)
points(jj, v.true[1:12], type="b",col="blue",pch=4)
par(xpd=2)
legend("topright", c("best estimate","true"), pch=c(1,4), col=c("red","blue"), lty=c(1,1))


###################################################
### code chunk number 22: differentThetas
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


