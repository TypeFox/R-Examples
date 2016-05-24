# -- Load datasets

data(tortoise)
data(whale)
data(teasel)

#--------------------------------------------------------------------------------------
# section 5.3.1 Age-specific survival  
# Survivorship plot like figure 5.1 in Caswell.  


x<-splitA(whale)
whaleT<-x$T
whaleF<-x$F

# Note example on page 120 uses matrix powers and not element by element 
# which is R default.  Matrix power is not part of base R, but for simple cases 
# this works to do A %*% A %*% A %*% A...

mp<-function(A,pow){
  if(pow==1){A}
  else{ x<-A
    for(i in (2:pow)){  A<-x %*% A}}
  A
}

## use colSums for sum of matrix columns e^T
surv<-matrix(numeric(150*4), ncol=4)
for(x in 1:150)
{
   surv[x,]<-colSums(mp(whaleT,x))
}
## Just plot first stage column?

plot(surv[,1]/surv[1,1], type="l", ylim=c(0,1), las=1, main="Fig 5.1",
  xlab="Age (years)", ylab=expression(paste("Survivorship ", italic(l(x)))))

#--------------------------------------------------------------------------------------
# section 5.3.2 Age-specific fertility  
#  equation 5.44
T20<- mp(whaleT,20)
whaleF %*% T20 %*% diag(1/colSums(T20))

##  Figure 5.2 in Caswell
fert<-numeric(200)
for(x in 1:200)
{
  T20<-mp(whaleT,x)
  phi<-whaleF %*% T20 %*% diag(1/colSums(T20))
  fert[x]<-phi[1,1]
}

op<-par(mar=c(5,5,4,1))
plot(fert, type="l", ylim=c(0,0.07), las=1, main="Fig 5.2",
  xlab="Age (years)", ylab="")
mtext("Age-specific fertility", 2,3.5)
par(op)


#--------------------------------------------------------------------------------------
## Sensitivity plot like figure 9.3 in Caswell 
## use text to add labels closer to x-axis

A<-tortoise[["med.high"]]
sens<-sensitivity(A)
def.par <- par(no.readonly = TRUE) 
par(mfrow=c(3,1), mar = c(1.5, 4.5, 1, 2) ,  oma=c(1.5,0,1.5,0), cex=1.2 )
## F in top row
ep<-barplot(sens[1,], ylim=c(0,.4), col="white", las=1,  
   ylab=expression(paste("to ", italic(F[i]))), names="")
box()
 text(ep, -.05,  1:8, xpd = TRUE)
## P on diagonal
ep<-barplot(diag(sens), ylim=c(0,.4), col="white", las=1, 
   ylab=expression(paste("to ", italic(P[i]))), names="")
box()
 text(ep, -.05,  1:8, xpd = TRUE)
# G on subdiagonal
barplot(c(sens[row(sens)==col(sens)+1],0), ylim=c(0,.4), col="white", las=1, 
   ylab=expression(paste("to ", italic(G[i]))) )
box()
 text(ep, -.05,  c(1:7, NA), xpd = TRUE)
mtext(expression(paste("Fig 9.3. Sensitivity of ", lambda, "...")), 3, outer=TRUE, cex=1.4)
mtext(expression(paste("Size class ", italic(i))), 1, outer=TRUE, cex=1.2)
par(def.par)



sens<-sensitivity(teasel)

### IMAGE plot in 9.4b requires lattice 

# z<-cbind( expand.grid(fate=1:6,stage=1:6),x=as.vector(sens))
## re-order fate
# z$fate <- ordered(z$fate, levels = 6:1)
# n<-range(log10(z$x))
# at<-seq(n[1],n[2], length.out=48)
# levelplot(log10(x)~ stage*fate, data=z, at=at, col.regions=grey(1:48/48))
# title("Fig 9.4b Levelplot", line=2.5)



# Stair step plot like  figure 9.4 in Caswell 
plot(log10(c(sens)), type="s", las=1, ylim=c(-5, 2),
xlab="Stage at time t", xaxt="n",
ylab= expression(paste(Log[10], " sensitivity of ",lambda)), 
main="Fig 9.4c Stair step plot")
axis(1, seq(1,36,6), 1:6)
text(log10(c(sens)), cex=.7, adj=c(0,-.3), 
     labels=paste(" ", 1:6, rep(1:6,each=6), sep=""))





#--------------------------------------------------------------------------------------
## Triangle plot like figure 9.11 in Caswell but for tortoise ( not sea turtle).  
## Nicer triangle plots are found in many different packages.

elas<-elasticity(tortoise[[3]])
el<-c(F=sum(elas[1,]), P=sum(diag(elas)), G=sum(elas[row(elas)==col(elas)+1]))


plot(c(0, 1, 2, 0), c(0, sqrt(3), 0, 0), type = "l", lwd = 2, main="Fig 9.11 Triangle plot",
 xlab = "Tortoise summed elasticities", ylab = "", axes = FALSE)
text(c(0, 2, 1), c(0,0, sqrt(3)), names(el), cex = 1.5, pos=c(1,1,3), xpd=TRUE)
points(2 - 2 * el[1] - el[3], el[3] * sqrt(3), cex=1.5, pch=16, col='blue')


#--------------------------------------------------------------------------------------
# Random design and variance decomposition Fig 10.10

# whale Pods

pods<-c(
0,0.0067,0.1632,0,0.9535,0.8827,0,0,0,0.0802,0.9586,0,0,0,0.0414,0.9752,
0,0.0062,0.1737,0,1,0.9020,0,0,0,0.0694,0.9582,0,0,0,0.0418,0.9855,
0,0.0037,0.0988,0,0.9562,0.9030,0,0,0,0.0722,0.9530,0,0,0,0.0406,0.9798,
0,0.0043,0.1148,0,1,0.9015,0,0,0,0.0727,0.9515,0,0,0,0.0485,0.9667,
0,0.0042,0.1054,0,0.8165,0.8903,0,0,0,0.0774,0.9515,0,0,0,0.0485,0.9810,
0,0.0027,0.0732,0,1,0.9123,0,0,0,0.0730,0.9515,0,0,0,0.0485,0.9545,
0,0.0025,0.0651,0,1,0.9254,0,0,0,0.0746,0.9515,0,0,0,0.0485,0.9810,
0,0.0047,0.1159,0,1,0.9200,0,0,0,0.0800,0.9706,0,0,0,0.0294,0.9608,
0,0.0068,0.1761,0,1,0.9241,0,0,0,0.0759,0.9562,0,0,0,0.0438,1,
0,0.0061,0.1418,0,1,0.9167,0,0,0,0.0833,0.9286,0,0,0,0.0714,1,
0,0.0050,0.1251,0,1,0.9216,0,0,0,0.0784,0.9515,0,0,0,0.0485,0.9810,
0,0.0021,0.0542,0,1,0.9254,0,0,0,0.0746,0.9515,0,0,0,0.0485,0.9810,
0,0.0027,0.0732,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
0,0.0045,0.1220,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,1,
0,0.0052,0.1428,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
0,0.0037,0.0998,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
0,0.0047,0.1273,0,1,0.9286,0,0,0,0.0714,0.9515,0,0,0,0.0485,0.9810,
0,0.0024,0.0797,0,1,0.8929,0,0,0,0.0595,0.9515,0,0,0,0.0485,1)



## Covariance matrix
p1<-matrix(pods, nrow=18, byrow=TRUE)

# addcolumn names
colnames(p1)<- paste("a", rep(1:4,each=4), 1:4, sep="")
## re-order columns to plot matrix by columns 
x<- order(paste("a", 1:4, rep(1:4,each=4), sep=""))
p1<-p1[,x]
covmat<- cov(p1)

## plots matching figure 3 in Brault & Caswell 1993 or figure 10.10 in Caswell 2001 (p 272)
persp(covmat, theta=45, phi=15, box=FALSE,
 main=expression("Fig 10.10a Killer Whale Covariances"))
w1<- matrix(apply(p1, 2, mean), nrow=4)
wS<-sensitivity(w1)
## V matrix of contributions
contmat <- covmat * c(wS) %*% t(c(wS))
persp(contmat, theta=45, phi=15, box=FALSE,
main=expression(paste("Fig 10.10b Contributions to V(", lambda, ")")) )
## contributions of V associated with aij  (matrix on page 271 in Caswell)
A<-matrix(apply(contmat, 2, mean), nrow=4)
dimnames(A)<-dimnames(whale)
# matrix on page 271
round(A/sum(A),3)


