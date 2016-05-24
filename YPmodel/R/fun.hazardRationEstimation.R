fun.hazardRationEstimation <-
function(Data,GroupData,ru,p,pl,bt,deni,kall,sm,gama,b,h,...){

## loading data
#-----------------------------------------------------------------#
#Temp parameters
n <- Data$length
Z <- Data$Z
X <- Data$X
Delta <- Data$Delta
#-----------------------------------------------------------------#

r <- ru
rl <- t(cbind(0,t(r[1:n-1])))
po <- p
plo <- pl
dr <- r-rl
rl <- r
denil <- deni

den <- bt[1]+bt[2]*r
hr <- (1+r) / den
denl <- bt[1]+bt[2]*rl
ff <- cbind(bt[1]/den,bt[2]*r/den)
ffl <- cbind(bt[1]/den,bt[2]*rl/denl)

inr1 <- ffl[,1]*(kall$Num1/sm)*kall$Num2*hr*(bt[2]*hr-1)/(1+r)^2/po*dr
inr2 <- ffl[,2]*(kall$Num1/sm)*kall$Num2*hr*(bt[2]*hr-1)/(1+r)^2/po*dr
inr1 <- inr1+sum(inr1)-fun.cumsum(inr1)
inr2 <- inr2+sum(inr2)-fun.cumsum(inr2)
rmul <- plo/sm
inr1 <- inr1*rmul
inr2 <- inr2*rmul

di <- deni
dil <- denil
xi1 <- Z*gama$Num1/dil-di/den*ffl[,1]*kall$Num2/sm+di*inr1
xi2 <- Z*gama$Num2*rl/dil-di/den*ffl[,2]*kall$Num2/sm+di*inr2
xi1 <- xi1*Delta
xi2 <- xi2* Delta
eti <- rmul* di* Delta

m <- 4
b <- cbind(b+h*matrix(c(1,0),nrow  <-  2, ncol  <-  1),b-h*matrix(c(1,0),nrow  <-  2, ncol  <-  1),b+h*matrix(c(0,1),nrow  <-  2, ncol  <-  1),b-h*matrix(c(0,1),nrow  <-  2, ncol  <-  1))

data11 <- fun.oldp2(b,m,Data)
s <- data11$s
ru <- data11$ru 
u <- data11$u
gama <- data11$gama 
p <- data11$p 
pl <- data11$pl
deni <- data11$deni 
sm <- data11$sm 

qf <- rbind(u[1,],u[2,])
pq <- cbind(qf[,1]-qf[,2],qf[,3]-qf[,4])/2/h
pr <- cbind((ru[,1]-ru[,2])/2/h,(ru[,3]-ru[,4])/2/h)

a <- matrix(1,nrow  <-  n, ncol  <-  2)
a[,1] <- hr*ff[,1]+(b[1]-b[2])/den^2*pr[,1]
a[,2] <- hr*ff[,2]+(b[1]-b[2])/den^2*pr[,2]
br <- (b[1]-b[2])/den^2/po

data12 <- svd(-pq/n)
u <- data12$u
s <- data12$d
v <- data12$v

s <- t(s) + 1.0e-8 * t(runif(2, min <- 0, max <- 1))
s <- 1 / s
inq <- v%*%diag(c(s))%*%t(u)

x1 <- matrix(1,nrow  <-  n, ncol  <-  2)
xn <- matrix(1,nrow  <-  n, ncol  <-  2)
x1[,1] <- -kall$Num2/sm*ffl[,1]*hr + inr1*(1+r)
x1[,2] <- -kall$Num2/sm*ffl[,2]*hr + inr2*(1+r)
xn[,1] <- kall$Num1/sm*ffl[,1] + inr1*den
xn[,2] <- kall$Num1/sm*ffl[,2] + inr2*den

nu1 <- rmul*(1+r)
nun <- rmul*den

d11 <- kall$Num1*dr/(1+r)
dnn <- kall$Num2*dr/den
sig1 <- a%*%inq%*%(t(x1)%*%((d11%*%matrix(1,nrow  <-  1, ncol  <-  2))*x1) + t(xn)%*%((dnn%*%matrix(1,nrow  <-  1, ncol  <-  2))*xn))%*%t(inq)

#############################
sig1 <- colSums(t(sig1*a))
#############################

sig1 <- t(sig1)/n
sig2 <- n*br^2*fun.cumsum(nu1^2*d11 + nun^2*dnn)
sig12 <- fun.cumsum(((d11*nu1)%*%matrix(1,nrow  <-  1, ncol  <-  2))*x1 + ((dnn*nun)%*%matrix(1,nrow  <-  1, ncol  <-  2))*xn)
sig12 <- sig12*(br%*%matrix(1,nrow  <-  1, ncol  <-  2))
sig12 <- colSums(t(sig12*(a%*%inq)))

sig <- t(sig1)+sig2+2*sig12

wt1 <- 1 + sig
sig <- sqrt(sig)
wt2 <- sig
wt3 <- hr

temp1 <- matrix(c(1:n), nrow  <-  n, ncol  <-  1)
yt <- temp1[Delta == 1]
k <- length(yt)
ow1 <- wt1[yt]
ow2 <- wt2[yt]
ow3 <- wt3[yt]

oow1 <- sort(ow1,decreasing=FALSE)

oow2 <- sort(ow2,decreasing=FALSE)

oow3 <- sort(ow3,decreasing=FALSE)


temp2 <- matrix(c(1:n), nrow  <-  n, ncol  <-  1)
a1 <- temp2[wt1<oow1[floor(k*.9)]]
a2 <- temp2[wt2<oow2[floor(k*.9)]]
a3 <- temp2[wt3<oow3[floor(k*.9)]]
ld1 <- min(a1)
ud1 <- max(a1)
ld2 <- min(a2)
ud2 <- max(a2)
ld3 <- min(a3)
ud3 <- max(a3)

repnum <- 1000
mb <- matrix(c(1:repnum), nrow  <-  repnum, ncol  <-  1)
mb21 <- mb
mb22 <- mb
mb23 <- mb


for (irep in 1:repnum){
	g <-  (1 + 1/sqrt(min(c(GroupData$length$Num1,GroupData$length$Num2)))) * rnorm(n, mean=0, sd=1)
	intdr <- cumsum(eti*g)
	wtild1 <- hr*(ff%*%(inq%*%rbind(t(xi1)%*%g,t(xi2)%*%g)))/sqrt(n)
	wtild2 <- (b[1]-b[2])/den^2*(pr%*%(inq%*%rbind(t(xi1)%*%g,t(xi2)%*%g)))/sqrt(n)
	wtild3 <- sqrt(n)*br*intdr
	wtild <- wtild1+wtild2+wtild3

	bn1 <- wtild/wt1
	bn2 <- wtild/wt2
	bn3 <- wtild/wt3

	mb21[irep] <- max(abs(bn1[ld1:ud1]))
	mb22[irep] <- max(abs(bn2[ld2:ud2]))
	mb23[irep] <- max(abs(bn3[ld3:ud3]))
}

m22 <- sort(mb22,decreasing=FALSE)

ca22 <- m22[repnum*0.95]
upp22 <- exp(ca22/sqrt(n)*wt2/hr)*hr
low22 <- exp(-ca22/sqrt(n)*wt2/hr)*hr

ca3 <- 1.96;
upp3 <- exp(ca3/sqrt(n)*wt2/hr)*hr
low3 <- exp(-ca3/sqrt(n)*wt2/hr)*hr
ca90 <- m22[repnum*0.9]
upp90 <- exp(ca90/sqrt(n)*wt2/hr)*hr
low90 <- exp(-ca90/sqrt(n)*wt2/hr)*hr

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    output<- list(hr=hr,ld2=ld2,ud2=ud2,upp22=upp22,low22=low22,upp3=upp3,low3=low3,upp90=upp90,low90=low90)
    return(output)
#-----------------------------------------------------------------#

}
