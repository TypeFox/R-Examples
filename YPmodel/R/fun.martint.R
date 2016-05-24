fun.martint <-
function(data,best,u0,fb,kall,kk,RandomData,repnum,jh,...){


#set.seed(0)
#dev.new()

n <- data$length
oz <- data$Z
Z <- data$Z
od <- data$Delta
oy <- data$X
fb2 <- fb$Num2
fb2 <- matrix(fb2, nrow=n, ncol=1)
k1all <- kall$Num1
k2all <- kall$Num2


#-------------------------------------#
# Random Data Set
#-------------------------------------#
uniformRandomData1 <- RandomData$uniformRandomData1
uniformRandomData2 <- RandomData$uniformRandomData2
gaussianRandomDataAll <- RandomData$gaussianRandomDataAll

m <- 4
b <- cbind(t(best)+jh*matrix(c(1,0),nrow  <-  2, ncol  <-  1),t(best)-jh*matrix(c(1,0),nrow  <-  2, ncol  <-  1),t(best)+jh*matrix(c(0,1),nrow  <-  2, ncol  <-  1),t(best)-jh*matrix(c(0,1),nrow  <-  2, ncol  <-  1))

## ntitr0.m
data15 <- fun.ntitr0(u0,m,b,oz,od,Z,n)
u <- data15$u
s <- data15$s
ru <- data15$ru

qf <- u
pq <- cbind(qf[,1]-qf[,2],qf[,3]-qf[,4])/2/jh

data16 <- svd(-pq/n)
u <- data16$u
s <- data16$d
v <- data16$v
#s <- t(s) + 1.0e-8 * t(runif(2, min <- 0, max <- 1))
s <- t(s) + 1.0e-8 * t(uniformRandomData1)
s <- 1 / s
inq <- v%*%diag(c(s))%*%t(u)

bt <- exp(-best)
r <- u0
rl <- t(cbind(0,t(r[1:n-1])))
dr <- r-rl
den <- bt[1]+bt[2]*r
hr <- (1+r) / den
#denl <- bt[1]+bt[2]*rl
ff <- cbind(bt[1]/den,bt[2]*r/den)

ffl <- -ff

lam2 <- log(1+bt[2]/bt[1]*r)/bt[2]

fb2l <- t(cbind(1,t(fb2[1:(n-1),1])))
wt1 <- n /k2all
wt2 <- matrix(1,nrow  <-  n, ncol  <-  1)
wt3 <- -n*fb2l/k2all

tem <- k2all/den
tem2 <- tem*(bt[2]*hr-1)
tem3 <- k2all/k1all*hr
inr1 <- ffl[,1]*tem2*dr
inr2 <- ffl[,2]*tem2*dr
inr1 <- inr1+sum(inr1)-fun.cumsum(inr1)
inr2 <- inr2+sum(inr2)-fun.cumsum(inr2)
inr1 <- inr1/k1all
inr2 <- inr2/k1all

ind1 <- od*(1-oz)
ind2 <- od*oz
mu1 <- inr1-ffl[,1]*tem3
mu2 <- inr2-ffl[,2]*tem3

wk <- matrix(c(1:n), nrow=n, ncol=1)/n
wk <- 1+4*wk*(1-wk)
wk <- wk*k1all/(k1all+k2all)

#set.seed(0)
obs <- fun.cumsum(ind2*wk- wk*k2all*dr/den)/sqrt(n)

inr <- wk*tem2*dr
inr <- inr+sum(inr)-fun.cumsum(inr)
inr <- inr/k1all
wtt1 <- inr-tem3*wk
qe1 <- fun.cumsum((tem*wk/n)*dr*ff[,1])
qe2 <- fun.cumsum((tem*wk/n)*dr*ff[,2])
a <- cbind(qe1,qe2)

#rdsd <- ceiling(100*runif(100, min<-0, max<-1))
rdsd <- ceiling(100*uniformRandomData2)
rdsd[1] <- 1

mb1 <- matrix(1, nrow=repnum, ncol=1)
mb2 <-mb1

lineCount <- 0
wtildCount <- c()
for(irep in 1:repnum){
    #set.seed(0)
    #g <- rnorm(n, mean = 0, sd = 1)
    gTemp <- gaussianRandomDataAll[(1+(irep-1)*n):(irep*n)]
    g <- matrix(gTemp, nrow=n, ncol=1)
    dm1 <- ind1*g/sqrt(n)
    dm2 <- ind2*g/sqrt(n)

    q1 <- t(mu1)%*%dm1+t(ffl[,1])%*%dm2
    q2 <- t(mu2)%*%dm1+t(ffl[,2])%*%dm2

    wtild1 <- fun.cumsum(wk*dm2)
    wtild2 <- fun.cumsum(wtt1*dm1)
    wtild3 <- -a%*%inq%*%t(cbind(q1, q2))
    wtild <- wtild1+wtild2+wtild3
    mb1[irep] <- max(abs(wtild[1:kk]))

    if(min(abs(irep-rdsd))==0){
        lineCount <- lineCount + 1
        wtildCount <- cbind(wtildCount,wtild[1:kk])
        #if(irep==1){
            #dev.hold
        #}
    }

}

#plot(oy[1:kk]*365,obs[1:kk],"l",col="blue",xlab="Days", ylab="p-value")
#title(main="Process of the martingale residual-based test (p-value)")
#for(i in 1:lineCount){
#    lines(oy[1:kk],wtildCount[,lineCount],"l",col="red")
#}
#dev.hold

mobs1 <- max(abs(obs[1:kk]))
pvalu1 <- mean(mb1>mobs1)

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    output<- list(mobs1=mobs1,pvalu1=pvalu1,obs=obs,wtildCount=wtildCount,lineCount=lineCount)
    return(output)
#-----------------------------------------------------------------#


}
