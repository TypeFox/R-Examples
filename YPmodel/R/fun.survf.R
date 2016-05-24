fun.survf <-
function(data,best,u0,fb,kall,kk,l,dfb1,RandomData,repnum,jh,...){

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
l2 <- l$Num2

#repnum <- 1000
#-------------------------------------#
# Random Data Set
#-------------------------------------#
uniformRandomData1 <- RandomData$uniformRandomData1
uniformRandomData2 <- RandomData$uniformRandomData2
gaussianRandomDataAll <- RandomData$gaussianRandomDataAll

#if(repnum*n>1000000){
#    tempNum <- repnum*n
#    discrepancy <- repnum*n - 1000000
#    set.seed(0)
#    guassTemp <- rnorm(discrepancy, mean = 0, sd = 1)
#    gaussianRandomDataAll <- t(cbind(gaussianRandomDataAll,guassTemp))
#}
#-------------------------------------#


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

ff <- cbind(bt[1]/den,bt[2]*r/den)

ffl <- -ff
dff2 <- fun.Rjudge(fb2,0,'>')
dfb2 <- fb2*dff2+(1-dff2)

lam2 <- log(1+bt[2]/bt[1]*r)/bt[2]
obs1 <- sqrt(n)*(l2-lam2)
obs21 <- sqrt(n)*fb2
obs21 <- sqrt(n)*exp(-l2)
fbm2 <- exp(-lam2)
obs22 <- sqrt(n)*fbm2
obs2 <- obs21-obs22

tem <- k2all/den*(bt[2]*hr-1)*dr

tem3 <- k2all/k1all*hr
inr1 <- ffl[,1]*tem
inr2 <- ffl[,2]*tem
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

#rdsd <- ceiling(100*runif(100, min<-0, max<-1))
rdsd <- ceiling(100*uniformRandomData2)
rdsd[1] <- 1
#repnum <- 1000
mb1 <- matrix(1, nrow=repnum, ncol=1)
mb2 <-mb1

lineCount <- 0
wtildCount <- c()
for(irep in 1:repnum){
    #g <- rnorm(n, mean = 0, sd = 1)
    #gTemp <- gaussRandomData[,irep]
    gTemp <- t(gaussianRandomDataAll)[(1+(irep-1)*n):(irep*n)]
    g <- (1+1/sqrt(n))*matrix(gTemp, nrow=n, ncol=1)
    #g <- (1+1/sqrt(n))*rnorm(n, mean = 0, sd = 1)
    dm1 <- ind1*g/sqrt(n)
    dm2 <- ind2*g/sqrt(n)

    q1 <- t(mu1)%*%dm1+t(ffl[,1])%*%dm2
    q2 <- t(mu2)%*%dm1+t(ffl[,2])%*%dm2

    c1 <- 1/den/dfb1

    wtild1 <- -c1*fun.cumsum(dm1/(k1all/n))
    wtild2 <- fun.cumsum(dm2/(k2all/n))
    a <- cbind(r/den,lam2-r/den)
    wtild3 <- -a%*%inq%*%t(cbind(q1, q2))
    wtild <- wtild1+wtild2+wtild3
    
    ce1 <- c1*fbm2
    ce2 <- -fb2

    wtild1 <- -wtild1*fbm2
    wtild2 <- -wtild2*fb2
    ce3 <- -fbm2
    wtild3 <- -wtild3*fbm2
    wtild <- wtild1+wtild2+wtild3

    #ckmdvar2;%
    #------------------------------------------------------------------------------------------------------#
    ttmm <- matrix(c(sum(mu1^2 *ind1+ffl[,1]^2 *ind2),sum(mu1*mu2 *ind1+ffl[,1]*ffl[,2] *ind2),sum(mu2^2*ind1+ffl[,2]^2 *ind2)), nrow=1, ncol=3)
    varq <- matrix(c(ttmm[1],ttmm[2],ttmm[2],ttmm[3]), nrow=2, ncol=2)/n
    aa <- (fbm2%*%matrix(1,nrow=1, ncol=2))*a
    sig1 <- aa%*%inq%*%varq%*%t(inq)
    sig1 <- colSums(t(sig1*aa))
    sig1 <- matrix(sig1, nrow=n,ncol=1)
    sig2 <- ce2^2 *fun.cumsum(ind2 /((k2all/n)^2))/n+ce1^2 *fun.cumsum(ind1 /((k1all/n)^2))/n
    sig12 <- (ce1%*%matrix(1,nrow=1, ncol=2))*fun.cumsum(cbind(mu1*ind1/(k1all/n),mu2*ind1/(k1all/n)))/n
    sig12 <- colSums(t(sig12*(aa%*%inq)))
    sig12 <- matrix(sig12, nrow=n,ncol=1)
    ttm <- (ce2%*%matrix(1,nrow=1, ncol=2))*fun.cumsum(cbind(ffl[,1]*ind2/(k2all/n),ffl[,2]*ind2/(k2all/n)))/n
    ttm <- colSums(t(ttm *(aa%*%inq)))
    #ttm <- t(ttm)
    sig12 <- sig12+ttm
    sig <- sig1+sig2+2*sig12
    #------------------------------------------------------------------------------------------------------#

    rts1 <- sqrt(sig)
    drr1 <- fun.Rjudge(rts1,0,'>')
    drts1 <- rts1*drr1+(1-drr1)

    wtild <-wtild/drts1
    wtild <-wtild*drr1

    mb2[irep] <- max(wk[1:kk]*abs(wtild[1:kk]))

#    if  min(abs(irep-rdsd))==0
#        plot(oy(1:kk),wtild(1:kk),'r:')
#        if irep==1
#        hold
#        end
#    end

    if(min(abs(irep-rdsd))==0){
        lineCount <- lineCount + 1
        wtildCount <- cbind(wtildCount,wtild[1:kk])
        #if(irep==1){
            #dev.hold
        #}
    }

}

obs2 <- obs2/drts1
obs2 <- obs2*drr1

#plot(oy[1:kk]*365,obs2[1:kk],"l",col="blue",xlab="Days", ylab="p-value")
#title(main="Process of the contrast-based test (p-value)")
#for(i in 1:lineCount){
#    lines(oy[1:kk]*365,wtildCount[,lineCount],"l",col="red")
#}

mobs2 <- max(wk[1:kk]*abs(obs2[1:kk]))

pvalu2 <- mean(mb2>mobs2)

#-----------------------------------------------------------------#
## Output Resuts 
#-----------------------------------------------------------------#
    output<- list(mobs2=mobs2,pvalu2=pvalu2,obs2=obs2,wtildCount=wtildCount,lineCount=lineCount)
    return(output)
#-----------------------------------------------------------------#


}
