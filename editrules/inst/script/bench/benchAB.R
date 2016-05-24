# Some benchmarks on elimination Algorithm A and B
# mvdl, 05.01.2011
#

cat('-------------------------------------------------------------------\n')
cat('This benchmark requires packages editrules and Hmisc\n')
require(editrules)
if (!require(Hmisc)) stop('install Hmisc')

fdir <- system.file('script/bench',package='editrules')
edits <- file.path(fdir,'edits.R')
files <- file.path(fdir,c('eliminator.R','randomEdits.R'))

for ( f in files ) source(f)

##  The double conversion makes sure that the empty character ('') 
##  is interpreted as empty value.
E <- editarray(as.character(editfile(edits,type='cat')))

vars <- getVars(E)

cat('timing algorithms A and B...\n')

t1 <- system.time(e1 <- elimInfo(E, algorithmA.follow))
t2 <- system.time(e2 <- elimInfo(E, algorithmB.follow))

for ( i in 1:99 ){
    t1 <- t1 + system.time(e1 <- elimInfo(E, algorithmA.follow))
    t2 <- t2 + system.time(e2 <- elimInfo(E, algorithmB.follow))
}


## eliminate one variable from E, do this for all variables and follow progress of elimination
tA <- tB <- numeric(length(vars))
names(tA) <- names(tB) <- vars
for ( v in vars ){ 
    tA[v] <- system.time({for ( i in 1:50 ) u <- algorithmA(E,v)})['user.self']/50*1000
    tB[v] <- system.time({for ( i in 1:50 ) u <- algorithmB(E,v)})['user.self']/50*1000
}
times <- cbind(tA,tB,diff=tB-tA,ratio=tB/tA)
o <- order(times[,3],decreasing=TRUE)
times <- cbind(times,e1$Created > e2$Created)

totA <- sum(times[,'tA'])
totB <- sum(times[,'tB'])
cat('Ratio tB/tA over all variables :', totB/totA,'\n')

totA <- sum(times[rev(o)[1:13],'tA'])
totB <- sum(times[rev(o)[1:13],'tB'])
cat('Ratio tB/tA over variables where iterations are necessary :', totB/totA,'\n\n')


cat('Creating plots\n')

## Create dotchart of timing difference and scatterplot of timings. 

par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(4,2,2,1.1))

lb <- rep("",36)
lb[seq(1,36,5)] <- seq(1,36,5)
pch <- 1+numeric(nrow(times))
I <- e1$Created > e2$Created
pch[I] <- 16
dotchart(times[o,'diff'], 
    pch=pch[o],
    main=expression(paste('Timing ',t[B]-t[A])),
    xlab=expression(t[B]-t[A]),
    labels=lb
)
abline(v=0,lty=1,col='grey')

pc1 <- 2 + numeric(36)
pc1[I] <- 17
plot(e1$contained,times[,'tA'],
    xlab='Nr of edits containing variable',
    main='Timing A and B (ms)',
    xlim=c(0,20),
    ylim=c(0,40),
    pch=pc1,
    oma=c(1,1,1,1))
pc2 <- 22 + numeric(36)
pc2[I] <- 15
points(e2$contained,times[,'tB'],pch=pc2)
legend('topright',
    legend=c('A','B'),
    pch=c(2,22),
    bty='n')


cat('-------------------------------------------------------------------\n')

cat('-------------------------------------------------------------------\n')
cat('Determine complexity as function of number of edits.\n')
M <- c(seq(10,100,10),150,200,250,500)
timing <- array(0,dim=c(length(M),2))
colnames(timing) <- c('A','B') 
i <- 0

timeAlg <- function(E,alg,n=10){
    vars <- getVars(E)
    t = 0
    t = t + system.time({ for(v in vars ) E <- alg(E,v)})['user.self']
    t/n
}

nAverage = 100
TA <- array(0,dim=c(length(M),8))
colnames(TA) <-c(names(summary(1:10)), 
        names(quantile(1:10,probs=c(0.05,0.95))  ))
TB <- TAB <- TA 

i <- 0
ms=1/1000
for ( m in M ){
    cat('\rworking on m=',m)
    tA <- numeric(nAverage)
    tB <- numeric(nAverage)
    tBA <- numeric(nAverage)
    i <- i + 1
    for ( j in 1:nAverage ){
        E <- genedits(m,group='A')
        vars <- getVars(E)
        tA[j] <- timeAlg(E,algorithmA,n=3)
        tB[j] <- timeAlg(E,algorithmB,n=3)
    }
    TA[i,] <- c(summary(tA), quantile(tA,probs=c(0.05,0.95)))/ms
    TB[i,] <- c(summary(tB), quantile(tB,probs=c(0.05,0.95)))/ms
    TAB[i,] <- c(summary(tA-tB),quantile(tA-tB,probs=c(0.05,0.95)))/ms
}
    cat('\n')


dev.new()
par(mfrow=c(1,2),oma=c(0,0,0,0),mar=c(4,2,2,1.1))
errbar(M,sqrt(TA[,'Mean']),sqrt(TA[,'5%']),sqrt(TA[,'95%']),
    xlab='nr of edits',
    ylab=expression(paste(sqrt(time),' (',sqrt(ms),')')),
    )
errbar(M,sqrt(TB[,'Mean']),sqrt(TB[,'5%']),sqrt(TB[,'95%']),add=TRUE,pch=2)
legend('topleft',
    legend=c('A','B'),
    pch=c(16,2),
    bty='n')
title(main=expression(paste('Timing (',sqrt(ms),')')))

errbar(M,TAB[,'Mean'],TAB[,'5%'],TAB[,'95%'],
    xlab='nr of edits',
    ylab=expression(paste(t[A]-t[B],'  (ms)'))
    )
abline(h=0)
title(main=expression(paste(t[A]-t[B],' ','(ms)')))

cat('finished\n')








