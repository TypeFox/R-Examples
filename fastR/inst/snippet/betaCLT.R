sampleMeans05 <- replicate(1000,mean(rbeta(5,0.5,0.5)))
sampleMeans10 <- replicate(1000,mean(rbeta(10,0.5,0.5)))
sampleMeans20 <- replicate(1000,mean(rbeta(20,0.5,0.5)))
sampleMeans40 <- replicate(1000,mean(rbeta(40,0.5,0.5)))
betaSim <- data.frame(
    mean=c(sampleMeans05,sampleMeans10,sampleMeans20,sampleMeans40),
    size=rep(c(5,10,20,40),each=1000))
myplot <- qqmath(~mean|factor(size),betaSim,scales=list(relation='free'))
