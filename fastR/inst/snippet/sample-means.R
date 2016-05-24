# 1000 sample means of samples of size 16 from N(100,12):
sampleMeans <- replicate(5000,mean(rnorm(16,100,12)))
mean(sampleMeans)
sd(sampleMeans)
myplot<-xhistogram(~sampleMeans,n=20,v=100,
    density=TRUE,args=list(mean=100,sd=3))
