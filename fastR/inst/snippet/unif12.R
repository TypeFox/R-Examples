sampleSums <- replicate(2000,sum(runif(12,-0.5,0.5)))
myplot1<-qqmath(~sampleSums)
myplot2<-histogram(~sampleSums)
