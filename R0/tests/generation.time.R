#Loading package
library(R0)

# GT for children at house(from Cauchemez PNAS 2011)

GT.chld.hsld1<-generation.time("empirical", c(0,0.25,0.2,0.15,0.1,0.09,0.05,0.01))
plot(GT.chld.hsld1, col="green")
GT.chld.hsld1
# Discretized Generation Time distribution
# mean: 2.729412 , sd: 1.611636 
# [1] 0.00000000 0.29411765 0.23529412 0.17647059 0.11764706 0.10588235 0.05882353
# [8] 0.01176471

GT.chld.hsld2<-generation.time("gamma", c(2.45, 1.38))
GT.chld.hsld2
# Discretized Generation Time distribution
# mean: 2.504038 , sd: 1.372760
# [1] 0.0000000000 0.2553188589 0.3247178420 0.2199060781 0.1144367560
# [6] 0.0515687896 0.0212246257 0.0082077973 0.0030329325 0.0010825594
#[11] 0.0003760069 0.0001277537


# GT for school & community
GTs1<-generation.time("empirical", c(0,0.95,0.05))
plot(GTs1, col='blue')


plot(GT.chld.hsld1, ylim=c(0,0.5), col="red")
par(new=TRUE)
plot(GT.chld.hsld2, xlim=c(0,7), ylim=c(0,0.5), col="black")
