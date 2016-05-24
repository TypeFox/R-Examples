set.seed(1)                      # use fixed "random" seed
someData <- data.frame(x=runif(300),group=factor(rep(1:3,each=100)))
p <- qqmath(~x|group, data=someData, distribution=qunif)
print(p)
