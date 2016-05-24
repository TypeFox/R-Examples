dat10 <- data.frame(
            x = rnorm(8*10),                     # 8 samples of size 10
            size=rep(10,8*10),                   # record sample size
            sample=rep(1:8,each=10)              # record sample number
            )
dat25 <- data.frame(
            x = rnorm(8*25),                     # 8 samples of size 25
            size=rep(25,8*25),                   # record sample size
            sample=rep(1:8,each=25)              # record sample number
            )
dat100 <- data.frame(
            x = rnorm(8*100),                    # 8 samples of size 100
            size=rep(100,8*100),                 # record sample size
            sample=rep(1:8,each=100)             # record sample number
            )
simdata <- rbind(dat10,dat25,dat100)

# generate the normal-quantile plots for each of the 30 samples
myplot <- qqmath(~x|sample*factor(size),data=simdata,
    layout=c(8,3), as.table=TRUE,cex=0.5)     
