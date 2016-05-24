results<-as.data.frame(t(replicate(1000,beta.mom(rbeta(50,2,5)))))
plot1 <- xhistogram(~shape1, results, type='density',v=2)
plot2 <- xhistogram(~shape2, results, type='density',v=5)
