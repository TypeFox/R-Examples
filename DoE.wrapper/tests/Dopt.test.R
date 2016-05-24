require(DoE.wrapper)
## try out all available designs with and without factor names
set.seed(1234)

## design with factor.names and constraint
   plan <- Dopt.design(36,factor.names=list(eins=c(100,250),zwei=c(10,30),drei=c(-25,25)),
                          nlevels=c(4,3,6), 
                          formula=~quad(.), 
                          constraint="!(eins>=200 & zwei==30 & drei==25)")
   summary(plan)
   design.info(plan)
   run.order(plan)
   cor(plan)
   y <- rnorm(36)
   r.plan <- add.response(plan, y)
   summary(r.plan)
   summary(lm(r.plan))
   #plan2 <- Dopt.augment(r.plan, m=10)
   #cor(plan2)

## design with candidates and constraint
candplan <- expand.grid(eins=c(100,150,200,250),zwei=c(10,20,30),drei=c(-25,-15,-5,5,15,25))
planc <- Dopt.design(36, candplan, formula=~quad(.), 
                          constraint="!(eins>=200 & zwei==30 & drei==25)", center=TRUE)
planc
cor(desnum(planc)[,-1])

## design with blocking without wholeBlockData (i.e. blocked, not splitplot)
planc <- Dopt.design(36, candplan, formula=~quad(.), 
                          constraint="!(eins>=200 & zwei==30 & drei==25)", center=TRUE, 
                          blocks=3)
summary(planc)
cor(desnum(planc)[,-1])

## design with blocking without wholeBlockData (i.e. blocked, not splitplot)
## variable block sizes
planc <- Dopt.design(36, candplan, formula=~quad(.), 
                          constraint="!(eins>=200 & zwei==30 & drei==25)", center=TRUE, 
                          blocks=c(6,6,12,12))
summary(planc)
cor(desnum(planc)[,-1])

## design with blocking with wholeBlockData (splitplot)
within<-expand.grid(A=c(-1,0,1),B=c(-1,0,1),C=c(-1,0,1))
whole<-expand.grid(D=factor(1:3),E=factor(1:3))

planc <- Dopt.design(54, within, formula=~D+E*(quad(A,B,C)), 
                          center=TRUE, 
                          blocks=rep(6,9), wholeBlockData=whole)
summary(planc)

whole <- data.frame(semester=1:3,reader=c(1,2,1))
planc <- Dopt.design(36, candplan, formula=~semester+reader+(eins+zwei+drei)^2,center=TRUE, 
                          constraint="!(eins>=200 & zwei==30 & drei==25)",
                          blocks=c(12,12,12), wholeBlockData=whole)
summary(planc)
cor(desnum(planc)[,-1])
r.planc <- add.response(planc, rnorm(36))
summary(lm(r.planc))