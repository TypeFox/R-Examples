library(OptInterim)

### Example 1 ###
### single group and two group 2-stage without pause
B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)  ### p0=.4, p1=.6 at x=1
x <- 1


set.seed(12357)
object12 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="ETSL",
sf="futility",num.arm=1,num.stage=2,control=OptimDesControl(n.int=c(1,5)),pause=0)
simout12<-SimDes(object12,sim.n=10000)
simout12adj<-SimDes(object12,sim.n=10000,CMadj=TRUE)
simout12_2 <- SimDes(object12,sim.n=10000,m.init = c(5, 5, 25, 25, 25))
simout12_3 <- SimDes(object12,sim.n=10000,interimRule = "t1",attainI = 0.8)
save.image('examples.RData')

set.seed(12357)
object12_np <- np.OptimDes(B.init,m.init,alpha,beta,param,x,pn=1.1,target="ETSL",
sf="futility",num.arm=1,num.stage=2,control=OptimDesControl(n.int=c(1,5)))

### Example 2 ###
### single group and two group 2-stage with pause=3
B.init <- 1:72
m.init <- rep(3,72)
alpha <- 0.10
beta <- 0.2
x <- 6
pnull<-.45
palt<-.6
param <- c(1, weibPmatch(x,pnull,shape=1), 
           1, weibPmatch(x,palt,shape=1))  ### p0=.45, p1=.6 at x=1

set.seed(12357)
object12P3 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="ETSL",
sf="futility",num.arm=1,num.stage=2,control=OptimDesControl(n.int=c(1,5)),pause=3)
#simout12P3<-SimDes(object12P3,sim.n=10000)
#simout12P3adj<-SimDes(object12P3,sim.n=10000,CMadj=TRUE)





### Example 3 ###
###2-group three stage designs with pause=0.3
B.init <- c(1, 2, 3, 4, 5)
m.init <- 4*c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
x<-1
#p0=.2, p1=.35 at x=1
param <- c(1.5, 0.7281438, 1.75, 0.9725991)

set.seed(12357)
object23P3 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="ES",sf="OF",
num.arm=2,num.stage=3,control=OptimDesControl(aboveMin=c(1.05,1.10)),pause=0.3)



save.image('examples.RData')
q('no')



