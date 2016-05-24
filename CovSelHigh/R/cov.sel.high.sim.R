cov.sel.high.sim <- function(N, Setting, rep, Models){

  x111213 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x11 <- as.numeric(x111213[,1] > 0)
  x12 <- x111213[,2]
  x13 <- as.numeric(x111213[,3] > 0)
  x1415 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x14 <- as.numeric(x1415[,1] > 0)
  x15 <- x1415[,2]
  x1617 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x16 <- as.numeric(x1617[,1] > 0)
  x17 <- x1617[,2]
  x1819 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x18 <- x1819[,1]
  x19 <- x1819[,2]
  x20 <- rbinom(N, 1, prob = 0.5)
  x212223 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x21 <- as.numeric(x212223[,1] > 0)
  x22 <- x212223[,2]
  x23 <- as.numeric(x212223[,3] > 0)
  x2425 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x24 <- as.numeric(x2425[,1] > 0)
  x25 <- x2425[,2]
  x2627 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x26 <- as.numeric(x2627[,1] > 0)
  x27 <- x2627[,2]
  x2829 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x28 <- x2829[,1]
  x29 <- x2829[,2]
  x30 <- rbinom(N, 1, prob = 0.5)
  x313233 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x31 <- as.numeric(x313233[,1] > 0)
  x32 <- x313233[,2]
  x33 <- as.numeric(x313233[,3] > 0)
  x3435 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x34 <- as.numeric(x3435[,1] > 0)
  x35 <- x3435[,2]
  x3637 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x36 <- as.numeric(x3637[,1] > 0)
  x37 <- x3637[,2]
  x3839 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x38 <- x3839[,1]
  x39 <- x3839[,2]
  x40 <- rbinom(N, 1, prob = 0.5)
  x414243 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x41 <- as.numeric(x414243[,1] > 0)
  x42 <- x414243[,2]
  x43 <- as.numeric(x414243[,3] > 0)
  x4445 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x44 <- as.numeric(x4445[,1] > 0)
  x45 <- x4445[,2]
  x4647 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x46 <- as.numeric(x4647[,1] > 0)
  x47 <- x4647[,2]
  x4849 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x48 <- x4849[,1]
  x49 <- x4849[,2]
  x50 <- rbinom(N, 1, prob = 0.5)
  x515253 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x51 <- as.numeric(x515253[,1] > 0)
  x52 <- x515253[,2]
  x53 <- as.numeric(x515253[,3] > 0)
  x5455 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x54 <- as.numeric(x5455[,1] > 0)
  x55 <- x5455[,2]
  x5657 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x56 <- as.numeric(x5657[,1] > 0)
  x57 <- x5657[,2]
  x5859 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x58 <- x5859[,1]
  x59 <- x5859[,2]
  x60 <- rbinom(N, 1, prob = 0.5)
  x616263 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x61 <- as.numeric(x616263[,1] > 0)
  x62 <- x616263[,2]
  x63 <- as.numeric(x616263[,3] > 0)
  x6465 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x64 <- as.numeric(x6465[,1] > 0)
  x65 <- x6465[,2]
  x6667 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x66 <- as.numeric(x6667[,1] > 0)
  x67 <- x6667[,2]
  x6869 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x68 <- x6869[,1]
  x69 <- x6869[,2]
  x70 <- rbinom(N, 1, prob = 0.5)
  x717273 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x71 <- as.numeric(x717273[,1] > 0)
  x72 <- x717273[,2]
  x73 <- as.numeric(x717273[,3] > 0)
  x7475 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x74 <- as.numeric(x7475[,1] > 0)
  x75 <- x7475[,2]
  x7677 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x76 <- as.numeric(x7677[,1] > 0)
  x77 <- x7677[,2]
  x7879 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x78 <- x7879[,1]
  x79 <- x7879[,2]
  x80 <- rbinom(N, 1, prob = 0.5)
  x818283 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x81 <- as.numeric(x818283[,1] > 0)
  x82 <- x818283[,2]
  x83 <- as.numeric(x818283[,3] > 0)
  x8485 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x84 <- as.numeric(x8485[,1] > 0)
  x85 <- x8485[,2]
  x8687 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x86 <- as.numeric(x8687[,1] > 0)
  x87 <- x8687[,2]
  x8889 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x88 <- x8889[,1]
  x89 <- x8889[,2]
  x90 <- rbinom(N, 1, prob = 0.5)
  x919293 <- mvrnorm(n=N, mu=c(0,0,0), Sigma=cbind(c(1,0,0.25),c(0,1,0.25),c(0.25,0.25,1)) )
  x91 <- as.numeric(x919293[,1] > 0)
  x92 <- x919293[,2]
  x93 <- as.numeric(x919293[,3] > 0)
  x9495 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x94 <- as.numeric(x9495[,1] > 0)
  x95 <- x9495[,2]
  x9697 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
  x96 <- as.numeric(x9697[,1] > 0)
  x97 <- x9697[,2]
  x9899 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
  x98 <- x9899[,1]
  x99 <- x9899[,2]
  x100 <- rbinom(N, 1, prob = 0.5)

##Unconfoundedness holds given X 
if(Setting==1){
      if(Models=="Linear"){
                x_12 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x1 <- as.numeric(x_12[,1] > 0)
                x2 <- x_12[,2]     
                x3 <- rbinom(N, 1, prob = 0.5)
                x4 <- rnorm(N, mean = 0, sd = 1)
                x_56 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x5 <- x_56[,2]
                x6 <- as.numeric(x_56[,1] > 0)
                x_78 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
                x7 <- x_78[,1]
                x8 <- x_78[,2]
                x9 <- rnorm(N, mean = 0, sd = 1)
                x10 <- rbinom(N, 1, prob = 0.5)
                e0<-rnorm(N)
                e1<-rnorm(N)
                p <- 1/(1+exp(3 -2*x1 -1*x2 -2*x3 -1*x4 -2*x7))
                treat <- rbinom(N, 1, prob = p)
                y0 <- 2 + 4*x1 + 2*x2 + 2*x5  + 4*x6 + 4*x8 +e0
                y1 <- 4 + 4*x1 + 2*x2 + 2*x5  + 4*x6 + 4*x8 +e1
                Y <-ifelse(treat == 1, y1, y0)
                T<-treat
                }
     if(Models=="Nonlinear"){
                x_12 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x1 <- as.numeric(x_12[,1] > 0)
                x2 <- x_12[,2]     
                x3 <- rbinom(N, 1, prob = 0.5)
                x4 <- rnorm(N, mean = 0, sd = 1)
                x_56 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x5 <- x_56[,2]
                x6 <- as.numeric(x_56[,1] > 0)
                x_78 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
                x7 <- x_78[,1]
                x8 <- x_78[,2]
                x9 <- rnorm(N, mean = 0, sd = 1)
                x10 <- rbinom(N, 1, prob = 0.5)
                e0<-rnorm(N)
                e1<-rnorm(N)
                p <- 1/(1+exp(3 -2*x1 -1*x2 -2*x3 -1*x4 -2*x7))
                treat <- rbinom(N, 1, prob = p)
                #y0 <- 2 - 6*x6/log((x1 + 1.4)^2) + 7*x2 + 2*x5  + 4*x8  + e0
                #y1 <- 4 - 9*x6/log((x1 + 1.4)^3) + 2*x2 + 2*x5  + 4*x8  + e1
                y0 <- 2 - 6*x6/(0.5+(x2 + 1.4)^2) + 7*x1 + 2*x5^2  + 4*x8  + e0
                y1 <- 6.4 - 9*x6/(0.5+(x2 + 1.4)^4) + 4*x1 + 2*x5^2  + 4*x8  + e1
                
                
                Y <-ifelse(treat == 1, y1, y0)
                T<-treat
     }
  if(Models=="Binary"){
    x_12 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
    x1 <- as.numeric(x_12[,1] > 0)
    x2 <- x_12[,2]     
    x3 <- rbinom(N, 1, prob = 0.5)
    x4 <- rnorm(N, mean = 0, sd = 1)
    x_56 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
    x5 <- x_56[,2]
    x6 <- as.numeric(x_56[,1] > 0)
    x_78 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
    x7 <- x_78[,1]
    x8 <- x_78[,2]
    x9 <- rnorm(N, mean = 0, sd = 1)
    x10 <- rbinom(N, 1, prob = 0.5)
    e0<-rnorm(N)
    e1<-rnorm(N)
    p <- 1/(1+exp(3 -2*x1 -1*x2 -2*x3 -1*x4 -2*x7))
    treat <- rbinom(N, 1, prob = p)
    f<-4*x1 + 2*x2 + 2*x5  + 4*x6 + 4*x8
    py0 <- 1/(1+exp(-2 + f))
    py1 <- 1/(1+exp(-4 + f))
    y0<-rbinom(N, 1, py0)
    y1<-rbinom(N, 1, py1)
    Y <-ifelse(treat == 1, y1, y0)
    T<-treat
  }
            }


  
#M-bias given X
if(Setting==2){
      if(Models=="Linear"){
                U1<-rnorm(N, mean = 0, sd = 1)
                U2<-rnorm(N, mean = 0, sd = 1)
                U3<-rnorm(N, mean = 0, sd = 1)
                ex9<-rnorm(N,mean=0,sd=0.5)
                ex4<-rnorm(N,mean=0,sd=0.5)
                x_12 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x1 <- as.numeric(x_12[,1] > 0)
                x2 <- x_12[,2]      
                x3 <- rbinom(N, 1, prob = 0.5)
                x4 <- 0.2+0.8*U3+ex4
                x_56 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x5 <- x_56[,2]
                x6 <- as.numeric(x_56[,1] > 0)
                x_78 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
                x7 <- x_78[,1]
                x8 <- x_78[,2]
                x9 <- 1+2*U1+3*U2+ex9
                x10 <- rbinom(N, 1, prob = 0.5)
                e0<-rnorm(N)
                e1<-rnorm(N)
                p <- 1/(1+exp(3 -2*x1 -1*x2 -2*x3 -1*x4 -2*x7-1*U1))
                treat <- rbinom(N, 1, prob = p)
                y0 <- 2 + 4*x1 + 2*x2 + 2*x5  + 4*x6 + 4*x8+7*U2+2*U3 +e0
                y1 <- 4 + 4*x1 + 2*x2 + 2*x5  + 4*x6 + 4*x8 +7*U2+2*U3+e1
                Y <-ifelse(treat == 1, y1, y0)
                T<-treat
                }
      if(Models=="Nonlinear"){
                U1<-rnorm(N, mean = 0, sd = 1)
                U2<-rnorm(N, mean = 0, sd = 1)
                U3<-rnorm(N, mean = 0, sd = 1)
                ex9<-rnorm(N,mean=0,sd=0.5)
                ex4<-rnorm(N,mean=0,sd=0.5)
                x_12 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x1 <- as.numeric(x_12[,1] > 0)
                x2 <- x_12[,2]      
                x3 <- rbinom(N, 1, prob = 0.5)
                x4 <- 0.2+0.8*U3+ex4
                x_56 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
                x5 <- x_56[,2]
                x6 <- as.numeric(x_56[,1] > 0)
                x_78 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
                x7 <- x_78[,1]
                x8 <- x_78[,2]
                x9 <- 1+2*U1+3*U2+ex9
                x10 <- rbinom(N, 1, prob = 0.5)
                e0<-rnorm(N)
                e1<-rnorm(N)
                p <- 1/(1+exp(3 -2*x1 -1*x2 -2*x3 -1*x4 -2*x7-1*U1))
                treat <- rbinom(N, 1, prob = p)
                #y0 <- 2 - 6*x6/log((x1 + 1.4)^2) + 7*x2 + 2*x5  + 4*x8 +7*U2+2*U3 + e0
                #y1 <- 4 - 9*x6/log((x1 + 1.4)^3) + 2*x2 + 2*x5  + 4*x8 +7*U2+2*U3 + e1
                y0 <- 2 - 6*x6/(0.5+(x2 + 1.4)^2) + 7*x1 + 2*x5^2  + 4*x8 +7*U2+2*U3 + e0
                y1 <- 6.4 - 9*x6/(0.5+(x2 + 1.4)^4) + 4*x1 + 2*x5^2  + 4*x8 +7*U2+2*U3 + e1
                
                
                Y <-ifelse(treat == 1, y1, y0)
                T<-treat
      }
  if(Models=="Binary"){
    U1<-rnorm(N, mean = 0, sd = 1)
    U2<-rnorm(N, mean = 0, sd = 1)
    U3<-rnorm(N, mean = 0, sd = 1)
    ex9<-rnorm(N,mean=0,sd=0.5)
    ex4<-rnorm(N,mean=0,sd=0.5)
    x_12 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
    x1 <- as.numeric(x_12[,1] > 0)
    x2 <- x_12[,2]      
    x3 <- rbinom(N, 1, prob = 0.5)
    x4 <- 0.2+0.8*U3+ex4
    x_56 <- mvrnorm(n=N, mu=c(0,0), Sigma=cbind(c(1,0.5),c(0.5,1)) )
    x5 <- x_56[,2]
    x6 <- as.numeric(x_56[,1] > 0)
    x_78 <- rmvbin(N, bincorr=cbind(c(1,0.7),c(0.7,1)), margprob=c(0.5,0.5))
    x7 <- x_78[,1]
    x8 <- x_78[,2]
    x9 <- 1+2*U1+3*U2+ex9
    x10 <- rbinom(N, 1, prob = 0.5)
    e0<-rnorm(N)
    e1<-rnorm(N)
    p <- 1/(1+exp(3 -2*x1 -1*x2 -2*x3 -1*x4 -2*x7-1*U1))
    treat <- rbinom(N, 1, prob = p)
    f<-4*x1 + 2*x2 + 2*x5  + 4*x6 + 4*x8+7*U2+2*U3
    py0 <- 1/(1+exp(-2 + f))
    py1 <- 1/(1+exp(-4 + f))
    y0<-rbinom(N, 1, py0)
    y1<-rbinom(N, 1, py1)
    Y <-ifelse(treat == 1, y1, y0)
    T<-treat
  }
}
  
  
  
a<-data.frame(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,
                x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50,x51,x52,x53,x54,x55,x56,x57,x58,x59,x60,
                x61,x62,x63,x64,x65,x66,x67,x68,x69,x70,x71,x72,x73,x74,x75,x76,x77,x78,x79,x80,x81,x82,x83,x84,x85,x86,x87,x88,x89,x90,
                x91,x92,x93,x94,x95,x96,x97,x98,x99,x100,Y ,T)

binind<-c(1,3,6,7,8,10,11,13,14,16,18,19,20,21,23,24,26,28,29,30,31,33,34,36,38,39,40,41,43,44,46,48,49,50,
          51,53,54,56,58,59,60,61,63,64,66,68,69,70,71,73,74,76,78,79,80,81,83,84,86,88,89,
          90,91,93,94,96,98,99,100)

a[,binind]<-lapply(a[,binind],as.factor)
a[,101]<-as.numeric(a[,101])
a[,102]<-as.numeric(a[,102])
da <- a  

ut <- list(dat=da)
ut
  
}

