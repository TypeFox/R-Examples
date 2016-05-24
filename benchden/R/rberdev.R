gen01 <- function(n) {runif(n, min=0, max=1)}                   # U(0,1)
gen02 <- function(n) {rexp(n)}                                  # Exp
gen03 <- function(n) {tmp <- runif(n,min=0,max=1)               # Maxwell
            outvek <- sqrt(-2*log(1-tmp))
            outvek
            }
gen04 <- function(n) {                                          # Double Exponential 
    vek<-rexp(n,rate=1)
    signvek<-sample(x=c(-1,1),size=n,replace=T)
    vek*signvek
    } 
gen05 <- function(n) {rlogis(n,location=0,scale=1)}             # logistic
gen06 <- function(n) {rcauchy(n, location=0, scale=1)}                                #Cauchy
gen07 <- function(n) { tmp <- runif(n,min=0,max=1)
            outvek <- -log(-log(tmp))
            outvek
            }
gen08 <- function(n) {tmp <- runif(n,min=0,max=1)               # 
            outvek <- tmp^2
            outvek
            }
gen09 <- function(n) {tmp <- runif(n,min=0,max=1)               # 
            outvek <- 1/((1-tmp)^2)
            outvek
            }
gen10 <- function(n) {tmp <- runif(n,min=0,max=1)               # 
            signvek<-sample(x=c(-1,1),size=n,replace=T)
            outvek <- 1/((1-tmp)^2)-1
            outvek <- outvek*signvek
            outvek
            }
    
gen11 <- function(n) {rnorm(n, mean=0, sd=1)}                                #N(0,1)
gen12 <- function(n) {rlnorm(n, meanlog = 0, sdlog = 1)}                         #Standard Lognormal
gen13 <- function(n) {outvek<-rep(0,n)                                   #BD 26 (trimodal uniform)
            component<-sample(2,size=n,replace=T,prob=c(0.5,0.5))
            outvek[component==1] <- runif(sum(component==1),min=-0.5, max=0.5 )
            outvek[component==2] <- runif(sum(component==2),min=-5, max=5 )
            outvek
}
gen14 <- function(n) {tmp <- runif(n,min=0,max=1)                
            outvek <- rep(0,n)
            outvek[tmp!=0.5]<-sign(tmp-0.5)*(exp(sign(tmp-0.5)/(0.5-tmp)))[tmp!=0.5]
            outvek
            }
      
gen15 <- function(n) {runif(n,min=0,max=1)*runif(n,min=0,max=1)}
gen16 <- function(n) {runif(n,min=-1,max=0)+runif(n,min=0,max=1)}
gen17 <- function(n) {rbeta(n, shape1=2, shape2=2)}                     # Beta
gen18 <- function(n) {rchisq(n,df=1)}
gen19 <- function(n) {rnorm(n,mean=0,sd=1)^3}

gen20 <- function(n) {
            1/rexp(n)^2
}

gen21 <- function(n) { outvek<-rep(0,n)                                  #BD 21 (Marronite)
            component<-sample(2,size=n,replace=T,prob=c(1/3,2/3))
            outvek[component==1] <- rnorm(sum(component==1),mean=-20,sd=0.25 )
            outvek[component==2] <- rnorm(sum(component==2),mean=0,sd=1 )
            outvek
}
gen22 <- function(n) { outvek<-rep(0,n)                                  #BD 22 (Skewed Bimodal)
            component<-sample(2,size=n,replace=T,prob=c(3/4,1/4))
            outvek[component==1] <- rnorm(sum(component==1),mean=0,sd=1)
            outvek[component==2] <- rnorm(sum(component==2),mean=1.5,sd=1/3)
            outvek
}

gen23 <- function(n) { outvek<-rep(0,n)                                  #BD 23 (Claw)
            component<-sample(6,size=n,replace=T,prob=c(1/2,1/10,1/10,1/10,1/10,1/10))
            outvek[component==1] <-rnorm(sum(component==1),mean=0,sd=1)
            outvek[component==2] <-rnorm(sum(component==2),mean=-1,sd=0.1)
            outvek[component==3] <-rnorm(sum(component==3),mean=-0.5,sd=0.1)
            outvek[component==4] <-rnorm(sum(component==4),mean=0,sd=0.1)
            outvek[component==5] <-rnorm(sum(component==5),mean=0.5,sd=0.1)
            outvek[component==6] <-rnorm(sum(component==6),mean=1,sd=0.1)
            outvek
}

gen24 <- function(n) { outvek<-rep(0,n)                                  #BD 24 (Smooth Comb)
            component<-sample(6,size=n,replace=T,prob=c(32/63,16/63,8/63,4/63,2/63,1/63))
            outvek[component==1] <-rnorm(sum(component==1),mean=-31/21,sd=32/63)
            outvek[component==2] <-rnorm(sum(component==2),mean=17/21,sd=16/63)
            outvek[component==3] <-rnorm(sum(component==3),mean=41/21,sd=8/63)
            outvek[component==4] <-rnorm(sum(component==4),mean=53/21,sd=4/63)
            outvek[component==5] <-rnorm(sum(component==5),mean=59/21,sd=2/63)
            outvek[component==6] <-rnorm(sum(component==6),mean=62/21,sd=1/63)
            outvek
}
gen25 <- function(n) { outvek<-rep(0,n)                                  #BD 25 (caliper)
		n2<-n
		outvek<-NULL
		while (n2 > 0) {
			x<-runif(n2)
			u<-runif(n2)
			ind<-(u<=1-x^(1/3))
			outvek<-c(outvek,x[ind==T])
			n2<-n2-sum(ind)
			}
                s<-sample(c(1,-1),replace=T,size=n)
                outvek<-s*(outvek+0.1)
            outvek
}
gen26 <- function(n) {outvek<-rep(0,n)                                   #BD 26 (trimodal uniform)
            component<-sample(3,size=n,replace=T,prob=c(0.5,0.25,0.25))
            outvek[component==1] <- runif(sum(component==1),min=-1, max=1 )
            outvek[component==2] <- runif(sum(component==2),min=20, max=20.1 )
            outvek[component==3] <- runif(sum(component==3),min=-20.1, max=-20 )
            outvek
}
gen27 <- function(n) {sample(c(-9,-7,-5,-3,-1,1,3,5,7,9),size=n,replace=T)+gen16(n)}    #BD27 (sawtooth)
gen28 <- function(n) {outvek<-runif(n,min=0,max=1)*runif(n,min=0,max=1)                  #BD28 (bilogarithmic peak)
            signvek<-sample(2,size=n,replace=T)         
            outvek[signvek==2]<-1-outvek[signvek==2]
            outvek
}
rberdev <- function(n=1,dnum=1) {
        if (is.nan(dnum) || ! dnum %in% 1:28)
            stop("dnum must be between 1 and 28")
        return (
            eval(               
                parse(text = sprintf("gen%02d(n)", dnum)) # evaluate "gen[dnum](n)"-string
            )
        )
        # gen01(n)   "uniform"                 
        # gen02(n)   "exponential"             
        # gen03(n)   "Maxwell"                 
        # gen04(n)   "double exponential"      
        # gen05(n)   "logistic"                
        # gen06(n)   "Cauchy"                  
        # gen07(n)   "extreme value"           
        # gen08(n)   "infinite peak"           
        # gen09(n)   "Pareto"                  
        # gen10(n)   "symmetric Pareto"       
        # gen11(n)   "normal"                 
        # gen12(n)   "lognormal"              
        # gen13(n)   "uniform scale mixture"  
        # gen14(n)   "Matterhorn"             
        # gen15(n)   "logarithmic peak"       
        # gen16(n)   "isosceles triangle"     
        # gen17(n)   "beta (2,2)"             
        # gen18(n)   "chi-square (1)"         
        # gen19(n)   "normal cubed"           
        # gen20(n)   "inverse exponential"    
        # gen21(n)   "Marronite"              
        # gen22(n)   "skewed bimodal"         
        # gen23(n)   "claw"                       
        # gen24(n)   "smooth comb"                   
        # gen25(n)   "caliper"                              
        # gen26(n)   "trimodal uniform"             
        # gen27(n)   "sawtooth"                     
        # gen28(n)   "bilogarithmic peak"                  
}

