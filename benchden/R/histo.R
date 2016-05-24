rhisto<-function(n,dnum=1) {
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")

	outvek<-rep(0,n)                                  
	if (dnum==1) {
	component<-sample(5,size=n,replace=T,prob=c(0.15,0.35,0.2,0.1,0.2))
		outvek[component==1] <- runif(sum(component==1),min=0.0, max=0.2 )
            outvek[component==2] <- runif(sum(component==2),min=0.2, max=0.4 )
            outvek[component==3] <- runif(sum(component==3),min=0.4, max=0.6 )
            outvek[component==4] <- runif(sum(component==4),min=0.6, max=0.8 )
            outvek[component==5] <- runif(sum(component==5),min=0.8, max=1.0 )
          	}
	if (dnum==2) {
	component<-sample(5,size=n,replace=T,prob=c(0.15,0.35,0.2,0.1,0.2))
		outvek[component==1] <- runif(sum(component==1),min=0.0, max=0.13)
            outvek[component==2] <- runif(sum(component==2),min=0.13, max=0.34 )
            outvek[component==3] <- runif(sum(component==3),min=0.34, max=0.61 )
            outvek[component==4] <- runif(sum(component==4),min=0.61, max=0.65 )
            outvek[component==5] <- runif(sum(component==5),min=0.65, max=1.0 )
          	}
	if (dnum==3) {
	component<-sample(10,size=n,replace=T,prob=c(0.01, 0.18, 0.16, 0.07, 0.06, 0.01, 0.06 ,0.37, 0.06, 0.02))
		outvek[component==1] <- runif(sum(component==1),min=0.0, max=0.1 )
            outvek[component==2] <- runif(sum(component==2),min=0.1, max=0.2 )
            outvek[component==3] <- runif(sum(component==3),min=0.2, max=0.3 )
            outvek[component==4] <- runif(sum(component==4),min=0.3, max=0.4 )
            outvek[component==5] <- runif(sum(component==5),min=0.4, max=0.5 )
          	outvek[component==6] <- runif(sum(component==6),min=0.5, max=0.6 )
            outvek[component==7] <- runif(sum(component==7),min=0.6, max=0.7 )
            outvek[component==8] <- runif(sum(component==8),min=0.7, max=0.8 )
            outvek[component==9] <- runif(sum(component==9),min=0.8, max=0.9 )
            outvek[component==10] <- runif(sum(component==10),min=0.9, max=1.0 )
          	}
	if (dnum==4) {
	component<-sample(10,size=n,replace=T,prob=c(0.01, 0.18, 0.16, 0.07, 0.06, 0.01, 0.06 ,0.37, 0.06, 0.02))
		outvek[component==1] <- runif(sum(component==1),min=0.0, max=0.02 )
            outvek[component==2] <- runif(sum(component==2),min=0.02, max=0.07 )
            outvek[component==3] <- runif(sum(component==3),min=0.07, max=0.14 )
            outvek[component==4] <- runif(sum(component==4),min=0.14, max=0.44 )
            outvek[component==5] <- runif(sum(component==5),min=0.44, max=0.53 )
          	outvek[component==6] <- runif(sum(component==6),min=0.53, max=0.56 )
            outvek[component==7] <- runif(sum(component==7),min=0.56, max=0.67 )
            outvek[component==8] <- runif(sum(component==8),min=0.67, max=0.77 )
            outvek[component==9] <- runif(sum(component==9),min=0.77, max=0.91 )
            outvek[component==10] <- runif(sum(component==10),min=0.91, max=1.0 )
          	}
	outvek
}

dhisto<-function(x,dnum=1) {
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")

	outvek<-rep(0,length(x))
	if (dnum==1) {
		outvek[(x>=0)&(x<=0.2)] <- (0.15/0.2)
		outvek[(x>0.2)&(x<=0.4)] <- (0.35/0.2)
		outvek[(x>0.4)&(x<=0.6)] <- (0.2/0.2)
		outvek[(x>0.6)&(x<=0.8)] <- (0.1/0.2)
		outvek[(x>0.8)&(x<=1)] <- (0.2/0.2)
	}
	if (dnum==2) {
		outvek[(x>=0)&(x<=0.13)] <- (0.15/0.13)
		outvek[(x>0.13)&(x<=0.34)] <- (0.35/0.21)
		outvek[(x>0.34)&(x<=0.61)] <- (0.2/0.27)
		outvek[(x>0.61)&(x<=0.65)] <- (0.1/0.04)
		outvek[(x>0.65)&(x<=1)] <- (0.2/0.35)
	}
	if (dnum==3) {
		outvek[(x>=0)&(x<=0.1)] <- (0.01/0.1)
		outvek[(x>0.1)&(x<=0.2)] <- (0.18/0.1)
		outvek[(x>0.2)&(x<=0.3)] <- (0.16/0.1)
		outvek[(x>0.3)&(x<=0.4)] <- (0.07/0.1)
		outvek[(x>0.4)&(x<=0.5)] <- (0.06/0.1)
		outvek[(x>0.5)&(x<=0.6)] <- (0.01/0.1)
		outvek[(x>0.6)&(x<=0.7)] <- (0.06/0.1)
		outvek[(x>0.7)&(x<=0.8)] <- (0.37/0.1)
		outvek[(x>0.8)&(x<=0.9)] <- (0.06/0.1)
		outvek[(x>0.9)&(x<=1.0)] <- (0.02/0.1)
	}
	if (dnum==4) {
		outvek[(x>=0)&(x<=0.02)] <- (0.01/0.02)
		outvek[(x>0.02)&(x<=0.07)] <- (0.18/0.05)
		outvek[(x>0.07)&(x<=0.14)] <- (0.16/0.07)
		outvek[(x>0.14)&(x<=0.44)] <- (0.07/0.3)
		outvek[(x>0.44)&(x<=0.53)] <- (0.06/0.09)
		outvek[(x>0.53)&(x<=0.56)] <- (0.01/0.03)
		outvek[(x>0.56)&(x<=0.67)] <- (0.06/0.11)
		outvek[(x>0.67)&(x<=0.77)] <- (0.37/0.10)
		outvek[(x>0.77)&(x<=0.91)] <- (0.06/0.14)
		outvek[(x>0.91)&(x<=1.0)] <- (0.02/0.09)
	}
	outvek
}

phisto<-function(q,dnum=1) {
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")
	if (dnum==1) out<-0.15*punif(q,min=0,max=0.2)+0.35*punif(q,min=0.2,max=0.4)+0.2*punif(q,min=0.4,max=0.6)+
		0.1*punif(q,min=0.6,max=0.8)+0.2*punif(q,min=0.8,max=1)
	if (dnum==2) out<-0.15*punif(q,min=0,max=0.13)+0.35*punif(q,min=0.13,max=0.34)+0.2*punif(q,min=0.34,max=0.61)+
		0.1*punif(q,min=0.61,max=0.65)+0.2*punif(q,min=0.65,max=1)
	if (dnum==3) out<-0.01*punif(q,min=0,max=0.1)+0.18*punif(q,min=0.1,max=0.2)+0.16*punif(q,min=0.2,max=0.3)+
		0.07*punif(q,min=0.3,max=0.4)+0.06*punif(q,min=0.4,max=0.5)+0.01*punif(q,min=0.5,max=0.6)+
		0.06*punif(q,min=0.6,max=0.7)+0.37*punif(q,min=0.7,max=0.8)+0.06*punif(q,min=0.8,max=0.9)+
		0.02*punif(q,min=0.9,max=1.0)
	if (dnum==4) out<-0.01*punif(q,min=0,max=0.02)+0.18*punif(q,min=0.02,max=0.07)+0.16*punif(q,min=0.07,max=0.14)+
		0.07*punif(q,min=0.14,max=0.44)+0.06*punif(q,min=0.44,max=0.53)+0.01*punif(q,min=0.53,max=0.56)+
		0.06*punif(q,min=0.56,max=0.67)+0.37*punif(q,min=0.67,max=0.77)+0.06*punif(q,min=0.77,max=0.91)+
		0.02*punif(q,min=0.91,max=1.0)
	out
}

qhisto<-function(p,dnum=1) {
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")
    p[p>1]<-1
    p[p<0]<-0
	if (dnum==1) out<-(p>=0)*(p<=0.15)*(p/0.75)+(p>0.15)*(p<=0.5)*((p-0.15)/1.75+0.2)+
    (p>0.5)*(p<=0.7)*((p-0.5)+0.4)+(p>0.7)*(p<=0.8)*(2*(p-0.7)+0.6)+(p>0.8)*(p<=1)*p
	if (dnum==2) out<-(p>=0)*(p<=0.15)*(p*13/15)+(p>0.15)*(p<=0.5)*((p-0.15)*0.6+0.13)+
    (p>0.5)*(p<=0.7)*((p-0.5)*1.35+0.34 )+(p>0.7)*(p<=0.8)*((p-0.7)*0.4+0.61)+(p>0.8)*(p<=1)*((p-0.8)*1.75+0.65)
	if (dnum==3) out<-(p>=0)*(p<=0.01)*(p*10)+(p>0.01)*(p<=0.19)*((p-0.01)/1.8+0.1)+
    (p>0.19)*(p<=0.35)*((p-0.19)*0.625+0.2)+(p>0.35)*(p<=0.42)*((p-0.35)/0.7+0.3)+ 
    (p>0.42)*(p<=0.48)*((p-0.42)/0.6+0.4)+(p>0.48)*(p<=0.49)*((p-0.48)*10+0.5)+
    (p>0.49)*(p<=0.55)*((p-0.49)/0.6+0.6)+(p>0.55)*(p<=0.92)*((p-0.55)/3.7+0.7)+
    (p>0.92)*(p<=0.98)*((p-0.92)/0.6+0.8)+(p>0.98)*(p<=1)*((p-0.98)*5+0.9)
	if (dnum==4) out<-(p>=0)*(p<=0.01)*(p*2)+(p>0.01)*(p<=0.19)*((p-0.01)/3.6+0.02)+
    (p>0.19)*(p<=0.35)*((p-0.19)*0.4375+0.07)+(p>0.35)*(p<=0.42)*((p-0.35)*30/7+0.14)+ 
    (p>0.42)*(p<=0.48)*((p-0.42)*1.5+0.44)+(p>0.48)*(p<=0.49)*((p-0.48)*3+0.53)+
    (p>0.49)*(p<=0.55)*((p-0.49)*11/6+0.56)+(p>0.55)*(p<=0.92)*((p-0.55)*10/37+0.67)+
    (p>0.92)*(p<=0.98)*((p-0.92)*7/3+0.77)+(p>0.98)*(p<=1)*((p-0.98)*4.5+0.91)
  out
}



nhisto<-function(dnum=1){
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")
	out<-NULL
	if (dnum==1) out<-"5 bin regular histogram"
	if (dnum==2) out<-"5 bin irregular histogram"
	if (dnum==3) out<-"10 bin regular histogram"
	if (dnum==4) out<-"10 bin irregular histogram"
	out
}

### breaks for histo testbeds

bhisto<-function(dnum=1) {
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")

	if (dnum==1) out<-c(0,0.2,0.4,0.6,0.8,1)
	if (dnum==2) out<-c(0,0.13,0.34,0.61,0.65,1)
	if (dnum==3) out<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
	if (dnum==4) out<-c(0,0.02,0.07,0.14,0.44,0.53,0.56,0.67,0.77,0.91,1)
  out
}

histo<-function(dnum=1) {
  if (is.nan(dnum) || ! dnum %in% 1:4)
    stop("dnum must be between 1 and 4")

	if (dnum==1) out<-list(
   name = nhisto(1)  ,
   peaks = c(0.3,0.9),
   support = matrix(c(0,1), nrow = 1),
   breaks = bhisto(1)
  )                     
	if (dnum==2) out<-list(
   name = nhisto(2)  ,
   peaks = c(0.235,0.63),
   support = matrix(c(0,1), nrow = 1),
   breaks = bhisto(2)
  )                     
	if (dnum==3) out<-list(
   name = nhisto(3)  ,
   peaks = c(0.15,0.75),
   support = matrix(c(0,1), nrow = 1),
   breaks = bhisto(3)
  )                     
	if (dnum==4) out<-list(
   name = nhisto(4)  ,
   peaks = c(0.045,0.485,0.72),
   support = matrix(c(0,1), nrow = 1),
   breaks = bhisto(4)
  )                     
  out

}

	
            
            