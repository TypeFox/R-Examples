#functions for calculating Jm, recursively, up to m=16
  J2.fun<- function(Psi.d2,list.J,list.Psi){
    	2*Psi.d2*list.J[[2]]-list.Psi[[1]]^2*list.J[[1]]
    	} 
    	
  J3.fun<- function(Psi.d3,list.J,list.Psi){
    	coef <- sapply(seq(1,3,by=1),function(x){(-1)^(x-1)*choose(3,3-x)})
    	coef[1]*Psi.d3*list.J[[3]]+coef[2]*list.Psi[[2]]^2*list.J[[2]]+coef[3]*list.Psi[[1]]^3*list.J[[1]] 
    	}
    	
  J4.fun<- function(Psi.d4,list.J,list.Psi){
    	coef <- sapply(seq(1,4,by=1),function(x){(-1)^(x-1)*choose(4,4-x)})
    	coef[1]*Psi.d4*list.J[[4]]+coef[2]*list.Psi[[3]]^2*list.J[[3]]+coef[3]*list.Psi[[2]]^3*list.J[[2]]+coef[4]*list.Psi[[1]]^4*list.J[[1]]
    	}
    	
  J5.fun<- function(Psi.d5,list.J,list.Psi){
    	coef <- sapply(seq(1,5,by=1),function(x){(-1)^(x-1)*choose(5,5-x)})
    	coef[1]*Psi.d5*list.J[[5]]+ coef[2]*list.Psi[[4]]^2*list.J[[4]]+coef[3]*list.Psi[[3]]^3*list.J[[3]]+coef[4]*list.Psi[[2]]^4*list.J[[2]]+coef[5]*list.Psi[[1]]^5*list.J[[1]]
    	}
  J6.fun<- function(Psi.d6,list.J,list.Psi){
    	coef <- sapply(seq(1,6,by=1),function(x){(-1)^(x-1)*choose(6,6-x)})
    	coef[1]*Psi.d6*list.J[[6]]+ coef[2]*list.Psi[[5]]^2*list.J[[5]]+coef[3]*list.Psi[[4]]^3*list.J[[4]]+coef[4]*list.Psi[[3]]^4*list.J[[3]]+coef[5]*list.Psi[[2]]^5*list.J[[2]]+ coef[6]*list.Psi[[1]]^6*list.J[[1]]
    	}	
  J7.fun<- function(Psi.d7,list.J,list.Psi){
    	coef <- sapply(seq(1,7,by=1),function(x){(-1)^(x-1)*choose(7,7-x)})
    	coef[1]*Psi.d7*list.J[[7]]+ coef[2]*list.Psi[[6]]^2*list.J[[6]]+coef[3]*list.Psi[[5]]^3*list.J[[5]]+coef[4]*list.Psi[[4]]^4*list.J[[4]]+coef[5]*list.Psi[[3]]^5*list.J[[3]]+ coef[6]*list.Psi[[2]]^6*list.J[[2]]+coef[7]*list.Psi[[1]]^7*list.J[[1]]
    	}
  J8.fun<- function(Psi.d8,list.J,list.Psi){
    	coef <- sapply(seq(1,8,by=1),function(x){(-1)^(x-1)*choose(8,8-x)})
    	coef[1]*Psi.d8*list.J[[8]]+ coef[2]*list.Psi[[7]]^2*list.J[[7]]+coef[3]*list.Psi[[6]]^3*list.J[[6]]+coef[4]*list.Psi[[5]]^4*list.J[[5]]+coef[5]*list.Psi[[4]]^5*list.J[[4]]+ coef[6]*list.Psi[[3]]^6*list.J[[3]]+coef[7]*list.Psi[[2]]^7*list.J[[2]]+coef[8]*list.Psi[[1]]^8*list.J[[1]]
    	}
  J9.fun<- function(Psi.d9,list.J,list.Psi){
    	coef <- sapply(seq(1,9,by=1),function(x){(-1)^(x-1)*choose(9,9-x)})
    	coef[1]*Psi.d9*list.J[[9]]+ coef[2]*list.Psi[[8]]^2*list.J[[8]]+coef[3]*list.Psi[[7]]^3*list.J[[7]]+coef[4]*list.Psi[[6]]^4*list.J[[6]]+coef[5]*list.Psi[[5]]^5*list.J[[5]]+ coef[6]*list.Psi[[4]]^6*list.J[[4]]+coef[7]*list.Psi[[3]]^7*list.J[[3]]+coef[8]*list.Psi[[2]]^8*list.J[[2]]+coef[9]*list.Psi[[1]]^9*list.J[[1]]
    	}
    	
  J10.fun<- function(Psi.d10,list.J,list.Psi){
    	coef <- sapply(seq(1,10,by=1),function(x){(-1)^(x-1)*choose(10,10-x)})
    	coef[1]*Psi.d10*list.J[[10]]+coef[2]*list.Psi[[9]]^2*list.J[[9]]+coef[3]*list.Psi[[8]]^3*list.J[[8]]+coef[4]*list.Psi[[7]]^4*list.J[[7]]+coef[5]*list.Psi[[6]]^5*list.J[[6]]+coef[6]*list.Psi[[5]]^6*list.J[[5]]+coef[7]*list.Psi[[4]]^7*list.J[[4]]+coef[8]*list.Psi[[3]]^8*list.J[[3]]+coef[9]*list.Psi[[2]]^9*list.J[[2]]+coef[10]*list.Psi[[1]]^10*list.J[[1]]
    	}
  
  J11.fun<- function(Psi.d11,list.J,list.Psi){
    	coef <- sapply(seq(1,11,by=1),function(x){(-1)^(x-1)*choose(11,11-x)})
    	coef[1]*Psi.d11*list.J[[11]]+coef[2]*list.Psi[[10]]^2*list.J[[10]]+coef[3]*list.Psi[[9]]^3*list.J[[9]]+coef[4]*list.Psi[[8]]^4*list.J[[8]]+coef[5]*list.Psi[[7]]^5*list.J[[7]]+coef[6]*list.Psi[[6]]^6*list.J[[6]]+coef[7]*list.Psi[[5]]^7*list.J[[5]]+coef[8]*list.Psi[[4]]^8*list.J[[4]]+coef[9]*list.Psi[[3]]^9*list.J[[3]]+coef[10]*list.Psi[[2]]^10*list.J[[2]]+coef[11]*list.Psi[[1]]^11*list.J[[1]]
    	}
    	
  J12.fun<- function(Psi.d12,list.J,list.Psi){
    	coef <- sapply(seq(1,12,by=1),function(x){(-1)^(x-1)*choose(12,12-x)})
    	coef[1]*Psi.d12*list.J[[12]]+coef[2]*list.Psi[[11]]^2*list.J[[11]]+coef[3]*list.Psi[[10]]^3*list.J[[10]]+coef[4]*list.Psi[[9]]^4*list.J[[9]]+coef[5]*list.Psi[[8]]^5*list.J[[8]]+coef[6]*list.Psi[[7]]^6*list.J[[7]]+coef[7]*list.Psi[[6]]^7*list.J[[6]]+coef[8]*list.Psi[[5]]^8*list.J[[5]]+coef[9]*list.Psi[[4]]^9*list.J[[4]]+coef[10]*list.Psi[[3]]^10*list.J[[3]]+coef[11]*list.Psi[[2]]^11*list.J[[2]]+coef[12]*list.Psi[[1]]^12*list.J[[1]]
    	}
    	
  J13.fun<- function(Psi.d13,list.J,list.Psi){
    	coef <- sapply(seq(1,13,by=1),function(x){(-1)^(x-1)*choose(13,13-x)})
    	coef[1]*Psi.d13*list.J[[13]]+coef[2]*list.Psi[[12]]^2*list.J[[12]]+coef[3]*list.Psi[[11]]^3*list.J[[11]]+coef[4]*list.Psi[[10]]^4*list.J[[10]]+coef[5]*list.Psi[[9]]^5*list.J[[9]]+coef[6]*list.Psi[[8]]^6*list.J[[8]]+coef[7]*list.Psi[[7]]^7*list.J[[7]] +coef[8]*list.Psi[[6]]^8*list.J[[6]]+coef[9]*list.Psi[[5]]^9*list.J[[5]]+coef[10]*list.Psi[[4]]^10*list.J[[4]]+coef[11]*list.Psi[[3]]^11*list.J[[3]]+coef[12]*list.Psi[[2]]^12*list.J[[2]]+coef[13]*list.Psi[[1]]^13*list.J[[1]]
    	}
    	
  J14.fun<- function(Psi.d14,list.J,list.Psi){
    	coef <- sapply(seq(1,14,by=1),function(x){(-1)^(x-1)*choose(14,14-x)})
    	coef[1]*Psi.d14*list.J[[14]]+ coef[2]*list.Psi[[13]]^2*list.J[[13]]+coef[3]*list.Psi[[12]]^3*list.J[[12]]+coef[4]*list.Psi[[11]]^4*list.J[[11]]+coef[5]*list.Psi[[10]]^5*list.J[[10]]+ coef[6]*list.Psi[[9]]^6*list.J[[9]]+coef[7]*list.Psi[[8]]^7*list.J[[8]]+coef[8]*list.Psi[[7]]^8*list.J[[7]] + coef[9]*list.Psi[[6]]^9*list.J[[6]]+coef[10]*list.Psi[[5]]^10*list.J[[5]] + coef[11]*list.Psi[[4]]^11*list.J[[4]] + coef[12]*list.Psi[[3]]^12*list.J[[3]]+coef[13]*list.Psi[[2]]^13*list.J[[2]]+coef[14]*list.Psi[[1]]^14*list.J[[1]]
    	}
    	
  J15.fun<- function(Psi.d15,list.J,list.Psi){
    	coef <- sapply(seq(1,15,by=1),function(x){(-1)^(x-1)*choose(15,15-x)})
    	coef[1]*Psi.d15*list.J[[15]]+coef[2]*list.Psi[[14]]^2*list.J[[14]]+coef[3]*list.Psi[[13]]^3*list.J[[13]]+coef[4]*list.Psi[[12]]^4*list.J[[12]]+coef[5]*list.Psi[[11]]^5*list.J[[11]]+coef[6]*list.Psi[[10]]^6*list.J[[10]]+coef[7]*list.Psi[[9]]^7*list.J[[9]]+coef[8]*list.Psi[[8]]^8*list.J[[8]]+coef[9]*list.Psi[[7]]^9*list.J[[7]]+coef[10]*list.Psi[[6]]^10*list.J[[6]]+coef[11]*list.Psi[[5]]^11*list.J[[5]]+coef[12]*list.Psi[[4]]^12*list.J[[4]]+coef[13]*list.Psi[[3]]^13*list.J[[3]]+coef[14]*list.Psi[[2]]^14*list.J[[2]]+coef[15]*list.Psi[[1]]^15*list.J[[1]]
    	}
    	
  J16.fun<- function(Psi.d16,list.J,list.Psi){
    	coef <- sapply(seq(1,16,by=1),function(x){(-1)^(x-1)*choose(16,16-x)})
    	coef[1]*Psi.d16*list.J[[16]]+coef[2]*list.Psi[[15]]^2*list.J[[15]]+coef[3]*list.Psi[[14]]^3*list.J[[14]]+coef[4]*list.Psi[[13]]^4*list.J[[13]]+coef[5]*list.Psi[[12]]^5*list.J[[12]]+coef[6]*list.Psi[[11]]^6*list.J[[11]]+coef[7]*list.Psi[[10]]^7*list.J[[10]]+coef[8]*list.Psi[[9]]^8*list.J[[9]]+coef[9]*list.Psi[[8]]^9*list.J[[8]]+coef[10]*list.Psi[[7]]^10*list.J[[7]]+coef[11]*list.Psi[[6]]^11*list.J[[6]]+coef[12]*list.Psi[[5]]^12*list.J[[5]]+coef[13]*list.Psi[[4]]^13*list.J[[4]]+coef[14]*list.Psi[[3]]^14*list.J[[3]]+coef[15]*list.Psi[[2]]^15*list.J[[2]]+coef[16]*list.Psi[[1]]^16*list.J[[1]]
    	}    	
    	
    	
   
