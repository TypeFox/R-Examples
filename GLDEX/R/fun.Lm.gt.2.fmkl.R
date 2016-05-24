`fun.Lm.gt.2.fmkl` <-
function(L3, L4, r) {

k<-0:(r-1)

if(L3==0 & L4==0){

result<-sum((-1)^(r-k-1)*choose(r-1,k)*choose(r+k-1,
k)*(-1/(k^2+2*k+1)-((-psigamma(k+2,0)+psigamma(1,0))/(k+1))))

return(result)

}

if(L3!=0 & L4==0){

result<-sum((-1)^(r-k-1)*choose(r-1,k)*choose(r+k-1,
k)*(1/L3*(1/(L3+k+1)-1/(k+1))-((-psigamma(k+2,0)+psigamma(1,0))/(k+1))))

return(result)

}

if(L3==0 & L4!=0){

result<-sum((-1)^(r-k-1)*choose(r-1,k)*choose(r+k-1,
k)*(-1/(k^2+2*k+1)-(((k+1)*fun.beta(k+1,L4+1)-1)/((k+1)*L4))))

return(result)

}

if(L3!=0 & L4!=0){

result<-sum((-1)^(r-k-1)*choose(r-1,k)*choose(r+k-1,
k)*(1/L3*(1/(L3+k+1)-1/(k+1))-(((k+1)*fun.beta(k+1,L4+1)-1)/((k+1)*L4))))

return(result)

}
}

