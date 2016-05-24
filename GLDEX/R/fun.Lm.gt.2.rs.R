`fun.Lm.gt.2.rs` <-
function(L3, L4, r) {

k<-0:(r-1)

result<-sum((-1)^(r-k-1)*choose(r-1,k)*choose(r+k-1,
k)*(1/(k+1+L3)+(-1)^r/(k+1+L4)))

return(result)

}

