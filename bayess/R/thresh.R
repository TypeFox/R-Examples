thresh=function(k,n1,c2,c3,r2,q1){
# ACCEPT-REJECT BOUND

choose(n1-c2,k)*0.9^k*choose(n1-k,c3+1)/(choose(n1,k)*choose(n1,c3+1))
}
