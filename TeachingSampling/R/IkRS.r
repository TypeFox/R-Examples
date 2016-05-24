IkRS <- function(N){
sam <- matrix(0, ncol=N, nrow=1)
for(k in 1:N){
sam<-rbind(sam, Ik(N,k))
}
sam
}