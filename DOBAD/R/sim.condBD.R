#Charlesdoss
##### The functions here are related to simulating from the distribution
##### of the BD process conditional on the endpoints.

sim.condBD <- function(bd.PO=list(states=c(5,7,3), times=c(0,.4,1)), N=1, 
                       L=.5, m=.7, nu=.4, n.fft=1024, prevSims=NULL){
  return(c(prevSims,
           replicate(n=N, sim.condBD.main(bd.PO=bd.PO, L=L,m=m,nu=nu, n.fft=n.fft),
            simplify=FALSE)));
}
