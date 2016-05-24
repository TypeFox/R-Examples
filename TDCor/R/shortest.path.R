shortest.path <-
function(bootstrap,BS_thr)

{ 

 dim_net=dim(bootstrap)[1]



# Transform the bootstrap matrix (weighted, oriented and "typed" (+ or -) edges to a simple network with oriented edges)

# The matrix of shortest paths will be built by iterative modification of sp.



 sp=sign(floor(abs(bootstrap/BS_thr)))*sign(bootstrap/BS_thr)



# Extract the bootstrap value and transform into log (help calculate the geometric mean)



 bb=abs(bootstrap)

 bb[bb==0]=1

 bb=log(bb)



### Iterative building of three matrices

# sp contains the shortest path (number of edges) linking each pair of gene

# bb contains the geometric mean of the bootstrap values of the single edges of the best shortest pathway

# mi contains the value of the index of directness of the edge showing the highest index in the best shortest pathway



 bool=T ; n=1 ; bb0=bb ;sp0=sp ; sp1=matrix(0,dim_net,dim_net)



while (prod(sp1==sp)!=1)

{  sp1=sp ; n=n+1 ; bb1=bb;



   for (j in 1:dim_net)

  {  



# Adding the targets of j of rank n to the list of targets of j.

# The bootstrap value of the (gene j -> target of rank n) interaction at that stage 

# is set as the maximum sum of the bootstrap of the (gene j -> target of rank n-1) interaction 

# and the bootstrap of the (target of rank n-1 -> target of rank n) interaction.



b=t(matrix(rep(bb1[,j],dim_net),dim_net,dim_net))

      bb[,j]=bb1[,j]+(1-sign(bb1[,j]))*apply((b+bb0)*(sign(b)*sign(bb0)),1,max)



# Building a matrix -mind- that contains 1 at each position in the matrix 

# corresponding to a new node added downstream of j with best bootstrap edge and 0 everywhere else.



mind0= t(matrix(rep(1:dim_net,dim_net),dim_net,dim_net))

mind1= matrix(rep(apply((b+bb0)*(sign(b)*sign(bb0)),1,which.max),dim_net),dim_net,dim_net)

mind= (1-sign(abs(mind1- mind0)))*rep(apply(sign(b)*sign(bb0),1,max))



# update of the rank n targets of gene j

      sp[,j]=sp1[,j]+(1-sign(abs(sp1[,j])))*n*sign(apply(t(matrix(rep(sp1[,j],dim_net),dim_net,dim_net))*sp0,1,sum))

  }

}



# Getting the geometric mean of the best bootstrap values

bb2=bb ; sp2=abs(sp) ; sp2[sp2==0]=1

bb=round(exp(bb2/abs(sp2))*sign(abs(sp)))



return(list(SP=sp,BS=bb))}
