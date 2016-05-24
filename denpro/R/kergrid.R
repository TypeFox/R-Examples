kergrid<-function(dendat,h,N){
#
#dendat is n*d- matrix of observations, 
#h is vector of positive smoothing parameters
#N is d-vector of the (dyadic) number of grid points for each direction
#
#dendat<-matrix(rnorm(20),10)
#h<-c(0.8,1,1.2)
#N<-c(4,4)
#
hnum<-length(h)
n<-dim(dendat)[1]
d<-dim(dendat)[2]
depth<-log(N,base=2)   
depoftree<-sum(depth)+1
#
minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])  
  maxim[i]<-max(dendat[,i])
}
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)
#
mindelta<-min(delta)
maxpositive<-n*(2*hmax/mindelta)^d
bigd<-sum(log(N,base=2))
maxnode<-bigd*ceiling(maxpositive)
#
numnode<-1
left<-matrix(0,maxnode,1)
right<-matrix(0,maxnode,1)
parent<-matrix(0,maxnode,1)
infopointer<-matrix(0,maxnode,1)
low<-matrix(0,maxnode,1)
low[1]<-1
upp<-matrix(0,maxnode,1)
upp[1]<-N[1]
#
numpositive<-0
value<-matrix(0,maxpositive,hnum)
index<-matrix(0,maxpositive,d)
nodefinder<-matrix(0,maxpositive,1)
#
gridlow<-matrix(0,d,1)
gridupp<-matrix(0,d,1)
#
for (i in 1:n){
 for (hrun in 1:hnum){ 
   #find the grid points in the support 
   for (j in 1:d){  
      gridlow[j]<-floor(((dendat[i,j]-minim[j])/delta[j])+1)
      gridupp[j]<-ceiling(((dendat[i,j]-minim[j]+2*h[hrun])/delta[j])-1)
   }
   base<-gridupp-gridlow+1
   gridcard<-prod(base)
   k<-0
   while (k<=(gridcard-1)){
      if (d>1){  
          inde<-digit(k,base)   #inde is d-vector
          inde<-inde+gridlow
      }
      else{
          inde<-gridlow+k
      }
      point<-minim-h[hrun]+delta*inde     #point is d-vector  
      val<-epane(point-dendat[i,],h[hrun])
      #find whether gridpoint is already in tree
      fe<-findend(inde,left,right,low,upp,N)
      if (fe$exists){
           pointer<-infopointer[fe$location]
           curval<-value[pointer,hrun]
           value[pointer,hrun]<-curval+val/n
      }
      else{  #gridpoint was not yet in the tree
         curre<-fe$location
         curdep<-fe$dep
         #
         ad<-addnode(inde,curre,curdep,left,right,parent,low,upp,N,numnode)
         numnode<-ad$numnode
         left<-ad$left
         right<-ad$right
         parent<-ad$parent
         low<-ad$low
         upp<-ad$upp
         nodeloc<-ad$nodeloc
         #
         numpositive<-numpositive+1
         infopointer[numnode]<-numpositive
         value[numpositive,hrun]<-val/n
         index[numpositive,]<-inde
         nodefinder[numpositive]<-nodeloc
      }
      k<-k+1 
   }
 }
}
left<-left[1:numnode]
right<-right[1:numnode]
parent<-parent[1:numnode]
infopointer<-infopointer[1:numnode]
#deplink<-deplink[1:numnode]
low<-low[1:numnode]
upp<-upp[1:numnode]
#
value<-value[1:numpositive,]
index<-index[1:numpositive,]
nodefinder<-nodefinder[1:numpositive]
return(list(left=left,right=right,parent=parent,infopointer=infopointer,
low=low,upp=upp,value=value,index=index,nodefinder=nodefinder))
}                              








