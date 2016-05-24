simukern<-function(n,d,seed,newvalue,index,delta,minim,h){
#
set.seed(seed)
dendat<-matrix(0,n,d)
volofatom<-prod(delta)
newvalue<-newvalue/sum(volofatom*newvalue)
#
for (i in 1:n){
   uni<-runif(1)
   cumyla<-newvalue[1]*volofatom
   run<-1
   while (cumyla<uni){
      run<-run+1 
      cumyla<-cumyla+newvalue[run]*volofatom
   }
   #run leads to the right bin
   inde<-index[run,]
   uni2<-runif(1)                #add blur in bin
   obse<-minim-h+delta*inde+delta*(uni2-1/2)
   dendat[i,]<-obse
}
return(dendat)
}
