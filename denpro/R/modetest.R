modetest<-function(pk,pknum,
h=NULL,N=NULL,Q=NULL,bootnum=NULL,delta=NULL,nsimu=NULL,minim=NULL,
type="boots",kernel="gauss",
n=NULL)
{

#pk is a list of level set trees
#h is vector of smoothing parameter values
#M is the number of bootstrap samples to be generated

run<-1
while (run<=pknum){
   curlst<-pk[[run]]

   if (type=="boots"){
       curh<-h[run]
       curmotes<-modetestydin(curlst,curh,N,Q,bootnum,delta,nsimu,minim,kernel)
   }
   else{
       curmotes<-modetestgauss(curlst,n)
   }

   if (run==1){
      if (pknum==1){
          moteslist<-curmotes
      }
      else{
          moteslist=list(curmotes)
      }
   }
   else{
      moteslist=c(moteslist,list(curmotes))
   }
   run<-run+1
}
#
return(moteslist)
}



