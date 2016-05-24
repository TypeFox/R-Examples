

.getpval=function(k,null){

 i=1
 maxn=length(null)
 maxv=null[maxn]

   if(k<=maxv){
     
     while(k>null[i] & i<maxn){
       i=i+i
     }

     if(i>maxn){
       r=maxn
       l=i/2
     } else{
       r=i
       l=i/2
     }

     m=.binsearch(k,l,r,null)
     pval=1-(m/maxn)
   } else{
     pval=1/maxn   
   }
 
return(pval)
}

