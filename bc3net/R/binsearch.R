
.binsearch=function(k,l,r,null){
   m=round((l+r)/2)
   if(k==null[m]){     
      return(m)
   }else if(k>null[m] & k<null[m+1]){
      return(m)
   }else if(k<null[m]){         
      .binsearch(k,l,m-1,null)
   }else if(k>null[m]){
      .binsearch(k,m+1,r,null)
   }
}


