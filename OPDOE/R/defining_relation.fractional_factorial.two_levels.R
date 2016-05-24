defining_relation.fractional_factorial.two_levels=function(M)
{
 l=dim(M)[1]
 if(is.null(l))
   l=1
 M=as.matrix(M,nrow=1,ncol=length(M))
 char=function(k)
    {
   a=rep("0",length(k))
   for (i in (1:length(k)))
      {if (k[i]==1) a[i]=LETTERS[i]}
   return(a)
    }
 z=rep("",l)
 b=vector("list",l)
 for (i in 1:l)
   {
    b[[i]]=grep("[A-Z]",char(M[i,]));
    z[i]=toString(char(M[i,])[b[[i]]]);
    z[i]=gsub(", ","",z[i])
   }
 return(z)
}

# finds minimal complete sets of confounded effects
# from a list of defining relations (all the products
# of given generators) 
