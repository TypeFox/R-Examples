conf.matrix=function(v,p)
{l=length(v)
 a=matrix(0,l,p)
 for (i in 1:l) 
     {for (j in 1:p) 
            {if (gregexpr(LETTERS[j],v[[i]])[[1]][1]>0) a[i,j]=1
            }
     }
 colnames(a)=LETTERS[1:p]
 return(a)
}
 
# converts the list of defining relations into
# a matrix ready for input to "conf.design"
