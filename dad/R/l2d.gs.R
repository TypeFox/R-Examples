l2d.gs <-
function(x1, x2, check=FALSE)
{  
   var1<-var(x1);
   var2<-var(x2);
   if(check)
    {if(abs(det(var1))<.Machine$double.eps | abs(det(var2))<.Machine$double.eps )
      {stop("One of the sample variances is degenerate")
      }
    }  
   p<-ncol(x1);
   m1<-colMeans(x1);
   m2<-colMeans(x2);
   return(as.numeric(l2d.gp(m1,var1,m2,var2)))
}

