nancorr <-
function(x,meth){

    nx=nrow(x)
    np=ncol(x)
      nancols=which(sum(is.nan(x))==nx)   
      goodcols=setdiff(seq(1,np),nancols)
    C=NaN*matrix(1,np,np)
    N=NaN*matrix(1,np,np)
      ngood=length(goodcols) 
      Cgood=matrix(0,ngood,ngood)
      Ngood=matrix(0,ngood,ngood)

      if (meth==1)
      { for (j in 1:ngood)  {                                            
            c1=goodcols[j]
            x1=x[,c1]                               
            for (k in j:ngood)
                { c2=goodcols[k]
                  x2=x[,c2]                                    
                  prod=x1*x2  
                                                    
                  goodprod=!is.nan(prod)   
                  N[c1,c2]=sum(goodprod)
                  Ngood[j,k]=N[c1,c2]             
                  C[c1,c2]=sum(prod[goodprod])/nx         
                  Cgood[j,k]=C[c1,c2]
                  C[c2,c1]=C[c1,c2] }
           }
       } else {
        goodrows=which(sum(c(is.nan(t(x[,goodcols])),1))==0)
        Cgood=cov(x[goodrows,goodcols])                  
        C[goodcols,goodcols]=Cgood                      
        Ngood=length(goodrows)*matrix(1,ngood,ngood)
        N[goodcols,goodcols]=Ngood
       }
    result = list(C=C)
    class(result) = "nancorr"
    result
}

