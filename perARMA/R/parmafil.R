parmafil <-
function(b,a,x){

      Ta=nrow(a)
      na=ncol(a)
      Tb=nrow(b)
      nb=ncol(b)
      nrx=nrow(x)
      ncx=ncol(x)


   if (ncx!=1)   {stop("input data has more than 1 column") }
   
   if (Ta!=Tb)   
           {stop("Dimension error: Ta!=Tb")
             } else {
           if (sum(a[,1]!=1))
                {for (j in 1:na)                      
                     { a[,j]=a[,j]/a[,1]}             
                }
        T=Ta
        xpad=c(matrix(0,nb-1,1),x)                                          
        yold=matrix(0,na-1,1)      
        ntimes=nrx                     
        y=999*matrix(1,nrx,1)

       for (i in 1:ntimes)
           {index=matlab::mod((i-1),T)+1
             xtmp=matlab::flipud(xpad[i:(i+nb-1)])
             xma=b[index,]%*%xtmp  

            if (na > 1)
            { ytmp=-a[index,2:na]%*%yold+xma  
              if (na>2)
                 {yold=rbind(ytmp, yold[1:(na-2),1])
                    } else {
                    yold=ytmp}
               }  else   {                           
            ytmp=xma}    
      y[i]=ytmp 
    }
}
      result = list(y=y) 
      class(result) = "parmafil"
      result
}