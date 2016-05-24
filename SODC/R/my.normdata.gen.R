my.normdata.gen <-
function(u,p,nobs,nvar){

   vec25=as.vector(rep(1,nvar/2))
   vec50=as.vector(rep(1,nvar))
   idenmat50=diag(nvar)



   mean.y1=c(-u*vec25,u*vec25) 
   mean.y2=c(u*vec50)
   mean.y3=c(u*vec25,-u*vec25)
   #mean.y4=c(-u*vec50)



   x=matrix(0, nrow=nobs, ncol=p)

   
   y = sample(c(1,2,3), nobs, replace = TRUE)

   for(i in 1:nobs){

        
         if(y[i]== 1)
         x[i,1:nvar]= mvrnorm(n=1, mean.y1, idenmat50)
         else if(y[i]== 2)
         x[i,1:nvar]= mvrnorm(n=1, mean.y2, idenmat50)
         else if(y[i]== 3)
         x[i,1:nvar]= mvrnorm(n=1, mean.y3, idenmat50)
        
         if(p>nvar)
            x[i,(nvar+1):p] = rnorm(n=(p-nvar), mean = 0, sd = 1)

      
    }
    return (list(x=x,y=y))
    

}
