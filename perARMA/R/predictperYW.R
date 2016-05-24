predictperYW <-function(x,T,p,missval,start,...)
{
   predictperYW_full <- function(x, T, p, missval, start, end ,realcol, predcol)
   {    estimators <- perYW(x, T, p, missval) 	
        phi = estimators$phi			
        phi = as.matrix(phi)

        old=x  
        xp<-matrix(length(x),1,0)
        lim=length(x)+1

if (start<lim)                                
    {for(i in start:end)
      {  if (i%%T==0) {xp[i]=phi[T,]%*%x[(i-1):(i-p)]
                       x[i]=xp[i]    }
                       else
                       {xp[i]=phi[i%%T,]%*%x[(i-1):(i-p)]
                       x[i]=xp[i]    }
    }
    plot(old[start:end], type="l", col="red",xlab="time",ylab="series after removing periodic mean")  
    lines(x[start:end], type="l")
    legend("bottomright", c(expression(real), expression(prediction)), 
     fill = c(realcol, predcol), ncol = 2, title = "legend")
     title(main = "Prediction of the series", 
            sub = paste("forecast from", start,"to",end))
   } else {
   for(i in lim:start)
      { 
       if (i%%T==0) {xp[i]=phi[T,]%*%x[(i-1):(i-p)]
                     x=c(x,xp[i])}
                     else 
                 {   xp[i]=phi[i%%T,]%*%x[(i-1):(i-p)]
                     x=c(x,xp[i])
                  }
      }
   new=x[lim:start]
   plot(seq(1, (lim-1)),x[1:(lim-1)], type="l", col="red",xlab="time",ylab="series after removing periodic mean")
   lines(seq(lim, start),new, type="l", col="blue")
   legend("bottomright", c(expression(real), expression(prediction)), 
   fill = c(realcol, predcol), ncol = 2, title = "legend")
   title(main = "Prediction of the series", 
            sub = paste("forecast form",length(old),"to", start))
}
   result = list(x = x, new=new)
   class(result) = "predictperYW"
   result
  }
    L <- modifyList(list(end=end,realcol = "blue", predcol = "red"), list(x = x, T = T, p = p, missval=missval, start=start, ...))
    do.call(predictperYW_full, L)
}