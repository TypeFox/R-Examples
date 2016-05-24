postest<-function(data, y, pik, NG, description=FALSE)
{
if (missing(data) | missing(y) | missing(pik) | missing(NG)) 
    stop("incomplete input")
str <- function(st, h, n) .C("str", as.double(st), as.integer(h), 
        as.integer(n), s = double(n), PACKAGE = "sampling")$s
data=as.data.frame(data)
sample.size=nrow(data)
t_post=0
if(!is.null(colnames(data)))
{m = match("Stratum", colnames(data))
if(!any(is.na(m))) 
        {
          m = match("poststratum", colnames(data))
          if (any(is.na(m))) 
            stop("the column 'poststratum' is missing")
        h=unique(data$Stratum)
        g=unique(data$poststratum)
        for(j in 1:length(g))   
         {p=str(data$poststratum, g[j], sample.size)
          Ng=sum(NG[,j])
          t1=t2=0
          for (i in 1:length(h))         
            {s = str(data$Stratum, h[i], sample.size)
             shg=s*p
             if(!all(shg==0))  
                {
                nhg=length(shg[shg==1])
                t1=t1+sum(y[shg==1]/pik[shg==1])
                t2=t2+sum(1/pik[shg==1])               
                if(description)
                     {cat("Stratum ",j,", postratum ", i," \n")
                      cat("the postratified estimator is:",Ng*t1/t2,"\n")
                      }
                t_post=t_post+Ng*t1/t2
               }
             else if(description) 
                cat("Stratum ",j,", postratum ", i," empty intersection set \n")
             }}
    }
else {
      g=unique(data$poststratum)
      for(j in 1:length(g))   
         {p=str(data$poststratum, g[j], sample.size)
          Ng=NG[j]
          t1=Ng*sum(y[p==1]/pik[p==1])/sum(1/pik[p==1])
          t_post=t_post+t1
          if(description)
                     {cat("postratum ", j," \n")
                      cat("the postratified estimator is:",t1,"\n")
                     }
          }
       }
}
else  stop("the column names in data are missing")  
t_post
}

