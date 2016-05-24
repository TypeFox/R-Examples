poststrata<-function(data, postnames = NULL) 
{
 if (missing(data) | missing(postnames)) stop("incomplete input")
 data = data.frame(data)
 if(is.null(colnames(data))) stop("the column names in data are missing")  
 index = 1:nrow(data)
 m = match(postnames, colnames(data))
 if (any(is.na(m))) 
            stop("the names of the poststrata are wrong")
 data2 = cbind.data.frame(data[, m])
 x1 = data.frame(unique(data[, m]))
 colnames(x1) = postnames
 nr_post=0
 post=numeric(nrow(data))
 nh=numeric(nrow(x1))
 for(i in 1:nrow(x1)) 
    { expr=rep(FALSE, nrow(data2))
      for(j in 1:nrow(data2)) expr[j]=all(data2[j, ]==x1[i, ])
      y=index[expr]
      if(is.matrix(y)) 
        nh[i]=nrow(y)
      else nh[i]=length(y)
      post[expr]=i
      }
result=cbind.data.frame(data,post)
names(result)=c(names(data),"poststratum")
list(data=result, npost=nrow(x1))
}

