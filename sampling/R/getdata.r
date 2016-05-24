getdata<-function(data,m)
{
if(!is.data.frame(data)) data=as.data.frame(data)
if(is.null(names(data))) stop("the column names are missing")
if(is.vector(m) & !is.list(m))
   {res=NULL 
    if(is.null(names(m)))
      if(all(m %in% c(0,1))) 
      {
      res=NULL
      if(!("ID_unit" %in% names(data))) 
      {res=cbind.data.frame((1:length(m))[m==1],data[m==1,])
      names(res)=c("ID_unit",names(data))
      }
      else res=data[m==1,]
        }
      else  
   res=data[m,]
    }
else
if(is.data.frame(m))
     {res=NULL 
      if(!is.null(names(m)))
		{   
      mm=match(names(data),names(m),nomatch=0)
      index=(1:ncol(data))[mm==0]
      if(length(index)>0)
	      {res=cbind.data.frame(data[m$"ID_unit",index],m)
             names(res)[1:length(index)]=names(data)[index] 
 		}
       else res=m
          }
     }
else
if(is.list(m))
 {res=list()
  if(length(m)>=1)
  for(j in 1:length(m)) 
	{
     mm=match(names(data),names(m[[j]]),nomatch=0)
     index=(1:ncol(data))[mm==0]
     if(length(index)>0)
     {
     res[[j]]=cbind.data.frame(data[m[[j]]$"ID_unit",index],m[[j]])
     names(res[[j]])[1:length(index)]=names(data)[index] 
      }}
  else res=m
} 
res
}




