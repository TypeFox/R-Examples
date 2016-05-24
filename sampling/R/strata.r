strata<-function(data, stratanames=NULL, size, method=c("srswor","srswr","poisson","systematic"),pik,description=FALSE)
{
if(missing(method)) {warning("the method is not specified; by default, the method is srswor")
                     method="srswor"
                    }
if(!(method %in% c("srswor","srswr","poisson","systematic"))) 
  stop("the method name is not in the list")
if(method %in% c("poisson","systematic") & missing(pik)) stop("the vector of probabilities is missing")
if(missing(stratanames)|is.null(stratanames))
   {
   if(length(size)>1) stop("the argument giving stratification variable is missing. The argument size should be a value.")
   if(method=="srswor")
	result=data.frame((1:nrow(data))[srswor(size,nrow(data))==1],rep(size/nrow(data),size))
   if(method=="srswr")
     {
      s=srswr(size,nrow(data))
      st=s[s!=0]
      l=length(st)
      result=data.frame((1:nrow(data))[s!=0])
      result=cbind.data.frame(result,st,prob=rep(1-(1-1/nrow(data))^size,l))
      colnames(result)=c("ID_unit","Replicates","Prob")
     }
   if(method=="poisson")
    {
     pikk=inclusionprobabilities(pik,size)
     s=(UPpoisson(pikk)==1)
     if(length(s)>0)
     result=data.frame((1:nrow(data))[s],pikk[s])
     if(description) 
      cat("\nPopulation total and number of selected units:",nrow(data),sum(s),"\n") 
     }
    if(method=="systematic")
    {
     pikk=inclusionprobabilities(pik,size)
     s=(UPsystematic(pikk)==1)
     result=data.frame((1:nrow(data))[s],pikk[s])
     }
   if(method!="srswr")
          colnames(result)=c("ID_unit","Prob")
   if(description & method!="poisson") cat("\nPopulation total and number of selected units:",nrow(data),sum(size),"\n")
   }
else
{
data=data.frame(data)
index=1:nrow(data)
m=match(stratanames,colnames(data))
if(any(is.na(m))) stop("the names of the strata are wrong")
data2=cbind.data.frame(data[,m],index)
colnames(data2)=c(stratanames,"index")
x1=data.frame(unique(data[,m]))
colnames(x1)=stratanames
result=NULL
for(i in 1:nrow(x1))
	{
if(is.vector(x1[i,]))
data3=data2[data2[,1]==x1[i,],]
else
{as=data.frame(x1[i,])
 names(as)=names(x1)
data3=merge(data2, as, by = intersect(names(data2), names(as)))
}
y=sort(data3$index)
if(description & method!="poisson") 
 {cat("Stratum" ,i,"\n")
  cat("\nPopulation total and number of selected units:",length(y),size[i],"\n")
  }
if(method!="srswr" & length(y)<size[i]) 
	{stop("not enough obervations in the stratum ", i, "\n") 
	 st=c(st,NULL)
       }
else
	{if(method=="srswor")
		{st=y[srswor(size[i],length(y))==1]
             r=cbind.data.frame(data2[st,],rep(size[i]/length(y),size[i]))
             }
       if(method=="systematic")
		{pikk=inclusionprobabilities(pik[y],size[i])
             s=(UPsystematic(pikk)==1)
             st=y[s]
             r=cbind.data.frame(data2[st,],pikk[s])
            }
 
	 if(method=="srswr") 
	   {s=srswr(size[i],length(y))
          st=rep(y[s!=0],s[s!=0]) 
          l=length(st)
          r=cbind.data.frame(data2[st,],prob=rep(1-(1-1/length(y))^size[i],l))
          }
       if(method=="poisson") 
		{pikk=inclusionprobabilities(pik[y],size[i])
             s=(UPpoisson(pikk)==1)
             if(any(s))
             { 
             st=y[s]
             r=cbind.data.frame(data2[st,],pikk[s])
             if(description) 
		 {cat("Stratum" ,i,"\n")
		  cat("\nPopulation total and number of selected units:",length(y),length(st),"\n")
		  }
               }
             else 
		{if(description) 
		 {cat("Stratum" ,i,"\n")
		  cat("\nPopulation total and number of selected units:",length(y),0,"\n")
              }
             r=NULL 
     	       }
            }
  
	}
if(!is.null(r))
	{
r=cbind(r,i)
result=rbind.data.frame(result,r)
        } 

 }  
if(method=="srswr")
colnames(result)=c(stratanames,"ID_unit","Prob","Stratum")
else
colnames(result)=c(stratanames,"ID_unit","Prob","Stratum")
if(description)  {cat("Number of strata ",nrow(x1),"\n")
                  if(method=="poisson")
   			 cat("Total number of selected units", nrow(result),"\n")
                  else       
                  cat("Total number of selected units", sum(size),"\n")
                  }
}
result
}

