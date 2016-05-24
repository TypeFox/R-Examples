Remove.intrastate <-
function (Bdata)
{
# Check whether intrastate transitions have been removed. 
  if (is.null(attr(Bdata,"param")))   # attribute trans is missing
    {    stop ("Attribute 'param' is missing from Biograph object. First run Parameters and add attribute.",quote=FALSE)	
    }
  z<- !is.na(diag(attr(Bdata,"param")$tmat))
  if (!TRUE %in%z) # all diagonal elements are NA
 {   print ("No intra-state transitions present in the data.",quote=FALSE)
 	 zz <- transitions (Bdata)
     return (Bdata)
  }  else
 { #  REMOVE INTRASTATE TRANSITIONS 
 	print (". . . . . .   Removes intrastate transitions . . . .",quote=FALSE) 
 	 locpat <- locpath(Bdata) 
 	 nsample <- nrow(Bdata)
     nn<- ncol(Bdata)-locpat
     z<- array(NA,dim=c(nsample,nn))
     for (k in 2:(nn+1))
     {z[,k-1]<- ifelse (substr(Bdata$path,k,k)!="" & substr(Bdata$path,k,k)==substr(Bdata$path,k-1,k-1),k-1,NA)
      }
      pp <- Bdata$path
      for(k in 2:(nn+1))
      { substr(pp,k,k) <- ifelse (substr(Bdata$path,k,k)!="" & substr(Bdata$path,k,k)==substr(Bdata$path,k-1,k-1)," ",substr(Bdata$path,k,k))
      }
      for (i in 1:nsample)
      { pp[i]<- string.blank.omit(pp[i])
      }
      dates <- Bdata[,(locpat+1):ncol(Bdata)]
      for (j in 1:nn)
      { dates[!is.na(z[,j]),j]<- NA}
      dd <- apply(dates,1,function(x) sort(x,na.last=TRUE))
      dates <- t(dd)
 print (". . . . Intrastate transitions removed. Recalculating transitions . . . . ",quote=FALSE)
 Bdata2 <- Bdata
 Bdata2$path <- pp
 Bdata2[,(locpat+1):ncol(Bdata2)] <- dates  
 z <- transitions (Bdata2)
 attr(Bdata2,"param") <- Parameters (Bdata2)
 attr(Bdata2,"format.date") <- attr(Bdata,"format.date")
 print ("A new Biograph object without intrastate transitions is returned.",quote=FALSE)
 return (Bdata2)
 }
 }
