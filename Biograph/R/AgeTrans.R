AgeTrans <-
function (Bdata) 
{   z<- check.par (Bdata)
	format.in <- attr(Bdata,"format.date") 
	if (format.in=="age"|format.in=="ages"|format.in=="day"|format.in=="days") Bdata2 <- Bdata else 
	  	Bdata2 <- date_b(Bdata=Bdata,format.out="age")
    nsample <- nrow(Bdata2)
	maxtrans <- (ncol(Bdata2)-locpath(Bdata2))
   ages <- array(0,c(nrow(Bdata2),maxtrans))
   for (j in 1:maxtrans) 
      { j22 <- locpath(Bdata2)+j ; ages[,j] <- Bdata2[,j22]  }  
   dimnames (ages) <- list (ID=Bdata2$ID,paste("tr",1:maxtrans,sep=""))
   ageentry <- Bdata2$start
   agecens <- Bdata2$end
   ns = nchar (Bdata2$path)
   return (list (ages =ages,
                 ageentry = ageentry,
                 agecens = agecens,
                 st_entry = substr(Bdata2$path,1,1),
                 st_censoring = substr(Bdata2$path,ns,ns)))
}
