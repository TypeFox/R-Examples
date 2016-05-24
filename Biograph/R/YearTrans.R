YearTrans <-
function (Bdata) 
{   z <- check.par (Bdata)
	format.in <- attr(Bdata,"format.date")
	selectday <- 1
	 if (format.in=="year"|format.in=="years") Bdata2 <- Bdata else if (format.in=="age"|format.in=="ages") Bdata2 <- Bdata else Bdata2 <- date_b(Bdata,format.in,selectday,format.out="year")
   nsample <- nrow(Bdata2)
   maxtrans <- ncol(Bdata2)-locpath(Bdata2)
   yeartrans <- array(0,c(nsample,maxtrans))
   for (j in 1:maxtrans) 
     { j22 <- locpath(Bdata2)+j 
       yeartrans[,j] <- Bdata2[,j22]  }
   yearborn <-Bdata2$born
   yearentry <- Bdata2$start
   yearcens <- Bdata2$end
   year_trans <- array(0,c(nsample,18))
    year_trans <- cbind(Bdata2$ID,yearborn,yearentry,yearcens,yeartrans)
   dimnames (year_trans) <- list (ID=Bdata2$ID,transition=c("ID","born","entry","censored",paste("tr",1:maxtrans,sep="")))
  return (year_trans)
}
