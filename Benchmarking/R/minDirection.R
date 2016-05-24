# $Id: minDirection.R 140 2015-05-15 21:48:02Z B002961 $

# Function to calculate the min step for each input or max step for
# each output to the frontier, retninger i MEA. A series of LP problems

minDirection <- function(lps, m, n, ORIENTATION, LP=FALSE)  {
   # 'md' antal elementer/varer i direction
   if ( ORIENTATION=="in" )  md <- m
   if ( ORIENTATION=="out" )  md <- n
   if ( ORIENTATION=="in-out" )  md <- m+n

   # Saet taeller for fooerste vare
   mn0 <- switch(ORIENTATION, "in"=0, "out"=m, "in-out"=0)
  
   Direct <- rep(NA,md)
   for ( h in 1:md )  {
      if (LP) print(paste(" -->  Vare",h),quote=FALSE)
      # Saet 1 i 0'te raekke svarende til kriterifunktion -- default er 0 -- og
      # -1 i raekken for den relevante vare/element
      set.column(lps, 1, c(1,-1), c(0,mn0+h))
      if (LP) print(lps)
      set.basis(lps, default=TRUE)
      status <- solve(lps)
      if (LP) print(paste("Status =",status),quote=FALSE)
      if (LP) print(get.objective(lps))
      Direct[h] <- get.objective(lps)
   }
   lpcontr <- lp.control(lps)
   eps <- sqrt(lpcontr$epsilon["epsint"])
   if (LP) print(paste("eps =",eps))
   ## Direct[ abs(Direct) < eps ] <- 0
   if (LP) { print("Min direction:"); print(Direct) }
   return(Direct)
}



# mea er en wrapper for dea med special vaerdi af DIRECT

mea <-  function(X,Y, RTS="vrs", ORIENTATION="in", XREF=NULL, YREF=NULL,
         FRONT.IDX=NULL, param=NULL,
         TRANSPOSE=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)  {

   e <- dea(X,Y, RTS, ORIENTATION, XREF, YREF,
            FRONT.IDX, SLACK=FALSE, DUAL=FALSE, DIRECT="min", param=param, 
            TRANSPOSE=FALSE, LP=LP, CONTROL=CONTROL, LPK=LPK)

   return(e)

}



# Tegn linjer for MEA

mea.lines <- function(N, X, Y, ORIENTATION="in")  {
orientation <- c("in-out","in","out")
if ( is.numeric(ORIENTATION) )  {
	ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
	ORIENTATION <- ORIENTATION_
}
ORIENTATION <- tolower(ORIENTATION)
if ( !(ORIENTATION %in% orientation) ) {
	stop(paste("Unknown value for ORIENTATION:",ORIENTATION))
}

for ( n in N )  {
   if (ORIENTATION=="in")  {
      abline(h=X[n,2],lty="dotted")
      abline(v=X[n,1],lty="dotted")
      dir1 <- c(X[n,1],0)
      dir2 <- c(0,X[n,2])
   } else if (ORIENTATION=="out")  {
      abline(h=Y[n,2],lty="dotted")
      abline(v=Y[n,1],lty="dotted")
      dir1 <- c(Y[n,1],0)
      dir2 <- c(0,Y[n,2])
   } else if (ORIENTATION=="in-out")  {
      abline(h=Y[n,1],lty="dotted")
      abline(v=X[n,1],lty="dotted")
      dir1 <- c(X[n,1],0)
      dir2 <- c(0,Y[n,1])
   } else stop("Only directins in, out or in-out allowed in mea.lines")
   vn <- dea(X[n,,drop=FALSE],Y[n,,drop=FALSE], ORIENTATION=ORIENTATION,
             XREF=X, YREF=Y,FAST=TRUE,DIRECT=dir1)
   hn <- dea(X[n,,drop=FALSE],Y[n,,drop=FALSE], ORIENTATION=ORIENTATION,
             XREF=X, YREF=Y,FAST=TRUE,DIRECT=dir2)
   print(paste("Nr",n,":  vn =",vn,";  hn =",hn), quote=FALSE)
   print(paste("dir = (",dir1,",",dir2,")"), quote=FALSE)
   if (ORIENTATION=="in")  {
      abline(h=(1-hn)*X[n,2],lty="dotted")
      abline(v=(1-vn)*X[n,1],lty="dotted")
      # lines(c((1-vn)*X[n,1], X[n,1]), c(h=(1-hn)*X[n,2], X[n,2]),lw=2)
      arrows(X[n,1], X[n,2], (1-vn)*X[n,1], (1-hn)*X[n,2], lwd=2 )
	   dir <- c(vn*X[n,1], hn*X[n,2])
		print(paste("dir =", dir), quote=FALSE)
		mm <- dea(X[n,,drop=FALSE],Y[n,,drop=FALSE], ORIENTATION=ORIENTATION,
		          XREF=X, YREF=Y,FAST=TRUE,DIRECT=dir)
	   print(paste("mm =",mm) , quote=FALSE)
	   points(X[n,1]-mm*dir[1], X[n,2]-mm*dir[2], pch=16, col="green")
  
      abline(0, X[n,2]/X[n,1],lty="dashed", col="red")
      abline(X[n,2] - X[n,1]*(X[n,2]-(1-hn)*X[n,2])/(X[n,1]-(1-vn)*X[n,1]), 
          (X[n,2] - (1-hn)*X[n,2])/(X[n,1] - (1-vn)*X[n,1]) , 
          lty="dashed", col="blue")
   } else if (ORIENTATION=="out")  {
      # print(paste("(1+vn)*Y =",(1+hn)*Y[n,2],";  (1+hn)*Y =",(1+vn)*Y[n,1]))
      abline(h=(1+hn)*Y[n,2],lty="dotted")
	   abline(v=(1+vn)*Y[n,1],lty="dotted")
	   arrows(Y[n,1], Y[n,2], (1+vn)*Y[n,1], (1+hn)*Y[n,2], lwd=2 )
	   dir <- c(vn*Y[n,1], hn*Y[n,2])
		print(paste("dir =", dir), quote=FALSE)
		mm <- dea(X[n,,drop=FALSE],Y[n,,drop=FALSE], ORIENTATION=ORIENTATION,
		          XREF=X, YREF=Y,FAST=TRUE,DIRECT=dir)
	   print(paste("mm =",mm) , quote=FALSE)
	   points(Y[n,1]+mm*dir[1], Y[n,2]+mm*dir[2], pch=16, col="green")

      abline(0, Y[n,2]/Y[n,1],lty="dashed", col="red")
      abline(Y[n,2] - Y[n,1]*(Y[n,2]-(1-hn)*Y[n,2])/(Y[n,1]-(1-vn)*Y[n,1]), 
          (Y[n,2] - (1-hn)*Y[n,2])/(Y[n,1] - (1-vn)*Y[n,1]) , 
          lty="dashed", col="blue")
   } else if (ORIENTATION=="in-out")  {
      print(paste("(1+vn)*Y =",(1+hn)*Y[n,1],";  (1+hn)*X =",(1+vn)*X[n,1]), quote=FALSE)
      abline(h=(1+hn)*Y[n,1],lty="dotted")
	   abline(v=(1-vn)*X[n,1],lty="dotted")
	   arrows(X[n,1], Y[n,1], (1-vn)*X[n,1], (1+hn)*Y[n,1], lwd=2 )
	   
	   dir <- c(vn*X[n,], hn*Y[n,])
	   print(paste("dir =", dir), quote=FALSE)
	   mm <- dea(X[n,,drop=FALSE],Y[n,,drop=FALSE], ORIENTATION=ORIENTATION,
             XREF=X, YREF=Y,FAST=TRUE,DIRECT=dir)
      print(mm)
      print(paste("mm =",mm) , quote=FALSE)
      print(paste("x0 =",X[n,]+mm*dir[1], ";  y0 =",Y[n,]+mm*dir[2]), quote=FALSE)
      points(X[n,]-mm*dir[1], Y[n,]+mm*dir[2], pch=16, col="green")
   }

} 
}


# 
# smea <-  function(X,Y, RTS="vrs", ORIENTATION="in", XREF=NULL, YREF=NULL,
#          FRONT.IDX=NULL,
#          TRANSPOSE=FALSE, LP=FALSE, CONTROL=NULL, LPK=NULL)  {
# 
#    e <- sdea(X,Y, RTS, ORIENTATION, DIRECT="min", 
#              TRANSPOSE=FALSE, LP=LP)
# 
#    return(e)
# 
# }
# 
