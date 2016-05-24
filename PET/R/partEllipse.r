partEllipse <- function(mod="hcirc1", x=0, y=0, intensity=0.1, n=257, re1=0.2,                                  re2=0.2,ring.wide=0.05,
               no.xax1=-1, no.xax2=1, no.yax1=-1, no.yax2=1, in.r=0.01,
               DebugLevel="HardCore")
{
#=============================================================================
#                           Dokumentation fehlt noch !!!!
# vom mittelpunkt verrueckte punkte
# x <- 0.2
# y <- -0.2
# radius der ellipsen
# re1 <- re2 <- .2
# interessierende werte
#no.xax1 <- 0.2
#no.xax2 <- 1
#no.yax1 <- -0.2
#no.yax2 <- 0.5
# radius der neuen kreise
#in.r <- 0.01
# ringbreite
#ring.wide <- 0.05
#=============================================================================

 DL1 <- logDebug(DebugLevel)
 DebugLevel <- DL1[[3]]
 DL2 <- DL1[[2]]
 DL1 <- DL1[[1]]

 if ( (mod != "hcirc1") && (mod != "hcirc2") && (mod != "hcirc3") ){
   cat("The mode mod =", mod, "is unknown. The default 'hcirc1' is used. \n")
 }
 if (intensity == 0) stop("intensity=", intensity, " don't be allowed.")
 if (ring.wide <= 0 || ring.wide > 1){
   cat("WARNING: 'ring.wide' of", ring.wide, "not supported. Default is use. \n")
   ring.wide=0.05
 }
 
 if(DL1) cat("Start of function halfCircle. \n")
 data <- matrix(0,n,n)
 Ax <- matrix(seq(-1,1,length.out=n), nrow=n, ncol=n, byrow=TRUE)
 # Matrix Ax 90° gedreht
 rot90Ax <- matrix(seq(1,-1,length.out=n), nrow=n, ncol=n, byrow=FALSE)
 Axnew <- Ax-x
 Aynew <- rot90Ax - y
 index1 <- which(  (Axnew > no.xax1 & Axnew < no.xax2 & 
                    Aynew > no.yax1 & Aynew < no.yax2  )&
                 (((Axnew)^2)/re1 + ((Aynew)^2)/re2  <= 1)  & 
                 (((Axnew)^2)/re1 + ((Aynew)^2)/re2  >= 1-ring.wide)   )

 if (mod == "hcirc2"){
   if(DL1) cat("The calculation need a few of minutes --> ")
   lengInd <- length(index1)
   for (i in 1:lengInd){
       data[which( ((Axnew - Axnew[index1[i]])^2 + (Aynew - Aynew[index1[i]])^2)/in.r  <= 1  )] <- intensity
   }
   if(DL1) cat("complete. \n")
 } else if (mod == "hcirc3"){
   if(DL1) cat("The calculation need a few of minutes --> ")
   lengInd <- length(index1)
   for (i in 1:lengInd){
       Axneww <- Axnew - Axnew[index1[i]]
       Ayneww <- Aynew - Aynew[index1[i]]
 	 index2 <- which( ((Axneww)^2 + (Ayneww)^2)/in.r  <= 1  )
	 data[index2] <- data[index2] + intensity
   }
   # resize the inensity of data, but is grow up
   data <- (data/abs(max(data)))*abs(intensity)
   if(DL1) cat("complete. \n") 
 } else data[index1] <- data[index1] + intensity

return(data)

}