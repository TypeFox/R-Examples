# $Id: malmq.R 144 2015-06-24 19:49:51Z B002961 $

# Beregner Malmquist indeks og dekomponering af samme ud fra to perioder

# Det er ikke forudsat at der er samme antal units i hver periode, men
# indekset bliver kun beregnet for foreningsmaengden af units

# Metoden kan bruges alene, men er ogsaa taenkt som kaldt fra en anden
# metode for haandtering af flere perioder.

 malmq <- function(X0, Y0, ID0=NULL,  X1,Y1, ID1=NULL, 
          RTS="vrs", ORIENTATION="in",
          SLACK=FALSE, DUAL=FALSE, DIRECT=NULL, param=NULL,
          TRANSPOSE=FALSE, FAST=TRUE, LP=FALSE, CONTROL=NULL, LPK=NULL)  
{

   # Det er underforstaaet at X'er og Y'er er matricer 'units x var'

   rts <- c("fdh","vrs","drs","crs","irs","irs2","add","fdh+","fdh++","fdh0")
   if ( missing(RTS) ) RTS <- "vrs" 
   if ( is.numeric(RTS) )  {
      if (LP) print(paste("Number '",RTS,"'",sep=""),quote=F)
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      if (LP) print(paste("' is '",RTS,"'\n",sep=""),quote=F)
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) )  stop(paste("Unknown scale of returns:", RTS))

   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      stop(paste("Unknown value for ORIENTATION:",ORIENTATION))
   }

   m0<- dim(X0)[2]  # number of inputs
   n0<- dim(Y0)[2]  # number of outputs
   K0 <- dim(X0)[1]  # number of units, firms, DMUs
   m1 <- dim(X1)[2]  # number of inputs
   n1 <- dim(Y1)[2]  # number of outputs
   K1 <- dim(X1)[1]  # number of units, firms, DMUs

   # Hvis ID'er mangler bliver de sat til 1,...,K
   if ( is.null(ID0) )  ID0 <- seq(1,K0)
   if ( is.null(ID1) )  ID1 <- seq(1,K1)

#print(X0)
#print(class(X0))
#print(K0)
#print(ID0)
#print(length(ID0))
   # Laengden af ID skal svare til antal units i X og Y
   if ( K0 != length(ID0) ||  K0 != length(ID0) )
      stop("Number of units in X0 and Y0 must correspont to length of ID0")
   if ( K1 != length(ID1) ||  K1 != length(ID1) )
      stop("Number of units in X1 and Y1 must correspont to length of ID1")


   # Find faellesmaengden af units og tilhoerende indeks
   idlab <- intersect(ID0, ID1)
   id0 <- ID0 %in% idlab   # seq(1,length(id0))[id0 %in% idlab]
   id1 <- ID1 %in% idlab

   # Input og output for faellesmaengden af units
   x0 <- X0[id0,]   
   y0 <- Y0[id0,]
   x1 <- X1[id1,]   
   y1 <- Y1[id1,]

   # Skal teknologien i en periode bestemmes af alle dem der i
   # perioden uafhaengigt af om de er i den anden periode?  Som det er
   # gjort nu er teknologien bestemt af dem der er i begge perioder.

   e00 <- dea(x0, y0, RTS=RTS, ORIENTATION=ORIENTATION,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)
   e10 <- dea(x1, y1, RTS=RTS, ORIENTATION=ORIENTATION, XREF=x0,YREF=y0,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)
   e11 <- dea(x1, y1, RTS=RTS, ORIENTATION=ORIENTATION,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)
   e01 <- dea(x0, y0, RTS=RTS, ORIENTATION=ORIENTATION,, XREF=x1,YREF=y1,
         SLACK=SLACK, DUAL=DUAL, DIRECT=DIRECT, param=param,
         TRANSPOSE=TRANSPOSE, FAST=TRUE, LP=LP, CONTROL=CONTROL, LPK=LPK)

   tc <- sqrt(e10/e11 * e00/e01)
   ec <- e11/e00
   m  <- tc * ec
   mq <- sqrt(e10/e00 * e11/e01)

   
   return( list(m=m, tc=tc, ec=ec, mq=mq, id=idlab, id0=id0, id1=id1, 
                e00=e00, e10=e10, e11=e11, e01=e01) )

}