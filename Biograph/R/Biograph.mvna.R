Biograph.mvna <-
function (Bdata)
# !!!!!!!!!!  REMOVE RECORDS WITH TRANSITION (5 to Cens)   !!!!!!!!!!!!!!!!!!!
#  because it generates "There is undefined transitions in the data set"

{
# ========== Bdata TO MVNA FORMAT ======================
# need to remove intrastate transitions
  if (is.null(attr(Bdata,"param"))) print ("Biograph.mvna: Parameters missing. Run Parameters . . . . ",quote=FALSE)
  z <- check.par(Bdata)
  namstates <-   attr(Bdata,"param")$namstates
  removed <- Remove.intrastate (Bdata) 
  Bdata2<- removed
  tmat <- attr(Bdata2,"param")$tmat
  namstates2 <- vector (mode="numeric",length=length(namstates))
  for (i in 1:length(namstates))
    { namstates2[i] <- grep(namstates[i],namstates)} 

  print ("      Biograph.mvna: Calls function Biograph.long . . . ")
  zm <- Biograph.long (Bdata2)
  Dmvna <- zm$Depisode
  colnames(Dmvna)[1:3] <- c("id","from","to")  # c("id","from","to","entry","exit")
  Dmvna <- subset(Dmvna,!is.na(Dmvna[,1]))
  Dmvna  <- data.frame (Dmvna)
   
  #Dmvna$from <- apply (Dmvna,1,function(x) grep(x[2],namstates))
  #Dmvna$to[Dmvna$DES=="cens"] <- "cens"
  # to <- ifelse (Dmvna$to=="cens","cens",apply(Dmvna,1,function (x) grep(x[3],namstates)))
 #  Dmvna$to <- unlist (to)
  if (is.null(attr(Bdata,"format.date"))) stop ('Function Parameters: date format (attribute format.in) missing from Biograph object (data). Please add, e.g.: <attr(GLHS,"format.date") <- "CMC">')
  format.in <- attr(Bdata,"format.date") 
  format.born <- attr(Bdata,"format.born")
 # print ("Biograph.mvna ")
 # print (attr(Bdata,"format.date")) # is 'age' if called from Cumrates
 # print (c(format.in,Dmvna$Tstart[1:20])) # format.in ="age"
  if (format.in=="age")
   {   Dmvna$entry <- Dmvna$Tstart 
       Dmvna$exit <- Dmvna$Tstop
   	
   }  else
   {  ystart <- date_convert (Dmvna$Tstart,format.in=format.in,format.out="year",format.born=format.born)
      ystop <- date_convert (Dmvna$Tstop,format.in=format.in,format.out="year",format.born=format.born)
      yborn <- date_convert (Dmvna$born,format.in=format.in,format.out="year",format.born=format.born)
      Dmvna$entry <- ystart  -yborn  # = Dmvna$Tstarta
      Dmvna$exit <- ystop  -yborn    # = Dmvna$Tstopa
    }
  # colnames(Dmvna)[1] <- "id"
  D2 <- data.frame(id=Dmvna$id,from=Dmvna$from,to=as.character(Dmvna$to),
    entry=as.numeric(Dmvna$entry),exit=as.numeric(Dmvna$exit))
  D2$exit <- ifelse (D2$exit <= D2$entry,D2$exit+1,D2$exit) # CORRECT if CMC exit = CMC entry
  attr(D2, "param") <- attr(Bdata2,"param")   
  attr(Dmvna, "param") <- attr(Bdata2,"param")   
  attr(Dmvna, "format.date") <- "age"
  attr(D2, "format.date") <- "age"
      
  print ("     Biograph.mvna completed: Object produced by Biograph.mvna: .$mvna ",quote=FALSE)
  return (list (D=D2,  
                D.cov=Dmvna,
                cens = "cens"))
}
