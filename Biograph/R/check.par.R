check.par <-
function (Bdata)
{ # print ('check.par: the function checks whether "format.date" is present (attribute of Biograph object) and all parameters are defined (by function Parameters) ')
	#  Check whether Parameters was called
  # Parameters: attr (Bdata,"param")
    format.in <- attr(Bdata,"format.date") 
  if (is.null(format.in))
    { print (" ",quote=FALSE)
      if(max(Bdata$born,na.rm=TRUE) > 500 & max(Bdata$born,na.rm=TRUE) < 1800) format.in <- "CMC" else format.in <- "year"
      stop (paste ('Function Parameters: date format (attribute format.in) missing from Biograph object (data). Please add the attribute. Biograph expects date format to be: ',format.in,sep=""))
    } 

  if (is.null(attributes(Bdata)$param)) z<- "Run Parameters first" else z <- "OK"
  # Check whether state space differs from state space in memory (maybe from previous data set)
   # NOTE: namstates (sequence) may be different after calling Parameters! 
#  t<- namstates3%in%namstates
#  if (FALSE%in%t) 
#    {  print ("check.par: NOTE: Your data have a different state space. Function Parameters is called.",quote=FALSE)
#       z<- Parameters(Bdata)
#      print ("statespace: ")
#       print (namstates)  }
  # check whether path is character variable (should not be factor)
  if (!is.character(Bdata$path)) stop("Path variable is not character. Please check the data.")
  if (!is.null(Bdata$param))
 { # Check iagelow and iagehigh
    iagelow <- attr (Bdata,"param")$iagelow
    iagehigh <- attr (Bdata,"param")$iagehigh
   if (iagelow>100) stop(paste("Lowest age is ",iagelow,", which seems too high. Please run Parameters and check lowest and highest ages.",sep=""),quote=FALSE)
  nage <- iagehigh - iagelow + 1 
if (nage > 150) 
   { warning("check.par: Please check ages. The number of age groups exceeds 150. Details preceding this message.")
   	print (paste("Lowest age = ",iagelow,sep=""),quote=FALSE)
   	print (paste("Highest age = ",iagehigh,sep=""),quote=FALSE) 
   	print (paste("Number of age groups = ",nage,sep=""),quote=FALSE) }
}
  return(z)
}
