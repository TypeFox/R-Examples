Parameters <-
function (Bdata,newnamstates) {
if (missing(newnamstates)) newnamstates <- NULL
print (". . . . Running function Parameters . . . . ")
#-------------   state space   ----------
# ---- Parameters from Bdata: state space, sequence of states occupied -----
nsample <- nrow (Bdata)

statespace <- StateSpace (Bdata,newnamstates)
namstates <- statespace$namstates
numstates <- length(namstates)

#  ---------  Determine date format ----------
  format.in <- attr(Bdata,"format.date") 
  if (is.null(format.in))
    { print (" ",quote=FALSE)
      if(max(Bdata$born,na.rm=TRUE) > 500 & max(Bdata$born,na.rm=TRUE) < 1800) format.in <- "CMC" else format.in <- "year"
      stop (paste ('Function Parameters: date format (attribute format.in) missing from Biograph object (data). Please add, e.g.: <attr(GLHS,"format.date") <- "CMC">. Biograph expects date format to be: ',format.in,sep=""))
    } 

# --------- Determine minimum and highest age --------------------
if (attr(Bdata,"format.date") == "age")
  { iagelow <- min(Bdata$start)
    iagehigh <- max(Bdata$end)  }
if (attr(Bdata,"format.date")=="CMC"|attr(Bdata,"format.date")=="cmc") 
        { iagelow <- min(Bdata$start-Bdata$born)
          iagehigh <- max (Bdata$end-Bdata$born)
          iagelow <- iagelow/12
          iagehigh <- iagehigh/12}  
if (attr(Bdata,"format.date")=="year"|attr(Bdata,"format.date")=="YEAR" |attr(Bdata,"format.date")=="Year" ) 
    	{  zx <- Bdata$start-Bdata$born
    	   iagelow <- as.numeric(min(zx,na.rm=TRUE))
    	   zy <- Bdata$end-Bdata$born 
           iagehigh <- as.numeric(max(zy,na.rm=TRUE))  }
if (substr(attr(Bdata,"format.date"),1,1)=="%") 
    	{  zx <- as.Date(Bdata$start,format=attr(Bdata,"format.date"))-as.Date(Bdata$born,format=attr(Bdata,"format.born"))
    	   iagelow <- as.numeric(min(zx,na.rm=TRUE))/365.24
    	   zy <- as.Date(Bdata$end,format=attr(Bdata,"format.date"))-as.Date(Bdata$born,format=attr(Bdata,"format.born"))
           iagehigh <- as.numeric(max(zy,na.rm=TRUE))/365.25  }
if (attr(Bdata,"format.date")=="day"|attr(Bdata,"format.date")=="DAY")     
        {  iagelow <- 0
           iagehigh <- as.numeric(max(Bdata$end),na.rm=TRUE)
        }       
iagelow <- trunc(iagelow)
iagehigh <- trunc(iagehigh)+1

nage <- iagehigh-iagelow+1
namage <- iagelow:iagehigh

# ---------- covariates -------------
locpat <- locpath(Bdata)
ncovariates  <- locpat - 5
covariates <- colnames(Bdata)[5:(locpat-1)]
# ncmc_tr <-  max(nchar(Bdata$path)) - 1
maxtrans <- max(nchar(Bdata$path)) - 1


# ---------  Flow table of transitions  -----------
 print ("Exploring types of transitions")
 zt <- transitions (Bdata,newnamstates=namstates)  

#assign("format.in",format.in,envir=.GlobalEnv)
#assign("iagelow",iagelow,envir=.GlobalEnv)
#assign("iagehigh",iagehigh,envir=.GlobalEnv)
#assign("nage",nage,envir=.GlobalEnv)
#assign("namage",namage,envir=.GlobalEnv)
#assign ("nsample",nrow(Bdata),envir=.GlobalEnv)  
#assign("ncovariates",ncovariates,envir=.GlobalEnv)
#assign("ncmc_tr",ncmc_tr,envir=.GlobalEnv)
#assign("maxtrans",maxtrans,envir=.GlobalEnv)

fff <- list (nsample = nsample,
              numstates=numstates,
              namstates = namstates,
              absorbstates=statespace$absorbstates,
              iagelow=iagelow,
              iagehigh=iagehigh,
              namage = namage,
              nage=nage, # number of age groups
              maxtrans = maxtrans, # maximum number of transitions by an individual during observation window
              ntrans = zt$ntrans,  # number_of_transitions
              trans_possible = zt$trans_possible,
              tmat = zt$tmat,
              transitions = zt$transitions,
              nntrans = zt$nntrans,
              locpat = locpat, # location of path (state sequence)
              ncovariates=ncovariates,  # number_of_covariates
              covariates=covariates,
              format.date = attr(Bdata,"format.date"),
              format.born=attr(Bdata,"format.born"))
return (fff)
}
