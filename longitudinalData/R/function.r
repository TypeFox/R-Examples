cat("\n####################################################################
############################## Function ############################
####################################################################")

cat("\n### Functions accepting NA ###")
meanNA <- function(x){mean(x,na.rm=TRUE)}
medianNA <- function(x){median(x,na.rm=TRUE)}

sdNA   <- function(x){sd(c(x),na.rm=TRUE)}
sdcNA <- function(x){leng<-length(x);sd(c(x),na.rm=TRUE)*sqrt((leng-1)/leng)}
varNA   <- function(x){var(c(x),na.rm=TRUE)}
rangeNA   <- function(x){range(x,na.rm=TRUE)}

which.minNA <- function(x){
  y <- which.min(x)
  if(length(y)==0){y<-NA}
  return(y)
}

### TRUE for Truly NA : false for NaN
is.tna <- function(x){
    if(length(x)==0){
        return(TRUE)
    }else{
        if(is.list(x)){x <- unlist(x)}else{}
        return(is.na(x)&!is.nan(x))
    }
}


### Printing long line shortening them
catShort <- function(x,nToCat=10){
    if(length(x)<=nToCat){
        cat(x)
    }else{
        cat(x[1:nToCat],"...")
    }
}

printMatrixShort <- function(mat,nColToPrint=10,nRowToPrint=5){
   if(ncol(mat)!=0){
        if(ncol(mat)>nColToPrint){
            trajToShow <- as.data.frame(mat[,1:nColToPrint])
            trajToShow$more <- "..."
        }else{
            trajToShow <- as.data.frame(mat)
        }
        if(nrow(mat)>nRowToPrint){
            print(trajToShow[1:nRowToPrint,])
            cat("... ...\n")
        }else{
            print(trajToShow)
        }
    }else{cat("   <no trajectories>\n")}
}


printShort <- function(x)  {
    if (length(x) <= 10) {
        print(x)
    }else {
        x <- cbind(x[1:10],"...") 
        names(x)[11] <- ""
        print(x)
    }
}
 

printOneTraj <- function(name,oneTraj){
   value <- data.frame(t(oneTraj[,2]))
   names(value) <- oneTraj[,1]
   row.names(value) <- paste(name,":",sep="")
   printShort(value)
}

printTrajLong <- function(trajLong,nRowToPrint=5){
   id <- unique(trajLong[,1])
   if(length(id)>nRowToPrint){id <- id[1:nRowToPrint]}else{}
   printOneTraj(id[1],trajLong[trajLong[,1]==id[1],2:3])
   for(i in id[-1]){
      cat("-------------------------------------------\n")
      printOneTraj(i,trajLong[trajLong[,1]==i,2:3])
   }
}


#NAtrunc <- function(x) x[1:max(which(!is.na(x)))]


reshapeLongToWide <- longToWide <- function(trajLong){
    if(ncol(trajLong)!=3){stop("[reduceNbTimesLong] The data.frame 'trajLong' has to be (no choice) in the following format:
    - first column should be the individual indentifiant;
    - the second should be the times at which the measurement are made;
    - the third one should be the measurement.")}else{}

    namesCol <- names(trajLong)
    trajLong <- trajLong[order(trajLong[,2]),]
    return(reshape(trajLong,idvar=namesCol[1],timevar=namesCol[2],v.names=namesCol[3],direction="wide"))
}


reshapeWideToLong <- wideToLong <- function(trajWide,times=1:(ncol(trajWide)-1)){
    id <- trajWide[,1]
    nbId <- nrow(trajWide)
    nbTimes <- ncol(trajWide)-1
    trajLong <- data.frame(id=rep(id,each=nbTimes),times=rep(times,nbId),values=as.numeric(t(as.matrix(trajWide[,-1]))))
    return(trajLong[!is.na(trajLong[,3]),])
}


## reshapeWide <- longToWide <- function(trajLong,idCol,timesCol,varyingCol){
##     toDrop <- names(trajLong)[!names(trajLong)%in%c(idCol,varyingCol,timesCol)]
##     trajLong <- trajLong[order(trajLong[,timesCol]),]
##     result <- reshape(trajLong,idvar=idCol,v.names=varyingCol,timevar=timesCol,drop=toDrop,direction="wide")
##     return(result)
## }

## reshapeLong <- wideToLong <- function(trajWide,idCol,varyingCol,times){
##     if(is.numeric(varyingCol)){varyingCol <- list(varyingCol)}
##     if(missing(times)){times <- 1:length(varyingCol[[1]])}
##     toDrop <- names(trajWide)[-unlist(varyingCol)]
##     toDrop <- toDrop[!(toDrop%in%idCol)]
##     return(reshape(trajWide,idvar=idCol,varying=varyingCol,times=times,drop=toDrop,direction="long"))
## }








cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++ Fin Function ++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")



