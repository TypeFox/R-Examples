demographics<-function (dem = "P0010001", state, statefips = FALSE, level = c("tract", 
    "blk", "blkgrp", "cdp", "msa", "county"), msaname = NULL) 
{
    demographics.aux <- function(dem = "P0010001", state, statefips = FALSE, 
        level = c("tract", "blk", "blkgrp", "cdp", "msa", "county"), 
        msaname = NULL) {
        	
      if(!any(level%in%c("tract","blk", "blkgrp", "cdp", "msa", "county"))){
      	stop("Not appropriate choice!")
      }  	
        	
	 state <- check.state(state, statefips)
        if (is.null(state)) {
            stop("Not a State! \n")
        }
state2 <- state
dem.fun <- function(dem, state, level) {
            require(paste("UScensus2010", level, sep = ""), character.only = TRUE)
            x <- paste(state,".", level,"10", sep = "")
            data(list = x, envir = parent.frame())
            temp <- get(x)
            out <- temp@data[,dem]
            out<-matrix(out,nrow=NROW(temp@data),ncol=length(dem))
            ifelse(level=="cdp",rownames(out)<-temp$nameAlevel,rownames(out)<-temp$fips)
            colnames(out)<-dem
            out
        }

if(level!="msa"){
     require(paste("UScensus2010", level, sep = ""), character.only = TRUE) 
     data(list=paste(state,".",level,"10",sep=""))
     temp<-get(paste(state,".",level,"10",sep=""))
}


     if(level=="county"){
		out<-as.matrix(temp@data[,dem],nrow=NROW(temp@data),ncol=length(dem))
		rownames(out)<-temp$NAMELSAD10
		colnames(out)<-dem	
     }else if (level == "msa") {
            temp <- MSA(msaname = msaname, state = toupper(state2),level = "county")
            out<-as.matrix(temp@data[,dem],nrow=NROW(temp@data),ncol=length(dem))
			rownames(out)<-temp$NAMELSAD10
			colnames(out)<-dem				
        }else if (level == "tract") {
            out <- dem.fun(dem, state, level)
        }else if (level == "blkgrp") {
            out <- dem.fun(dem, state, level)
        }
        else if (level == "blk") {
            out <- dem.fun(dem, state, level)
        }
        else if(level == "cdp") {
            out <- dem.fun(dem, state, level)
        }else{
        	stop("Not appropriate choice!")
        	}
       
   out
}
    demographics.aux(dem, state, statefips, level, msaname)
}