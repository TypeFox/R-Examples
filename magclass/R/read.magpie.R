read.magpie <- function(file_name,file_folder="",file_type=NULL,as.array=FALSE,old_format=FALSE,comment.char="*",check.names=TRUE) {
  
  file_name <- paste(file_folder,file_name,sep="")  

	if(length(Sys.glob(file_name))==0) {
	  stop(paste("file",file_name,"does not exist"))
	}

  #expand wildcards
  file_name_unexpanded <- file_name	
  file_name <- Sys.glob(file_name)
  if(length(file_name)>1) {
    file_name <- file_name[1]
    warning(paste("file name",file_name_unexpanded,"is ambiguous, only first alternative is used!"))
  } else if(length(file)==0) {
    stop("File ",file_name_unexpanded," could not be found!")
  }

  #if file-type is not mentioned file-ending is used as file-type
  if(is.null(file_type)) {
    file_type <- tail(strsplit(file_name,'\\.')[[1]],1)
  }
  if(!(file_type %in% c('m','mz','csv','cs2','cs3','cs4','csvr','cs2r','cs3r','cs4r','put',"asc","nc"))) stop(paste("Unkown file type:",file_type))

  .readComment <- function(file_name,comment.char="*") {
    comment <- NULL
    if(!is.null(comment.char)) {
      if(comment.char!="") {
        zz <- file(file_name)
        open(zz)
        read_repeat <- TRUE
        while(read_repeat){
          tmp <- readLines(zz,1)
          if(length(grep(paste("^",.escapeRegex(comment.char),sep=""),tmp))) {
            comment <- c(comment,tmp)
          } else {
            read_repeat <- FALSE
          }
        }
        close(zz)
      }
    }
    return(substring(comment,2))
  }

	if(file.exists(file_name)) {
    if(file_type=="m" | file_type=="mz") {

      if(file_type=="mz") {
        zz <- gzfile(file_name,"rb")
      } else {
        zz <- file(file_name,"rb")
      }
 
      if(!old_format) {
        fformat_version <- readBin(zz,integer(),1,size=2)
        nchar_comment <- readBin(zz,integer(),1,size=4)
        empty <- 94
        if(fformat_version > 1) {
          nchar_sets <- readBin(zz,integer(),1,size=2)
          empty <- empty - 2
        }
        readBin(zz,integer(),empty,size=1) #Bytes reserved for later file format improvements
      } else {
        fformat_version <- 0
      }
      nyears <- readBin(zz,integer(),1,size=2)
      year_list <- readBin(zz,integer(),nyears,size=2)
      nregions <- readBin(zz,integer(),1,size=2)
      nchar_regions <- readBin(zz,integer(),1,size=2)
      
      regions <- strsplit(readChar(zz,nchar_regions),"\n")[[1]]
      
      cpr <- readBin(zz,integer(),nregions,size=4)
      nelem <- readBin(zz,integer(),1,size=4)
      nchar_data <- readBin(zz,integer(),1,size=4)
      
      datanames <- strsplit(readChar(zz,nchar_data),"\n")[[1]]
      
      if(old_format) readBin(zz,integer(),100,size=1) #100 Byte reserved for later file format improvements

      output <- array(readBin(zz,numeric(),nelem,size=4),c(sum(cpr),nyears,nelem/sum(cpr)/nyears))
      output[is.nan(output)] <- NA
      if(any(cpr!=1)) {
        cellnames <- paste(rep(regions,cpr),1:sum(cpr),sep=".")
      } else {
        cellnames <- regions
      }
      if(length(cellnames)==1) cellnames <- list(cellnames)
      dimnames(output)[[1]] <- cellnames
      if(year_list[1]>0) dimnames(output)[[2]] <- paste("y",year_list,sep="")
      if(length(datanames)>0) dimnames(output)[[3]] <- datanames
      
      if(fformat_version > 0) {
        if(nchar_comment>0) attr(output,"comment") <- strsplit(readChar(zz,nchar_comment),"\n")[[1]]  
      }
      if(fformat_version > 1) {
        if(nchar_sets > 0) names(dimnames(output)) <- strsplit(readChar(zz,nchar_sets),"\n")[[1]]
      }
      close(zz)     
      attr(output,"FileFormatVersion") <- fformat_version
      read.magpie <- new("magpie",output)
    
    } else if(file_type=="cs3" | file_type=="cs3r") {
      x <- read.csv(file_name,comment.char=comment.char, check.names=check.names)
      datacols <- grep("^dummy\\.?[0-9]*$",colnames(x))
      xdimnames <- apply(x[datacols],2,unique)
      if(!is.list(xdimnames)) xdimnames <- list(xdimnames)
      xdimnames[[length(xdimnames)+1]] <- colnames(x)[-datacols]
      names(xdimnames) <- NULL
      tmparr <- array(NA,dim=sapply(xdimnames,length),dimnames=xdimnames)
      for(i in xdimnames[[length(xdimnames)]]) {
        j <- as.matrix(cbind(x[datacols],i))
        .duplicates_check(j)
        tmparr[j] <- x[,i]
      }
      read.magpie <- as.magpie(tmparr)   
      attr(read.magpie,"comment") <- .readComment(file_name,comment.char=comment.char)
    } else if(file_type=="cs4" | file_type=="cs4r") {
      x <- read.csv(file_name,comment.char=comment.char,header=FALSE)
      read.magpie <- as.magpie(x,tidy=TRUE)
      attr(read.magpie,"comment") <- .readComment(file_name,comment.char=comment.char)
    } else if(file_type=="asc"){
      grid<-suppressWarnings(try(maptools::readAsciiGrid(file_name,dec="."),silent=T))
      if(is(grid,"try-error")){
        grid<-try(maptools::readAsciiGrid(file_name,dec=","))
        if(is(grid,"try-error")) stop("File cannot be read. Make sure the file is in AsciiGrid format with either '.' or ',' as decimal point character.")
      }
      if(!all(grid@grid@cellsize==0.5)) stop("Only 0.5 degree data supported. Input data is in (",paste(grid@grid@cellsize,collapse=","),") degree (x,y).")
      #Convert to SpatialPixelsDataFrame
      sp::fullgrid(grid)<-FALSE
      magpie_coords<-as.matrix(magclassdata$half_deg[,c("lon","lat")])
      rowmatch<- function(A,B) {
        # Rows in A that match the rows in B
        f <- function(...) paste(..., sep=":")
        if(!is.matrix(B)) B <- matrix(B, 1, length(B))
        a <- do.call("f", as.data.frame(A))
        b <- do.call("f", as.data.frame(B))
        match(b, a)
      }
      mp_rows<-rowmatch(grid@coords,magpie_coords)
      names(mp_rows)<-1:59199
      if(any(is.na(mp_rows)))warning(sum(is.na(mp_rows))," magpie cells are missing in the grid file. They will be set to NA.")
      omitted_cells<-which(!(1:length(grid@data[[1]]))%in%mp_rows)
      if(length(omitted_cells)>0){
        omitted_fraction<-sum(grid@data[[1]][omitted_cells]/sum(grid@data[[1]]))
        warning(length(omitted_cells)," of ",length(grid@data[[1]])," cells in the file that contain data are discarded because they do not correspond to magpie cells.\n  Those cells contain ",omitted_fraction*100," percent of the global sum of the input file.")
      }
      read.magpie<-rep(-1001,59199)
      names(read.magpie)<-1:59199
      goodcells<-names(mp_rows)[which(!is.na(mp_rows))]
      read.magpie[goodcells]<-grid@data[[1]][mp_rows[goodcells]]
      read.magpie[is.na(mp_rows)]<-NA
      names(read.magpie)<-paste(magclassdata$half_deg$region,1:59199,sep=".")
      read.magpie<-as.magpie(read.magpie)      
    } else if(file_type=="nc") {
      nc_file <- ncdf4::nc_open(file_name)
      options("magclass.verbosity" = 1)
      if(is.null(nc_file$dim$time$len)) nc_file$dim$time$len <- 1
      if(is.null(nc_file$dim$time$vals)) nc_file$dim$time$vals <- 1995
      #create a single array of all ncdf variables
      nc_data <- array(NA,dim=c(nc_file$dim$lon$len,nc_file$dim$lat$len,nc_file$dim$time$len,nc_file$nvars))
      for (i in 1:nc_file$nvars) {
        nc_data[,,,i] <- ncdf4::ncvar_get(nc_file,varid=names(nc_file$var)[i])
      }
      #reorder ncdf array into magpie cellular format (still as array)
      coord <- magclassdata$half_deg[,c("lon","lat")]
      if(nc_file$dim$lat$vals[1] > 0) coord[,"lat"] <- coord[,"lat"]*-1
      coord[,"lon"]<-(coord[,"lon"]+180)/0.5+0.5
      coord[,"lat"]<-(coord[,"lat"]+90)/0.5+0.5
      read.magpie <- array(NA,dim=c(59199,nc_file$dim$time$len,nc_file$nvars),dimnames=list(paste(magclassdata$half_deg$region,1:59199,sep="."),paste("y",nc_file$dim$time$vals,sep=""),names(nc_file$var)))
      for(i in 1:ncells(read.magpie)){
        read.magpie[i,,] <- nc_data[coord[i,1],coord[i,2],,]
      }
      #convert array to magpie object
      read.magpie <- as.magpie(read.magpie)
    } else {
      #check for header
      if(file_type=="put") {
  		  temp <- read.csv(file_name,nrow=1,header=FALSE,sep="\t",comment.char=comment.char)      
      } else {
  		  temp <- read.csv(file_name,nrow=1,header=FALSE,comment.char=comment.char)
      }      

      #check for numeric elements in first row, which means a missing header
  		header <- TRUE
      for(temp_elem in temp) {
        if(is.numeric(temp_elem)) header <- FALSE
      }  		

      if(file_type=="put") {
  		  temp <- read.csv(file_name,header=header,sep="\t",comment.char=comment.char)      
      } else {
  		  temp <- read.csv(file_name,header=header,comment.char=comment.char)
      }
  		
  		#analyse column content
  		coltypes <- rep(0,dim(temp)[2])
  		for(column in 1:dim(temp)[2]) {
  		  if(sum(coltypes=="year")==0 & length(grep("^(y[0-9]{4}|[0-2][0-9]{3})$",temp[,column]))==dim(temp)[1]) {
  		    coltypes[column] <- "year"
        } else if(sum(coltypes=="region")==0 & sum(coltypes=="regiospatial")==0 & length(grep("^[A-Z]{3}$",temp[,column]))==dim(temp)[1]) {   
  		    coltypes[column] <- "region"
        } else if(sum(coltypes=="regiospatial")==0 & sum(coltypes=="region")==0 & length(grep("^[A-Z]{3}_[0-9]+$",temp[,column]))==dim(temp)[1]) {   
  		    coltypes[column] <- "regiospatial"  
		    } else if(!is.numeric(temp[1,column])) {
		      coltypes[column] <- "other"
        } else if(sum(coltypes=="cell")==0 & length(temp[,column])%%max(temp[,column])==0 & suppressWarnings(try(all(unique(temp[,column])==1:max(temp[,column])),silent=TRUE)==TRUE)) {
          coltypes[column] <- "cell"
        } else {
          coltypes[column] <- "data"
        }
  		}
  		  		
  		if(any(coltypes=="year")) {
  		  temp <- temp[order(temp[,which(coltypes=="year")]),]
  		  if(length(grep("y",temp[,which(coltypes=="year")]))==0) {
  		    temp[,which(coltypes=="year")] <- as.factor(paste("y",temp[,which(coltypes=="year")],sep=""))
  		  }
  		}
  		
  		#backup check if cell column is really a cell column
  		if(any(coltypes=="cell")){
  		  if(dimnames(temp)[[2]][which(coltypes=="cell")]=="iteration") {
          temp[,which(coltypes=="cell")] <- paste("iter",format(temp[,which(coltypes=="cell")]),sep="")
  		    coltypes[which(coltypes=="cell")] <- "other"
  		  } else if(header & !(dimnames(temp)[[2]][which(coltypes=="cell")]%in%c("dummy","dummy.1","dummy.2","dummy.3",""," ","cell","cells","Cell","Cells"))){
  		    coltypes[which(coltypes=="cell")] <- "data"
  		  } 
  		}
  		
  		if(any(coltypes=="cell")) {
  		  ncells <- dim(temp)[1]
  		  if(any(coltypes=="year")) ncells <- ncells/length(unique(temp[,which(coltypes=="year")]))
        if(any(coltypes=="other")) ncells <- ncells/length(unique(temp[,which(coltypes=="other")]))         
  		  if(!all(temp[1:ncells,which(coltypes=="cell")]==1:ncells)) coltypes[which(coltypes=="cell")] <- "data"
  		}
  		
  		#set all coltypes after the first occurrence of "data" to "data"
  		if(any(coltypes=="data")) coltypes[min(which(coltypes=="data")):length(coltypes)] <- "data"
     
      #set first columntype from "cell" to "data" if it seems that the data set is just a vector of numbers
      if(all(coltypes == c("cell",rep("data",length(coltypes)-1))) & dim(temp)[1]==1) coltypes[1] <- "data"
  		
  		#check coltypes for consistency
  		if(length(which(coltypes=="data"))==0) {
  		  print(coltypes)
  		  stop(paste("Inconsistency in data columns! No data column found in",file_name))  		
  		}
  		
      if(sum(coltypes=="data")!=(length(coltypes)-min(which(coltypes=="data"))+1)){
  		  print(coltypes)
  		  stop("Inconsistency in data columns!")
  		}
  		if(!all(which(coltypes=="data")==min(which(coltypes=="data")):length(coltypes))){
  		  print(coltypes)
  		  stop("Inconsistency in data columns!")
  		} 		
  		if(sum(coltypes=="data")==0){
  		  print(coltypes)
  		  stop("No data column found!")  		
  		}
  		if(sum(coltypes=="other")>1){
  		  print(coltypes)
  		  stop("Invalid format. More than one \"other\" column is not allowed!")
  		}		
      
      if(header) {
        if(length(grep("^y+[0-9]{4}$",dimnames(temp)[[2]][which(coltypes=="data")[1]]))==1) {
          headertype <- "year"
        } else if(length(grep("^[A-Z]{3}$",dimnames(temp)[[2]][which(coltypes=="data")[1]]))==1) {
          headertype <- "region"
        } else {
          headertype <- "other"
        }
      } else {
        headertype <- "none"
      }
  		
  		if(any(coltypes=="other")){
  		  othernames <- levels(as.factor(temp[,which(coltypes=="other")]))
  		  nother <- length(othernames)
  		  if(header) {
  		    if(headertype=="other") {
  		      elemnames <- dimnames(temp)[[2]][which(coltypes=="data")]
            elemnames <- paste(rep(othernames,each=length(elemnames)),elemnames,sep=".")
          } else {
            elemnames <- othernames
          }
        } else {
          if(sum(coltypes=="data")==1) {
            elemnames <- othernames
          } else {
  		      elemnames <- 1:sum(coltypes=="data")
            elemnames <- paste(rep(othernames,each=length(elemnames)),elemnames,sep=".")
          }
  		  }
  		  ncols <- length(elemnames)
		  } else {
		    nother <- 1
		    if(header) {
		      if(headertype=="other") {
            elemnames <- dimnames(temp)[[2]][which(coltypes=="data")]
            ncols <- length(elemnames)
          } else {
            elemnames <- NULL
            ncols <- 1
          }
		    } else {
		      elemnames <- NULL
		      ncols <- sum(coltypes=="data")
        }
		  }
      
      
      if(any(coltypes=="year")){
  		  yearnames <- levels(temp[,which(coltypes=="year")])
  		  nyears <- length(yearnames)
		  } else if(headertype=="year") {
		    yearnames <- dimnames(temp)[[2]][which(coltypes=="data")]
		    nyears <- length(yearnames)
		  } else {
		    yearnames <- NULL
		    nyears <- 1
      }
      
      if(any(coltypes=="cell")){
  		  ncells <- max(temp[,which(coltypes=="cell")])
		  } else {
		    if(headertype!="year"){
		      ncells <- dim(temp)[1]/(nyears*nother)
        } else {
          ncells <- dim(temp)[1]/nother  
        }
      }
      
      if(any(coltypes=="regiospatial")) {
        cellnames <- gsub("_",".",temp[1:ncells,which(coltypes=="regiospatial")],fixed=TRUE)
      } else {
        if(any(coltypes=="region")){
    		  regionnames <- levels(temp[,which(coltypes=="region")])
  		  } else if(headertype=="region") {
          regionnames <- dimnames(temp)[[2]][which(coltypes=="data")]       		  
          ncells <- ncells*length(regionnames)
  		  } else {
  		    regionnames <- "GLO"
        }
    		if(length(unique(regionnames)) < length(regionnames)) {
          cellnames <- paste(regionnames,1:ncells,sep=".")
    		} else {
          cellnames <- regionnames
    		}
  		}
  		if(length(cellnames)==1) cellnames <- list(cellnames)

  		
  		if(any(coltypes=="other") & (headertype=="other" | headertype=="none")) {
  		  output <- array(NA,c(ncells,nyears,ncols))
  		  dimnames(output)[[1]] <- cellnames
  		  dimnames(output)[[2]] <- yearnames
  		  dimnames(output)[[3]] <- elemnames 
  		  counter <- 0
        for(other.elem in othernames){
          output[,,(1:sum(coltypes=="data"))+counter] <- array(as.vector(
               as.matrix(temp[which(temp[,which(coltypes=="other")]==other.elem),
               which(coltypes=="data")])),c(ncells,nyears,sum(coltypes=="data")))
          counter <- counter + sum(coltypes=="data")
        } 
      } else if(!any(coltypes=="other") & headertype=="region") {
  		  output <- array(NA,c(ncells,nyears,ncols))
  		  dimnames(output)[[1]] <- cellnames
  		  dimnames(output)[[2]] <- yearnames
  		  dimnames(output)[[3]] <- elemnames 
     	  for(i in 1:length(cellnames)) {
      	  output[i,,1] <- temp[,which(coltypes=="data")[i]]
      	} 
      } else if(!any(coltypes=="other") & headertype=="year"){
        output <- array(NA,c(ncells,nyears,ncols))
  		  dimnames(output)[[1]] <- cellnames
  		  dimnames(output)[[2]] <- yearnames
  		  dimnames(output)[[3]] <- elemnames
      	for(year in yearnames) {
      	  output[,year,1] <- temp[,year]
      	}
      } else if(any(coltypes=="other") & headertype=="region") {
  		  output <- array(NA,c(ncells,nyears,ncols))
  		  dimnames(output)[[1]] <- cellnames
  		  dimnames(output)[[2]] <- yearnames
  		  dimnames(output)[[3]] <- elemnames 
     	  for(i in 1:length(cellnames)) {
     	    for(elem in elemnames) {
      	    output[i,,elem] <- temp[which(temp[,which(coltypes=="other")]==elem),which(coltypes=="data")[i]]
   	      }
      	} 
      } else if(any(coltypes=="other") & headertype=="year"){
        output <- array(NA,c(ncells,nyears,ncols))
  		  dimnames(output)[[1]] <- cellnames
  		  dimnames(output)[[2]] <- yearnames
  		  dimnames(output)[[3]] <- elemnames
      	for(year in yearnames) {
      	  for(elem in elemnames) {
      	    output[,year,elem] <- temp[which(temp[,which(coltypes=="other")]==elem),year]
    	    }
      	}        	    		
  		} else {
  		  output <- array(as.vector(as.matrix(temp[,which(coltypes=="data")])),c(ncells,nyears,ncols))
		    dimnames(output)[[1]] <- cellnames
  		  dimnames(output)[[2]] <- yearnames
  		  dimnames(output)[[3]] <- elemnames 		
      }
  		read.magpie <- output
      attr(read.magpie,"comment") <- .readComment(file_name,comment.char=comment.char)
    }
	} else {
	  warning(paste("File",file_name,"does not exist"))
		read.magpie <- NULL
	}
	if(as.array){
	  read.magpie <- as.array(as.magpie(read.magpie))[,,]
	} else {
	  read.magpie <- as.magpie(read.magpie)
	}
	return(read.magpie)
}
