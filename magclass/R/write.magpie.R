write.magpie <- function(x,file_name,file_folder="",file_type=NULL,append=FALSE,comment=NULL,comment.char="*",mode=NULL,nc_compression=9) {
  if(!is.null(mode)) {
    umask <- Sys.umask()
    umask_mode <- as.character(777-as.integer(mode))
    Sys.umask(umask_mode)
  }
  if(is.null(x)) x <- as.magpie(numeric(0))
  if(is.magpie(x)) {
    years <- !(is.null(dimnames(x)[[2]]))
    
    #if file-type is not mentioned file-ending is used as file-type
    if(is.null(file_type)) {
      file_type <- tail(strsplit(file_name,'\\.')[[1]],1)
    }
    if(!file_folder==""){
      file_path <- paste(file_folder,file_name,sep="/")
    }
    else{
      file_path <- file_name
    }  
    
    #look for comment/addtitional information
    if(is.null(comment) & !is.null(attr(x,"comment"))) comment <- attr(x,"comment")
    if(is.null(comment)) comment <- "" 
    
    #expand wildcards
    file_path <- paste(Sys.glob(dirname(file_path)),basename(file_path),sep="/")
    if(length(file_path)>1) {
      file_path <- file_path[1]
      warning("file name is ambiguous, only first alternative is used!")
    }
    
    if(append & file.exists(file_path)) {
      x2 <- read.magpie(file_path)
      x <- mbind(x2,x)
    }
    
    if(file_type=="m" | file_type=="mz") {
      fformat_version <- "2"  #File format version 1 (older data has version 0)
      comment <- paste(comment,collapse="\n")
      ncells <- dim(x)[1]
      nyears <- dim(x)[2]
      ndata  <- dim(x)[3]    
      rle <- rle(gsub("\\..*$","",dimnames(x)[[1]]))
      regions <- rle$values
      cpr <- rle$lengths
      nregions <- length(regions)
      regions_collapsed <- paste(regions,collapse='\n')
      datanames <- dimnames(x)[[3]]
      datanames_collapsed <- paste(datanames,collapse='\n')    
      sets_collapsed <- paste(getSets(x,fulldim = FALSE), collapse = '\n')
      
      if(years) {
        year_list <- as.integer(substr(dimnames(x)[[2]],2,5))
      } else {
        year_list <- 0
      }
      
      if(file_type=="mz") {
        zz <- gzfile(file_path,"wb")
      } else {
        zz <- file(file_path,"wb")
      }
      
      writeBin(as.integer(fformat_version),zz,size=2)
      writeBin(as.integer(nchar(comment)),zz,size=4)
      writeBin(as.integer(nchar(sets_collapsed)),zz,size=2)
      writeBin(as.integer(rep(0,92)),zz,size=1) #92 Byte reserved for later file format improvements
      writeBin(as.integer(c(nyears,year_list,nregions,nchar(regions_collapsed))),zz,size=2)
      writeChar(regions_collapsed,zz,eos=NULL)
      writeBin(as.integer(c(cpr,ndata*ncells*nyears,nchar(datanames_collapsed))),zz,size=4)
      if(datanames_collapsed!="") writeChar(datanames_collapsed,zz,eos=NULL)
      writeBin(as.numeric(as.vector(x)),zz,size=4)
      if(comment!="") writeChar(comment,zz,eos=NULL)
      if(nchar(sets_collapsed)>0) writeChar(sets_collapsed,zz,eos=NULL)
      close(zz)     
    } else if(file_type=="asc") {
      coord <- magclassdata$half_deg[,c("lon","lat")]
      if(dim(coord)[1]!=dim(x)[1]) stop("Wrong format! Only 0.5deg data can be written as ascii grid!")
      if(any(comment!="")) warning("asc format does not support comments!")
      for(y in 1:nyears(x)) {
        tmp_file <- ifelse(nyears(x)>1,sub("(\\.[^\\.]*)$",paste("_",getYears(x)[y],"\\1",sep=""),file_path),file_path)
        for(d in 1:ndata(x)) {
          tmp2_file <- ifelse(ndata(x)>1,sub("(\\.[^\\.]*)$",paste("_",getNames(x)[d],"\\1",sep=""),tmp_file),tmp_file)
          data <- as.data.frame(as.vector((x[,y,d])))
          grid <- suppressWarnings(sp::SpatialPixelsDataFrame(points = coord[c("lon", "lat")], data = data))
          sp::write.asciigrid(grid,tmp2_file)          
        }
      }
    } else if(file_type=="nc") {
      if (is.null(getNames(x)) | is.null(getYears(x))) stop("Year and Data name are necessary for saving to NetCDF format")
      mag <- as.array(x)
      
      #grid
      coord <- magclassdata$half_deg[,c("lon","lat")]
      coord[,"lon"]<-(coord[,"lon"]+180)/0.5+0.5
      coord[,"lat"]<-(coord[,"lat"]+90)/0.5+0.5
      # netcdf generation ####
      NODATA <- -9999
      
      # 4D array: lon, lat, time, data
      lon <- seq(-179.75,179.75,by=0.5)
      lat <- seq(-89.75,89.75,by=0.5)
      time <- as.numeric(unlist(lapply(strsplit(dimnames(mag)[[2]],"y"),function(mag) mag[2])))
      data <- dimnames(mag)[[3]]
      #data <- gsub("\\.","-",data)
      
      #Convert magpie data to array
      cat("Converting MAgPIE Data to 720 x 360 array")
      netcdf <- array(NODATA,dim=c(720,360,dim(mag)[2],dim(mag)[3]),dimnames=list(lon,lat,time,data))
      pb <- txtProgressBar(min = 0, max = dim(mag)[1], style = 3)
      for(i in 1:dim(mag)[1]){
        netcdf[coord[i,1],coord[i,2],,] <- mag[i,,,drop=FALSE]
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # NC file dimensions
      dim_lon <- ncdf4::ncdim_def("lon","degrees_east",lon)
      dim_lat <- ncdf4::ncdim_def("lat","degrees_north",lat)
      dim_time <- ncdf4::ncdim_def("time","years",time,calendar = "standard")
      
      #Define variables
      ncv <- list()
      for (i in dimnames(netcdf)[[4]]) ncv[[i]] <- ncdf4::ncvar_def(i, comment, list(dim_lon,dim_lat,dim_time), NODATA, prec="double",compression=nc_compression)
      
      #Create file
      if (file.exists(file_path)) file.remove(file_path)
      ncf <- ncdf4::nc_create(file_path, ncv)
      
      #Put data into file
      cat("Saving to NetCDF format")
      pb <- txtProgressBar(min = 0, max = dim(netcdf)[4], style = 3)
      for (i in dimnames(netcdf)[[4]]) {
        ncdf4::ncvar_put(ncf, ncv[[i]], netcdf[,,,i])
        setTxtProgressBar(pb, which(dimnames(netcdf)[[4]] == i))
      }
      close(pb)
      ncdf4::nc_close(ncf)
      } else if(file_type=="cs3" | file_type=="cs3r") {
      if(file_type=="cs3r") dimnames(x)[[2]] <- sub("y","",dimnames(x)[[2]])
      x <- unwrap(x)
      if(dim(x)[1]==1 & length(grep("GLO",dimnames(x)[[1]]))==1) {
        dimnames(x)[[1]] <- "TODELETE"  
      } else {
        if(nregions(x) == dim(x)[1]) {
          dimnames(x)[[1]] <- sub("\\..*$","",dimnames(x)[[1]])  
        } else {
          dimnames(x)[[1]] <- sub("\\.","_",dimnames(x)[[1]])
        }
      }
      x <- wrap(x,map=list(NA,length(dim(x))))
      dimnames(x)[[1]] <- sub("^([^\\.]*)\\.([^\\.]*)","\\2\\.\\1",dimnames(x)[[1]])
      
      dimnames(x)[[1]] <- gsub("TODELETE","",dimnames(x)[[1]])  
      dimnames(x)[[1]] <- gsub("^\\.","",dimnames(x)[[1]])  
      dimnames(x)[[1]] <- gsub("\\.$","",dimnames(x)[[1]])  
      dimnames(x)[[1]] <- gsub("\\.\\.","\\.",dimnames(x)[[1]])
      dimnames(x)[[1]] <- gsub("\\.",",",dimnames(x)[[1]])
      
      
      header <- dimnames(x)[[2]]
      x <- cbind(dimnames(x)[[1]],x)
      dimnames(x)[[2]] <- c(gsub("[^,]*(,|$)","dummy\\1",x[1,1]),header)
      zz <- file(file_path,open="w")
      if(any(comment!="")) writeLines(paste(comment.char,comment,sep=""),zz)
      write.csv(x,file=zz,quote=FALSE,row.names=FALSE)
      close(zz)
    } else if(file_type=="cs4" | file_type=="cs4r") {
      print_cells <- nregions(x)<ncells(x)
      print_regions <- getRegions(x)[1]!="GLO"
      print_data <- ((ndata(x)>1) | !is.null(getNames(x)))
      
      output <- as.data.frame(x)
      output <- output[c("Year","Region","Cell",names(output)[-c(1:3)])]
      
      if(!print_cells) output["Cell"] <- NULL
      if(!print_regions) output["Region"] <- NULL
      if(!print_data) output["Data1"] <- NULL
      if(!years) {
        output["Year"] <- NULL
      } else {
        if(file_type=="cs4") levels(output[["Year"]]) <- paste0("y",levels(output[["Year"]]))
      }
      zz <- file(file_path,open="w")
      if(any(comment!="")) writeLines(paste(comment.char,comment,sep=""),zz)
      write.table(output,file=zz,quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
      close(zz)
      
    } else {
      print_cells <- nregions(x)<ncells(x)
      print_regions <- getRegions(x)[1]!="GLO"
      print_data <- ((ndata(x)>1) | !is.null(getNames(x)))
      
      #non-cellular data
      if(!print_cells & (!print_data | !years | !print_regions )) {
        if(file_type=="csvr" | file_type=="cs2r") dimnames(x)[[2]] <- sub("y","",dimnames(x)[[2]])
        if(!print_data) {
          output <-  array(x,dim=dim(x)[1:2],dimnames=list(dimnames(x)[[1]],dimnames(x)[[2]]))
          output <- aperm(output)
          if(print_regions) {
            output <- rbind(substring(dimnames(x)[[1]],1,3),output)
            if(years) output <- cbind(c("dummy",dimnames(x)[[2]]),output)
          } else {
            if(years) output <- cbind(dimnames(x)[[2]],output)          
          }
          header <- FALSE
        } else if(!years) {
          output <-  array(x,dim=dim(x)[c(1,3)],dimnames=list(dimnames(x)[[1]],dimnames(x)[[3]]))
          header <- !is.null(dimnames(output)[[2]])
          if(print_regions) output <- cbind(substring(dimnames(x)[[1]],1,3),output)
          if(header & !print_regions) {
            output <- t(output)
            header <- FALSE
            output <- cbind(dimnames(x)[[3]],output)
          }
        } else {
          output <-  array(x,dim=dim(x)[2:3],dimnames=list(dimnames(x)[[2]],dimnames(x)[[3]]))
          header <- !is.null(dimnames(output)[[2]])
          output <- cbind(dimnames(x)[[2]],output)
          dimnames(output)[[2]][1] <- "dummy"
        }
        if(header & print_regions) dimnames(output)[[2]][1] <- "dummy"
        zz <- file(file_path,open="w")
        if(any(comment!="")) writeLines(paste(comment.char,comment,sep=""),zz)
        write.table(output,zz,sep=",",col.names=header,row.names=FALSE,quote=FALSE)  
        close(zz)
      } else {      
        if(file_type=="csvr" | file_type=="cs2r") dimnames(x)[[2]] <- sub("y","",dimnames(x)[[2]])
        if(file_type=="cs2" | file_type=="cs2r") print_regions <- FALSE
        output <- array(NA,c(dim(x)[1]*dim(x)[2],dim(x)[3]+print_regions+print_cells+years))
      	output[,(1+print_regions+print_cells+years):dim(output)[2]] <- as.vector(as.matrix(x))
        if(years) {
          yearvec <- c()
          for(year in dimnames(x)[[2]]) yearvec <- c(yearvec,rep(year,dim(x)[1]))
          output[,1] <- yearvec
        }
      	if(print_regions) output[,1+years] <- substring(rep(dimnames(x)[[1]],dim(x)[2]),1,3)
      	if(print_cells) {
          if(file_type=="cs2"  | file_type=="cs2r") {
            output[,1+print_regions+years] <- rep(gsub(".","_",dimnames(x)[[1]],fixed=TRUE),dim(x)[2])
          } else {
            output[,1+print_regions+years] <- rep(1:dim(x)[1],dim(x)[2])
          }
       }
      	if(!is.null(dimnames(x)[[3]])) {
      	  dimnames(output)[[2]] <- c(rep("dummy",print_regions+print_cells+years),dimnames(x)[[3]])
      	  header <- TRUE
      	} else {
      	  header <- FALSE
      	}
        zz <- file(file_path,open="w")
        if(any(comment!="")) writeLines(paste(comment.char,comment,sep=""),zz)
        write.table(output,zz,sep=",",col.names=header,row.names=FALSE,quote=FALSE)
        close(zz)
      }
    }
  } else {
    stop("Input is not in MAgPIE-format!")
  }
  if(!is.null(mode)) Sys.umask(umask)
}
