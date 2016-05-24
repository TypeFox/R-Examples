# Title:  ModisDownload 
# Version: 4.2, Augues 2015
# Author: Babak Naimi (naimi.b@gmail.com)

# Major changes have been made on this version comparing to the 2.x. Since the FTP is not supported anymore,
# the functions have been adjusted to support HTTP!

# (Thanks Tomislav Hengl as his script [spatial-analyst.net] was used as the main core of the first version)
# Description: This is the source for two functions (ModisDownload and MRTproj) in R language that assist you to download, mosaic, subset and reproject the MODIS image products.
# Dependent library: RCurl
# Reference: http://www.r-gis.net

# MRTproj function is a wrapper to MRT tools, can be used to resample and subset an HDF image and convert to GeoTiff

# Licence GPL v3


modisProducts <- function() {
  .ModisLPxxx <- NULL
  load(system.file("external/ModisLP.RData", package="rts"))
  return(.ModisLPxxx)
  rm(.ModisLPxxx)
}

#-----------------------

.modisHTTP <- function(x,v='005') {
  if (!requireNamespace("RCurl",quietly = TRUE)) stop("Package RCurl is not installed")
  mp <- modisProducts()
  if (is.numeric(x)) {
    if (x > nrow(mp)) stop("The product code is out of subscription!")
    x <- as.character(mp[x,1])
  }
  x <- trim(x)
  if ('e4ftl01.cr.usgs.gov' %in% strsplit(x,'/')[[1]]) {
    if (strsplit(x,'/')[[1]][1] != 'http:') x <- paste('http://',x,sep='')
    if (strsplit(x,'')[[1]][length(strsplit(x,'')[[1]])] != "/") x <- paste(x,"/",sep="")
    if (!RCurl::url.exists(x)) stop("the http address does not exist OR Server is down!")
  } else {
    w <- which(mp[,1] == x)
    if (length(w) != 1) stop("The Name does not exist in MODIS Land produnct list!")
    if (as.character(mp[w,2]) == "Terra") ad <- "MOLT"
    else if (as.character(mp[w,2]) == "Aqua") ad <- "MOLA"
    else ad <- "MOTA"
    xx <- paste("http://e4ftl01.cr.usgs.gov/",ad,"/",x,".",v,"/",sep="")
    if (!RCurl::url.exists(xx)) {
      if (!RCurl::url.exists(paste("http://e4ftl01.cr.usgs.gov/",ad,"/",sep=""))) stop("the http address does not exist! Version may be incorrect OR Server is down!")
      else {
        items <- try(strsplit(RCurl::getURL(paste("http://e4ftl01.cr.usgs.gov/",ad,"/",sep="")), "\r*\n")[[1]],silent=TRUE)
        dirs <- unlist(lapply(strsplit(unlist(lapply(strsplit(items[-c(1:19)],'href'),function(x){strsplit(x[2],'/')[[1]][1]})),'"'),function(x) {x[2]}))
        dirs <- na.omit(dirs)
        w <- which(unlist(lapply(strsplit(dirs,'\\.'),function(x) x[[1]])) == x)
        if (length(w) > 0) v <- unlist(lapply(strsplit(dirs,'\\.'),function(x) x[[2]]))[w]
        x <- paste("http://e4ftl01.cr.usgs.gov/",ad,"/",x,".",v,"/",sep="")
      }
    } else x <- xx
  }
  x
}

#-----------------
.getModisList <- function(x,h,v,dates) {
  if (!requireNamespace("RCurl",quietly = TRUE)) stop("Package RCurl is not installed")
  if (inherits(dates,"character")) dates <- as.Date(dates,format='%Y.%m.%d')
  dates <- na.omit(as.Date(dates))
  if (length(dates) == 0) stop("dates is not appropriately selected!")
  dates <- sort(dates)
  try.nr <- 3
  items <- 0
  class(items) <- "try-error"
  ce <- 0
  while(class(items) == "try-error") { 
    items <- try(strsplit(RCurl::getURL(x), "\r*\n")[[1]],silent=TRUE)
    if (class(items) == "try-error") {
      Sys.sleep(5)
      ce <- ce + 1
      if (ce == (try.nr+1)) stop("Download error: Server does not response!")
    }
  }
  items <- items[-c(1:19)]
  # get the directory names (available dates)
  dirs <- unlist(lapply(strsplit(unlist(lapply(strsplit(items,'href'),function(x){strsplit(x[2],'/')[[1]][1]})),'"'),function(x) {x[2]}))
  dirs <- na.omit(dirs)
  d <- as.Date(dirs,format='%Y.%m.%d')
  
  # extract the selected dates
  if (length(dates) == 2) {
    dirs <- dirs[which(d >= dates[1] & d <= dates[2])]
  } else dirs <- dirs[unlist(lapply(dates,function(x){which(d == x)}))]
  
  if (length(dirs) < 1) stop("No available data for the selected dates")
  if (length(dates) > 2 & length(dirs) < length(dates)) warning("The images are not available for some specified dates!")
  
  # creating a vector of available Modis images in selected tiles and dates
  Modislist <- list()
  
  for (i in 1:length(dirs)) {
    getlist <- 0
    class(getlist) <- "try-error"
    ce <- 0
    while(class(getlist) == "try-error") {
      getlist <- try(strsplit(RCurl::getURL(paste(x,dirs[i], "/", sep="")), "\r*\n")[[1]],silent=TRUE)
      if (class(getlist) == "try-error") {
        Sys.sleep(5)
        ce <- ce + 1
        if (ce == (try.nr+1)) stop("Download error: Server does not response!")
      }
    }
    getlist <- getlist[-c(1:19)]
    getlist <- unlist(lapply(strsplit(getlist,"href"),function(x){strsplit(x[2],'"')[[1]][2]}))
    w <- which(is.na(getlist))
    if (length(w) > 0) getlist <- getlist[-w]
    m <- c()
    for (vv in v) {
      for (hh in h) {
        if (vv < 10) vc <- paste('0',as.character(vv),sep='')
        else vc <- as.character(vv)
        if (hh < 10) hc <- paste('0',as.character(hh),sep='')
        else hc <- as.character(hh)
        ModisName <- grep(".hdf$",grep(paste('h',hc,'v',vc,sep=''),getlist,value=TRUE),value=TRUE)
        #if (length(ModisName) == 1) {
          m <- c(m,paste(x,dirs[i], "/",ModisName,sep='')[length(ModisName)])
          Modislist[[dirs[i]]] <- m
        #}
      }
    }
  }
  Modislist
}
#--------

.downloadHTTP <- function(x,filename) {
  if (!requireNamespace("RCurl",quietly = TRUE)) stop("Package RCurl is not installed")
  success <- FALSE
  er <- try(writeBin(RCurl::getBinaryURL(x),con=filename),silent=TRUE)
  if (class(er) == "try-error") print("Download Error: Server does not response!!")
  else success <- TRUE
  return(success)
}

.getMODIS <- function(x, h, v, dates, version='005') {
  xx <- .modisHTTP(x,v=version)
  Modislist <- .getModisList(xx,h=h,v=v,dates=dates)
  if (length(Modislist) == 0) stop("There is NO available images for the specified product!")
  out <- data.frame(matrix(nrow=0,ncol=2))
  names(out) <- c("Date","Name")
  dirs <- names(Modislist)
  
  for (d in dirs) {
    cnt <- 1
    for (ModisName in Modislist[[d]]) {
      n <- strsplit(ModisName,"/")[[1]]
      n <- n[length(n)]
      if (.downloadHTTP(ModisName,n)) {
        out <- rbind(out,data.frame(Date=d,Name=n))
      }
      cnt <- cnt + 1
    }
    
  }
  out
}

#---------


#------------



if (!isGeneric("getMODIS")) {
  setGeneric("getMODIS", function(x,h,v,dates,version='005')
    standardGeneric("getMODIS"))
}

if (!isGeneric("mosaicHDF")) {
  setGeneric("mosaicHDF", function(hdfNames,filename,MRTpath, bands_subset,delete=FALSE)
    standardGeneric("mosaicHDF"))
}
if (!isGeneric("reprojectHDF")) {
  setGeneric("reprojectHDF", function(hdfName,filename,MRTpath, ...)
    standardGeneric("reprojectHDF"))
}

if (!isGeneric("ModisDownload")) {
  setGeneric("ModisDownload", function(x,h,v,dates, ...)
    standardGeneric("ModisDownload"))
}

setMethod("mosaicHDF", "character",
          function(hdfNames,filename,MRTpath,bands_subset,delete=FALSE) {
            if (missing(MRTpath)) stop("MRTpath argument should be specified...")
            if (length(hdfNames) < 2) stop("mosaic cannot be called for ONE image!")
            if (missing(bands_subset))  bands_subset <- ''
            if (missing(delete)) delete <- FALSE
            
            mosaicname = file(paste(MRTpath, "/TmpMosaic.prm", sep=""), open="wt")
            write(paste(getwd(),"/",hdfNames[1], sep=""), mosaicname)
            for (j in 2:length(hdfNames)) write(paste(getwd(),"/",hdfNames[j], sep=""),mosaicname,append=T)
            close(mosaicname)
            # generate mosaic:
            
            if (bands_subset != '') {
              e <- system(paste(MRTpath, '/mrtmosaic -i ', MRTpath, '/TmpMosaic.prm -s "',bands_subset,'" -o ',getwd(), '/',filename, sep=""))
              if (e != 0) warning ("Mosaic failed! 'bands_subset' may has incorrect structure!")
            } else {
              e <- system(paste(MRTpath, '/mrtmosaic -i ', MRTpath, '/TmpMosaic.prm -o ',getwd(), '/',filename, sep=""))
              if (e != 0) warning ("Mosaic failed!")
            }
            if (delete & e == 0) for (ModisName in hdfNames) unlink(paste(getwd(), '/', ModisName, sep=""))
            if (e == 0) return (TRUE)
            else return (FALSE)
          }
)


setMethod("reprojectHDF", "character",
          function(hdfName,filename,MRTpath,UL="",LR="",resample_type='NEAREST_NEIGHBOR',proj_type='UTM',
                   bands_subset='',proj_params='0 0 0 0 0 0 0 0 0 0 0 0',datum='WGS84',utm_zone=NA,pixel_size=1000) {
            
            fname = file('tmp.prm', open="wt")
            write(paste('INPUT_FILENAME = ', getwd(), '/',hdfName, sep=""), fname) 
            if (bands_subset != '') {
              write(paste('SPECTRAL_SUBSET = ( ',bands_subset,' )',sep=''),fname,append=TRUE)
            }
            if (UL[1] != '' & LR[1] != '') {
              write('SPATIAL_SUBSET_TYPE = OUTPUT_PROJ_COORDS', fname, append=TRUE)
              write(paste('SPATIAL_SUBSET_UL_CORNER = ( ', as.character(UL[1]),' ',as.character(UL[2]),' )',sep=''), fname, append=TRUE)
              write(paste('SPATIAL_SUBSET_LR_CORNER = ( ', as.character(LR[1]),' ',as.character(LR[2]),' )',sep=''), fname, append=TRUE)
            }
            write(paste('OUTPUT_FILENAME = ', filename, sep=""), fname, append=TRUE)
            write(paste('RESAMPLING_TYPE = ',resample_type,sep=''), fname, append=TRUE)
            write(paste('OUTPUT_PROJECTION_TYPE = ',proj_type,sep=''), fname, append=TRUE)
            write(paste('OUTPUT_PROJECTION_PARAMETERS = ( ',proj_params,' )',sep=''), fname, append=TRUE)
            write(paste('DATUM = ',datum,sep=''), fname, append=TRUE)
            if (proj_type == 'UTM') write(paste('UTM_ZONE = ',utm_zone,sep=''), fname, append=TRUE)
            write(paste('OUTPUT_PIXEL_SIZE = ',as.character(pixel_size),sep=''), fname, append=TRUE)
            close(fname)
            e <- system(paste(MRTpath, '/resample -p ',getwd(),'/','tmp.prm', sep=''))
            if (e == 0) return (TRUE)
            else return(FALSE)
          }
          
)


setMethod("getMODIS", "character",
          function(x,h,v,dates,version='005') {
            xx <- .modisHTTP(x,v=version)
            Modislist <- .getModisList(xx,h=h,v=v,dates=dates)
            if (length(Modislist) == 0) stop("There is NO available images for the specified product!")
            
            dirs <- names(Modislist)
            out <- data.frame(matrix(nrow=length(dirs),ncol=3))
            out[,1] <- dirs
            dc <- 1
            for (d in dirs) {
              dwnld <- rep(FALSE,length(Modislist[[d]]))
              cnt <- 1
              out[dc,2] <- length(dwnld)
              for (ModisName in Modislist[[d]]) {
                n <- strsplit(ModisName,"/")[[1]]
                n <- n[length(n)]
                if (.downloadHTTP(ModisName,n)) dwnld[cnt] <- TRUE
                cnt <- cnt + 1
              }
              out[dc,3] <- length(which(dwnld))
              dc <- dc+1
            }
            if (sum(out[,3]) > 0) {
              cat(paste('from ', sum(out[,2]),' available images, ',sum(out[,3]),' images are successfully downloaded.',sep=''))
            } else cat('Download is failed!')
            
          }
)


setMethod("getMODIS", "numeric",
          function(x,h,v,dates,version='005') {
            xx <- .modisHTTP(x,v=version)
            Modislist <- .getModisList(xx,h=h,v=v,dates=dates)
            if (length(Modislist) == 0) stop("There is NO available images for the specified product!")
            
            dirs <- names(Modislist)
            out <- data.frame(matrix(nrow=length(dirs),ncol=3))
            out[,1] <- dirs
            dc <- 1
            for (d in dirs) {
              dwnld <- rep(FALSE,length(Modislist[[d]]))
              cnt <- 1
              out[dc,2] <- length(dwnld)
              for (ModisName in Modislist[[d]]) {
                n <- strsplit(ModisName,"/")[[1]]
                n <- n[length(n)]
                if (.downloadHTTP(ModisName,n)) dwnld[cnt] <- TRUE
                cnt <- cnt + 1
              }
              out[dc,3] <- length(which(dwnld))
              dc <- dc+1
            }
            if (sum(out[,3]) > 0) {
              cat(paste('from ', sum(out[,2]),' available images, ',sum(out[,3]),' images are successfully downloaded.',sep=''))
            } else cat('Download is failed!')
            
          }
)


setMethod("ModisDownload", "character",
          function(x,h,v,dates,version='005',MRTpath,mosaic=FALSE,bands_subset='',delete=FALSE,proj=FALSE,UL="",LR="",resample_type='NEAREST_NEIGHBOR',proj_type='UTM', proj_params='0 0 0 0 0 0 0 0 0 0 0 0',datum='WGS84',utm_zone=NA,pixel_size) {
            dHDF <- .getMODIS(x,h,v,dates,version)
            dHDF$Date <- as.character(dHDF$Date)
            dHDF$Name <- as.character(dHDF$Name)
            if (nrow(dHDF) < 2) mosaic <- FALSE
            
            if (mosaic) {
              
              du <- unique(dHDF$Date)
              
              for (d in du) {
                dw <- dHDF[which(dHDF$Date == d),]
                if (nrow(dw) > 1){
                  date_name <- sub(sub(pattern="\\.", replacement="-", d), pattern="\\.", replacement="-", d)
                  name <- paste("Mosaic_",date_name,".hdf",sep='')
                  Mosaic.success <- mosaicHDF(dw[,2],name,MRTpath=MRTpath,bands_subset=bands_subset,delete=delete)
                  if (Mosaic.success) {
                    if (delete) for (ModisName in dw[,2]) unlink(paste(getwd(), '/', ModisName, sep=""))
                    if (proj) {
                      pref <- strsplit(dw[1,2],'\\.')[[1]][1]
                      e <- reprojectHDF(name,filename=paste(pref,'_',date_name,'.tif',sep=''),MRTpath=MRTpath,UL=UL,LR=LR,proj_type=proj_type,proj_params=proj_params,utm_zone=utm_zone,pixel_size=pixel_size)
                      if (e & delete) unlink(paste(name))
                      if (!e) warning (paste("The procedure has failed to REPROJECT the mosaic image for date ",d,"!",sep=""))
                    }
                  } else {
                    warning(paste("The procedure has failed to MOSAIC the images for date ",d,"!",sep=""))
                    if (proj) {
                      warning ("Since the mosaic is failed, the individual hdf images are reprojected...")
                      pref <- strsplit(dw[1,2],'\\.')[[1]]
                      pref <- paste(pref[1],"_",pref[3],sep="")
                      for (ModisName in dw[,2]) {
                        e <- reprojectHDF(ModisName,filename=paste(pref,'_',date_name,'.tif',sep=''),MRTpath=MRTpath,UL=UL,LR=LR,bands_subset=bands_subset,proj_type=proj_type,proj_params=proj_params,utm_zone=utm_zone,pixel_size=pixel_size)
                        if (e & delete) unlink(paste(ModisName))
                        if (!e) warning (paste("The procedure has failed to REPROJECT the individual HDF image ",ModisName,"!",sep=""))
                      }
                    }
                  } 
                }
              }   
            } else {
              if (proj) {
                for (i in 1:nrow(dHDF)) {
                  ModisName <- dHDF[i,2]
                  pref <- strsplit(ModisName,'\\.')[[1]]
                  pref <- paste(pref[1],"_",pref[3],sep="")
                  d <- dHDF[i,1]
                  date_name <- sub(sub(pattern="\\.", replacement="-", d), pattern="\\.", replacement="-", d)
                  e <- reprojectHDF(ModisName,filename=paste(pref,'_',date_name,'.tif',sep=''),MRTpath=MRTpath,UL=UL,LR=LR,bands_subset=bands_subset,proj_type=proj_type,proj_params=proj_params,utm_zone=utm_zone,pixel_size=pixel_size)
                  if (e & delete) unlink(paste(ModisName))
                  if (!e) warning (paste("The procedure has failed to REPROJECT the individual HDF image ",ModisName,"!",sep=""))
                }
              }
            }
            
          }
)


setMethod("ModisDownload", "numeric",
          function(x,h,v,dates,version='005',MRTpath,mosaic=FALSE,bands_subset='',delete=FALSE,proj=FALSE,UL="",LR="",resample_type='NEAREST_NEIGHBOR',proj_type='UTM', proj_params='0 0 0 0 0 0 0 0 0 0 0 0',datum='WGS84',utm_zone=NA,pixel_size) {
            dHDF <- .getMODIS(x,h,v,dates,version)
            dHDF$Date <- as.character(dHDF$Date)
            dHDF$Name <- as.character(dHDF$Name)
            if (nrow(dHDF) < 2) mosaic <- FALSE
            if (mosaic) {
              
              du <- unique(dHDF$Date)
              
              for (d in du) {
                dw <- dHDF[which(dHDF$Date == d),]
                if (nrow(dw) > 1){
                  date_name <- sub(sub(pattern="\\.", replacement="-", d), pattern="\\.", replacement="-", d)
                  name <- paste("Mosaic_",date_name,".hdf",sep='')
                  Mosaic.success <- mosaicHDF(dw[,2],name,MRTpath=MRTpath,bands_subset=bands_subset,delete=delete)
                  if (Mosaic.success) {
                    if (delete) for (ModisName in dw[,2]) unlink(paste(getwd(), '/', ModisName, sep=""))
                    if (proj) {
                      pref <- strsplit(dw[1,2],'\\.')[[1]][1]
                      e <- reprojectHDF(name,filename=paste(pref,'_',date_name,'.tif',sep=''),MRTpath=MRTpath,UL=UL,LR=LR,proj_type=proj_type,proj_params=proj_params,utm_zone=utm_zone,pixel_size=pixel_size)
                      if (e & delete) unlink(paste(name))
                      if (!e) warning (paste("The procedure has failed to REPROJECT the mosaic image for date ",d,"!",sep=""))
                    }
                  } else {
                    warning(paste("The procedure has failed to MOSAIC the images for date ",d,"!",sep=""))
                    if (proj) {
                      warning ("Since the mosaic is failed, the individual hdf images are reprojected...")
                      pref <- strsplit(dw[1,2],'\\.')[[1]]
                      pref <- paste(pref[1],"_",pref[3],sep="")
                      for (ModisName in dw[,2]) {
                        e <- reprojectHDF(ModisName,filename=paste(pref,'_',date_name,'.tif',sep=''),MRTpath=MRTpath,UL=UL,LR=LR,bands_subset=bands_subset,proj_type=proj_type,proj_params=proj_params,utm_zone=utm_zone,pixel_size=pixel_size)
                        if (e & delete) unlink(paste(ModisName))
                        if (!e) warning (paste("The procedure has failed to REPROJECT the individual HDF image ",ModisName,"!",sep=""))
                      }
                    }
                  } 
                }
              }   
            } else {
              if (proj) {
                for (i in 1:nrow(dHDF)) {
                  ModisName <- dHDF[i,2]
                  pref <- strsplit(ModisName,'\\.')[[1]]
                  pref <- paste(pref[1],"_",pref[3],sep="")
                  d <- dHDF[i,1]
                  date_name <- sub(sub(pattern="\\.", replacement="-", d), pattern="\\.", replacement="-", d)
                  e <- reprojectHDF(ModisName,filename=paste(pref,'_',date_name,'.tif',sep=''),MRTpath=MRTpath,UL=UL,LR=LR,bands_subset=bands_subset,proj_type=proj_type,proj_params=proj_params,utm_zone=utm_zone,pixel_size=pixel_size)
                  if (e & delete) unlink(paste(ModisName))
                  if (!e) warning (paste("The procedure has failed to REPROJECT the individual HDF image ",ModisName,"!",sep=""))
                }
              }
            }
            
          }
)
