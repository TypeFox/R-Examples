read.grid <- function(file,record_size=4, endian="little"){

  if (!is.na(grep(".asc$", file, ignore.case = TRUE) || grep(".ascii$", file, ignore.case = TRUE) || grep(".txt$", file, ignore.case = TRUE)))
  { #read ASCII
    myhead.orig <- read.table(file, nrows=6)
    myhead <- myhead.orig[,2]
    names(myhead) <- as.character(myhead.orig[,1] )
    tab <- read.table(file, na.strings=-9999, skip=6)
    nodata_value = myhead[6]
    for(i in 1:NROW(tab)){
        tab[i,]<-rev(tab[i,])
    }
  }
  else
  {   #read binary grid
    zz <- file(file, "rb")
    tab=readBin(zz, numeric(),size=record_size,n= 12, endian=endian) #read header
    
    nncols       = tab[1] #retrieve header information
    nnrows       = tab[2]
    xllcorner    = tab[3]
    yllcorner    = tab[4]
    cellsize     = tab[5]
    nodata_value = tab[6]
    myhead=tab[1:6] 
    names(myhead)= c("ncols","nrows","xllcorner","yllcorner","cellsize","nodata_value")
    
    tab=readBin(zz, numeric(),size=record_size,n= nncols*nnrows, endian=endian)  #read actual content
    close(zz)
    dim(tab)=c(nncols,nnrows) #rearrange to matrix
  } 

  tab[tab==nodata_value]=NA #replace na-values

  return(list(head=myhead, tab=tab[,NCOL(tab):1]))
}
