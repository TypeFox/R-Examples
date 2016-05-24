#' @title Write a SaTScan geo file
#' @description Write a SaTScan geo file
#' @details Writes the input data frame to a file in the OS, using the .geo extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .geo extension will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.geo = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  if (dim(x)[2] < 2) stop("Need a data frame with 2 or more columns")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".geo"), 
              row.names=userownames, col.names=FALSE)
} 
# check to see if data.frame is really needed


#' @title Write a SaTScan cas (case) file
#' @description Write a SaTScan cas (case) file
#' @details Writes the input data frame to the OS, using the .cas extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .cas will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.cas = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  if (dim(x)[2] < 2) stop("Need a data frame with 2 or more columns")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".cas"), 
              row.names=userownames, col.names=FALSE)
} 


#' @title Write a SaTScan ctl (control) file
#' @description Write a SaTScan ctl (control) file
#' @details Writes the input data frame to the OS, using the .geo extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .ctl will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export

write.ctl = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  if (dim(x)[2] > 3) stop("Need a data frame with 3 or fewer columns")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".ctl"), 
              row.names=userownames, col.names=FALSE)
} 


#' @title Write a SaTScan nbr (neighbor) file
#' @description Write a SaTScan nbr (neighbor) file
#' @details Writes the input data frame to the OS, using the .nbr extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .nbr will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.nbr = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".nbr"), 
              row.names=userownames, col.names=FALSE)
} 
# This one needs a test, since rows likely have different n of adjacent centroids

#' @title Write a SaTScan met file
#' @description Write a SaTScan met file
#' @details Writes the input data frame to the OS, using the .met extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .met will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.met = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".met"), 
              row.names=userownames, col.names=FALSE)
} 

#' @title Write a SaTScan max file
#' @description Write a SaTScan max file
#' @details Writes the input data frame to the OS, using the .max extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .max will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.max = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".max"), 
              row.names=userownames, col.names=FALSE)
} 
# This one needs a test, since rows likely have different n of adjacent centroids

#' @title Write a SaTScan adj file
#' @description Write a SaTScan adj file
#' @details Writes the input data frame to the OS, using the .adj extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .adj will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.adj = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".adj"), 
              row.names=userownames, col.names=FALSE)
} 

#' @title Write a SaTScan ha (alternative hypothesis) file
#' @description Write a SaTScan ha (alternatove hypothesis) file
#' @details Writes the input data frame to the OS, using the .ha extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .ha will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.ha = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".ha"), 
              row.names=userownames, col.names=FALSE)
} 

#' @title Write a SaTScan pop (population) file
#' @description Write a SaTScan pop (population) file
#' @details Writes the input data frame to the OS, using the .pop extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .pop will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.pop = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".pop"), 
              row.names=userownames, col.names=FALSE)
} 

#' @title Write a SaTScan grd (grid) file
#' @description Write a SaTScan grd (grid) file
#' @details Writes the input data frame to the OS, using the .grd extension.  Contents of the data 
#' frame should be only what you want SaTScan to see.  
#' This is a simple function that calls write.table, since SaTScan just needs ASCII files.
#' @param x Your data frame.
#' @param location Directory location where the file should be written
#' @param filename Name for the output file in the OS; .grd will be added.
#' @param userownames If TRUE, will write the row names into the file.
#' @export
write.grd = function(x, location, filename, userownames = FALSE){
  if (class(x) != "data.frame") stop("Need a data frame")
  utils::write.table(x, quote=F, file = paste0(location,"/",filename,".grd"), 
              row.names=userownames, col.names=FALSE)
} 



