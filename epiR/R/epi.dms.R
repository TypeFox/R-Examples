"epi.dms" <- function(dat){

   # If matrix is comprised of one column, assume conversion FROM decimal degrees TO degree, minutes, seconds:
   if(dim(dat)[2] == 1){
      dat. <- abs(dat)
      deg <- floor(dat.)
      ms <- (dat. - deg) * 60
      min <- floor(ms)
      sec <- (ms - min) * 60
      rval <- as.matrix(cbind(deg, min, sec), dimnames = NULL)
      id <- dat[,1] < 0
      id <- ifelse(id == TRUE, -1, 1)
      rval[,1] <- rval[,1] * id
      # names(rval) <- c("deg", "min", "sec")
      }

   # If matrix is comprised of two columns, assume conversion is FROM degrees and decimal minutes TO decimal degrees:
   else if(dim(dat)[2] == 2){
      deg <- abs(dat[,1])
      min <- dat[,2] / 60
      rval <- as.matrix(deg + min, dimnames = NULL)
      id <- dat[,1] < 0
      id <- ifelse(id == TRUE, -1, 1)
      rval <- rval * id
      # names(rval) <- "ddeg"
      }

   # If matrix is comprised of three columns, assume conversion FROM degrees, minutes, seconds TO decimal degrees:
   else if(dim(dat)[2] == 3){
      deg <- abs(dat[,1])
      min <- dat[,2] / 60
      sec <- dat[,3] / (60 * 60)
      rval <- as.matrix(deg + min + sec, dimnames = NULL)
      id <- dat[,1] < 0
      id <- ifelse(id == TRUE, -1, 1)
      rval <- rval * id
      # names(rval) <- "ddeg"
      }
   return(rval)    
}