setOldClass("zoo")

setGeneric(
  name = 'phenoAmp',
  def = function(x, ...)
    standardGeneric("phenoAmp")
)

setMethod(
  f = "phenoAmp",
  signature = "ts",
  definition = function(x, mon.range = c(1, 12)) {

    d1 <- data.frame(yr = floor(time(x)), mon = cycle(x), val =
    	as.numeric(x))
    mons <- mon.range[1]:mon.range[2]
    d2 <- d1[d1$mon %in% mons, ]
    yrs <- unique(d2$yr)
    yrs.ok <- table(d2$yr, is.na(d2$val))[, 1] == length(mons)
    
    # range
    diff.range <- function(x) {
      if (sum(!is.na(x)) == 0) {
        return(NA)
      } else {
        diff(range(x, na.rm=TRUE))
      }
    }  
    a1 <- aggregate(d2$val, list(d2$yr), diff.range)
    range <- ifelse(yrs.ok, a1$x, NA)
    
    # range/mean
    a2 <- aggregate(d2$val, list(d2$yr), mean, na.rm=TRUE)
    range.mean <- ifelse(yrs.ok, a1$x/a2$x, NA)
  
    # cv
    a3 <- aggregate(d2$val, list(d2$yr), sd, na.rm=TRUE)
    cv <- ifelse(yrs.ok, a3$x/a2$x, NA)
    
    as.data.frame(cbind(year=yrs, range, range.mean, cv), row.names =
    	length(yrs))
  }
)


setMethod(
  f = "phenoAmp",
  signature = "zoo",
  definition = function(x, mon.range = c(1, 12)) {

    # validate args
    if (!is(index(x), "DateTime"))
      stop('time index must be a DateTime object')
    indexx <- as.Date(index(x))  
      
    d1 <- data.frame(time = indexx, yr = years(indexx), mon =
    	monthNum(indexx), val = as.numeric(x))
    mons <- mon.range[1]:mon.range[2]
    d2 <- d1[d1$mon %in% mons, ]
    yrs <- unique(d2$yr)
    n <- table(d2$yr, is.na(d2$val))[, 1]
    
    # range
    diff.range <- function(x) {
      if (sum(!is.na(x)) == 0) {
        return(NA)
      } else {
        diff(range(x, na.rm=TRUE))
      }
    }  
    a1 <- aggregate(d2$val, list(d2$yr), diff.range)
    range <- a1$x
    
    # range/mean
    a2 <- aggregate(d2$val, list(d2$yr), mean, na.rm=TRUE)
    range.mean <- a1$x/a2$x
  
    # cv
    a3 <- aggregate(d2$val, list(d2$yr), sd, na.rm=TRUE)
    cv <- a3$x/a2$x
    
    as.data.frame(cbind(year=yrs, range, range.mean, cv, n), row.names =
    	length(yrs))
  }
)
