
## Specials <- function(f,spec,split2="+",...) {
##   tt <- terms(f,spec)
##   pos <- attributes(tt)$specials[[spec]]
##   if (is.null(pos)) return(NULL)
##   x <- rownames(attributes(tt)$factors)[pos]
##   st <- gsub(" ","",x)
##   res <- unlist(strsplit(st,"[()]"))[2]
##   if (is.null(split2)) return(res)
##   unlist(strsplit(res,"+",fixed=TRUE))
## }

## f <- Surv(lefttime,time,status)~x1+id(~1+z,cluster)
## spec <- "id"
## split1=","
## split2="+"

Specials <- function(f,spec,split1=",",split2=NULL,...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x) ## trim
##  res <- unlist(strsplit(st,"[()]"))
  spec <- unlist(strsplit(st,"[()]"))[[1]]
  res <- substr(st,nchar(spec)+2,nchar(st)-1) 
  if (!is.null(split1))    
    res <- unlist(strsplit(res,split1))
  res <- as.list(res)
  for (i in seq(length(res))) {
    if (length(grep("~",res[[i]]))>0)
      res[[i]] <- as.formula(res[[i]])
  }
  return(res)
##  if (is.null(split2)) return(res)
##  strsplit(res,"+",fixed=TRUE)
}

decomp.specials <- function (x, pattern = "[()]", sep = ",", ...) 
  {
    st <- gsub(" ", "", x)
    if (!is.null(pattern)) 
      st <- rev(unlist(strsplit(st, pattern, ...)))[1]
    unlist(strsplit(st, sep, ...))
  }



## myspecials <- c("id","strata","f")
## f <- Event(leftime,time,cause) ~ id(~1+z+z2,cluster) + strata(~s1+s2) + f(a) + z*x
## ff <- Specials(f,"id",split2=",")

