## rename to split.zoo to make generic method for zoo objects

##' split a zoo object
##'
##' @param x data
##' @param f a factor
##' @return the split data
##'
##' @export
"Split.zoo" <-
  function(x,f) {
    splitUnivariate = function(x,f) {
      tmp = x[f == values[1]]
      for(i in 2:length(values)) {
        tmp = merge(tmp,x[f == values[i]])
      }
      colnames(tmp) = values
      return(tmp)
    }
    size = dim(as.matrix(zoo::coredata(x)))
    values = unique(f)
    if(length(zoo::index(x)) != length(f)) {
      warning("Length of grouping variable modified to make commensurate with zoo object")
      f = rep(f, length.out=length(zoo::index(x)))
    }
    
    if(size[2] == 1) {
      return(splitUnivariate(x,f))
    } else {
      tmp = list()
      theNames =colnames(x)
      for(i in 1:size[2]) {
        tmp[[i]] = splitUnivariate(x[,i],f)
      }
      names(tmp) = theNames
      return(tmp)
    }
  }

