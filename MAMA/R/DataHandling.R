 mergedata <- function (x, varname) 
{
   
    N <- length(GEDM(x))
    merge <- list()
    for (i in c(1:N)) {
      class.col <- grep(varname,names(clinical(x)[[i]]))
        if (i == 1) {
            merge$dat = GEDM(x)[[i]]
            merge$cl = as.numeric(clinical(x)[[i]][, class.col])
            merge$origin = c(rep(i, length(as.numeric(clinical(x)[[i]][, 
                class.col]))))
        }
        else {
            merge$dat = cbind(merge$dat, GEDM(x)[[i]])
            merge$cl = c(merge$cl, as.numeric(clinical(x)[[i]][, 
                class.col]))
            merge$origin = c(merge$origin, rep(i, length(as.numeric(clinical(x)[[i]][, 
                class.col]))))
        }
    }
    return(merge)
}

selectClass <- function (x, varname, type) {
classList <- list ()
if (type == "factor")
{
for (i in 1:length(clinical(x))) {
  class.col <- grep(varname,names(clinical(x)[[i]]))
  classList[[i]] <- clinical(x)[[i]][, class.col]
}
}
if (type == "binary")
{
for (i in 1:length(clinical(x))) {
  class.col <- grep(varname,names(clinical(x)[[i]]))
  classList[[i]] <- clinical(x)[[i]][, class.col]
  levels(classList[[i]])<-c(1,0)
}
}

names(classList) <- datanames(x)
return(classList)
}


