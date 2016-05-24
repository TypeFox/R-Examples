sim = function(
  n=10, m = 250,
  statistic = "mean",
  family = "norm",...) {

  

  ## fix args, make a string of name=value
  args = list(...)
  if(length(args) >0) {
    nms = names(args)
    args = paste(names(args),"=",args,sep="")
    args = paste(args,sep="",collapse=",")
  } else {
    args = NULL
  }

  myPaste = function(...) paste(...,sep="",collapse="")
  
  f = function(n,statistic, family,...) {
    ## a single sample
    doThis = myPaste("r",family,"(",n,",",args,")")
    if(is.null(args)) {
      doThis = myPaste("r",family,"(",n,")")
    }
    x = eval(parse(text=doThis))

    return(do.call(statistic,list(x)))
  }

  sapply(1:m,function(x) f(n,statistic,family,...))
}
