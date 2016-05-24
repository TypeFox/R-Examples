#<<BEGIN>>
plot.mccut <- function(x, stat = c("median","mean"), lim = c(0.025, 0.25, 0.75, 0.975), griddim = NULL, xlab = names(x), ylab = "Fn(x)", main = "", draw=TRUE, ...)
#ISALIAS plot.mc
#--------------------------------------------
{
  summ <- function(x) c(mean=mean(x,na.rm=TRUE),quantile(x,probs=c(0.5,lim),na.rm=TRUE),Nas=sum(is.na(x)))[c(2,1,3:(length(lim)+3))]


  quel <- which(sapply(x,inherits,what="plot.mccut"))
  lquel <- length(quel)
  if(lquel == 0) stop("plot.mc was not evaluated in evalmccut : no summary to produce")
  if(lquel > 1) stop("More than one plot.mc was evaluated in evalmccut : impossible to produce a summary")

  x <- x[[quel]]

  typen <- sapply(x,attr,"type")

  LEQUANT <- function(x,type){
    if(is.list(x)) return(lapply(x,LEQUANT,type))

    if(type=="0") {
      x <- as.matrix(x[1])
      dimnames(x) = list("NoInc","NoVar")
      return(x)}

    if(type=="V") {
      x <- t(as.matrix(x[,1,]))
      rownames(x) <- "NoUnc"
      return(x)}

    x <- drop(apply(x,c(1,3),summ))

    if(type=="U") {
      x <- as.matrix(x)
      colnames(x) <- "NoVar"}
      
      return(x)
  }

  x <- mapply(LEQUANT,x,typen)
  
  class(x) <- "plotmc"
  plot.mc(x,stat = stat, lim=lim, griddim = griddim, xlab = xlab, ylab = ylab, main = main, draw=draw, ...)
  return(invisible(x))
}

