`is.quantitative` <-
function(formula, data){
  if (!missing(data))
   { 
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    ans<-length(table(unique(mf[,1])))>2
   }
  else
   {
    ff<-as.character(formula[[2]])
    ans<-length(unique(get(ff)))>2
   }
  ans
}

