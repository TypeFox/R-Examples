
mcalc <- function(x,f,dim=NULL,append=FALSE) {
  x <- clean_magpie(x)
  f <- as.formula(f)
  vars <- all.vars(f[[3]])
  
  if(is.null(dim)) dim <- .getDim(vars,x)
  
  for(v in vars) {
    l <- list()
    l[dim] <- v
    tmp <- mselect(x,l,collapseNames=FALSE)
    getNames(tmp,dim=dim) <- all.vars(f[[2]])
    assign(v,tmp) 
  }
  if(append) {
    assign(as.character(as.list(match.call())$x),mbind(x,eval(f[[3]])),envir =  parent.frame()) 
  } else {
    return(eval(f[[3]]))
  }
}
