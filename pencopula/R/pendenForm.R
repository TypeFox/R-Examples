 
pendenForm <- function(penden.env) {
  pf <- get("frame",penden.env)
  dim.data <- dim(get("Y",penden.env))
  p <- dim.data[2]
  if(p<=1) stop("At least two variables are needed!")
  n <- dim.data[1]
  if(!is.null(colnames(get("Y",penden.env)))) {
    names <- colnames(get("Y",penden.env))
    assign("names",names,penden.env)
    for(j in 1:p) assign(names[j],get("Y",penden.env)[,p],penden.env)
  }
  else {
    names <- c()
    for(j in 1:p) {
      names[j] <- paste("y",j,sep="")
      assign(names[j],get("Y",penden.env)[,p],penden.env)
    }
    assign("names",names,penden.env)
  }
  assign("n",n,penden.env)
  assign("p",p,penden.env)
}
