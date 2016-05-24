"padding" <- function(str, space, with, to=c("left","right","center"))
{
  mto <- match.arg(to)
  free  <- space - nchar(str)
  fill   <- substring(paste(rep(with, ceiling(free / nchar(with))), collapse = ""),1,free)
## cat("	free=",free,",  fill=",fill,",  mto=",mto,"\n")
  if(free <= 0)
	  invisible("")
  else if  (mto == "left") paste(str,fill,sep = "")
  else if  (mto == "right") paste(fill,str,sep = "")
  else  paste(substring(fill,1,free %/% 2),str,substring(fill,1+free %/% 2,free), sep = "")
}

