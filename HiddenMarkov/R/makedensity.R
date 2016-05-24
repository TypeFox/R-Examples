makedensity <- function(distn){
    dname <- paste("d", distn[1], sep="")
    x <- paste("function(x, pm, pn=NULL, log=FALSE)
         do.call(\"", dname, "\", c(list(x=x), pm, pn,", sep="")
    if (distn[1]=="glm") x <- paste(x, " list(family=\"", distn[2],
         "\", link=\"", distn[3], "\"),", sep="")
    eval(parse(text=paste(x, " list(log=log)))", sep="")))
}

