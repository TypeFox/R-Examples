makedistn <- function(distn){
    pname <- paste("p", distn[1], sep="")
    x <- paste("function(x, pm, pn=NULL, log=FALSE)
         do.call(\"", pname, "\", c(list(q=x), pm, pn,", sep="")
    if (distn[1]=="glm") x <- paste(x, " list(family=\"", distn[2],
         "\", link=\"", distn[3], "\"),", sep="")
    eval(parse(text=paste(x, " list(log.p=log)))", sep="")))
}
