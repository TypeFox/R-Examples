law.plot2d <- function(data,density=FALSE,densfunc=dnorm,probfunc=pnorm,tinf=-5,tsup=5) {
  
  nmax <- ncol(data)
  M <- nrow(data)
  
  tt <- tktoplevel()
  nn <- 1
  img <- tkrplot(tt, function() {
    if (density==TRUE) {
      curve(densfunc,xlim=c(tinf,tsup),xlab="",ylab="",main="Convergence in law?",ylim=c(0,0.5))
      par(new=TRUE)
      hist(data[,nn],xlim=c(tinf,tsup),col="red",xlab="t",ylab="F(t)",main="",ylim=c(0,0.5),
           freq=FALSE,breaks=ceiling(nn/10)+1)
    }
    else {
      plot.ecdf(ecdf(data[,nn]),do.points=FALSE,xlim=c(tinf,tsup),col.h="red",xlab="t",ylab=expression(F~"(t) and " ~ hat(F)[n]~"(t)"),main="Convergence in law?")
      curve(probfunc,xlim=c(tinf,tsup),xlab="",ylab="",main="",add=TRUE)
#      par(new=TRUE)
#      legend(x=tinf, y=1, legend=c("F(t)",expression(hat(F)[n](t))), col = c("black","red"),lty=1)
    }				}
                 )
  f <- function(...) {
    n <- as.numeric(tclvalue("nn"))
    if (n != nn) {
      nn <<- n
      tkrreplot(img)
    }
  }
  s <- tkscale(tt, command=f, from=1, to=nmax, variable="nn",
               showvalue=TRUE, resolution=1, orient="horiz")
  mytext <- tklabel(tt,text="n=")
  tkgrid(img,columnspan=2)
  tkgrid(mytext,s,columnspan=1)
  tkgrid.configure(s,sticky="w")
  tkgrid.configure(mytext,sticky="es")

  return(tt)
}
