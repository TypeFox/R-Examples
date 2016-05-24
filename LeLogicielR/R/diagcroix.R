diagcroix <- function(x,titre,col,family="HersheyScript") {
  if (missing(titre)) {
    titre <- paste("Diagramme en croix",paste("de la variable",
                                              deparse(substitute(x))),sep="\n")
  }
  lx <- levels(x)
  n <- if (is.null(lx)) length(unique(x)) else length(lx)
  if (missing(col)) col <- 1:n
  ux <- if (is.factor(x)) mapLevels(x) else unique(x) # vient du package gdata
  plot(rep(as.numeric(ux),table(x)),sequence(table(x)),pch="-",axes=F,
       xlab="Modalites",xlim=c(0,n+1),ylab="Effectifs",family="HersheyScript",
       col=rep(col,table(x)))
  title(titre,family=family)
  at <- if(is.null(lx)) as.numeric(unique(x)) else 1:n
  for (i in 1:n) axis(side=1,at=at[i],labels=lx[i],col.axis=col[i])
  axis(side=1,at=at,labels=rep("",n))
  axis(side=2)
  fleches(x=F,y=T)
}

