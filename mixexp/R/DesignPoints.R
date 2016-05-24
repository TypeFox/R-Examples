DesignPoints = function (des = NULL,nmxcmp=3, x = NULL, y = NULL, z = NULL,x1lower=0,x1upper=0,
                         x2lower=0, x2upper=0,x3lower=0,x3upper=0,
                         cornerlabs = c("x3","x2","x1"),
                         axislabs=c("x1","x2","x3"),pseudo=FALSE)
{
  if (nmxcmp > 3 ){ stop("DesignPonts function only works for designs with three mixture components")}
  
  for (i in 1:3) {
    axislabs[i]<-paste("Fraction ",axislabs[i])
  }  
  check1 <- is.null(des)
  check2 <- is.null(x)
  if (check1 & check2) {
    design = FALSE
  }
  else {
    design = TRUE
  }
  if (check2) {
    x <- c(0, 0, 1)
    y <- c(0, 1, 0)
    z <- c(1, 0, 0)
  }
  if (check1) {
  }
  else {
    if (ncol(des) !=3) {cat("Warning: the design matrix has more than three columns; the DesignPoints function ","\n")
                        cat("only plots design points for designs with three mixture components. Component x1 is","\n")
                        cat("assumed to to be the first column of the design, x2 the second and x3 the third. Other","\n")
                        cat("columns are ignored. Use cornerlabs and axislabs to change variable names in the plot.","\n")}
    x<-des[ ,3]
    y<-des[ ,2]
    z<-des[ ,1]
  }
  w <- runif(length(x))

cls<-c(rep(0,6))
cls[1]<-min(z)
cls[2]<-max(z)
cls[3]<-min(y)
cls[4]<-max(y)
cls[5]<-min(x)
cls[6]<-max(x)
if (max(cls[1],cls[3],cls[5])>0 | min(cls[2],cls[4],cls[6])<1) {constraints= TRUE} else {constraints=FALSE}

  MixturePlot(x, y, z, w, x3lab = axislabs[3], x2lab = axislabs[2], 
              x1lab = axislabs[1], corner.labs = cornerlabs, 
              lims = cls, constrts = constraints, contrs = FALSE, cols = FALSE, 
              mod = 1, n.breaks = 9, despts = design,pseudo=pseudo)
}