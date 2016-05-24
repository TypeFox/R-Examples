pairs2 <- function(x,y, mars=c(4,4,0.1,0.1),...)
  {
  nx <- ncol(x)
  ny <- ncol(y)
  namesX <- if(is.null(colnames(x))) paste("x[,",1:nx,"]", sep="") else colnames(x)
  namesY <- if(is.null(colnames(y))) paste("y[,",1:ny,"]", sep="") else colnames(y)
  
  dots <- list(...)
    nmdots <- names(dots)
  oma <- if ("oma" %in% nmdots) 
        dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
        dots$main
    else NULL
    if (is.null(oma)) {
        oma <- c(2, 2, 2, 2)
        if (!is.null(main)) 
            oma[3L] <- 4}
   opar <- par(mfcol = c(ny, nx), mar = mars, oma = oma)
   on.exit(par(opar))
  
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)

  for(i in 1:nx)
     {
     for(j in 1:ny)
        {
         if (i==1)  Ylab=namesY[j] else Ylab <-""  
         if (j==ny) Xlab=namesX[i] else Xlab <-""           
         localPlot(x[, i], y[, j], xlab = Xlab, ylab = Ylab, ...)
        }
     
     
     }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, 1, TRUE, 0.5, cex = cex.main, font = font.main)}
}
