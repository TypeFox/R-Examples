"plotFLIMresid" <- function (multimodel, multitheta, plotoptions, dind) 
{
  if(dev.cur() != 1)
    dev.new()
  model <- multimodel@modellist[[1]]
  if (!identical(model@title, "")) 
    tit <- c(0, 0, 1, 0)
  else tit <- c(0, 0, 0, 0)
  plotrow <- 2
  plotcol <- 2
  par(plotoptions@paropt)
  par(mfrow = c(plotcol, plotrow))
  
  nt <- model@nt
  nl <- model@nl
  x <- model@x
  x2 <- model@x2
  increasing_x2 <- model@x2[2] > model@x2[1]
    
  m <- multimodel@modellist
  t <- multitheta
  resultlist <- multimodel@fit@resultlist
  #require(fields)
  #require(gclus)
  tracemat <- matrix(0, nrow(model@inten), ncol(model@inten))
  residlist <- svdresidlist <- list()
  for (i in 1:length(m)) {
        residuals <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
        for (j in 1:length(resultlist[[i]]@resid)) {
            residuals[, j] <- resultlist[[i]]@resid[[j]]
        }
        svdresidlist[[length(svdresidlist) + 1]] <- doSVD(residuals, 
            2, 2)
        residlist[[length(residlist) + 1]] <- residuals
    }
    ## matplot function with "log" option is not compatible with 
    ## neg. x values; do the below to avoid warning
    xpos<- x
    xpos[which(x<=0)]<-NA
    if (increasing_x2) { 
        limd<- max(  max(residlist[[1]]), abs(min(residlist[[1]]))) 
        image.plot(xpos, x2, residlist[[1]], xlab = plotoptions@xlab, ylab = plotoptions@ylab, 
            main = "Residuals Dataset 1", zlim = c(-limd, limd),
	    col = diverge_hcl(40, h = c(0, 120), c = 60, l = c(45, 90), 
	    power= 1.2))
    }
    if (nt > 1 && nl > 1) {
        matplot(x, svdresidlist[[1]]$left[,i], type = "l", main = "Left sing. vec. residuals ", 
            ylab = "", log = "x", xlab = plotoptions@xlab, col = 1)
        if (length(m) > 1) {
            for (i in 2:length(m)) {
                matlines(m[[i]]@x, svdresidlist[[i]]$left, log = "x", 
                  type = "l", col = i)
            }
        }
        
    }
### get selected indices 
    csave <- rep(-1, nrow(tracemat)*ncol(tracemat))
    csave[model@x2] <- 1 
    dim(csave) <- dim(tracemat)
     c1 <- c2 <- TRUE 
   for(i in 1:ncol(csave)) { 
	 if(c1)
	    if(1 %in% csave[,i]) {	
	        rowstart <- i
		c1 <- FALSE 
            }
	 if(c2)
	    if(1 %in% csave[,ncol(csave)-(i-1)]) {	
		rowend <- ncol(csave)-(i-1) 
		c2 <- FALSE 
	    }
   }
   c1 <- c2 <- TRUE 
   for(i in 1:nrow(csave)) { 
	  if(c1) 
	      if(1 %in% csave[i,]) {	
	          colstart <- i
		  c1 <- FALSE
	      }
	   if(c2)     
	      if(1 %in% csave[nrow(csave)-(i-1),] ) {	
	        colend <- nrow(csave)-(i-1)
	        c2 <- FALSE
	      }
   }

    ############# plot 1 right singular vector as image
   
        
            resmat <- as.vector(tracemat)
            resmat[model@x2] <- as.vector(svdresidlist[[1]]$right[1,]) 
            dim(resmat) <- dim(tracemat)
            image.plot(resmat[colstart:colend, rowstart:rowend], xlab = "", 
	    axes = FALSE, ylab = "", 
	    col = diverge_hcl(40, h = c(0, 120), c = 60, l = c(45, 90), 
	    power = 1.2), 
            main = "Right sing. vec. residuals")
        
    svddatalist <- list()
    for (i in 1:length(m)) {
        svddatalist[[length(svddatalist) + 1]] <- doSVD(multimodel@data[[i]]@psi.df, 
            2, 2)
    }
    if (nt > 1 && nl > 1) {
        plot(1:length(svddatalist[[1]]$values), log10(svddatalist[[1]]$values), 
            ylab = "", xlab = "", main = "Sing. values data", 
            type = "b")
        if (length(m) > 1) {
            for (i in 2:length(m)) {
                lines(1:length(svddatalist[[i]]$values), log10(svddatalist[[i]]$values), 
                  type = "b", col = i)
            }
        }
    }
    if(length(plotoptions@title) != 0){
			tit <- plotoptions@title
			if(plotoptions@addfilename) tit <- paste(tit,m[[i]]@datafile)
    }
    else {
                        tit <- ""
		        if(plotoptions@addfilename) tit <- paste(tit, m[[i]]@datafile)
    }
    mtext(tit, side = 3, outer = TRUE, line = 1)
    par(las = 2)
    
    if (dev.interactive() && length(plotoptions@makeps) != 0) {
      if(plotoptions@output == "pdf")
        pdev <- pdf 
      else  pdev <- postscript
        dev.print(device = pdev, file = paste(plotoptions@makeps,
	"_resids.", plotoptions@output, sep = ""))
    }
}
