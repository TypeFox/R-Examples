"plotFLIM" <- function (multimodel, multitheta, plotoptions) 
{
    m <- multimodel@modellist
    t <- multitheta
    resultlist <- multimodel@fit@resultlist
    if (plotoptions@residplot) 
        plotFLIMresid(multimodel, multitheta, plotoptions)
    if (plotoptions@noFLIMsummary) 
        return()
    if(plotoptions@writecon)
      {
        fn <- ifelse(plotoptions@makeps == "", "flimcon", plotoptions@makeps)
        xx<-getXList(result=list(currModel=multimodel, currTheta=multitheta),
                     group=vector(), 
                     file=fn) 
      }
    for (i in 1:length(m)) {
        k <- t[[i]]@kinpar
        model <- m[[i]]
        nt <- model@nt
        nl <- model@nl
        x <- model@x
        x2 <- model@x2
        irfmu <- vector()
        cohirfmu <- vector()
        irftau <- vector()
        if (dev.cur() != 1) 
            dev.new()
        par(plotoptions@paropt)
        par(mgp = c(2, 1, 0), mar = c(3, 3, 3, 2), oma = c(1, 
            0, 4, 0), cex.main = 0.95, mfrow = c(plotoptions@summaryplotrow, 
            plotoptions@summaryplotcol))

    	if (!(plotoptions@plotpulsefol) && model@cohcol!=0 ) 
        	spec <- as.matrix(getSpecList(multimodel, t)[[i]][, -model@cohcol])
    	else
        	spec <-  as.matrix(getSpecList(multimodel, t)[[i]])


        write.table(spec, file=paste(plotoptions@makeps,
                            "_spec_dataset_", i, ".txt", sep=""),
                    row.names = m[[i]]@x2,
                    quote=FALSE) 
##===========plotHistAmp==============

    	for (j in 1:ncol(spec)) {
        	hist(spec[, j], xlab = paste("tau=", signif(1/k[j], 5)), 
            	     main = paste("Comp.", j, "amplitude"))
    	}
##===========plotHistNormComp=============

	if (ncol(spec) > 1) {
        	sumAv <- matrix(0, ncol(spec), nrow(spec))
        	sumspec <- rowSums(spec)
        	for (j in 1:ncol(spec)) {
            		sumAv[j, ] <- spec[, j]/sumspec 
                	hist(sumAv[j, ], xlab = paste("mean", signif(mean(spec[, 
                  	     j])/mean(sumspec), 5)), main = paste("Norm. Component", j))
        	}
    	}

##===========plotIntenImage=============
        plotIntenImage(multimodel, multitheta, i)
##===========plotSelectedPix=============
        plotSelIntenImage(multimodel, multitheta, i)

##===============<tau>==============
	selpixmat <- matrix(0, nrow(model@inten), ncol(model@inten))
## find out the rows where the selected pixels begin and end
	selpixmat[x2]<-1
	j<-1
	while(!(1 %in% selpixmat[,j])) ##for safety reason need to add condition about j<ncol(selpixmat)
		j<-j+1
	colstart <- j
	j<-ncol(selpixmat)
	while(!(1 %in% selpixmat[,j])) ##for safety reason need to add condition about j>1
		j<-j-1
	colend<-j
	j<-1
	while(!(1 %in% selpixmat[j,])) ##for safety reason need to add condition about j<nrow(selpixmat)
		j<-j+1
	rowstart <- j
	j<-nrow(selpixmat)
	while(!(1 %in% selpixmat[j,])) ##for safety reason need to add condition about j>1
		j<-j-1
	rowend<-j
	resmat <- matrix(0, nrow(model@inten), ncol(model@inten))
	if (length(k) > 1) {
		xmat <- matrix(nrow = length(k), ncol = model@nl)
		vecres <- vector()
		for (j in length(k):1) xmat[j, ] <- sumAv[j, ] * (1/k[j])
			vecR <- colSums(xmat)
##============plotHistof<tau>==================
		hist(vecR, xlab = paste("mean", signif(mean(vecR), 5)), 
		     main = "Hist <tau>")
##============plotimage<tau>===================		
		resmat[x2] <- vecR

		zmin <- min(vecR) - (.10* min(vecR))
		zmax <- max(vecR) 

		if(length( plotoptions@ylimspec ) == 2)
			zlim <- plotoptions@ylimspec
		else
			zlim <- c(zmin, zmax)
 
		if(length(plotoptions@imagepal)==0)
			colp <- tim.colors()
		else
			colp <- plotoptions@imagepal
		image.plot(resmat[rowstart:rowend,colstart:colend],
		            xlab = "", axes = FALSE, ylab = "", 
		            main = "<tau>", zlim = zlim, col=colp)
	
	}
##===========plotAmplImages===================	
	if (ncol(spec) > 1) {
        	for (j in 1:ncol(spec)) {
	            	if (length(plotoptions@ylimcomp) == (ncol(spec)*2)) {
		    		zmin <- plotoptions@ylimcomp[j*2-1]
	                	zmax <- plotoptions@ylimcomp[j*2]
		    	}	
		    	else {	
		    		zmin <- 0
	                	zmax <- 1
		    	}	
			resmat[x2] <- sumAv[j, ]
			image.plot(resmat[rowstart:rowend,colstart:colend], 
		         	   ylab = "", axes = FALSE,xlab = paste("tau=",signif(1/k[j], 5)), 
		                   main = paste("Comp.", j, "norm. amp."),  zlim=c(zmin,zmax))
		}
	}	

        residlist <- svdresidlist <- list()
        residuals <- matrix(nrow = nt, ncol = nl)
        for (j in 1:length(resultlist[[i]]@resid)) {
            residuals[, j] <- resultlist[[i]]@resid[[j]]
        }
        svdresidlist[[length(svdresidlist) + 1]] <- doSVD(residuals, 2, 2)
        residlist[[length(residlist) + 1]] <- residuals
        limd <- max(max(residlist[[1]]), abs(min(residlist[[1]])))
        if (plotoptions@FLIMresidimag) {
            if (!(any(diff(x) <= 0) || any(diff(x2) <= 0))) 
                image.plot(x, x2, residlist[[1]], xlab = plotoptions@xlab, 
                           ylab = plotoptions@ylab, main = paste("Residuals Dataset",i), 
                           zlim = c(-limd, limd), col = two.colors())

##diverge_hcl(40,h = c(0, 120), c = 60, l = c(45, 90), power = 1.2))

        }

        xpos <- x
        xpos[which(x <= 0)] <- NA
        if (nt > 1 && nl > 1) {
            matplot(xpos, svdresidlist[[1]]$left[, 1], type = "l", 
                main = "Left sing. vec. residuals ", ylab = "", 
                log = "x", xlab = plotoptions@xlab, col = 1)
        }

        resmat[x2] <- as.vector(svdresidlist[[1]]$right[1,])

        limd <- max(max(svdresidlist[[1]]$right[1,]), abs(min(svdresidlist[[1]]$right[1,])))

        image.plot(resmat[rowstart:rowend, colstart:colend], 
            xlab = "",zlim = c(-limd, limd), col = two.colors(),
            axes = FALSE, ylab = "", 
            main = "Right sing. vec. residuals")

##diverge_hcl(40, h = c(0, 120), c = 60,l = c(45, 90), power = 1.2), 

        svddatalist <- list()
        svddatalist[[length(svddatalist) + 1]] <- doSVD(multimodel@data[[i]]@psi.df,2, 2)
        if (nt > 1 && nl > 1) {
            plot(1:length(svddatalist[[1]]$values), log10(svddatalist[[1]]$values), 
                ylab = "", xlab = "", main = "Sing. values data", 
                type = "b")
        }
        if (length(plotoptions@title) != 0) {
            if (length(m) > 1) 
                tit <- paste(plotoptions@title, ", dataset ", 
                  i, sep = "")
            else tit <- plotoptions@title
            if (plotoptions@addfilename) 
                tit <- paste(tit, m[[i]]@datafile)
            mtext(tit, side = 3, outer = TRUE, line = 1)
        }
        if (dev.interactive() && length(plotoptions@makeps) != 
            0) {
            if (plotoptions@output == "pdf") 
                pdev <- pdf
            else pdev <- postscript
            dev.print(device = pdev, file = paste(plotoptions@makeps, 
                "dataset_", i, "_summary.", plotoptions@output, 
                sep = ""))
        }
    }
    writeEst(multimodel, multitheta, plotoptions)
    # TODO: replace functionality provided by displayEst
    displayEst(plotoptions)
}

##===========new version=========================
"plotHistAmp" <- function (multimodel, t, i = 1) 
{
    k <- t[[i]]@kinpar
    spectralist <- getSpecList(multimodel, t)
    m <- multimodel@modellist

    if ((m[[i]]@cohcol != 0)[1]) 
        spec <- as.matrix(spectralist[[i]][, -m[[i]]@cohcol])
    else
        spec <- spectralist[[i]]

    for (j in 1:ncol(spec)) {
        hist(spec[, j], xlab = paste("tau=", signif(1/k[j], 5)), 
            main = paste("Comp.", j, "amplitude"))
    }
}
##===========new version=========================
"plotHistNormComp" <-function (multimodel, t, i = 1) 
{
    k <- t[[i]]@kinpar
    model <- multimodel@modellist[[i]]
    spectralist <- getSpecList(multimodel, t)

    if( (model@cohcol != 0)[1])
	spec <- as.matrix(spectralist[[i]][,-model@cohcol])
    else
        spec <- spectralist[[i]]
    
    if (ncol(spec) > 1) {
        sumAv <- matrix(0, ncol(spec), nrow(spec))
        sumspec <- rowSums(spec)
        for (j in 1:ncol(spec)) {
            sumAv[j, ] <- spec[, j]/sumspec 
                hist(sumAv[j, ], xlab = paste("mean", signif(mean(spec[, 
                  j])/mean(sumspec), 5)), main = paste("Norm. Component", 
                  j))
        }
    }
}
"plotIntenImage" <- function (multimodel, t, i=1,
                              tit=c("Intensity Image")) 
{
  model <- multimodel@modellist[[i]]
  colvec <- gray(seq(from = 0, to = 1, length = 100))
  image.plot(model@inten, col = colvec, xlab = "", axes = FALSE, 
             ylab = "", main = tit)
}
"plotSelIntenImage" <- function (multimodel, t, i=1,
                                 tit=c("Region of Interest"), cex=1) { 
  model <- multimodel@modellist[[i]]
  colvec <- gray(seq(from = 0, to = 1, length = 100))
  tt <- (model@inten - min(model@inten))/max(model@inten - min(model@inten))
  colmat <- gray(tt)
  colmat[model@x2] <- "#0000FF"
  dim(colmat) <- dim(model@inten)
  csave <- colmat
  colmat <- t(colmat)
  for (j in 1:(ncol(colmat))) 
    colmat[, j] <- rev(colmat[, j])
  plotcolors(colmat, main = tit, cex=cex)
}

"plotTau" <- function(multimodel, t, i=1, tit=" < tau > ", plotoptions=kinopt(), lifetimes = TRUE) {
	model <- multimodel@modellist[[i]]
        k <- t[[i]]@kinpar
        nt <- model@nt
        nl <- model@nl
        x <- model@x
        x2 <- model@x2
	spec <- getSpecList(multimodel, t)[[i]]

	selpixmat <- matrix(0, nrow(model@inten), ncol(model@inten))
## find out the rows where the selected pixels begin and end
	selpixmat[x2]<-1
	j<-1
	while(!(1 %in% selpixmat[,j])) ##for safety reason need to add condition about j<ncol(selpixmat)
		j<-j+1
	colstart <- j
	j<-ncol(selpixmat)
	while(!(1 %in% selpixmat[,j])) ##for safety reason need to add condition about j>1
		j<-j-1
	colend<-j
	j<-1
	while(!(1 %in% selpixmat[j,])) ##for safety reason need to add condition about j<nrow(selpixmat)
		j<-j+1
	rowstart <- j
	j<-nrow(selpixmat)
	while(!(1 %in% selpixmat[j,])) ##for safety reason need to add condition about j>1
		j<-j-1
	rowend<-j
        resmat <- matrix(0, nrow(model@inten), ncol(model@inten))
	if (length(k) > 1) {
		sumAv <- matrix(0, ncol(spec), nrow(spec))
		sumspec <- rowSums(spec) 
		for (j in ncol(spec):1) 
			sumAv[j, ] <- spec[,j]/sumspec

		xmat <- matrix(nrow = length(k), ncol = model@nl)
		vecres <- vector()
		for (j in length(k):1){
			kk<-if (lifetimes) (1/k[j]) else k[j]
 			xmat[j, ] <- sumAv[j, ] * kk
		}
		vecR <- colSums(xmat)
	
##============plotimage<tau>===================
		resmat[x2] <- vecR

		zmin <- min(vecR) - (.10* min(vecR))
		zmax <- max(vecR) 

		if(length( plotoptions@ylimspec ) == 2)
			zlim <- plotoptions@ylimspec
		else
			zlim <- c(zmin, zmax)
 
		if(length(plotoptions@imagepal)==0)
			colp <- tim.colors()
		else
			colp <- plotoptions@imagepal
		image.plot(resmat[rowstart:rowend,colstart:colend],
		            xlab = "", axes = FALSE, ylab = "", 
		            main = tit, zlim = zlim, col=colp)
	
	}
}
"plotNormComp" <- function(multimodel, t, i=1)  {

	model <- multimodel@modellist[[i]]
        k <- t[[i]]@kinpar
        nt <- model@nt
        nl <- model@nl
        x <- model@x
        x2 <- model@x2

    	if ((model@cohcol != 0)[1]) 
        	spec <- as.matrix(getSpecList(multimodel, t)[[i]][, -model@cohcol])
    	else
        	spec <-  getSpecList(multimodel, t)[[i]]

	selpixmat <- matrix(0, nrow(model@inten), ncol(model@inten))
## find out the rows where the selected pixels begin and end
	selpixmat[x2]<-1
	j<-1
	while(!(1 %in% selpixmat[,j])) ##for safety reason need to add condition about j<ncol(selpixmat)
		j<-j+1
	colstart <- j
	j<-ncol(selpixmat)
	while(!(1 %in% selpixmat[,j])) ##for safety reason need to add condition about j>1
		j<-j-1
	colend<-j
	j<-1
	while(!(1 %in% selpixmat[j,])) ##for safety reason need to add condition about j<nrow(selpixmat)
		j<-j+1
	rowstart <- j
	j<-nrow(selpixmat)
	while(!(1 %in% selpixmat[j,])) ##for safety reason need to add condition about j>1
		j<-j-1
	rowend<-j

	if (ncol(spec) > 1) {
		sumAv <- matrix(0, ncol(spec), nrow(spec))
		sumspec <- rowSums(spec) 
		for (j in ncol(spec):1) 
			sumAv[j, ] <- spec[,j]/sumspec
		resmat <- matrix(0, nrow(model@inten), ncol(model@inten))
		resmat[x2] <- sumAv[j, ]
        	for (j in 1:ncol(spec)) {
			resmat[x2] <- sumAv[j, ]
			image.plot(resmat[rowstart:rowend,colstart:colend], 
		         	   ylab = "", axes = FALSE,
                                   xlab = paste("tau=",signif(1/k[j], 5)), 
		                   main = paste("Comp.", j, "norm. amp."),
                                   zlim=c(0,1))
		}
	}	
}	

