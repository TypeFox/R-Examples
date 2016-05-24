plot.imageQT <-function (x, cires, unclassval = 0, class = FALSE, QT = TRUE, return = FALSE, qtl = 1, ...){

imsize <- x$imsize
statix <- x$indS

### index sorting function: imix ###
    
imix <- function(six, size) {
        
	sectix <- function(size, ixl, presix = NULL) {
		if (is.null(presix)) {
                	presix <- matrix(1:size^2, nrow = size)
            	}
            	nr <- nrow(presix)
            	ixm <- matrix(0, size/2, size/2)
            	switch(ixl, 
		`0` = {
                	ixm <- presix[1:(nr/2), 1:(nr/2)]
            	}, 
		`1` = {
                	ixm <- presix[(nr/2 + 1):nr, 1:(nr/2)]
            	}, 
		`2` = {
                	ixm <- presix[1:(nr/2), (nr/2 + 1):nr]
            	}, 
		`3` = {
                	ixm <- presix[(nr/2 + 1):nr, (nr/2 + 1):nr]
            	})
        return(ixm)
        }
        
l <- nchar(six)
six <- strsplit(six, "")[[1]]
outix <- NULL
for (i in 1:l) {
	outix <- sectix(size, six[i], presix = outix)
}

return(outix)

}

statmix <- lapply(statix, imix, size = imsize)
nonov <- (length(unique(unlist(statmix))) == length(unlist(statmix)))

if (!nonov) {
        stop("QT decomposition contains overlapping subimages!!\n")
}
immat <- matrix(unclassval, imsize, imsize)

if (class) {
        for (i in 1:length(statmix)) {
            immat[statmix[[i]]] <- cires[i]
        }
}

image(plotmtx(immat), ...)

if (QT) {
        m <- matrix(1:imsize^2, ncol = imsize)
        for (i in 1:length(statmix)) {
        	tmpim <- statmix[[i]]
        	ldim <- nrow(tmpim)
            	
		tl <- rev(which(m == tmpim[1, 1], arr.ind = T) - 1)/imsize
		br <- tl + c(1, 1) * ldim/imsize
		tr<-c(br[1],tl[2])
		bl<-c(tl[1],br[2])
	        x0 <- c(tl[1], tl[1],br[1],br[1])
                y0 <- c(tl[2], tl[2],br[2],br[2])
                x1 <- c(bl[1], tr[1],bl[1],tr[1])
                y1 <- c(bl[2], tr[2],bl[2],tr[2])

                segments(x0, 1-y0, x1,1 -y1, col = qtl)
        }
}

if (return) {
        return(immat)
}

}
