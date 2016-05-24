"readData" <-
function (filenm, typ = "", sep = "") 
{   
    inten <- matrix()
    if(typ == "plain") {
	   psisim <- as.matrix(read.table(filenm))
	   x <- as.numeric(rownames(psisim))
	   nt <- length(x) 
	   x2 <- as.numeric(scan(filenm, nlines=1))
	   nl <- length(x2)
    }
    else { 
    typ <- scan(filenm, sep = sep, skip = 2, what = "", nmax = 1)

    # taking out the NetCDF support since we have to work out
    # the dependency issue
    #filespl <- strsplit(filenm, split=".",fixed=TRUE)[[1]]
    #endn <- filespl[length(filespl)]
    #if(endn =="cdf" || endn == "CDF") {
    #  newdat <- readNetCDF_TIMP(filenm)
    #  return(newdat)
    #} 
    if (typ[1] == "FLIM") {
        file.str <- scan(filenm, sep = sep, skip = 3, what = "")
	origheight <- as.integer(file.str[1])
        origwidth <- as.integer(file.str[2])
        nt <- as.integer(file.str[3])
        x <- vector(mode = "numeric", length = nt)
        nl <- as.integer(file.str[4])
        x2 <- vector(mode = "numeric", length = nl)
        psisim <- matrix(nrow = nt, ncol = nl)
        for (i in 1:nt) {
            x[i] <- as.numeric(file.str[i + 4])
        }
        cnt <- nt + 5
        for (i in 1:nl) {
            for (j in 1:(nt + 1)) {
                if (j == 1) {
                  x2[i] <- as.numeric(file.str[cnt])
                  cnt <- cnt + 1
                }
                else {
                  psisim[(j - 1), i] <- as.numeric(file.str[cnt])
                  cnt <- cnt + 1
                }
            }
        }
        cnt <- cnt + 2
        inten <- matrix(nrow = origheight, ncol = origwidth)
        for (i in 1:origwidth) {
            for (j in 1:origheight) {
                inten[j, i] <- as.numeric(file.str[cnt])
                cnt <- cnt + 1
            }
        }
    }    
    if ((typ[1] == "time") || (typ[1] == "Time")) {
        file.str <- scan(filenm, sep = sep, skip = 3, what = "")
        nt <- as.integer(file.str[2])
        x <- as.numeric(file.str[3:(2 + nt)])
        ioff <- nt + 2
        nl <- as.integer((length(file.str) - ioff)/(nt + 1))
        x2 <- vector(mode = "numeric", length = nl)
        psisim <- matrix(nrow = nt, ncol = nl)
        for (i in 1:nl) {
            x2[i] <- as.numeric(file.str[ioff + 1])
            psisim[, i] <- as.numeric(file.str[(ioff + 2):(ioff + 
                nt + 1)])
            ioff <- ioff + nt + 1
        }
	if (x2[1] > x2[nl]) {
	 psi2 <- psisim 
	 for(j in 1:nrow(psi2)) 
	    psi2[j,] <- rev( psi2[j,])
	 psisim <- psi2 
	 x2 <- rev(x2)
	}
    }
    if ((typ[1] == "wavelength") || (typ[1] == "Wavelength")) {
        file.str <- scan(filenm, sep = sep, skip = 3, what = "")
        nl <- as.integer(file.str[2])
        x2 <- as.numeric(file.str[3:(2 + nl)])
        ioff <- nl + 2
        nt <- as.integer((length(file.str) - ioff)/(nl + 1))
        x <- vector(mode = "numeric", length = nt)
        psisim <- matrix(nrow = nt, ncol = nl)
        for (i in 1:nt) {
            x[i] <- as.numeric(file.str[ioff + 1])
            psisim[i, ] <- as.numeric(file.str[(ioff + 2):(ioff + 
                nl + 1)])
            ioff <- ioff + nl + 1
        }
	if (x2[1] > x2[nl]) {
	 psi2 <- psisim 
	 for(j in 1:nrow(psi2)) 
	    psi2[j,] <- rev( psi2[j,])
	 psisim <- psi2 
	 x2 <- rev(x2)
	}
    }
    if (typ[1] == "MG") {
	bd <- read.table(filenm, skip=4, header=FALSE, sep=";")  
	x <- na.omit(as.numeric(read.table(filenm, skip=3, header=FALSE, sep=";")[1,2:ncol(bd)]))
	attributes(x) <- NULL 
	x2 <- na.omit(bd[2:nrow(bd),1])
	attributes(x2) <- NULL
	nl <- as.integer(length(x2))
	nt <- as.integer(length(x))
	psisim <- t(as.matrix(bd[2:(nl+1),2:(nt+1)]))	
    }
   }
    dat(psi.df = psisim, x = x, nt = nt, x2 = x2, nl = nl, simdata = FALSE, 
        datCall = list(match.call()), inten = inten, datafile = filenm)
}
