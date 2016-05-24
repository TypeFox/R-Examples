#### Function supermat as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

supermat <-
function(rbcl = NULL, matk = NULL, trn = NULL){

	if(is.null(rbcl)& is.null(matk)){
	    stop("at least one phylip file for \"rbcl\" or \"matk\" should be provided")
	}
	if(!is.null(rbcl)){
	    rbcl1 <- read.phylip(file.path(getwd(),rbcl))
        rbcldat <- phy2dat(x = rbcl1)
	    rbcldat <- framsub(rbcldat) #Substitute the "-" to "?" on the edges of sequence
	}
	if(!is.null(matk)){
	    matk1 <- read.phylip(file.path(getwd(),matk))
        matkdat <- phy2dat(x = matk1)
	    matkdat <- framsub(matkdat)
	}
	if((!is.null(rbcl))&(!is.null(matk))){
	    name <- union(rbcldat[, 1], matkdat[, 1])
        value <- rep("", length(name))
	    datf <- data.frame(name, value)
        datf <- add.mat(datf, rbcldat)
        #datf <- appendchar(datf, "?")
	    datf <- add.mat(datf, matkdat)
        #datf <- appendchar(datf, "?")
	}
	if((!is.null(rbcl))&(is.null(matk))){
	    datf <- rbcldat
	}
	if((is.null(rbcl))&(!is.null(matk))){
	    datf <- matkdat
	}
	
	if(is.null(trn)){
	    return(datf)
	}
	if(!is.null(trn)){
	    trn2 <- list()
        for(i in 1:length(trn)){
            trn2[[i]] <- read.phylip(file.path(getwd(),trn[i]))
        }
        trndat <- lapply(trn2, phy2dat)
        trnhdat <- lapply(trndat, framsub)
	    nelements <- length(trn2)
        for(i in 1:nelements){
            datf <- add.mat(datf, trnhdat[[i]])
            #datf <- appendchar(datf, "?")
        }
	    return(datf)
	}
}

