"decomposeSurv" <- function
(
 formula, 
 data,
 sort=FALSE,
 offset=NULL,
 PTpreset=NA                            # are PT coefs preset?
 )
### decomposes complex survival formula
### with time-factor interactions
### and fp()-parts in formula
###  trans: I(), powM2, powM1, powM0.5, sqrt, pow2, pow3, log,  
### 2008-04
{
	orig.formula <- formula
        
	## expand formula if needed:
        repeat {
                terms <- terms(formula, "fp", data=data)
		fac <- attr(terms, "factors")
                
		needed <- rownames(fac)[!rownames(fac) %in% colnames(fac)][-1]
                
		if(length(needed) == 0) break
                
		formula <- as.formula(paste(as.character(formula)[2], "~", 
                                            as.character(formula)[3], "+",
                                            paste(needed, sep="+")))
	}
        
        ## determine position of <fp> terms (of all orders) in
        ## not expanded formula
        spec <- untangle.specials(terms, "fp", order=1:3)
        fpnames <- gsub("fp\\(([a-zA-Z0-9]*)\\)", "\\1", spec$vars) # c("a", "b")
        fpiden <- paste("(", fpnames, ")", sep="") # c("(a)", "(b)")

        formulaOrig <- formula
        ## replace fp() by the 8x2 possible terms (in formula[3])
        ## attention: fp-variables must consist of a-z,A-Z,0-9, otherwise
        ## must be made modifications here (gsub) ...
        ## FURTHER: 1st term is the linear one
        ## FURTHER: 1st half (=8) are the usual terms, 2nd half (9-16) are the repeated powers
        ## FURTHER: problems possible if one variable name is part of the other, e.g. ABC & ABCDE
        ##CODE <- paste("(I(\\1) + powM2(\\1) + powM1(\\1) + powM0.5(\\1) + log(\\1) + sqrt(\\1) + pow2(\\1) + pow3(\\1) + ",
        ##              "  RI(\\1) + RpowM2(\\1) + RpowM1(\\1) + RpowM0.5(\\1) + Rlog(\\1) + Rsqrt(\\1) + Rpow2(\\1) + Rpow3(\\1) )")
        CODE <- paste("(I(PT(\\1)) + powM2(PT(\\1)) + powM1(PT(\\1)) + powM0.5(PT(\\1)) + log(PT(\\1)) + ",
                      "sqrt(PT(\\1)) + pow2(PT(\\1)) + pow3(PT(\\1)) + ",
                      "RI(PT(\\1)) + RpowM2(PT(\\1)) + RpowM1(PT(\\1)) + RpowM0.5(PT(\\1)) + ",
                      "Rlog(PT(\\1)) + Rsqrt(PT(\\1)) + Rpow2(PT(\\1)) + Rpow3(PT(\\1)) )")                              
        sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", CODE, formulaOrig[3])
        formula <- as.formula(paste(as.character(formula)[2], "~", sub3))
        
        ## define simple transformations
        powM2 <- function(z) z^(-2)
        powM1 <- function(z) z^(-1)
        powM0.5 <- function(z) z^(-0.5)
        pow2 <- function(z) z^2
        pow3 <- function(z) z^3
        
        ## define repeated powers
        RI <- function(z) z * log(z)
        RpowM2 <- function(z) z^(-2) * log(z)
        RpowM1 <- function(z) z^(-1) * log(z)
        RpowM0.5 <- function(z) z^(-0.5) * log(z)
        Rlog <- function(z) log(z) * log(z)
        Rsqrt <- function(z) sqrt(z) * log(z)
        Rpow2 <- function(z) z^2 * log(z)
        Rpow3 <- function(z) z^3 * log(z)

        ## pretransformation function
        PT <- function(z) {
                if(!is.na(PTpreset[1])) {
                        ## *** use preset coefficients
                        varname <- eval(substitute(as.character(quote(z))))
                        colname <- paste("PT(", varname, ")", sep="")
                        (z + PTpreset["shift", colname]) / PTpreset["scale", colname]
                } else {
                        ## *** calc the coefficients
                        obj <- fp.scale(z)
                        (z + obj$shift) / obj$scale
                }
        }
        
	## construct 3-col response:
	resp <- model.extract(model.frame(formula, data = data), "response")
	if(ncol(resp) == 2)
          resp <- cbind(start=rep(0, nrow(resp)), resp)
        
	## sort by STOPtime and -Cens
	if(sort) {
                sort <- order(resp[, 2],  -resp[, 3])
                data <- data[sort, , drop=FALSE]
		resp <- resp[sort, ]
	}
        
	mm <- model.matrix(formula, data = data) ## Model-Matrix
        mm1 <- mm[, -1, drop=FALSE]	# w/o intercept

        ## *** RETREIVE PRETRANS COEFS
        sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", "PT(\\1)", formulaOrig[3])
        formulaPT <- as.formula(paste(as.character(formulaOrig)[2], "~", sub3))
        ## varied pretransformation function to get coefficients:        
        PT <- function(z) {
                obj <- fp.scale(z)
                c(shift=obj$shift, scale=obj$scale, rep(-999, length(z) - 2))
        }
        mmPT <- model.matrix(formulaPT, data = data)
        PTcoefs <- mmPT[1:2, !is.na(mmPT[1,]) & (mmPT[3,]==-999), drop=FALSE]
        rownames(PTcoefs) <- c("shift", "scale")
        
        ## offset ..
        if(length(offset) != 0)
          offset.values <- offset
        else
          offset.values <- NA
        
        terms <- terms(formula, "fp", data=data)
	fac <- attr(terms, "factors")
	labels <- attr(terms, "term.labels")
	
	## splits by special chars
	f <- function(str){
          for(chars in c("(", ")", ":", " ", ",", "*", "^"))
            str <- unlist(strsplit(str, split=chars, fixed=TRUE))
          str
          }
        
	rowSplit <- sapply(rownames(fac), f, simplify=FALSE)	# splitted effects
	stopName <- tail(rowSplit[[1]], 2)[1]	# name of stoptime
	rowInter <- unlist(lapply(rowSplit[-1], function(z) any(z == stopName)))
        ##  rowOffset <- unlist(lapply(rowSplit[-1], function(z) any(z == "offset")))
        ##  rowOffset
        
        
	fac <- fac[-1, , drop=FALSE]	# omit Surv
        
	colSplit <- lapply(colnames(fac), f)
	colInter <- unlist(lapply(colSplit, function(z) any(z == stopName)))
        
	nTimes <- colSums(fac[rowInter, , drop=FALSE])
	nFac   <- colSums(fac[!rowInter, , drop=FALSE])
        
	inters <- (nFac>0) & (nTimes>0)
	NTDE <- sum(inters)
        
	
	timedata <- matrix(0, nrow(data), 0)
	timeind <- c()
        
	## loop for (time x effect)
	for(i in which(inters)) {
		## search pure time:
		ind <- (colSums(fac[rowInter, i] != fac[rowInter, , drop=FALSE]) == 0) & (nFac==0)
		timedata <- cbind(timedata, mm1[, ind, drop=FALSE])
                
		## search pure effect:
		ind <- (colSums(fac[!rowInter, i] != fac[!rowInter, , drop=FALSE]) == 0) & (nTimes == 0)
		timeind <- c(timeind, which(ind[!colInter]))
	}
	mm1 <- mm1[, !colInter, drop=FALSE]
        
	covnames <- c(colnames(mm1),
                      paste(colnames(timedata), colnames(mm1)[timeind], sep=":")
                      )
        alter <- c(colnames(mm1),
                   paste(colnames(mm1)[timeind], colnames(timedata), sep=":")
                   )
        
	## indicator to identify the original formula:
	ind <- covnames %in% colnames(attr(terms(orig.formula, "fp", data=data), "factors")) |
               alter %in% colnames(attr(terms(orig.formula, "fp", data=data), "factors"))


        ## FP indicator matrix (1:16)
        NFP <- length(fpnames)
        fpind <- matrix(0, NFP, length(covnames), dimnames=
                        list(fpnames, covnames))
        for(i in seq(length=NFP)) {
                inds <- grep(fpiden[i], covnames)
                ## when interactions occur, 1:16,1:16 instead of 1:32 is needed
                fpind[i, inds] <- ((seq(along=inds) - 1) %% 16) + 1
        }
        
        ## return object
	list(
             fac=fac,                   # factor matrix ..
             resp=resp,                 # N x 3 - response matrix
             mm1=mm1,                   # model matrix without time effects

             NTDE=NTDE,                 # number time dep. effects
             timedata=timedata, 	# matrix with time functions as columns
             timeind=timeind, 		# indicator of time-dependend effects

             NFP=NFP,                   # number of frac.polys
             fpnames=fpnames,           # names of the variables used as fractional poly's
             fpind=fpind,               # matrix with frac.polyn. in each row, numbered by 1 to 5
             PTcoefs=PTcoefs,           # coefficients of pretransformations
             
             covnames=covnames,         # names of covariates
             ind=ind,			# indicator vector:
                                        # which terms are really part of the formula
             offset.values=offset.values # offset values
             )
}
