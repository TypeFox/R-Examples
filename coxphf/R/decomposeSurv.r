### old from version 1.05

"decomposeSurv" <-
function(
 formula,
 data,
 sort=FALSE
)
### decomposes complex survival formula
### with time-factor interactions
### 2006-10
{
	orig.formula <- formula

	## expand formula if needed:
	repeat {
		terms <- terms(formula, data=data)
		fac <- attr(terms, "factors")
		needed <- rownames(fac)[!rownames(fac) %in% colnames(fac)][-1]

		if(length(needed) == 0) break

		formula <- as.formula(paste(as.character(formula)[2], "~",
							as.character(formula)[3], "+",
							paste(needed, sep="+")))
	}

	## construct 3-col response:
	resp <- model.extract(model.frame(formula, data = data), "response")
	if(ncol(resp) == 2)
		resp <- cbind(start=rep(0, nrow(resp)), resp)

	## sortieren nach STOPzeit und -Cens
	if(sort) {
	      sort <- order(resp[, 2],  -resp[, 3])
      	data <- data[sort, , drop=FALSE]
		resp <- resp[sort, ]
	}

	mm <- model.matrix(formula, data = data) ## Model-Matrix
	mm1 <- mm[, -1, drop=FALSE]	# w/o intercept

	terms <- terms(formula, data=data)
	fac <- attr(terms, "factors")
	labels <- attr(terms, "term.labels")

	## splittes by special chars
#	f <- function(str)
#		for(chars in c("(", ")", ":", " ", ",", "*", "^"))
#			str <- unlist(strsplit(str, split=chars, fixed=TRUE))
  f<-function(str) {
          for(chars in c("(", ")", ":", " ", ",", "*", "^"))
            str <- unlist(strsplit(str, split=chars, fixed=TRUE))
          str
          }


	rowSplit <- sapply(rownames(fac), f)	# splitted effects
	stopName <- tail(rowSplit[[1]], 2)[1]	# name of stoptime
	rowInter <- unlist(lapply(rowSplit[-1], function(z) any(z == stopName)))

	fac <- fac[-1, , drop=FALSE]	# omit Surv

	colSplit <- sapply(colnames(fac), f)
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

	## indicator to identify the original formula:
	ind <- covnames %in% colnames(attr(terms(orig.formula, data=data), "factors"))

	list(NTDE=NTDE, 			# number time dep. effects
		fac=fac, 			# factor matrix ..
		resp=resp, 			# N x 3 - response matrix
		mm1=mm1, 			# model matrix without time effects
		timedata=timedata, 	# matrix with time functions as columns
		timeind=timeind, 		# indicator of time-dependend effect
		covnames=covnames,	# names of covariates
		ind=ind			# indicator if some terms of not in formula
	)
}
