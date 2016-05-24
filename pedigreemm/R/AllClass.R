#' Pedigree class
#'
#'
setClass("pedigree", representation =
	 list(sire = "integer", dam = "integer", label = "character"),
	 validity = function(object) {
	     n <- length(sire <- object@sire)
	     if (length(dam <- object@dam) != n)
		 return("sire and dam slots must be the same length")
	     if (length(object@label) != n)
		 return("'label' slot must have the same length as 'sire' and 'dam'")
	     if(n == 0) return(TRUE)
	     animal <- 1:n
	     snmiss <- !is.na(sire)
	     dnmiss <- !is.na(dam)
	     if (any(sire[snmiss] >= animal[snmiss]) ||
		       any(dam[dnmiss]  >= animal[dnmiss]))
		 return("the sire and dam must precede the offspring")
             if (any(sire[snmiss] < 1 | sire[snmiss] > n) |
                 any(dam[dnmiss] < 1 | dam[dnmiss] > n))
                 return(paste("Non-missing sire or dam must be in [1,",
                              n, "]", sep = ''))
	     TRUE
	 })

setClass("pedigreemm", representation = list(relfac = "list"),
         contains = "merMod")
setClass("glmerpedigreemm", representation = list(resp="glmResp"),
         contains = "pedigreemm")
setClass("lmerpedigreemm", representation = list(resp="lmerResp"),
         contains = "pedigreemm")




