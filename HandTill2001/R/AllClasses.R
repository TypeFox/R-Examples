setClass(Class = "cap"
	 , representation = representation(
					   response = "factor"
					   , "VIRTUAL"
					   )
	 , prototype = prototype(response = factor())
	 , validity = function(object){
	 }
	 )

setClass(Class = "bincap"
	 , contains = "cap"
	 , representation = representation(
					   predicted = "numeric"
					   , true = "character"
					   )
	 , prototype = prototype(predicted = numeric(), true = character() )
	 , validity = function(object){
	   if(length(object@response) != length(object@predicted))
	     return("response and predicted must have the same number of observations.")
	   if(any(object@predicted < 0, na.rm = TRUE) || any(object@predicted > 1, na.rm = TRUE))
	     return("probabilities should be in [0,1].")
	   if(length(object@true) > 1)
	     return("give a single character for the 'true'/'presence' class.")
	   if(length(levels(object@response)) > 2)
	     return("response has to be a two class factor.")
	 }
	 )
setClass(Class = "multcap"
	 , contains = "cap"
	 , representation = representation(
					   predicted = "matrix"
					   )
	 , prototype = prototype(predicted = matrix(nrow = 0, ncol = 0))
	 , validity = function(object){
	   p <- object@predicted
	   r <- object@response
	   if(! isTRUE(all.equal(sort(as.character(unique(r))), sort(levels(r)))))
	     warning(paste("found extraneous factor level(s) '"
			   , paste(setdiff(levels(r),as.character(unique(r))), collapse=', ')
			   ,"' of response.\n"
			   , "You may want to work around this by 'response <- factor(response)'."
			   , sep = ''))

	   if(! all(levels(r) %in% dimnames(p)[[2]]))
	     warning(paste("found factor level(s) '"
			   , paste(setdiff(levels(r), dimnames(p)[[2]]), collapse=', ')
			   ,"' of response unmatched by predicted."
			   , sep = ''))
	   if(! all(unique(r) %in% dimnames(p)[[2]], na.rm = TRUE))
	     return(paste("found value(s) '"
			  , setdiff(as.character(unique(r)), dimnames(p)[[2]])
			  , "' of response unmatched by predicted.\n"
			  , "You may want to add column(s) filled with NA to predicted."
			  , sep = '')
	   )
	   if(! all(dimnames(p)[[2]] %in% levels(r)))
	     return(paste("found column(s) '"
			  , paste(setdiff(dimnames(p)[[2]], levels(r)), collapse=', ')
			  ,"' of predicted unmatched by levels(response)."
			  , sep = '')
	   )
	   if(length(r) != nrow(p))
	     return("response and predicted must have the same number of observations.")
	   if(any(p < 0, na.rm = TRUE) || any(p > 1, na.rm = TRUE))
	     return("probabilities should be in [0,1].")
	   if(! isTRUE(all.equal(rep(1,nrow(p)) ,as.numeric(rowSums(p, na.rm = TRUE)))))
	     return("row sums of predicted must be 1.")
	 }
	 )


