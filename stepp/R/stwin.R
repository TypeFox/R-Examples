#################################################################
#
# stwin.R
#
################
# stepp window #
################
setClass("stwin",
	   representation(type   = "character",	# stepp window type
				r1     = "numeric",   	# largest number of patients in common
								# among consecutive subpopulations
				r2     = "numeric",
				nevent = "numeric"),
	   prototype(type="sliding", r1=5, r2=20, nevent=-1)
	   )

setValidity("stwin", 
	function(object){
	  status <- FALSE
	  err1   <- FALSE
	  err2   <- FALSE
	  if (is.na(match(object@type, c("sliding")))) {
		errmsg1 <- "invalid stepp window type:"
		errmsg1 <- paste(errmsg1, object@type)
		err1    <- TRUE
	  }
	  if (object@r1 >= object@r2){ 
		errmsg2 <- "minpatspop (r1) MUST be less than patspop (r2)"
		err2    <- TRUE
	  }
	  else status <- TRUE

	  if (err1) print(errmsg1)
	  if (err2) print(errmsg2)
	  return(status)
	}
)

setMethod("summary", signature="stwin",
	definition=function(object){
	  write(paste("        Number of patients per subpopulation (patspop r2):", object@r2),file="")
    	  write(paste("Largest number of patients in common among consecutive subpopulations(minpatspop r1):", object@r1),file="")
	}
)

# constructor function for stepp window
stepp.win <- function(type, r1, r2){
	sw <- new("stwin", type=type, r1=r1, r2=r2)
	return (sw)
}
