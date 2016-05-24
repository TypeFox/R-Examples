
debuglevel <- function(level=NULL) {
	if(is.null(level)){
		return(.Call(TMR_getTraMineRDebugLevel))
	}
	.Call(TMR_setTraMineRDebugLevel,as.integer(level))
	return(level)
	
}
