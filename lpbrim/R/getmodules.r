#' @title Get modules
#' @description This function takes a mod object, and returns a list of
#' matrices, each corresponding to a single module
#' @export
#' @param mod The output of findModules
getmodules = function(mod)
{
	Lev1 <- mod$S[c(1:NROW(mod$M)),]
	Lev2 <- mod$S[c((NROW(mod$M)+1):NROW(mod$S)),]
	LoW = list()
	for(i in 1:NCOL(mod$S))
	{
		LoW[[i]] <- mod$M[(Lev1[,i]==1),(Lev2[,i]==1)]
	}
	return(LoW)
}
