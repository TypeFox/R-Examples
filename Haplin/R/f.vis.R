"f.vis"<-
function(..., vis = T, fc = F)
{
# BESTE VERSJON?
	.ex.debug <- exists(".debug")
	.ex.vis <- !missing(vis)
	if(.ex.vis) .vis <- vis	# vis OVERSTYRER .debug
	if(!.ex.vis & .ex.debug) .vis <- get(".debug")
	# .debug SLAAR INN HVIS vis MANGLER
	if(!.ex.vis & !.ex.debug) .vis <- T	# DEFAULT BEHAVIOR
#
#
	.verdier <- list(...)	#
## return(substitute(...))
	.tullf <- function(x)
	deparse(substitute(x))	#
## return(do.call(".tullf", .verdier))
#
## .call <- sys.call()#
	.call <- match.call()	#
	.is.char <- sapply(.call, is.character)
	.navn <- names(.call)
	.pos.rem <- match(c("vis", "fc"), .navn, nomatch = 0)
	.arg <- as.character(.call)[ - c(1, .pos.rem)]
	.is.char <- .is.char[ - c(1, .pos.rem)]	## .arg <- deparse(.call)
#  
## cat("Advarsel: ogsaa problem med tekststreng med komma...\n")
## cat("Advarsel: maa ogsaa kunne ta x...\n")
#
#
	if(.vis) {
#
		if(fc) cat("\nFunction call:", deparse(sys.call(sys.parent())), 
				"\n")
	}
	for(i in seq(along = .arg)) {
		if(.vis) {
			if(.is.char[i])
				cat("\n", .verdier[[i]], "\n")
			else {
				cat("\n", .arg[i], ":\n", sep = "")
				print(.verdier[[i]])
			}
		}
##  else null()#
## eval(.verdier[[i]]) # SORGER FOR AT ARGUMENTENE ALLTID BLIR EVALUERT, UANSETT
	}
#
## return(.arg)
	if(.vis)
		cat("\n")
	invisible()
}
