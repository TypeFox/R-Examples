print.info <- function(info, args.only = F){
##
## ABBREVIATED PRINTING OF info OBJECT. IF args.only = T, SHOWS ONLY
## THOSE VALUES THAT ARE USABLE AS ARGUMENTS TO HAPLIN
##
## NOTE: ASSUMES filename IS THE FIRST ELEMENT, AND ALL OTHER ARGUMENTS
## TO HAPLIN ARE AT THE SECOND "LEVEL" OF LIST
#
## FOR SAFETY'S SAKE
.info <- unclass(info)
#
## NAMES OF ARGUMENTS TO haplin
if(args.only){
	#
	## IF ONLY TO SHOW ARGUMENTS TO A HAPLIN CALL, REMOVE THE OTHERS
	.names <- names(formals(haplin))
	for(i in 2:length(.info)){
		## REMOVE ELEMENTS IN EACH GROUP
		.info[[i]] <- .info[[i]][is.element(names(.info[[i]]), .names)]
	}
	## REMOVE GROUPS NOT RELATED TO ARGUMENTS
	.ind <- sapply(.info, function(x) length(x) > 0)
	.info <- .info[.ind]
}
#
## ABBREVIATED PRINTING
f.printlist(.info)
#
return(invisible(.info))

}
