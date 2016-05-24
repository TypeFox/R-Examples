 # # # # # # # # # # # # # # # # # #
# # #     Class  definition     # # #
 # # #          DGobj            # # #
  # # # # # # # # # # # # # # # # # 

.DGobj.validity <- function(object){
	if(ncol(object@demographic)!=3){
			stop("[DGobj:validity] Error: the demographic matrix should have 3 columns (2 columns for coordinates and 1 column for the growth variable).")
		return(FALSE)
	} else {}
	if(ncol(object@genetic)<4){
			stop("[DGobj:validity] Error: the genetic matrix should have at least 4 columns (2 columns for coordinates and at least 2 columns for the strain frequencies).")
		return(FALSE)
	} else {}
	if(nrow(object@demographic)<2){
			stop("[DGobj:validity] Error: the demographic matrix should have at least 2 lines.")
		return(FALSE)
	} else {}
	if(nrow(object@genetic)<2){
			stop("[DGobj:validity] Error: the genetic matrix should have at least 2 lines.")
		return(FALSE)
	} else {}
	return(TRUE)
}


setClass(
    Class="DGobj",
    slots=c(
        demographic="matrix",
        genetic="matrix"
   ),
    prototype=prototype(
        demographic=matrix(),
        genetic=matrix()
   ),
    validity=.DGobj.validity
)

### Show ###
.DGobj.show <- function(object){
    cat("   ~~~ Class:",class(object),"~~~ ")
    cat("\n ~ demographic: \n")
    print(object@demographic)
    cat("\n ~ genetic: \n")
    print(object@genetic)
    cat("\n")
    return(invisible())
}
setMethod(f="show",signature="DGobj",definition=.DGobj.show)


### Getteur ###
.DGobj.get <- function(x,i,j,drop){
    switch(EXPR=i,
       "demographic"={return(x@demographic)},
       "genetic"={return(x@genetic)},
       stop("[DGobj:get] ",i," is not a 'DGobj' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="DGobj",definition=.DGobj.get)


### Setteur ###
.DGobj.set <- function(x,i,j,value){
    switch(EXPR=i,
       "demographic"={x@demographic=value},
       "genetic"={x@genetic=value},
       stop("[DGobj:set] Error:",i," is not a 'DGobj' slot")
    )
    validObject(x)
    return(x)
}
setMethod(f="[<-",signature="DGobj",definition=.DGobj.set)


### Summary ###
.DGobj.summary=function(object){
	cat("DG object:\n")
	cat("Number of demographic sample points: ", nrow(object@demographic),"\n")
	cat("Summary of the growth variable: ","\n")
	print(summary(object@demographic[,3]))
	cat("Number of genetic sample points: ", nrow(object@genetic),"\n")
	cat("Number of genetic samples: ", sum(object@genetic[,3:ncol(object@genetic)]),"\n")
	cat("Number of strains: ", ncol(object@genetic)-2,"\n")
	return(invisible())
}
setMethod(f="summary",signature="DGobj",definition=.DGobj.summary)

.DGobj.names=function(x){
	slotNames(x)
}
setMethod(f="names",signature="DGobj",definition=.DGobj.names)
