
####################################################
scan.vector <- function( vec ){ 
    cat(vec , file="ex.data", sep="\n")
    vars <- scan( "ex.data" , what="character", quiet=TRUE)
    file.remove( "ex.data" )
    return(vars)
        }
scan.vec <- scan.vector
####################################################				
# scan function with default what = "character"
scan0 <- function( file="" , ...){
	scan( file=file , what="character" , ...)			
		}
#########################################################		