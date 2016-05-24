
###########################################
# Function for analyzing response patterns
md.pattern.sirt <- function(dat){
    dat <- as.matrix(dat)
	if (ncol(dat)>1000){
	   stop("Function only works for datasets with fewer than 1000 variables!\n")
					}
    res <- md_pattern_rcpp( dat_=dat )
    res$dat.ordered <- res$dat[ order( res$resp_patt ) , ]
    return(res)
            }
#*******************************************			
# calling the Rcpp function
md_pattern_rcpp <- function (dat_){ 
	.Call("md_pattern_csource", dat_ , PACKAGE = "sirt")
					}			
#*******************************************					