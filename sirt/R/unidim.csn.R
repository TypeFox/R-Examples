
#############################################################
unidim.test.csn <- function( dat , RR=400 , prop.perm=.75 ,
	progress=TRUE){

	dat <- as.matrix(dat)
	dat <- stats::na.omit(dat)
	#****
	score <- rowSums(dat)
	dat <- dat[ order(score) , ]
	score <- rowSums(dat)
	ds <- diff(score)
	ind <- c(1,which( ds > 0 )+1 )
	score_index <- cbind( ind , c( ind[-1]-1 , nrow(dat) ) )
	score_index <- cbind( score_index , score_index[,2] - score_index[,1] + 1)
	N <- nrow(dat)
	dat_perm <- matrix(NA , N , ncol(dat) )
	SS <- nrow(score_index)
	for (ss in 1:SS){
	for (ii in 1:ncol(dat)){
		iss <- score_index[ss,1]:score_index[ss,2]
		dat_perm[iss,ii] <- base::sample( dat[iss,ii] )
					}
			}
	score_index <- cbind( score_index , round( score_index[,3] * prop.perm ) )
	progress <- 1*progress
	progress_vec <- c(0,which( diff( floor( 10 * ( 1:RR ) / ( RR+1 ) ) ) == 1 ) )
	# perform test
#	res <- gooijer_csn_table( dat , dat_perm , RR=RR , NS=0 , progress , progress_vec ,
#				as.matrix(score_index) )
# extern "C" {
# SEXP gooijer_csn_table( SEXP dat_, SEXP dat_perm_, 
# SEXP RR_, SEXP NS_, SEXP progress_, SEXP progress_vec_, 
# SEXP score_index_) ;
# }
    if (progress==1){    cat("|**********|\n")  }
	res <- .Call( "gooijer_csn_table" , 
		        dat , dat_perm , RR=RR , NS=0 , progress , progress_vec ,
				as.matrix(score_index) ,
				package="sirt" )
	res$p <- mean(res$stat < res$stat_perm)
	res$H0_quantiles <- stats::quantile(res$stat_perm , c(.05 , .25 , .50 , .75 , .95 , .99 , .999 ) )
	if ( progress == 1 ){
		cat( paste0("CSN Statistic = " , round( res$stat , 5) ) )
		cat( " , p =" , round( res$p , 5 ) , "\n")		
		}
	return(res)
	}
#############################################################