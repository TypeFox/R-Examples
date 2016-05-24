#' Manage histogram tails
#' 
#' Manage histogram tails to remove very far outliers.
#' \code{histuple} is \code{list(hx = hist\$mids, hy = hist\$counts)},
#' which is an internal representation of histogram
#'
#' @param htu list, input histuple
#' @param merge_tails length-two numeric vector, points to be merged for left and right tails
#'
#' @return list, histuple
#'
#' @keywords data
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#' \dontrun{
#' htu2 <- ecd.manage_hist_tails(htu, c(1,2))
#' }
### <======================================================================>
ecd.manage_hist_tails <- function(htu, merge_tails=c(0,0))
{
	if(is.null(htu$hx) | length(htu$hx) == 0){
		stop("Hx in histuple is empty!\n")
	}
	if(is.null(htu$hy) | length(htu$hy) == 0){
		stop("Hy in histuple is empty!\n")
	}
	if(length(merge_tails) != 2){
		stop("Length of merge_tails must be two!\n")
	}

	# handle one tail at a time
	fn_merge_tail <- function(htu, allowed_merge )
	{
	  	hx <- htu$hx
	  	hy <- htu$hy

    	if(!(allowed_merge > 0)) return (htu) # nothing to do!
    	
    	if(is.null(hx) | length(hx) == 0){
        	stop("Hx in histuple is empty!\n")
    	}
    	if(is.null(hy) | length(hy) == 0){
        	stop("Hy in histuple is empty!\n")
    	}

		owed <- 0
		while ( allowed_merge > 0 ) {
			if ( hy[1] > 0 ) {
			  	hy <- hy[-1]
			  	hx <- hx[-1]
			  	allowed_merge <- allowed_merge - 1
			  	owed <- owed + 1
			  	while ( hy[1] == 0 ) {
					hy <- hy[-1]
					hx <- hx[-1]
			  	}
			}
			else break
		}
		hy[1] <- hy[1] + owed
		
		list(hx = hx, hy = hy)
	}
	
	revhtu <- function(htu) {
  		list(hx = rev(htu$hx), hy = rev(htu$hy))
	}
  
  	# manage left tail
    htu2 <- fn_merge_tail(htu, merge_tails[1])

  	# manage right tail
    htu3 <- fn_merge_tail(revhtu(htu2), merge_tails[2])
    revhtu(htu3)
}




### <---------------------------------------------------------------------->
