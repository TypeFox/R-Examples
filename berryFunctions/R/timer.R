#' Timer alarm
#' 
#' Beeps in a given interval and gives a progress bar in the console
#' 
#' @details defaults to practice useR lightning talks: 15 slides, each shown 20 secs, change automatically
#' 
#' @return none
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, June 2015
#' @seealso \code{\link{alarm}}, \code{\link{Sys.sleep}}, \code{\link{txtProgressBar}}
#' @references \url{http://user2015.math.aau.dk/lightning_talks}
#' @keywords utilities chron
#' @export
#' @examples
#' 
#' timer(interval=0.5, n=3)
#' timer(interval=0.2, n=8, write=TRUE) # a slight deviation occurs for a large n
#' # timer() # to practice lightning talks at useR! conferences
#'
#' @param interval \code{\link{alarm}} interval in seconds. DEFAULT: 20
#' @param n number of alarm signals to be given. DEFAULT: 15
#' @param write Should the actual estimated time be written for overhead computing time control purposes? DEFAULT: FALSE
#' 
timer <- function(
interval=20,
n=15,
write=FALSE
)
{
begin <- Sys.time()
pb <- txtProgressBar(max=n, style=3)
for(i in 1:n)
  {
  Sys.sleep(interval)
  setTxtProgressBar(pb, i)
  alarm()
  }
close(pb)
time_used <- round(difftime(Sys.time(), begin, units="secs"), 2)
if(write) message("Actual time passed by: ", time_used,
                  " secs. Deviance from target: ", round(time_used-interval*n, 2),
                  " (", round((time_used-interval*n)/(interval*n)*100, 2), "%).")
}

# Examples:
# timer(interval=1, n=3)
# timer(interval=0.2, n=15, write=TRUE)
# timer() # to practice lightning talks at useR! conferences
