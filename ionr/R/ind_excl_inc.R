#' Incrementally calls ind_excl_step
#'
#' @description See \code{\link{ind_excl}} for details.
#' @inheritParams ind_excl
#' @return Provides the results of a single step in indicator exclusion procedure. See example for details
#' @encoding utf-8

#' @examples

#' ## Create a scale-outcome set that violates ION. Only 2 last indicators out of 8
#' ## relate to the outcome, the others just relate to the 2 indicators
#' set.seed(466)
#' a<-scale_sim(n=2500, to_n=2, tn_n=6)
#' # run the exclusion procedure. Pcrit taken from Table 2 in Vainik et al., 2015,
#' # European Journal of Personality
#' res=ind_excl_inc(a[[1]],a[[2]], pcrit=0.0037)
#' # which indicators does the procedure exclude?
#' res
#'
#' @export


## for debugging

# set.seed(466) a<-scale_sim(n=2500, to_n=2, tn_n=6) indicators=a[[1]] outcome=a[[2]]
# #indicatornames=1:ncol(indicators) pcrit=0.0037 verbose=TRUE coruse='everything'

ind_excl_inc <- function(indicators, outcome, indicatornames = 1:ncol(indicators), pcrit = 0.05, verbose = F, coruse = "everything") {
    
    # run ind_excl_step function as long as something can be excluded, ie the smallest p-value is smaller than the
    # criterion
    
    tempcrit <- pcrit
    exclude <- vector()
    while (tempcrit <= pcrit) {
        # 16.01.08 the comparison value used to be 0.05. but as somewhone might want to operate with 0.07 etc, it is now
        # set to p crit
        iexc <- ind_excl_step(indicators, outcome, indicatornames, exclude, coruse)
        if (verbose == T) 
            print(iexc)
        flush.console()
        if (iexc[1, 2] > pcrit) {
            # break loop, if the first p-value is above the p-criterion. this means that no more exclusions will take place.
            break
        } else {
            # tempcrit=iexc[1,2] # his is done just that tempcrit would be smaller than pcrit. not really needed, commented
            # out on 16.01.08.
            exclude <- c(exclude, rownames(iexc)[1])
            
        }
    }
    
    return(exclude)
    
    
} 
