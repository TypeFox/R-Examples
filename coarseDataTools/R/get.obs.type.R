##' Tries to guess the observation types (SIC, DIC, or exact).
##'
##' @param dat a matrix of data, similar to what needs to be passed to \code{dic.fit()}.
##' @return vector of guessed types
##' @export
get.obs.type <- function(dat) {
        type <- rep(0, nrow(dat))
        ## get the single interval censored
        type[dat[,"EL"]==dat[,"ER"]]<-1
        type[dat[,"SL"]==dat[,"SR"]]<-1
        type[dat[,"ER"]>=dat[,"SL"]]<-1

        ## some of those are actually exact!
        type[(dat[,"EL"]==dat[,"ER"]) & (dat[,"SL"]==dat[,"SR"])]<- 2
        return(type)
}
