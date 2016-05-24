
#---------------------------------------------------------------------------
agree.sddm <- function(ratings, conf.level=0.95, method=c("vst", "delta"), NAaction=c("fail", "omit")){

    wscv <- agree.wscv(ratings, conf.level, method, NAaction)
    adj <- sqrt(2)*qnorm((1+conf.level)/2)

    list(value=wscv$value*adj, lbound=wscv$lbound*adj, ubound=wscv$ubound*adj)

    
}
