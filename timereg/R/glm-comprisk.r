
prep.glm.comprisk <- function(out,time="time",cause="cause",times,censmod=0,cens.code=0)
{ ## {{{
    ###
    out$id <- 1:nrow(out)
    mm <- c()
    for (h in times)
    {
        i2out  <- prep.comp.risk(out,time=time,cause=cause,times=h,cens.code=cens.code)
        Nt <- (i2out[,time] < h)*(i2out[,cause]==1)
        nocens <- (i2out[,time] < h)
        mm <- rbind(mm,cbind(i2out,Nt,h,nocens))
    }
    return(mm)
} ## }}} 

