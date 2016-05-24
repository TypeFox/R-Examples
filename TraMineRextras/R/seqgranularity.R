seqgranularity <- function(seqdata, tspan=3, method="last"){

    metlist <- c("first","last","mostfreq")
    if (!method %in% metlist) {
        stop(" [!] method must be one of: ", paste(metlist, collapse = " "),
            call. = FALSE)
    }

    n <- nrow(seqdata)
    lgth <- max(seqlength(seqdata))

    new.lgth <- ceiling(lgth/tspan)
    new.lgth.f <- floor(lgth/tspan)

    newseq <- seqdata[,1:new.lgth]
    cnames <- names(seqdata)
    newcnames <- cnames[seq(from=1, to=lgth, by=tspan)]

    prev <- 0
    if (method=="first") prev <- tspan - 1
    for (i in 1:new.lgth.f){
        newseq[,i] <- seqdata[,tspan*i - prev]
    }
    if (method=="mostfreq") {
        for(i in 1:new.lgth.f) {
            st.freq <- suppressMessages(seqistatd(seqdata[,(tspan*(i-1) + 1):(tspan*i)]))
            newseq[,i] <- apply(st.freq,1,function(x){names(which.max(x))})
        }
}
    if (new.lgth > new.lgth.f){
        if (method=="first"){
            newseq[,new.lgth] <- seqdata[,tspan*new.lgth.f + 1]
        } else if (method=="mostfreq") {
            st.freq <- suppressMessages(seqistatd(seqdata[,(tspan*new.lgth.f+1):lgth]))
            newseq[,new.lgth] <- apply(st.freq,1,function(x){names(which.max(x))})
        } else { # method=="last"
            newseq[,new.lgth] <- seqdata[,lgth]
        }
        newcnames[new.lgth] <- cnames[lgth]
    }
    colnames(newseq) <- newcnames
    attr(newseq,"xtstep") <- ceiling(attr(seqdata,"xtstep")/tspan)
	return(newseq)
}
