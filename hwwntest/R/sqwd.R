sqwd <-
function (x, filter.number=10, family="DaubExPhase", type = "station",
	m0=3) 
{
    ans1 <- wd(x, filter.number=filter.number,
	family=family, type = "station")
    ec <- sqcoefvec(m0 = m0, filter.number, family, stop.on.error = FALSE)
    if (ec$ecode != 0) {
        cat("sqcoefvec error\n")
        cat(ec$ians$message, "\n")
        return(0)
    }
    answd <- sqndwd(x = x, ec = ec)
    if (type == "wavelet") {
        answst <- convert(answd)
        nlev <- nlevelsWT(answst)
        answd <- wd(rep(0, 2^nlev), filter.number = answst$filter$filter.number, 
            family = answst$filter$family)
        for (i in 0:(nlev - 1)) {
            v <- getpacket(answst, level = i, index = 0)
            answd <- putD(answd, level = i, v = v)
            v <- getpacket(answst, level = i, index = 0, type = "C")
            answd <- putC(answd, level = i, v = v)
        }
        return(answd)
    }
    else if (type == "station") 
        return(answd)
    else stop(paste("Unknown type", type))
}
