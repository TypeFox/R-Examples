`spec` <-
function (x, taper = 0, detrend = FALSE, demean=TRUE,method = c("pgram", "ar"), 
ci.plot=FALSE,ylim=range(c(lower.conf.band,upper.conf.band)),
...) 
{
    x = ts(x, frequency = 1)
if(!ci.plot) {if(missing(ylim))   sp=switch(match.arg(method), pgram = spec.pgram(x, taper = taper, 
        detrend = detrend, demean = demean,...), ar = spec.ar(x,...))  else    
        sp=switch(match.arg(method), pgram = spec.pgram(x, taper = taper, 
        detrend = detrend, demean = demean,ylim=ylim,...), 
        ar = spec.ar(x,ylim=ylim,...))
return(sp)
}
if(ci.plot==TRUE && method=="ar") stop("Option ci.plot==TRUE is not implemented for the ar method")
if(ci.plot==TRUE && method=="pgram") {
sp = spec.pgram(x, taper = taper, detrend = detrend, 
            demean = demean, ylim = ylim, plot=FALSE,...)
v=df.kernel(sp$kernel)
lower.conf.band=sp$spec*v/qchisq(.025,v)
upper.conf.band=sp$spec*v/qchisq(.975,v)
sp = spec.pgram(x, taper = taper, 
        detrend = detrend, demean = demean,ylim=ylim, ...)
lines(sp$freq,sp$spec*v/qchisq(.025,v),lty='dotted')
lines(sp$freq,sp$spec*v/qchisq(.975,v),lty='dotted')
}
sp
}

