##### internal functions from Dr. Tibshirani's software package GSA
coxfunc <-
function(x, y, censoring.status)
{
        junk <- coxscor(x, y, censoring.status)
        scor<-junk$scor
        sd <- sqrt(coxvar(x, y, censoring.status, coxstuff.obj=junk$coxstuff.obj))
        tt <- scor/sd
        return(list(tt = tt, numer = scor, sd = sd ))
}
