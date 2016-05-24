"rsctrl" <-
function(burn_in=100, n_eff=100, step=16,
                 seed=0, tfixed=FALSE)
{   ier <- 0
    if(n_eff < 0 | n_eff > 8191)    ier <- ier + 4
    if(burn_in < 0)                 ier <- ier + 8
    if(step <= 0)                   ier <- ier + 16
    if(seed < 0 | seed > (2**31-2)) ier <- ier + 128
    if (ier>0) rserror(ier)
    RET<-list(burn_in=burn_in, n_eff=n_eff, step=step,
              seed=seed, tfixed=tfixed)
    class(RET)<-"RSctr"
    RET
}

