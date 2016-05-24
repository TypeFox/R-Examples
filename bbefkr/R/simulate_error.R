simulate_error <-
function(samplesize, errordensity)
{
    if(errordensity == "normal")
    {
        err = rnorm(samplesize)
    }
    if(errordensity == "asyclaw")
    {
        err = simasyclaw(samplesize) 
    }
    if(errordensity == "asydoubleclaw")
    {
        err = simasydoubleclaw(samplesize)
    }
    if(errordensity == "bimodal")
    {
        err = simbimodal(samplesize)
    }
    if(errordensity == "claw")
    {
        err = simclaw(samplesize)
    }
    if(errordensity == "discretecomb")
    {
        err = simdiscretecomb(samplesize)
    }
    if(errordensity == "doubleclaw")
    {
        err = simdoubleclaw(samplesize)
    }
    if(errordensity == "kurtotic")
    {
        err = simkurtotic(samplesize)
    }
    if(errordensity == "outlier")
    {
        err = simoutlier(samplesize)
    }
    if(errordensity == "sepbimodal")
    {
        err = simsepbimodal(samplesize)
    }
    if(errordensity == "skewbimodal")
    {
        err = simskewbimodal(samplesize)
    }
    if(errordensity == "skewunimodal")
    {
        err = simskewunimodal(samplesize)
    }
    if(errordensity == "smoothcomb")
    {
        err = simsmoothcomb(samplesize)
    }
    if(errordensity == "strongskew")
    {
        err = simstrongskew(samplesize)
    }
    if(errordensity == "trimodal")
    {
        err = simtrimodal(samplesize)
    }
    return(err)
}
