"confi.fun" <-
function(abw,kap.mei,method)
{
    # Using the already calculated derivation this function calculates
    # the upper and lower boundary of a confidence band by substracting
    # resp. adding the derivation "abw".
    
    if(method=="linear")
    {
        lower <- kap.mei-abw
        upper <- kap.mei+abw
    }
    if(method=="log")
    {
        lower <- kap.mei^(1/abw)
        upper <- kap.mei^abw
    }
    return(list(lower=lower,upper=upper))
}

