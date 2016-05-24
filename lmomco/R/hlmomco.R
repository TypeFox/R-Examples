"hlmomco" <-
function(x,para) {
    if(! are.par.valid(para)) return();
    the.pdf <- dlmomco(x,para);
    the.cdf <- plmomco(x,para);
    the.pdf[is.na(the.pdf)] <- 0 # should already be done, but not harm again
    return(the.pdf/(1 - the.cdf));
}
