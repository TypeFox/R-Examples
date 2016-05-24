"dlmomco" <-
function(x,para) {
    if(! are.par.valid(para)) return()
    f <- par2pdf(x,para,paracheck=FALSE)
    # Although all the pdfCCC() should be doing this, there
    # is little harm in trapping NAs again.
    f[is.na(f)] <- 0 # and reset NAs to zero
    return(f)
}
