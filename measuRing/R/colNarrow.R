colNarrow <- structure(
    function#Narrow rings
### This function can detect narrow rings in an object such as that
### produced by \code{\link{ringWidths}}.
    ##details<< Each ring is averaged with those rings on either side
    ##of it (t-1,t,t+1), and averages are divided by the highest
    ##computed average in the sample; such quotients are scaled from
    ##10 (the narrowest possible ring) to one (the broadest ring).
    (
        rwidths, ##<<a dataframe with the ring widths such as that
                 ##produced by \code{\link{ringWidths}}.
        marker = 5 ##<<a number from 1 to 10. Those rings with scaled
                   ##averages greater than or equal to this argument
                   ##will be identified as narrow rings.
        
    )
    {
    f.prop <- function(i,data){
        ifelse(i==1|i==nrow(data),NA,
               data[i,2]/mean(data[(i-1):(i+1),2]))}
    px.flag <- NA
        if(nrow(rwidths)!=0){
    frac <- sapply(1:nrow(rwidths),f.prop,rwidths)
    names(frac) <- rownames(rwidths)
    frac1 <- frac[frac<=1&!is.na(frac)]
    nw.dim <- round(10*((1-frac1)/max(1-frac1)))
    px.flag <- names(nw.dim[nw.dim>=marker])}
    return(px.flag)
    ###character vector with the columns in gray matrix corresponding
    ###to the narrow rings (see \code{\link{ringDetect}},
    ###\code{\link{multiDetect}}, and\code{\link{plotSegments}}).
    }
,
    ex=function(){
        ## (not run) Read one image section in package measuRing:
        image1 <- system.file("P105_a.png", package="measuRing")    
        ## (not run) compute a gray matrix from RGB in the image:
        gray <- imageTogray(image = image1,ppi=1000)
        ## (not run) Columns in gray matrix to be included/excluded:
        Toinc <- c(196,202,387,1564) 
        Toexc <- c(21,130,197,207,1444,1484)
        ## (not run) tree-ring widths:
        rwidths <- ringWidths(gray,inclu = Toinc,exclu = Toexc,last.yr=2012)
        ##(not run) narrow rings:
        narrows <- colNarrow(rwidths,marker = 8)
    }
)
