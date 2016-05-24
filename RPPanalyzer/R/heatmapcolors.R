`heatmapcolors` <-
function(dat,ncol,lowcol="blue",highcol="red",border=NULL) {
    # define color palette, set the 0 to 'black' if possible
    # if all values negative, only take blue values, if all positive,
      #only take red values
    # generates ncol color steps
    if(is.null(border)) {
    rng <- range(dat,na.rm=T)
    if(rng[1]<0){
        # all negative
        if(rng[2]<0) {
            border=ncol
        } else {
            border <- round(abs(rng[1])/abs(rng[2]-rng[1])*ncol)
        }
    } else { # all positive
        border <- 0
    }
    }
    colors <-c(colorpanel(border,low=lowcol,high="white"),
               colorpanel(ncol-border,low="white",high=highcol))
               
    return(colors)
}

