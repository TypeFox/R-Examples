ringSelect <- structure(
   function#Visual selection
###This function is used to visually include/exclude ring borders;
###consequently, graphics devices from \code{\link{ringDetect}} or
###\code{\link{plotSegments}} should be active.
    ##details<<Columns in gray matrix are either identified and stored
    ##by left-clicking the mouse over the central axis of gray image;
    ##pixel numbers of just added ring borders are highlighted on such
    ##a gray raster. The raphics devices are sequentially closed by
    ##right-clicking the mouse. After a graphics device has been
    ##closed, the graphics device of the following segment is
    ##activated, and visual selection on such a new segment can be
    ##performed. Closing the graphics devices with other procedures
    ##will stop the selection of ring borders.
    (
        rdetect,##<<a list containing data frames of ring widths and
                   ##ring borders such as that produced by
                   ##\code{\link{ringDetect}}.
        any.col = TRUE ##<<logical. If FALSE only those column numbers
                       ##in gray matrix previouly identified as ring
                       ##borders can be selected.
    )
    {
        pixtypes <- rdetect[['ringBorders']]
        dim. <- attributes(rdetect)[['gray.dim']]
        plot. <- attributes(rdetect)[['opar']]
        options(warn = -1)
        
        f.inclu <- function(gray,pixtypes,any.col){
           segs <- length(dev.list())
           if(segs==0)
               stop('missing image segments',call.=FALSE)
           f.rown <- function(x)as.numeric(rownames(x))

           if(any.col){
           rowpix <- f.rown(pixtypes)
           col. <- 'black'}
           else{
           rowpix <- f.rown(pixtypes[pixtypes[,"borders"],])
               col. <- 'red'}           
           range. <- function(n)c(0,n)
           yrange <- range.(dim.[1])
           xrange <- range(rowpix)
           
            fexcleca <- function(rowpix,col.){
                par(plot.)
                par(mfg=c(2,1))
               xy <- identify(rowpix,rep(yrange[2]/2,length(rowpix)),
                rowpix,tolerance=0.1,cex=0.55,col=col.,srt=90)
                return(xy)}
           
            xy <- list()
            for(i in 1:segs)
               repeat{
                dev.set(which=dev.list()[(segs+1)-i])
                if(length(dev.list())==1){
                    dev.set(which=dev.cur())}
                xy[[i]] <- fexcleca(rowpix,col.)
                print(xy[[i]])
                dev.off(which=dev.cur())
                break()}
            selected <- list()
            for(i in 1:segs)selected[[i]] <- xy[[i]]
            exclu <- unlist(selected)
            return(exclu)}
        
        return(f.inclu(gray,pixtypes,any.col))
        
        ### vector with column numbers in gray matrix of the
        ### identified ring borders.

    }
,
    ex=function(){
        ## Read one image in package folder:
        image1 <- system.file("P105_a.tif", package="measuRing")
        ## (not run) Initial diagnostic:
        detect1 <- ringDetect(image1,segs=2,marker=7)
        ##
        ## (not run) Choose other columns in gray matrix (see ringSelect);
        ## (not run) graphical devices from ringDetect should be active!
        ## (not run) Including columns:
        ##
        ## (uncomment and run):
        ## Toinc <- ringSelect(detect1)
        ## detect1 <- update(detect1, inclu = Toinc)
        ##
        ## (not run)  ring borders to be excluded:
        ## (uncomment and run):
        ## Toexc <- ringSelect(detect1,any.col = FALSE)
        ## detect1 <- update(detect1, exclu=Toexc)
        ## (not run) kill previous plot:
        graphics.off()
    }
)
