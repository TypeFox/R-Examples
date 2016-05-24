grayDarker <- structure(
    function#Gray extremes
    ###This function can detect the extremes of the smoothed gray.
    (
        smoothed,	##<<a data frame with the smoothed gray such
                        ##as that produced by
                        ##\code{\link{graySmoothed}}.
        origin = 0, 	##<<an origin to find the extremes.
        darker = TRUE 	##<<logical. If TRUE the function finds the
                        ##negative extremes. If FALSE the possitive
                        ##extremes are detected.
    )
    {

        f.rown <- function(x)as.numeric((rownames(x)))

        f.tit <- function(image){
            p <- '.tif'
            if(any(grepl('.png',image)))p <- '.png'
            bn <- basename(image)
            gsub(p,'',bn)}

        names. <- f.tit(attributes(smoothed)[['image']])
        smoothed[,'smooth'] <- smooth(na.omit(smoothed[,names.]),
        twiceit=FALSE)
        turnp <- turnpoints(smoothed[,'smooth'])[['tppos']] #CRAN::pastecs
        inidet0 <- smoothed[turnp,]
        inidet <- na.omit(inidet0[inidet0[,names.]>origin,])
        if(darker)
        inidet <- na.omit(inidet0[inidet0[,names.]<origin,])        
        inidet <- f.rown(inidet)
        return(inidet)
            ###vector with the columns in gray matrix corresponding to
            ###the extremes (see \code{\link{graySmoothed}} and
            ###\code{\link{linearDetect}}).
    }
,
    ex=function(){
        ## (not run) Read one image section:
        image1 <- system.file("P105_a.png", package="measuRing")    
        ## (not run) gray matrix from RGB in image:
        gray <- imageTogray(image = image1,ppi = 1000)
        ## (not run) smoothed gray:
        smoothed <- graySmoothed(gray)
        ## (not run) column numbers of possitive and negative
        ## extremes:
        posit <- grayDarker(smoothed,darker=FALSE)
        nega <- grayDarker(smoothed,darker=TRUE)
        str(nega)

   }
)
