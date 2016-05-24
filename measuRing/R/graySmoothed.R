graySmoothed <- structure(
    function#Smoothed gray
    ###Averaging, detrending, and smoothing of the columns in a gray matrix. 
    (
        image,##<<character or matrix. Either path of an image section
              ##or an array representing a gray matrix.
        all = FALSE, ##<< logical. If TRUE the column numbers and
                    ##moving averages are added to the output.
        ...##<< arguments to be passed to \code{\link{imageTogray}}.
    )
    {

        f.tit <- function(image){
            p <- '.tif'
            if(any(grepl('.png',image)))p <- '.png'
            bn <- basename(image)
            gsub(p,'',bn)}
        
        gray <- image
        if(is.character(gray))
            gray <- imageTogray(gray,...)
        names. <- f.tit(attributes(gray)[['image']])
        averageIngray <- function(gray){
        e <- apply(gray, 2, function(x) exp(mean(log(x))))
        smoothd <- data.frame(V_pixel=e)
        return(smoothd)}
        smoothd <- averageIngray(gray)
        lagsel <- lagIngray(gray)
        hanning <- function (x, n = 7) { ## from dplR!
            j <- 0:(n - 1)
            win <- 1 - cos(2 * pi/(n - 1) * j)
            win <- win/sum(win)
            as.vector(filter(x, win))}
        hann0 <- hanning(smoothd[,'V_pixel'],lagsel)
        names(hann0) <- rownames(smoothd)
        hann <- hann0[!is.na(hann0)]
        hanm <- hann[1]
        hanM <- hann[length(hann[!is.na(hann)])]
        hann0[as.numeric(names(hann0))<as.numeric(names(hanm))] <- hanm 
        hann0[as.numeric(names(hann0))>as.numeric(names(hanM))] <- hanM 
        smoothd[,'m.av'] <- hann0
        smoothd[,names.] <-smoothd[,'V_pixel']-smoothd[,'m.av']
        if(all==FALSE){smoothd <- data.frame(smoothd[,names.])
                   names(smoothd) <- names.}

          attrg <- attributes(gray)[-1L]         
         attributes(smoothd) <- c(attributes(smoothd),attrg)

        
        return(smoothd)
            ###data frame with the smoothed grays. If argument all is
            ###TRUE the output is extended the with the columns in
            ###gray matrix, and the moving averages.
    }
,
    ex=function(){
        ## (not run) Read one image section in package measuRing:
        image1 <- system.file("P105_a.png", package="measuRing")    
        ## (not run) the smoothed gray:
        smoothed <- graySmoothed(image1,ppi=1000)
        ## (not run) Plot of the smoothed gray:        
        Smooth <- ts(smoothed)
        main. <- 'Smoothed gray'
        plot(Smooth,xlab = 'Column', main=main.,
             ylab = 'Smoothed gray',col = 'gray')

    }
)
