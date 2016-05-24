linearDetect <- structure(
    function#  Linear detection
### Function for developing linear detection of ring borders.
## details<< This function assumes that negative extremes on smoothed
## grays (darker grays) correspond to late wood and that positive
## extremes on smoothed grays correspond to early wood; therefore, the
## intervals between late/early wood are used to detect the ring
## borders. Detection is centered on a constant controlled by the
## argument origin. Ring borders are detected within consecutive pairs
## of negative/positive extremes: those smoothed grays closest to
## origin in the itervals are accounted as ring borders.

    (
        smoothed,##<< a data frame with smoothed grays such as that
        ##produced by \code{\link{graySmoothed}}.
        origin = 0,##<<numeric. an origin in smoothed gray to find the
                   ##ring borders.
        darker = TRUE 	##<<logical. If TRUE the algorithm uses the
        ##negative extremes on smoothed grays to detect the ring
        ##borders. If FALSE the possitive extremes are used.
    )
    {
        f.rown <- function(x)as.numeric(rownames(x))

        f.tit <- function(image){
            p <- '.tif'
            if(any(grepl('.png',image)))p <- '.png'
            bn <- basename(image)
            gsub(p,'',bn)}

        names. <- f.tit(attributes(smoothed)[['image']])
        names.. <- paste(names.,'.',sep='')
        range. <- range(smoothed)
        if(darker){turneg <- grayDarker(smoothed,origin)}
        else{turneg <- grayDarker(smoothed,origin,darker=F)}
        smoothed[,names..] <- smoothed[,names.] - origin
        smoothed[,'cutZ'] <-
            abs(c(diff(ifelse(smoothed[,names..]<0,0,1)),0))*f.rown(smoothed)
        difcutZ <- diff(c(0,smoothed[,'cutZ'][smoothed[,'cutZ']!=0]))
        fracutZ <- list()
        for(i in 1:length(difcutZ))fracutZ[[i]] <- rep(i,difcutZ[i])
        codeZ <- unlist(fracutZ)
        smoothZ_data <- data.frame(smoothed[1:length(codeZ),],codeZ=codeZ)

        smoothZ_list <- split(smoothZ_data,smoothZ_data[,'codeZ'])
        det2ext <- unique(smoothZ_data[turneg,c('codeZ')])
        selaut <- smoothZ_list[as.character(det2ext)]
        newauto <- c(do.call(rbind,lapply(selaut,function(x)f.rown(x[1,]))))
        attri <- list(origin=origin,darker=darker)
        attributes(newauto) <- c(attributes(newauto),attri)
        return(newauto)
	###vector with column numbers in gray matrix of the detected
	###ring borders (see \code{\link{grayDarker}}, and
	###\code{\link{graySmoothed}}).
    }
,
    ex=function(){
        ## (not run) Read one image section in package measuRing:
        image1 <- system.file("P105_a.tif", package="measuRing")    
        ## (not run) smoothed gray:
        smoothed <- graySmoothed(image1)
        ## linear detection:
        borders <- linearDetect(smoothed)
        str(borders)
    }
)
