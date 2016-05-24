ringBorders <- structure(
    function# Ring borders
    ### This function can find the ring borders in a gray matrix.
    (
        image,##<<character or matrix. Either path of an image section
        ##or an array ##representing a gray matrix.
        auto.det = TRUE,##<<logical. If TRUE the linear detection is
                        ##implemented (see
                        ##\code{\link{linearDetect}}).
        darker = TRUE, 	##<<logical. If TRUE the algorithm uses the
                        ##negative extremes on smoothed grays to
                        ##detect the ring borders. If FALSE the
                        ##possitive extremes are used.
        origin = 0, 	##<<numeric. an origin in smoothed gray to
                        ##find the ring borders.
        inclu = NULL, ##<<NULL or vector with column numbers in gray
                      ##matrix, other than those automatically
                      ##detected, to be considered as ring borders.If
                      ##NULL no column numbers are included.
        exclu = NULL,##<<NULL or vector with column numbers in gray
        ...##<< arguments to be passed to \code{\link{imageTogray}}.


    )
    {

        gray <- image
        if(is.character(gray))
            gray <- imageTogray(gray,...)

        f.rown <- function(x)as.numeric((rownames(x)))

        f.tit <- function(image){
            p <- '.tif'
            if(any(grepl('.png',image)))p <- '.png'
            bn <- basename(image)
            gsub(p,'',bn)}
        
        names. <- f.tit(attributes(gray)[['image']])
        smoothed <- graySmoothed(gray,FALSE)
        if(auto.det){
           borders <- linearDetect(smoothed,origin,darker)}
        else{borders <- NA}
        ini1 <- borders
        inclu <- c(inclu,recursive=TRUE)
        exclu <- c(exclu,recursive=TRUE)
        confl <- exclu[exclu%in%inclu] 
        inclu <- setdiff(inclu,confl)
        ininoinc <- setdiff(ini1,inclu)
        ininoincexc <- setdiff(ininoinc,exclu)
        cases <- list(included=inclu,excluded=exclu,automatic=ininoincexc)
        f.type <- function(x,type){
            dtype <- data.frame(borders=rep(type,length(x)))
            dtype[,'N_pixel'] <- x
            return(dtype)}
        typelist <- list()
        for(i in 1:length(cases))
            typelist[[i]] <- f.type(cases[[i]],names(cases[i]))
        types1 <- do.call(rbind,typelist)
        noav <- c('included','automatic')
        types <- subset(types1,types1[,'borders']%in%noav)
        types[,'borders'] <- rep(TRUE,nrow(types)) 
        smoothed[,'N_pixel'] <- rownames(smoothed)
        avelum <- merge(smoothed,types,by.x='N_pixel',by.y='N_pixel',all.x=TRUE)
        rownames(avelum) <- avelum[,'N_pixel']
        avelum <- avelum[order(f.rown(avelum)),c(names.,'borders')]
          attrg <- attributes(gray)[-1L]         
        ## attributes(avelum)[names(attrg)] <- attrg
         attrc <- list(origin = origin, auto.det = auto.det,
           darker = darker, inclu = inclu, exclu = exclu)
         attributes(avelum) <- c(attributes(avelum),
           gray=c(attributes(gray)['dim']),attrg,attrc)

        return(avelum)
        ###a data frame with the smoothed grays and the identified
        ###ring borders (see \code{\link{grayDarker}},
        ###\code{\link{graySmoothed}}, and
        ###\code{\link{linearDetect}}).
    }
,
    ex=function(){
        ## (not run) Read one image sample in folder of package
        ## measuRing:
        image1 <- system.file("P105_a.tif", package="measuRing")        
        ## column numbers in gray matrix to be included/avoided:
        Toinc <- c(196,202,387,1564) 
        Toexc <- c(21,130,197,207,1444,1484)        
        ##(not run) the ring borders:
        borders <- ringBorders(image1,inclu = Toinc,exclu = Toexc)
        str(borders)
        ##(not run) Plot of smoothed grays with the ring borders:
        Smooth <- ts(borders[,1])
        inclupix <- subset(borders,borders%in%TRUE)
        inclucol <- as.numeric(rownames(inclupix))
        xyborders <- data.frame(column=inclucol,smooth=inclupix[,1])
        y.lim <- c(-0.05,0.05)
        main. <- 'Ring borders'
        {plot(Smooth,xlab = 'Column',ylab = 'Smoothed gray',
              main=main.,col = 'darkgoldenrod1')
         points(xyborders[,1],xyborders[,2],pch=19,cex=0.5,col='orangered')}
    }
)
