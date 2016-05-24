dataSegments <- structure(
    function #Data segments
    ###Segmented data sets required by function
    ###\code{\link{plotSegments}}.
    (
        image,##<<Either path of an image section or an array
              ##representing a gray matrix.
        segs = 1,##<<number of image segments.
        ...##<< arguments to be passed to \code{\link{ringWidths}}.

     )
    {
        rwidths <- ringWidths(image,...)
        attrw <- attributes(rwidths)
        ns <- names(attributes(rwidths))
        nm <- lapply(strsplit(ns,'rbord.'),function(x)x[length(x)])
        names(attrw) <- nm
        attrg. <- attrw[names(formals(imageTogray))]
        attrb1 <- setdiff(names(formals(ringBorders)),'...')

        attrb2 <- unique(attrb1,names(formals(imageTogray)))
        attrb. <- attrw[attrb2]
        gray <- do.call(imageTogray,attrg.)
        attrb. <- attrb.[!sapply(attrb.,is.null)]
        attrbg. <- c(attrb.,attrg.)[unique(names(c(attrb.,attrg.)))]
        
        coltypes <- do.call(ringBorders,attrbg.)
        ## attrc <- attributes(coltypes)
        
        f.chunk<- function(x,segs)
            {split(x,factor(sort(rank(x)%%segs)))}
        f.rown <- function(x)as.numeric(rownames(x))
        pixels. <- f.chunk(f.rown(coltypes),segs)
        minl <- min(unlist(lapply(pixels.,length)))
        f.split <- function(x,y){
         l <- list();for(i in 1:length(y))
         l[[i]] <- subset(x,f.rown(x)%in%y[[i]])
        return(l)}
        
        coltypes. <- f.split(coltypes,pixels.)
        rwidths. <- f.split(rwidths,pixels.)
        range. <- lapply(pixels.,range)
        gray.<- list()
        for(z in 1:segs)gray.[[z]] <- gray[,pixels.[[z]]]
        splitd <- list(gray.,coltypes.,rwidths.)
        attributes(splitd) <- c(attributes(splitd),attrw,segs=segs)
        names(splitd) <- c('imageTogray','ringBorders','ringWidths')
        return(splitd)
        ###a list with segmented sets of the gray matrix, the ring
        ###borders, and the ring widths (see \code{\link{ringWidths}},
        ###and \code{\link{plotSegments}}).
   }
,
   ex=function(){
        ## (not run) Read one image section in package measuRing:
        image1 <- system.file("P105_a.tif", package="measuRing")    
        ## (not run) compute a gray matrix from its RGB:
        gray <- imageTogray(image1)
        ## (not run) Columns in gray matrix to be included/excluded:
        Toinc <- c(196,202,387,1564) 
        Toexc <- c(21,130,197,207,1444,1484)
        ## (not run) segmented data:
        segm <- dataSegments(image1,segs = 3)
        lapply(segm,str)
        attributes(segm)
     }
 )
