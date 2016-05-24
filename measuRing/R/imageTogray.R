imageTogray <- structure(
    function#Gray matrix
    ### This function can compute a gray matrix from the RGB in an
    ### image section.Such image section can be compressed in either
    ### portable network graphics format (png) or tagged image file
    ### format (tif).
    
    (
        image,##<<character. path of an image section.
        ppi = NULL, ##<< NULL or integer. If NULL the image resolution
                    ##in points per inch is extracted from attributes
                    ##in \code{image}. If this attribute is not
                    ##embedded then users should provide it
        rgb = c(0.3,0.6,0.1), ##<<vector with three fractions, all of
        ##them adding to one, to combine RGB channels into gray
        ##matrix.
        p.row = 1 ##<<proportion of rows of gray matrix to be
                  ##processed.
     )
    {
        ftree <- function(image,native,...){
            if(is.array(image)){tree <- image}
            if(is.character(image)){
             if(grepl(".png",image)){
             tree <- do.call(readPNG,list(image,...))}
             if(grepl(".tif",image)){
             tree <- do.call(readTIFF,list(image,...))}}
            return(tree)}
        tree <- ftree(image,native=FALSE,info=TRUE)
        ppi.dfl <- function(tree){
            att. <- attributes(tree)
            if(is.null(att.))return(NULL)
            attnam. <- names(attributes(tree))
            resol. <- grepl('esolut',attnam.)|
                grepl('ESOLUT',attnam.)
            x. <- grepl('[xX]',attnam.[resol.])
            if(any(x.)){
                xX.resol <- attnam.[resol.][x.]
                attributes(tree)[[xX.resol]]}
            else{NULL}}
        if(is.null(ppi))ppi <- ppi.dfl(tree)
        if(is.null(ppi))
         {stop('missing image resolution, provide ppi',call.=FALSE)}
        gray <- t(apply(tree,1,function(x)x%*%rgb))
        if(p.row<1){
        width. <- min(dim(gray))/2 
        width.. <- round(p.row*width.)
        l1 <- width. + width..
        l2 <- width. - width..
        gray <- gray[l1:l2,]}
        attributes(gray) <- c(attributes(gray),list(image = image,
        dim = c(dim(gray),recursive=TRUE),
        ppi = ppi,rgb = rgb,p.row = p.row))
        return(gray)
        ###a gray matrix containing the image reflectances.
    }
,
    ex=function(){
        ## (not run) Read two image sections in package measuRing:
        image1 <- system.file("P105_a.tif", package="measuRing")
        image2 <- system.file("P105_a.png", package="measuRing")
        ## (not run) compute a gray matrix:
        gray <- imageTogray(image1)
        ## (not run) - the ppi is embedded in the image:
        attributes(gray)
        ## (not run) but, the ppi is not embedded in image2:
        ## - imageTogray will return an error:
        ## uncoment and run
        ## gray2 <- imageTogray(image2)
        ## attributes(gray2)
        ## - the ppi should be provided (i.e. ppi = 1200):
        gray3 <- imageTogray(image2,ppi = 1200)
        attributes(gray3)
        ##(not run) a plot of the gray matrix        
        xrange <- range(0:ncol(gray)) + c(-1,1)
        yrange <- range(0:nrow(gray)) + c(-1,1)    
        {plot(xrange,yrange,xlim=xrange,ylim=yrange,xlab='',
              ylab='',type='n',asp=0)
        rasterImage(gray,xrange[1],yrange[1],xrange[2],yrange[2])}    
    }
)
