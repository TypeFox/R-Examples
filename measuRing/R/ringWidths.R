ringWidths <- structure(
    function#Ring widths 
    ### This function can compute the ring widths from the ring
    ### borders detected on an image section.
    (
        image,##<<character or matrix. Either path of an image section
        ##or an array representing a gray matrix.
        last.yr = NULL,##<<year of formation of the newest ring. If
                       ##NULL then the rings are numbered from one
                       ##(right) to the number of detected rings
                       ##(left).
        ...##<< arguments to be passed to \code{\link{ringBorders}}.
    )
    {

        f.rown <- function(x)as.numeric((rownames(x)))

        f.tit <- function(image){
            p <- '.tif'
            if(any(grepl('.png',image)))p <- '.png'
            bn <- basename(image)
            gsub(p,'',bn)}
        
        pixtypes <- ringBorders(image,...)        
        attb <- attributes(pixtypes)
        ppi <- attb[['ppi']]
        names. <- f.tit(attb[['image']])
        scale <- 25.4/ppi ## (mm)
        pixtypes[,'distance'] <- f.rown(pixtypes)*scale
        finald <- pixtypes[pixtypes[,"borders"]%in%TRUE,]
        
        f.label <- function(finald,last.yr){
           finald[,'item'] <- c(1:nrow(finald))
           finald[,'growth'] <- with(finald,(max(distance) - distance))
           if(!is.null(last.yr))year1 <- last.yr + 1
           else{year1 <- nrow(finald)}
           finald[,'year'] <- with(finald,year1-item)
           finald[,'delta'] <- with(finald,c(rev(diff(rev(growth))),0))
           finald <- finald[1:(nrow(finald)-1),c('year','delta')]}
        
        if(nrow(finald)==0){
        trwd <- data.frame(year=vector(),delta=vector())}
        else{
         trwd <- f.label(finald,last.yr)
         last.yr <- max(trwd[,'year'])}
        names(trwd) <- c('year',names.)
        attributes(trwd) <- c(attributes(trwd),## attcol,
                              rbord=attb,last.yr=last.yr)
        return(trwd)
        ###data frame with the ring widths.


    }
,
    ex=function(){
        ## (not run) Read one image section:
        image1 <- system.file("P105_a.tif", package="measuRing")       
        ## (not run) columns in gray matrix to be included/excluded:
        Toinc <- c(196,202,387,1564) 
        Toexc <- c(21,130,197,207,1444,1484)
        ## (not run) tree-ring widths
        rwidths <- ringWidths(image1,inclu = Toinc,exclu = Toexc,last.yr=NULL)
        str(rwidths)
        ##plot of computed tree-ring widths:
        maint <- 'Hello ring widths!'
        plot(rwidths,type='l',col = 'red',main = maint,
             xlab = 'Year',ylab = 'Width (mm)')    
    }
)
