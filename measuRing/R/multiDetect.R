multiDetect <- structure(
   function# Multiple detection
### This function provides tools to recursively detect the ring
### borders in multiple image sections.
            ##details<<Users running R from IDEs and aiming to develop
            ##visual control on several image segments should be sure
            ##that such environments support multiple graphics devices
            ##(see \code{\link{ringDetect}}).
    (
        pattern,##<<character with a common pattern in the names of
                ##the image sections.
        is.png = FALSE, ##<< logical. If FALSE the tif images in
                        ##working directory are processed.
        from = 1,##<<character with a complementary pattern, or
                 ##position number in folder, of the initial image to
                 ##be processed.
        to = 'all',##<<character with a complementary pattern, or
                   ##position number in folder, of the final image to
                   ##be processed. If this argument is 'all' then all
                   ##the images matching the argument in pattern are
                   ##processed.
        inclu.dat = NULL, ##<< a data frame such as that contained in
                          ##the ouput of this same function with the
                          ##column numbers to be updated.
        ... ##<<arguments to be passed to \code{\link{ringDetect}}
            ##(ppi, last.yr, rgb, p.row, auto.det, darker, origin,
            ##inclu, exclu, and plot), or to
            ##\code{\link{plotSegments}}(segs, marker, col.marker,
            ##ratio, and tit).
        )
    {

        patt <- list.files(path=getwd(),pattern=pattern)
        patt.ext <- '.tif'        
        if(is.png)patt.ext <- '.png'
        PGen <- patt[grep(patt,pattern = patt.ext)]
        erimp <- c('missing','pattern','in path')
        if(length(PGen)==0){
           Er1. <- paste(erimp[1],patt.ext,erimp[2],erimp[3])
        stop(paste(Er1.),call.= FALSE)}

        if(to == 'all')to <- length(PGen)
        if(from == 'all')from <- length(PGen)

        fr_to <- function(PGen,to){
          patft <- grep(PGen,pattern = to)
          if(length(patft)==0){
          sft <- dQuote(to)
          Er2. <- paste(erimp[1],erimp[2],sft,erimp[3]) 
          stop(Er2.,call.=FALSE)}
          return(patft)}        

        if(is.character(to)){
           to <- fr_to(PGen,to)}
        if(is.character(from))from <- fr_to(PGen,from)

        frmp <- function(from){
            if(from>length(PGen)){
        stop(paste(erimp[1],patt.ext,'image nr:',from,erimp[3]),call.=FALSE)}
        return(from)}
        from <- frmp(from);to <- frmp(to)
        if(from > to){
        from.. <- to;to <- from;from <- from..}
        PGen.. <- PGen[from:to]
        indots <- as.list(substitute(list(...)))[-1L]
            options(warn = -1)

        fplot <- function(PGen,...){
            minc <- 'Including columns in'
            mex <- 'Excluding rings in'
            dots <- '...'
            autoGen <- do.call(ringDetect,list(PGen,   
               tit=paste(minc,PGen,dots),...))
            incGen <- ringSelect(autoGen)
            autoGen <- update(autoGen,inclu = incGen,
               tit=paste(mex,PGen,dots))
            excGen <- ringSelect(autoGen,FALSE)
            autoGen <- update(autoGen,
               inclu = incGen,exclu = excGen,plot=FALSE)}

        fnoplot <- function(PGen,...){
           autoGen <- do.call(ringDetect,list(PGen,...))}

        flinc <- function(x){
            x <- as.data.frame(x[,-1L])
            x. <- lapply(seq_len(ncol(x)),
            function(i)x[!is.na(x[i]),i])
            names(x.) <- names(x)
            return(x.)}
        
        if(!is.null(inclu.dat)){
            inc <- flinc(inclu.dat)
            fplotinc <- function(PGen,inc,...){
            minc <- 'Including columns'
            mex <- 'Excluding rings'
            dots <- '...'
            options(warn = -1)
            autoGen <- do.call(ringDetect,list(PGen,   
               tit=paste(minc,PGen,dots),inclu = inc,...))
            incGen <- ringSelect(autoGen)
            incGen <- unique(c(incGen,inc))
            autoGen <- update(autoGen,inclu = incGen,
               tit=paste(mex,PGen,dots),auto.det=FALSE)
            excGen <- ringSelect(autoGen,FALSE)
            autoGen <- update(autoGen,
               inclu = incGen,exclu = excGen,plot=FALSE)}        

            fnoplotinc <- function(PGen,inc,...){
           autoGen <- do.call(ringDetect,list(PGen,...))
           autoGen <- update(autoGen,inclu = inc,auto.det=FALSE)}
            
        if(indots$plot)
        autoGen <- lapply(PGen..,function(x)fplotinc(x,inc,...))
        else{
        autoGen <- lapply(seq_len(length(PGen..)),
                          function(i)fnoplotinc(PGen..[i],inc[i],...))}}

        else{        
        if(indots$plot)
        autoGen <- lapply(PGen..,function(x)fplot(x,...))
        else{
        autoGen <- lapply(PGen..,function(x)fnoplot(x,...))}}#}

        rwidths. <- lapply(seq_len(length(autoGen)),function(i)
          autoGen[[i]][['ringWidths']])
        rwidths. <- lapply(rwidths.,function(x)x[order(x$year),])
        
          rownam. <- lapply(rwidths.,function(x){
          nam. <- names(x)
          x[,2] <- as.numeric(rownames(x))
          names(x) <- nam.
          return(x)})
        rmarks. <- lapply(rwidths.,function(x){
          if(nrow(x)==0){
           x <- data.frame(year=vector(),delta=vector())}
          else{
           marker=indots$marker
           cols <- colNarrow(x,marker)
           x[cols,2] <- TRUE
           x[setdiff(rownames(x),cols),2] <- NA}
          return(x)})
        ldat <- list(rwidths.,rownam.,rmarks.)
        fmatch <- function(rownam.){
            Reduce(function(x,y){merge(x,y,all=TRUE)},rownam.)}
        alld <- lapply(seq_len(length(ldat)),function(i)fmatch(ldat[[i]]))
        names(alld) <- c('ringWidths','colNames','colNarrows')

        if(is.null(indots$marker))
        return(c(alld[1:2],call=sys.call()))
        return(c(alld,call=sys.call()))
        
        ###list with three data frames: the tree-ring widths, the
        ###column numbers of the detected ring borders, and the narrow
        ###rings (see \code{\link{ringBorders}},
        ###\code{\link{ringDetect}}, and \code{\link{plotSegments}}).
    }
,
    ex=function(){
        ## (not run) Set working directory:
        setwd(system.file(package="measuRing"))
        ## (not run) List the tif images the folder:
        list.files(path=getwd(),pattern='.tif')

        ## (not run) run multiDetect:
        ## -provide at least one argument of ringDetect
        tmp <- multiDetect(pattern = 'P105',last.yr=2012,plot = FALSE)
        ##
        ## (not run) Excluding/changing some column numbers in tmp:
        dd <- tmp$colNames
        ddtest <- dd #to be compared with final outputs
        dd[dd$year%in%1999:2012,] <- NA
        tail(dd,20)
        tmp1 <- update(tmp,inclu=dd,auto.det=FALSE)
        dm <- tmp1$colNames
        dmtest <- dm #to be compared with final outputs
        ##
        ## (not run)changing columns from tmp with visual control
        ## -choose five or six rings at the bark side and later
        ## exclude any one of them:
        tmp2 <- update(tmp,plot=TRUE,to='_a',segs = 1,auto.det=FALSE)
        dm2 <- tmp2$colNames
        newm <- merge(dm,dm2,by='year',all.x=TRUE)
        dm[,'P105_a'] <- newm[,ncol(newm)]
        tmp3 <- update(tmp,inclu=dm,plot=FALSE,auto.det=FALSE)
        dm3 <- tmp3$colNames
        ## compare initial and final columns in gray matrix
        tail(ddtest,15)
        tail(dmtest,15)
        tail(dm3,15)
    }
)
