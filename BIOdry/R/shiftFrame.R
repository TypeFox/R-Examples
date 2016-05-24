shiftFrame <- structure(function#Shifting of ring-data frames
### Ring-data frames are reshaped into multilevel data frames (or vice
### versa). SI units of the processed chronologies can be changed.
                        ##details<<Rows of ring-data frames should be
                        ##named with chronological years, and columns
                        ##of such data frames should be labeled with
                        ##level units defined by sample design. Level
                        ##units in labels are separated with dot (.),
                        ##beginning with names of higher level units
                        ##(i.e. ecorregion, climatic location, or
                        ##plot) and ending with the lowest level unit
                        ##(usually sample/replicate). For example, the
                        ##code name of core 'a' in tree '2' on plot
                        ##'P16001' on ecorregion 'M1' will have the
                        ##name: 'M1.P16001.2.a'.
(
    rd, ##<<\code{data.frame}. Ring-data frame (see details), or
        ##multilevel data frame (see value).
    lev.nm = c('plot','tree','sample'), ##<<for the case of ring-data
                                        ## frames, \code{character}
                                        ## vector with names of the
                                        ## factor-level columns in the
                                        ## final multilevel data
                                        ## frame, beginning with name
                                        ## of the highest level column
                                        ## and ending with the name of
                                        ## the lowest level column. If
                                        ## \code{rd} is a multilevel
                                        ## data frame then this
                                        ## argument is ignored.
    which.x = NULL, ##<<for the case of multilevel data frames,
                    ##\code{character} name of the column to be
                    ##reshaped into a ring-data frame. If NULL then
                    ##the first \code{numeric} column is processed. If
                    ##\code{rd} is a ring-data frame then this
                    ##argument is ignored.
    un = NULL ##<< \code{NULL}, one, or two \code{character} names of
              ##SI units to record/transform the processed
              ##variables. One character records metric system; two
              ##characters with the form c(initial, final) change SI
              ##units in the processed data. Defined SI units are
              ##micrometers 'mmm', milimeters 'mm', centimeters 'cm',
              ##decimeters 'dm', or meters 'm'. If NULL then no metric
              ##system is recorded.
) {
    fach <- c('factor','character')
    long <- any(sapply(rd,class)%in%fach)
    chun <- function(from,to){
        sm <- 10 ^ -c(6,3:0)
        un <- c('mmm','mm','cm','dm','m')
        names(sm) <- un
        eq <- sm[from]/sm[to]
        names(eq) <- to
        return(eq)}
    lu <- length(un)
    if(!long)
    {
        nm <- rep(names(rd),each=nrow(rd))
        splnm <- data.frame(
            do.call(rbind,
                    strsplit(nm,split = '\\.')))
        names(splnm) <- lev.nm
        splnm <- splnm[,rev(names(splnm))]
        yr <- as.numeric(rep(
            rownames(rd),ncol(rd)))
        dl <- unlist(c(rd),use.names = FALSE)
        dt <- na.omit(data.frame(x = dl,year = yr,splnm))
        rownames(dt) <- NULL
        
        if(lu == 2)
            dt[,'x'] <- with(dt, x * chun(un[1],un[2]))
        
    }
    else{
        nmx <- names(rd)
        n. <- sapply(rd,is.numeric)
        nux <- names(rd[,nmx[n.]])
        nfx <- names(rd[,!nmx%in%nux])
        nux1 <- nux[1]
        if(!is.null(which.x))nux1 <- which.x
        rd <- rd[,c(nux1,'year',nfx)]
        names(rd) <- c('x','year',nfx)
        ftosp <- lapply(rd[,rev(nfx)],as.factor)
        dsp <- split(rd,ftosp,drop = TRUE)
        dsp <- lapply(dsp,function(x)x[,c('year','x')])
        for(i in 1:length(dsp))
            names(dsp[[i]]) <- c('year',names(dsp[i]))
        fmatch <- function(tomatch.){
            Reduce(function(x,y){
                merge(x,y,all = TRUE)},tomatch.)}
        rP <- fmatch(dsp)
        rownames(rP) <- rP[,'year']
        dt <- rP[,names(rP)[-1L]]
        dt <- dt[,order(names(dt))]
        if(lu == 2)
            dt <- chun(un[1],un[2]) * dt
        attributes(dt)[['nmLong']] <- rev(nfx)
    }
    
    attributes(dt)[['un']] <- un[lu]
    return(dt)
### If \code{rd} is a ring-data frame then output is a multilevel data
### frame with the reshaped variable in the first column and years on
### the second one, followed by factor-level columns from the column
### with the lowest level units (sample/replicate) to the column with
### the higest possible level units. If \code{rd} is a multilevel data
### frame then the output is a ring-data frame (see details).
} , ex=function(){
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    
    ## Reshaping multilevel data into a ring-data frame:
    pwide <- shiftFrame(Prings05)
    str(pwide)
    ## Reshaping the ring-data frame into the initial multilevel data,
    ## and defining SI units:
    plong <- shiftFrame(pwide,un = 'mmm')
    str(plong)
})
