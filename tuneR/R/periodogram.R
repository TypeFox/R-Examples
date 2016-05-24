setGeneric("periodogram",
    function(object, ...) standardGeneric("periodogram"))

setMethod("periodogram", signature(object = "WaveGeneral"), 
function(object, width = length(object), overlap = 0,
    starts = NULL, ends = NULL, taper = 0, normalize = TRUE, 
    frqRange = c(-Inf, Inf), ...)
{
    validObject(object)

    if(width > length(object))
        stop("width must be less or equal length(object)")
    if(nchannel(object) > 1) 
        stop("Processing more than one channel is not yet implemented...")
    
    testwidth <- 2^ceiling(log(width, 2))
    if(width !=  testwidth) {
        width <- testwidth
        warning("'width' must be a potence of 2, hence using the ceiling: ", 
            width, "\n")
    }
        
    Wspec <- new("Wspec")
    Wspec@stereo <- FALSE
    Wspec@samp.rate <- object@samp.rate
    Wspec@taper <- taper
    temp <- width / 2
    Wspec@freq <- object@samp.rate * seq(1, width / 2) / width

    wo <- width - overlap
    lo <- length(object)
    lw <- lo - width
    n <- lw %/% wo
    add <- lw - n*wo
    lo <- lo + add
    dat <- c(if(is(object, "Wave")) object@left else object@.Data[,1], rep(0, add))
    dat <- dat - mean(dat)
    if(normalize) 
        dat <- dat / max(abs(dat))
        
    n <- n + 1
    if(is.null(starts) && is.null(ends)){
        starts <- seq(1, lo-width+1, by = wo)
        ends <- seq(width, lo, by = wo)
    }
    else{
        if(is.null(starts)) starts <- ends - width
        if(is.null(ends)) ends <- starts + width
    }

    temp <- spec.pgram(dat[starts[1]:ends[1]], taper = taper,
        pad = 0, fast = TRUE, demean = FALSE, detrend = FALSE,
        plot = FALSE, na.action = na.fail, ...)

    Wspec@kernel <- temp$kernel
    Wspec@df <- temp$df
    Wspec@starts <- starts
    Wspec@width <- width
    Wspec@overlap <- overlap
    Wspec@normalize <- normalize
    
    spec <- vector(n, mode = "list")
    spec[[1]] <- temp$spec
    for(i in (seq(along = starts[-1]) + 1)){
        spec[[i]] <- spec.pgram(dat[starts[i]:ends[i]], taper = taper, 
            pad = 0, fast = TRUE, demean = FALSE, detrend = FALSE,
            plot = FALSE, na.action = na.fail, ...)$spec
        
    }

    store <- (Wspec@freq >= min(frqRange)) & (Wspec@freq <= max(frqRange))

    Wspec@spec <- if(normalize)
        lapply(spec, function(x){
            sx <- sum(x)
            lx <- length(x)
            (if(!sx) rep(1 / lx, lx) else x / sum(x))[store]
        })
        else spec[store]

    Wspec@variance <- mapply(function(x,y) var(dat[x:y]), starts, ends)
    Wspec@energy <- 20 * log10(mapply(function(x,y) sum(abs(dat[x:y])), starts, ends))

    Wspec@freq <- Wspec@freq[store]
    return(Wspec)
})



setMethod("periodogram", signature(object = "character"), 
function(object, width, overlap = 0, from = 1, to = Inf, 
    units = c("samples", "seconds", "minutes", "hours"), 
    downsample = NA, channel = c("left", "right"), 
    pieces = 1, ...)
{
    channel <- match.arg(channel)
    units <- match.arg(units)
    ### the simple part for only one piece:
    if(pieces == 1){
        Wobj <- readWave(filename = object, 
            from = from, to = to, units = units)
        Wobj <- mono(Wobj, which = channel)
        if(!is.na(downsample))
            Wobj <- downsample(Wobj, downsample)
        if(missing(width)) 
            width <- length(Wobj)
        WSobj <- periodogram(Wobj, width = width, ...)
        return(WSobj) # finished
    }
    
    ### only if more than one piece:
    header <- readWave(filename = object, header = TRUE)
    fctr <- switch(units,
                   samples = 1,
                   seconds = header$sample.rate,
                   minutes = header$sample.rate * 60,
                   hours = header$sample.rate * 3600)
    if(fctr > 1) {
        from <- from * fctr + 1
        to <- to * fctr
    }
    to <- min(to, header$samples)
    slength <- to - from + 1
    
    dsfactor <- if(!is.na(downsample)) header$sample.rate / downsample else 1
    ovlength <- (pieces - 1) * overlap
    each <- floor((slength + ovlength) / pieces / width / dsfactor) * width * dsfactor
    last <- floor(((slength + ovlength) - each * (pieces - 1)) / width / dsfactor) * width * dsfactor
    to <- cumsum(c(rep(each, pieces - 1), last))
    from1 <- from
    from <- c(0, to[-length(to)]) + from
    to <- to + from1 - 1
    ovsub <- seq(from = 0, by = overlap * dsfactor, length.out = pieces)
    to <- to - ovsub
    from <- from - ovsub
    for(i in 1:pieces){
        Wobj <- readWave(filename = object, 
            from = from[i], to = to[i], units = "samples")
        Wobj <- mono(Wobj, which = channel)
        if(!is.na(downsample))
            Wobj <- downsample(Wobj, downsample)
        WSobj <- periodogram(Wobj, width = width, overlap = overlap, ...)
        if(i == 1){
            WspecObj <- WSobj
        }else{
            WspecObj@spec <- c(WspecObj@spec, WSobj@spec)
            WspecObj@starts <- c(WspecObj@starts, WSobj@starts + from[i] - 1)
            WspecObj@variance <- c(WspecObj@variance, WSobj@variance)
            WspecObj@energy <- c(WspecObj@energy, WSobj@energy)
        }
    }
    return(WspecObj)
})
