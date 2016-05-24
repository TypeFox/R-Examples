quantize <- function(notes, energy, parts){
    if(missing(notes))
        notes <- rep(-100, length(energy))
    if(missing(energy))
        energy <- rep(-100, length(notes))
    if(length(energy)!=length(notes))
        stop("energy and notes must be vectors of same length")
       
    starts <- round(seq(1, length(notes) + 1, length = parts + 1))
    ends <- starts[-1] - 1
    starts <- starts[-(parts + 1)]
    

    notes[is.na(notes)] <- -100
    notes <- mapply(function(x, y) notes[x:y], x = starts, y = ends, SIMPLIFY = FALSE)
    notes <- sapply(notes, function(x){
                    x <- table(x)
                    as.integer(names(x[which.max(x)]))})
    is.na(notes[notes == -100]) <- TRUE

    energy[is.na(energy)] <- -100
    energy <- mapply(function(x, y) energy[x:y], x = starts, y = ends, SIMPLIFY = FALSE)
    energy <- sapply(energy, mean)
    is.na(energy[energy == -100]) <- TRUE

    return(list(notes=notes, energy=energy))
}


quantMerge <- function(notes, minlength, barsize, bars){
    lengthunit <- bars * barsize
    notes <- quantize(notes, parts = lengthunit)$notes
    notes[is.na(notes)] <- -100
    duration <- rep(1, lengthunit)
    index <- punctuation <- slur <- logical(lengthunit)
    full <- 2^(0:6)
    punctated <- full * 1.5
    
    i <- 1
    barlength <- 0
    while(i <= lengthunit){
        index[i] <- TRUE
        restbarlength <- barsize - barlength
        temp <- diff(notes[i:lengthunit])
        notelength <- 
            if(length(temp) && temp[1] == 0) 
                max(cumsum(ifelse(temp == 0, 1, -Inf))) + 1
            else 1
        notelength <- min(restbarlength, notelength)
        while(!(notelength %in% c(full, punctated)))
            notelength <- notelength - 1 
        ## FixMe: Need to calculate slurs!!!
        barlength <- barlength + notelength
        j <- i + notelength
        if(notelength %in% punctated){
            punctuation[i] <- TRUE
            notelength <- notelength / 1.5
        }
        duration[i] <- minlength / notelength
        if(barlength == barsize) 
            barlength <- 0
        i <- j
    }
    
    is.na(notes[notes == -100]) <- TRUE
    notes <- data.frame(note = notes, duration = duration, 
                        punctuation = punctuation, slur = FALSE)
    notes[index,]
}
