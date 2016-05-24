FFpure <- 
function(object, peakheight = 0.01, diapason = 440, 
    notes = NULL, interest.frqs = seq(along = object@freq),
    search.par = c(0.8, 10, 1.3, 1.7))
{
# Interesting frequencies up to 1396 (--> 1421) Hz = 66. FF (Soprano) for 512 window-width
    if(is.character(interest.frqs)){
        if(interest.frqs == "bass") interest.frqs <- 4:17
        else if(interest.frqs == "tenor") interest.frqs <- 5:28
        else if(interest.frqs == "alto") interest.frqs <- 6:37
        else if(interest.frqs == "soprano") interest.frqs <- 11:66
        else stop("Value for 'interest.frqs' not valid")
    }

    if(!is.null(notes)){
        expected.frqs <- 2^(notes / 12) * diapason
        if(!is.null(interest.frqs)){
            warning("Argument 'interest.frqs' ignored. Only one of 'interest.frqs' or 'notes' can be specified.")
            interest.frqs <- NULL
        }
        interest.frqs <- sort(unique(unlist(
            lapply(expected.frqs, 
                function(x){ 
                    temp <- abs(x - object@freq)
                    which(temp %in% sort(temp)[1:2])
                })
        )))
    }

    N <- length(object@spec)
    interesse <- numeric(N)

    for(k in 1:N){
        spec <- object@spec[[k]]
        spec[!(seq(along = spec) %in% interest.frqs)] <- 0
        index <- (max(spec[interest.frqs]) * peakheight) < spec
        peak.1 <- interest.frqs[match(spec[index][1], spec[interest.frqs])]
        peak.1.weiter <- max(peak.1 + floor(peak.1 * search.par[1]), search.par[2])  # Heuristik!
        if(peak.1.weiter >= max(interest.frqs)) peak.1.weiter <- max(interest.frqs) - 1
        temp <- (peak.1:peak.1.weiter)[(peak.1:peak.1.weiter) %in% interest.frqs]
        peak.2 <- interest.frqs[match(max(spec[temp]), spec[interest.frqs])]
        temp <- na.omit(spec[round(peak.2*search.par[3]):round(peak.2*search.par[4])])
        if(length(temp) && any(temp > max(spec[interest.frqs]) * peakheight)){
            index <- (max(spec[interest.frqs]) * peakheight * 0.1) < spec
            peak.1 <- interest.frqs[match(spec[index][1], spec[interest.frqs])]
            peak.1.weiter <- max(peak.1 + floor(peak.1 * search.par[1]), search.par[2])  # Heuristik!
            if(peak.1.weiter >= max(interest.frqs)) peak.1.weiter <- max(interest.frqs) - 1
            temp <- (peak.1:peak.1.weiter)[(peak.1:peak.1.weiter) %in% interest.frqs]
            peak.2 <- interest.frqs[match(max(spec[temp]), spec[interest.frqs])]
        }
        if((peak.2 + 1) %in% interest.frqs){
            if((peak.2 - 1) %in% interest.frqs)
                peak.2.neben <- peak.2 + if(spec[peak.2 - 1] > spec[peak.2 + 1]) -1 else 1
            else
                peak.2.neben <- peak.2 + 1
        }
        else peak.2.neben <- peak.2 - 1
#        interesse[k] <- peak.2 + (((spec[peak.2.neben] / spec[peak.2])^exp(-1)) * (peak.2.neben - peak.2) / 2)
        interesse[k] <- peak.2 + (sqrt(spec[peak.2.neben] / spec[peak.2]) * (peak.2.neben - peak.2) / 2)
    }
    return(interesse * object@freq[1])
}    


FF <- 
function(object, peakheight = 0.01, silence = 0.2, minpeak = 9, 
    diapason = 440, notes = NULL, interest.frqs = seq(along = object@freq),
    search.par = c(0.8, 10, 1.3, 1.7))
{
    FFvalue <- FFpure(object, peakheight = peakheight, diapason = diapason, 
        notes = notes, interest.frqs = interest.frqs, search.par = search.par)
    N <- length(FFvalue)
    silence <- if(silence) sort(object@energy)[max(1, round(N * silence))] else -Inf
    energy <- object@energy < silence
    for(k in 1:N){
      if(energy[k] && (sum(object@spec[[k]] > peakheight) > minpeak)){
        is.na(FFvalue[k]) <- TRUE
      }
    }
    return(FFvalue)
}

noteFromFF <-
function(x, diapason = 440, roundshift = 0)
{
    round(12 * log(x / diapason, 2) + roundshift)
}
