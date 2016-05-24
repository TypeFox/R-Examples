PeptideSpectrum <- function(expt, theory, t = 0.4, b = 5, label = "", xlim = c(100, 1500), supress = FALSE) {

    # !!! Is there a way to check theory came from pepFrag? 

    if(!is.data.frame(expt)) stop("the expt argument must be a data frame")
    if(t < 0 | t > 100) stop("the t argument must be between 0 and 100")
    if(b < 0 | b > 100) stop("the b argument must be between 0 and 100")
    if(!is.character(label)) stop("the label argument must be a character string")


    ## format experimental data frame

    expt <- data.frame(mz = expt[,1], intensity = expt[,2])
    expt$normalized <- round((expt$intensity / max(expt$intensity)) * 100)
    expt <- subset(expt, expt$normalized >= b) 


    ## find experimental-theoretical matches within the specified m/z tolerance

    matches <- vector("list")   
    for(i in 1:nrow(expt)) {
        tmp_matches <- data.frame(NULL)
        tmp_matches <- theory[expt$mz[i] >= theory$ms2mz - t & expt$mz[i] <= theory$ms2mz + t,]
        num_tmp_matches <- nrow(tmp_matches)
        expt_mz <- rep(expt$mz[i], times = num_tmp_matches)
        expt_intensity <- rep(expt$normalized[i], times = num_tmp_matches)
        matches[[i]] <- data.frame(expt_mz, expt_intensity, tmp_matches)
    }
    identifications <- as.data.frame(do.call("rbind", matches))
    if(nrow(identifications) == 0) stop("no peaks were identified")
    row.names(identifications) <- 1:nrow(identifications)
    identifications$error <- round(identifications$expt_mz - identifications$ms2mz, digits = 3)
    num_identifications <- nrow(identifications)


    # get fragmentation locations, ion types, and set colors (b/c-ions red, y/z-ions blue)
                                        
    getLocation <- function(type) {
        tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
        tmp <- strsplit(tmp, split = "[[:alpha:]]")[[1]][2]
        return(as.numeric(tmp))
    }
    location <- sapply(identifications$ms2type, getLocation)

    getSeries <- function(type) {
        tmp <- strsplit(as.character(type), split = "[[:punct:]]")[[1]][2]
        tmp <- strsplit(tmp, split = "[[:digit:]]")[[1]][1]
        return(tmp)
    }
    series <- sapply(identifications$ms2type, getSeries)

    color <- sapply(series, function(x) ifelse(x == "b" | x == "c", "red", "blue"))
    

    ## plot spectrum
                                           
    plot(expt$mz,
         expt$normalized, 
         type = "h", 
         xlim, 
         ylim = c(0, 150), 
         xlab = "m/z", 
         ylab = "intensity (%)", 
         yaxs = "i", 
         yaxt = "n")
    axis(2, at = c(0, 50, 100), labels = c(0, 50, 100))


    ## annotate identified peaks

    if(supress == FALSE) {
        x_range <- xlim[2] - xlim[1]
        y_position <- vector("numeric") 
        for(i in 1:(num_identifications - 1)) {
            if((identifications$expt_mz[i+1] - identifications$expt_mz[i]) / x_range < 0.025 & 
                    all.equal(identifications$expt_intensity[i], 100) != TRUE) {
                y_position[i] <- identifications$expt_intensity[i] + 30
                lines(rep(identifications$expt_mz[i], 2), 
                      c(identifications$expt_intensity[i], 
                      identifications$expt_intensity[i] + 30), 
                      lty = "dotted", 
                      col = color[i])  
            } else y_position[i] <- identifications$expt_intensity[i]
        }    
        y_position[num_identifications] <- identifications$expt_intensity[num_identifications]
        text(identifications$expt_mz, y_position + 12, labels = identifications$ms2type, col = color, srt = 90, cex = 0.9)
    }


    ## draw sequence with fragmentation locations

    seq_vector <- strsplit(as.character(identifications$ms1seq)[1], split = "")[[1]]
    num_residues <- length(seq_vector)   
    plot.window(xlim = c(1, 20), ylim = c(0, 10))
    text(c(1:num_residues), 9, labels = seq_vector)
    for(i in 1:length(series)) {                                                       
        if(series[i] == "b" | series[i] == "c") lines(c(location[i] + 0.25, location[i] + 0.55), c(8.5, 9.5), col = "red") 
        else lines(c(num_residues - location[i] + 0.45, num_residues - location[i] + 0.75), c(8.5, 9.5), col = "blue")
    }

    text(18, 9, label)
    return(identifications)
}
