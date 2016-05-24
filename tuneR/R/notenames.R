notenames <- function(notes, language = c("english", "german")){
    language <- match.arg(language)
    if(!is.numeric(notes))
        stop("notes must be integer values")
    rg <- range(notes)
    
    ## How is the note called?
    name <- switch(language,
        english = c("C",  "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B",
            unlist(lapply(c("", "'", "''", "'''", "''''"), 
                function(x) paste(c("c", "c#", "d", "d#", "e", "f", "f#", "g", "g#", "a", "a#", "b"), x, sep="")))),
        german = c("C",  "Cis", "D", "Dis", "E", "F", "Fis", "G", "Gis", "A", "Ais", "H",
            unlist(lapply(c("", "'", "''", "'''", "''''"), 
                function(x) paste(c("c", "cis", "d", "dis", "e", "f", "fis", "g", "gis", "a", "ais", "h"), x, sep="")))),
        stop("currently only notenames in english and german are implemented")
    )
    ## Now we know notes between -33 and 38:
    low <- -33 - rg[1]
    high <- rg[2] - 38

    ## cutting off below and above:
    if(low < 0) name <- name[(-1):low]
    if(high < 0) name <- name[(-length(name)-high-1):(-length(name))]
    ## adding unknown names:
    name <- c(if(low > 0) rep(" ", low), name, if(high > 0) rep(" ", high))
    
    ## now sorting the stuff according to input vector
    notes <- notes - min(notes) + 1
    return(name[notes])
}
