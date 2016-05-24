mySMHandler <- function(c) {
    pkg <- startupPackage(c)
    endline <- startupEndline(c)
    linestarter <- paste(":",pkg,"> ", sep ="")
    ## effective number of chars / line:
    nceff <- getOption("width")-nchar(linestarter)-2

    ## make message a list with one char vector entry per line
    ## (split inbetween by "\n")
    wds0 <- unlist(conditionMessage(c))
    wdsA <- as.list(unlist(strsplit(wds0,"\n")))

    for( j in 1: length(wdsA)){
        wds <- unlist(strsplit(wdsA[[j]]," "))
        lwds <- length(wds)
        if(lwds){
            i <- 1
            tryc <- linestarter
            wc <- 0
            while(i<=lwds){
               tryold <- tryc
               tryc <- paste(tryc,wds[i])
               nt <- nchar(tryc)
               wc <- wc+1
               if(nt>= nceff){
                  if(wc == 1) {
                    ### one word already too long; print it nevertheless
                    writeLines(tryc,stderr())
                    i <- i+1
                  }else{
                    ### only print the string upto the current word
                    writeLines(tryold,stderr())
                  }
                  ### set word counter to 0
                  wc <- 0
                  tryc <- linestarter
               }else{
                  if(i==lwds)
                    writeLines(tryc,stderr())
                  i <- i+1
               }
            }
            if(j < length(wdsA) &&!endline) writeLines(linestarter, stderr())
        }
    }
    if(endline) writeLines("",stderr())
    else writeLines(linestarter,stderr())
}
