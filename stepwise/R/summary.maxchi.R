"summary.maxchi" <-
    function(object, ...) {
    sigSites <- 1
    siteWinlocs <- -1
    siteChisqs <- -1.0
    if (length(object$winlocs)<0) {
        cat("\nNo significant MaxChi statistics found.\n\n")
        return()
    }
    siteChisqs[1] <- object$chisqs[1];
    siteWinlocs[1] <- object$winlocs[1];
    pairs <- print.pair(object$pairmem1[1], object$pairmem2[1])
    for (i in 2:length(object$winlocs)) {
        if(object$winlocs[i] != siteWinlocs[sigSites]) { # found a new site 
            sigSites <- sigSites+1;       
            siteChisqs[sigSites]=object$chisqs[i];
            siteWinlocs[sigSites]=object$winlocs[i]; 
            pairs[sigSites] <- print.pair(object$pairmem1[i], object$pairmem2[i]);
        } 
        else {
            if(object$chisqs[i] > siteChisqs[sigSites]) {
                siteChisqs[sigSites] <- object$chisqs[i]
                pairs[sigSites] <- print.pair(object$pairmem1[i],object$pairmem2[i])
            } 
            else if(object$chisqs[i] == siteChisqs[sigSites]) {
            # concatenate this pair on end of current list of pairs 
                pairs[sigSites] <- paste(pairs[sigSites],"\n\t\t\t   ", 
                  print.pair(object$pairmem1[i],object$pairmem2[i]), sep="")
            }
        } # end of else
    } # end of for loop
    cat("--------------------------------------------------\n");
    cat("There were", sigSites,"site-specific MaxChi statistics significant at the\n");
    cat("10 percent level (90th percentile = ", sprintf("%5.3f", object$quants[1]), 
    ", 95th percentile = ",sprintf("%5.3f", object$quants[2]),"):\n\n", sep="");
    cat ("Number Location  MaxChi   pairs\n")
    for (i in 1:length(siteWinlocs)) {
        if(siteChisqs[i]>object$quants[1]) star<-"*" else star<-" "
        cat(sprintf("%6d  %7d  %5.3f%s   %s\n",i, as.integer(siteWinlocs[i]), siteChisqs[i], star,  pairs[i]));
    }
    cat("--------------------------------------------------\n");
    cat("Notes - \"Location\" is the polymorphic site just before the proposed breakpoint.\n");
    cat("      - MaxChi statistics significant at the 5 percent level indicated by a * \n\n");
    
    a <- list(siteWinlocs = siteWinlocs, siteChisqs = siteChisqs, pairs = pairs)
    invisible(a)
}
