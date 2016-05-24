### Cropping device ############################################################

##' @title Cropping Device
##' @param file An output file name (including extension) or NULL (no cropping)
##' @param warn A logical indicating whether warnings are given
##' @param ... additional arguments passed to dev.off()
##' @return invisible()
##' @author Marius Hofert
dev.off.crop <- function(file=NULL, warn=TRUE, ...)
{
    ## Close device
    dev.off(...)

    ## Check
    if(!is.null(file)) { # cropping should be done
        if(.Platform$OS.type != "unix") {
            if(warn) warning("dev.off.crop() only works on Unix(-like) systems; no cropping is done.")
        } else { # Unix-like systems

            ## Cropping
            f <- file.path(getwd(), file) # get whole path (incl. file name) to file
            file.ext <- file_ext(file) # could be pdf, Pdf, PDF...
            file.ext.low <- tolower(file.ext) # file ending in lower case
            switch(file.ext.low,
                   "ps" =, "eps" = {
                       if(!command.exists("epstool")) { # check if epstool is available
                           if(warn) warning("Command 'epstool' for cropping PS files not found; no cropping is done.")
                       } else { # crop
                           ## Note: epstool cannot directly output to file with same name
                           f. <- add_to_name(f, s="_crop")
                           system(paste("epstool --copy --bbox --quiet", f, f., "1>/dev/null 2>&1; mv", f., f))
                       }
                   },
                   "pdf" = {
                       if(!command.exists("pdfcrop")) { # check if pdfcrop is available
                           if(warn) warning("Command 'pdfcrop' for cropping PDF files not found; no cropping is done.")
                       } else { # crop
                           system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"))
                       }
                   },
                   "png" = {
                       if(!command.exists("mogrify")) { # check if mogrify is available
                           if(warn) warning("Command 'mogrify' for cropping PNG files not found; no cropping is done.")
                       } else { # crop
                           system(paste("mogrify -trim", f, f, "1>/dev/null 2>&1"))
                       }
                   },
                   stop("Non-supported file extension."))

        }
    }

    ## Return
    invisible()
}
