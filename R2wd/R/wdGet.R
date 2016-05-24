wdGet<-
    function (filename = NULL, path = getwd(), method="rcom",visible = TRUE)
{
    ## this first section selects the com client
    ## if rcom is selected
    if (method == "rcom") {
        if (!require(rcom)) {
            warning("The package rcom is unavailable.")
            if (require(RDCOMClient)) {
                warning("Using RDOMClient package instead of rcom.")
                client<-"RDCOMClient"
            }
            else {
                stop("Neither rcom or RDCOMClient packages are installed.")
                client<-"none"
            }
        }
        else {
            client<-"rcom"
            if ("package:RDCOMClient" %in% search()) {
                warning("\nUsing rcom package. Detaching RDCOMClient package to avoid conflicts.")
                try(detach("package:RDCOMClient"))
            }
        }

    }
    ## if RDCOMClient is selected
    if (method == "RDCOMClient") {
        if (!require(RDCOMClient)) {
            stop("The package RDCOMClient is unavailable. \n \n\t\tTo install RDCOMClient use:\n \n\t\tinstall.packages('RDCOMClient' repos = 'http://www.omegahat.org/R')")
        }
        client<-"RDCOMClient"
        if ("package:rcom" %in% search()) {
            warning("Using RDCOMClient package. Detaching rcom package to avoid conflicts.")
            try(detach("package:rcom"))
        }
    }

    switch(client,
           "rcom"={
               wdapp <- comGetObject("Word.Application")
               if (is.null(wdapp)) wdapp <- comCreateObject("Word.Application")
           },
           "RDCOMClient"={
               wdapp<-COMCreate("Word.Application")
           },
           none=stop("no client")
           )
    if (visible) wdapp[["visible"]] <- TRUE

    ## this section loads filename or adds an empty document

    ## if filename is NULL then add a new document otherwise open filename
    if (is.null(filename)) {
        ## if the word application is freshly opened, add a document
        if (wdapp[["Documents"]][["Count"]]==0) wdapp[["Documents"]]$Add()
        ## otherwise write to the document which is currently open
    } else {
        ## if the word application is already open see if the desired
        ## document is already open, if yes, activate it.
        wddocs<-wdapp[["Documents"]]
        found<-FALSE
        if (wddocs[["Count"]]>0) {
            for (i in 1:wddocs[["Count"]]){
                wddoc<-wddocs$Item(i)
                if (wddoc[["Name"]]==filename){
                    wddoc$Activate()
                    found<-TRUE
                    break
                    }
            }
        }
        ## if not found then try to open
        if (!found){
              wddoc<-try(wdapp[["Documents"]]$Open(paste(path, filename, sep = "/")))
              ## if file does not exist and there are no documents, close word
              if (class(wddoc)=="try-error" | is.null(wddoc)) {
                  if(wddocs[["Count"]]==0) wdapp$Quit()
                  print(paste("File",paste(path, filename, sep = "/"),"not found"))
              }
          }
    }
    .R2wd <<- wdapp
    invisible(wdapp)
}

