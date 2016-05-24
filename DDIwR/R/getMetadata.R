require(XML)
getMetadata <- function(xmlpath, OS = "windows", saveFile=FALSE, ...) {
    
    # TODO: detect DDI version or ask the version through a dedicated argument
    
    
    other.args <- list(...)
    enter <- getEnter(OS=OS)
    
    fromsetupfile <- FALSE
    if ("fromsetupfile" %in% names(other.args)) {
        fromsetupfile <- other.args$fromsetupfile
    }
    
    tp <- treatPath(xmlpath, type="XML")
    
    currdir <- getwd()
    # if (saveFile) {
        setwd(tp$completePath)
    # }
    
    singlefile <- length(tp$files) == 1
    
    if (!fromsetupfile) {
        cat("Processing:\n")
    }
    
    for (ff in seq(length(tp$files))) {
        if (!fromsetupfile) {
            cat(tp$files[ff], "\n")
        }
        
        if (saveFile) {
            sink(paste(tp$filenames[ff], "R", sep="."))
        }
        
        dd <- xmlTreeParse(tp$files[ff])$doc$children$codeBook
        
        
        #### !!! ####
        # NEVER use getNodeSet() it's toooooo slooooow!!!
        # use instead xmlElementsByTagName()
    
        dd <- xmlElementsByTagName(dd, "dataDscr")[[1]]
        dd <- xmlElementsByTagName(dd, "var")
              
        xmlVarNames <- as.vector(sapply(dd, xmlGetAttr, "name"))
        # return(drop(xmlVarNames))
        
        metadata <- list()
        metadata$varlab <- list()
        metadata$vallab <- list()
        
        if (saveFile) {
            cat("metadata <- list()", enter)
            cat("metadata$varlab <- list()", enter)
            cat("metadata$vallab <- list()", enter, enter)
        }
        
        for (i in seq(length(dd))) {
            
            # metadata$varlab[[xmlVarNames[i]]] <- xmlValue(getNodeSet(dd[[i]], "//labl[@level='variable']")[[1]])
            varlab <- xmlValue(xmlElementsByTagName(dd[[i]], "labl")[[1]])
            varlab <- gsub("\"", "'", varlab)
            varlab <- gsub("\\\\", "/", varlab)
            metadata$varlab[[xmlVarNames[i]]] <- varlab
            
            if (saveFile) {
                cat(paste("metadata$varlab$", xmlVarNames[i], " <- \"", varlab, "\"", enter, sep=""))
            }
            
            #vallabs <- unlist(lapply(getNodeSet(dd[[i]], "//labl[@level='category']"), xmlValue))
            vallabs <- xmlElementsByTagName(dd[[i]], "catgry")
            
            if (length(vallabs) > 0) {
                
                # metadata$vallab[[xmlVarNames[i]]] <- unlist(lapply(getNodeSet(dd[[i]], "//catValu"), xmlValue))
                values <- as.vector(unlist(lapply(lapply(vallabs, xmlElementsByTagName, "catValu"), function(x) {
                    return(xmlValue(x[[1]][[1]]))
                })))
                values <- gsub("\"", "'", values)
                values <- gsub("\\\\", "/", values)
                
                labl <- as.vector(lapply(vallabs, xmlElementsByTagName, "labl"))
                havelbls <- unlist(lapply(labl, function(x) length(x) > 0))
                
                values <- values[havelbls]
                labl <- labl[havelbls]
                
                if (length(values) > 0) {
                    metadata$vallab[[xmlVarNames[i]]] <- values
                    testNum <- tryCatch(as.numeric(values),
                                        warning = function(x) {
                                                     return("...string...!!!")
                                        })
                    
                    if (all(testNum != "...string...!!!")) {
                        metadata$vallab[[xmlVarNames[i]]] <- testNum
                        
                        if (saveFile) {
                            cat(paste("metadata$vallab$", xmlVarNames[i], " <- c(", 
                                paste(testNum, collapse=", "), ")", enter, sep=""))
                        }
                        
                        justlbls <- as.vector(unlist(lapply(labl, function(x) {
                            return(xmlValue(x[[1]][[1]]))
                        })))
                        
                        justlbls <- gsub("\"", "'", justlbls)
                        justlbls <- gsub("\\\\", "/", justlbls)
                        
                        names(metadata$vallab[[xmlVarNames[i]]]) <- justlbls
                        
                        if (saveFile) {
                            cat(paste("names(metadata$vallab$", xmlVarNames[i], ") <- c(\"",
                                    paste(justlbls, collapse="\", \""), "\")", enter, sep=""))
                        }
                    }
                    else {
                        
                        justlbls <- as.vector(unlist(lapply(lapply(vallabs, xmlElementsByTagName, "catValu"), function(x) {
                            return(xmlValue(x[[1]][[1]]))
                        })))
                        justlbls <- gsub("\"", "'", justlbls)
                        justlbls <- gsub("\\\\", "/", justlbls)
                        
                        if (saveFile) {
                            cat(paste("metadata$vallab$", xmlVarNames[i], " <- c(\"",
                                    paste(justlbls, collapse="\", \""), "\")", enter, sep=""))
                        }
                    }
                }
            }
            cat(enter)
        }
        if (saveFile) {
            sink()
        }
    }
    
    setwd(currdir)
    if (singlefile) {
        return(invisible(metadata))
    }
}

