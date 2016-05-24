sim = function(N, Freqs, rel = "UN", save = FALSE, strPath = "",
               strVer = "", BlockSize = N/100, fileName = NULL){

    rel = toupper(rel)
    if(!grepl("(UN|FS|PC)", rel, ignore.case = T, perl = TRUE)){
        stop("Unrecognized relationship. Must be one of 'U', 'FS', 'PC'")
    }


    f1 = NULL

    if(save){
        if(is.null(fileName)){
            if(nchar(strVer) > 0){
                fileName = paste(strPath, sprintf("results-sim-%s-%d-%s.csv.gz",
                                 rel, N, strVer),sep = "")
            }else{
                fileName = paste(strPath, sprintf("results-sim-%s-%d.csv.gz",
                                 rel, N),sep = "")
            }
        }

        f1 = gzfile(fileName, 'w')
        if(!isOpen(f1)){
            stop(paste("Couldn't open", fileName, "for writing"))
        }
    }

    ## make space for results
    lrsibs = rep(0, N)
    lrpc = lrsibs
    ibs = lrsibs

    ## progressBar

    pb = txtProgressBar(min = 0, max = 100, style = 3)
    step = N/100

    ## number of blocks (usually 100), but...
    numBlocks = N/BlockSize

    if(numBlocks - floor(numBlocks) > 0){
        stop("BlockSize must be a whole integer factor of sample size")
    }


    j = 0
    k = 0
    setTxtProgressBar(pb, k)

    ## predetermine some inputs
    nLoci = length(Freqs$loci)
    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)

    if(rel == "UN"){ ## put this out here so that it is faster
        for(block in 1:numBlocks){
            P = randomProfilePairs(Freqs, BlockSize)

            for(i in 1:BlockSize){
                i1 = (block-1)*BlockSize + i

                p1 = P[[i]]$prof1
                p2 = P[[i]]$prof2
                ibs[i1] = IBS(p1, p2, nLoci)
                lrsibs[i1] = lrSib(p1,p2, NULL, nLoci, f = f, n = n)
                lrpc[i1] = lrPC(p1, p2, NULL, nLoci, f = f, n = n)

                if(j == step){
                    k = k + 1
                    setTxtProgressBar(pb, k)
                    j = 0
                }

                j = j + 1

                if(save){ ## bit expensive to check every iteration
                    line = sprintf("%10.8E,%10.8E,%d", lrsibs[i1],
                                                       lrpc[i1], ibs[i1])
                    writeLines(line, f1)
                }
            }
         }
    }else if(rel == "FS"){
        for(block in 1:numBlocks){
            P = randomSibPairs(Freqs, BlockSize)

            for(i in 1:BlockSize){
                i1 = (block-1)*BlockSize + i

                p1 = P[[i]]$sib1
                p2 = P[[i]]$sib2
                ibs[i1] = IBS(p1, p2, nLoci)
                lrsibs[i1] = lrSib(p1,p2, NULL, nLoci, f = f, n = n)
                lrpc[i1] = lrPC(p1, p2, NULL, nLoci, f = f, n = n)

                if(j == step){
                    k = k + 1
                    setTxtProgressBar(pb, k)
                    j = 0
                }

                j = j + 1

                if(save){ ## bit expensive to check every iteration
                    line = sprintf("%10.8E,%10.8E,%d", lrsibs[i1],
                                                       lrpc[i1], ibs[i1])
                    writeLines(line, f1)
                }
            }
         }
    }else if(rel == "PC"){
        for(block in 1:numBlocks){
            P = randomPCPairs(Freqs, BlockSize)

            for(i in 1:BlockSize){
                i1 = (block-1)*BlockSize + i

                p1 = P[[i]]$parent
                p2 = P[[i]]$child
                ibs[i1] = IBS(p1, p2, nLoci)
                lrsibs[i1] = lrSib(p1,p2, NULL, nLoci, f = f, n = n)
                lrpc[i1] = lrPC(p1, p2, NULL, nLoci, f = f, n = n)

                if(j == step){
                    k = k + 1
                    setTxtProgressBar(pb, k)
                    j = 0
                }

                j = j + 1

                if(save){ ## bit expensive to check every iteration
                    line = sprintf("%10.8E,%10.8E,%d", lrsibs[i1],
                                                       lrpc[i1], ibs[i1])
                    writeLines(line, f1)
                }
            }
        }
    }

    if(save){
        close(f1)
        cat(paste("Results written to", fileName, "\n"))
    }

    invisible(data.frame(sib = lrsibs, pc = lrpc, ibs = ibs))
}



