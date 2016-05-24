## NOTE:
##   We do not support old-version DSSP any longer
##   Please update your DSSP program to the newest version
"dssp.pdb" <-
  function (pdb, exefile = "dssp", resno=TRUE, full=FALSE, verbose=FALSE, ...) {

    ## Log the call
    cl <- match.call()
    
    ## Check if the program is executable
    os1 <- .Platform$OS.type
    status <- system(paste(exefile, "--version"),
                     ignore.stderr = TRUE, ignore.stdout = TRUE)
    
###    if(!(status %in% c(0,1)))
###      stop(paste("Launching external program 'DSSP' failed\n",
###                 "  make sure '", exefile, "' is in your search path", sep=""))
    
    ## check atom composition - need backbone atoms to continue SSE analysis
    checkatoms <- TRUE
    if(checkatoms) {
      inds <- atom.select(pdb, "backbone", verbose=verbose)
      tmp <- trim.pdb(pdb, inds)
      
      resid <- paste(tmp$atom$resno, tmp$atom$chain, sep="-")
      musthave <- c("C", "CA", "N", "O")

      incomplete <- sapply(unique(resid), function(x) {
        inds <- which(resid==x)
        elety <- sort(tmp$atom$elety[inds])
        if(!all(musthave %in% elety))
          return(TRUE)
        else
          return(FALSE)
      })

      if(all(incomplete))
        stop("No residues found with a complete set of backbone atoms")
      if(any(incomplete))
        warning(paste("Residues with missing backbone atoms detected:",
                      paste(unique(resid)[incomplete], collapse=", "), 
                      collapse=" "))
    }
    
    infile <- tempfile()
    outfile <- tempfile()
    write.pdb(pdb, file = infile)
    cmd <- paste(exefile, infile, outfile)

    if(verbose)
      cat(paste("Running command:\n ", cmd , "\n"))
    
    if(os1 == "windows")
      success <- shell(cmd, 
                       ignore.stderr = !verbose, ignore.stdout = !verbose)
    else 
      success <- system(cmd,
                        ignore.stderr = !verbose, ignore.stdout = !verbose)
    
    if(success!=0)
      stop(paste("An error occurred while running command\n '",
                 cmd, "'", sep=""))
    
##
## For Debug (Tue Aug  3 18:22:11 PDT 2010)
##  -- Following multi chain error report from Heiko Strathmann
##    outfile <- "2jk2.dssp"
##    outfile <- "4q21.dssp"
##
    trim <- function(s) {
      s <- sub("^ +", "", s)
      s <- sub(" +$", "", s)
      s[(s == "")] <- NA
      s
    }

    split.line <- function(x, split=" ") {
      tmp <- unlist(strsplit(x, split=split))
      inds <- which(tmp!="")
      return(trim(tmp[inds]))
    }

    raw.lines <- readLines(outfile)
    unlink(c(infile, outfile))
    type <- substring(raw.lines, 1, 3)
    raw.lines <- raw.lines[-(1:which(type == "  #"))]

    ## delete chain breaking lines
    aa <- substring(raw.lines, 14, 14)
    if(any(aa == "!"))
       raw.lines <- raw.lines[-which(aa == "!")]

    cha      <- substring(raw.lines, 12, 12)
    sse      <- substring(raw.lines, 17, 17)
    res.name <- substring(raw.lines, 14, 14)
    
    res.id  <- as.numeric(substring(raw.lines, 1, 5))   ## dssp residue IDs
    res.num <- as.numeric(substring(raw.lines, 6, 10))  ## Residue numbers
    res.insert <- substring(raw.lines, 11, 11)          ## Insertion codes
    res.ind <- 1:length(res.num)                        ## Internal indices

    if(any(res.insert!=" ")) {
       if(resno) {
         warning("Insertions are found in PDB: Residue numbers may be incorrect.
                Try again with resno=FALSE")
       }
       else { 
         ii <- diff(res.num)
         ii[ii==0] <- 1     #Consecutive numbers at insertion residues
         ii[ii<0] <- 2      #Jumps at possible chain termination
         res.num <- res.num[1] + c(0, cumsum(ii))
       }
    }
    
    if(full) {
      ## Difference between sse res id and internal res indices
      diff        <- res.id - res.ind
      names(diff) <- res.id

      ## Beta bridge partner residue ids
      bp1         <- as.numeric(substring(raw.lines, 26, 29))
      bp2         <- as.numeric(substring(raw.lines, 30, 33))
      bp1[bp1==0] <- NA
      bp2[bp2==0] <- NA

      ## Convert from dssp SSE residue IDs to internal residue indices
      bp1[ !is.na(bp1) ] <- as.vector(bp1[ !is.na(bp1) ] - diff[as.character(bp1[!is.na(bp1)])])
      bp2[ !is.na(bp2) ] <- as.vector(bp2[ !is.na(bp2) ] - diff[as.character(bp2[!is.na(bp2)])])

      ## H-bond records
      hbonds <- split.line(split.line(substring(raw.lines, 40, 83), split=","), split=" ")
      hbonds <- matrix(as.numeric(hbonds), ncol=8, byrow=TRUE)
      hbonds <- as.data.frame(hbonds)
      
      for(i in seq(1,7,by=2)) {
        hbonds[[i]][ which(hbonds[[i]]==0) ] <- NA
        
        ## Convert from relative to absolute residue numbering
        hmm <- res.id + hbonds[[i]]
        
        ## Convert from dssp SSE residue IDs to internal residue indices
        hbonds[[i]][ !is.na(hmm) ] <- as.vector(hmm[ !is.na(hmm) ] -
                                           diff[as.character(hmm[!is.na(hmm)])])
      }
      
      ## Bind bridge pair and H-bond records to one matrix
      hbonds <- cbind(bp1, bp2, hbonds)
      cnames <- c("BP1", "BP2", "NH-O.1", "E1", "O-HN.1", "E2", "NH-O.2", "E3", "O-HN.2", "E4")
      colnames(hbonds) <- cnames

      if(resno) {
        ## 2 col matrix mapping the res.ind's to res.num and chain id
        tmp.map            <- cbind(res.num, cha)
        row.names(tmp.map) <- res.ind

        ## Add an additional matrix holding the Chain IDs
        hbonds           <- cbind(hbonds, data.frame(matrix(NA, ncol=6, nrow=nrow(tmp.map)),
                                                     stringsAsFactors=FALSE))
                                                     
        colnames(hbonds) <- c(cnames,
                              "ChainBP1", "ChainBP2", "Chain1", "Chain2", "Chain3", "Chain4")

        ## Add chain IDs for each entry
        tmp.inds                    <- which(!is.na(hbonds[,"BP1"]))
        tmp.names                   <- as.character(hbonds[tmp.inds,"BP1"])
        hbonds[tmp.inds,"BP1"]      <- as.numeric(tmp.map[tmp.names, "res.num"])
        hbonds[tmp.inds,"ChainBP1"] <- tmp.map[tmp.names, "cha"]

        tmp.inds                    <- which(!is.na(hbonds[,"BP2"]))
        tmp.names                   <- as.character(hbonds[tmp.inds,"BP2"])
        hbonds[tmp.inds,"BP2"]      <- as.numeric(tmp.map[tmp.names, "res.num"])
        hbonds[tmp.inds,"ChainBP2"] <- tmp.map[tmp.names, "cha"]

        tmp.inds                    <- which(!is.na(hbonds[,"NH-O.1"]))
        tmp.names                   <- as.character(hbonds[tmp.inds,"NH-O.1"])
        hbonds[tmp.inds,"NH-O.1"]   <- as.numeric(tmp.map[tmp.names, "res.num"])
        hbonds[tmp.inds,"Chain1"]   <- tmp.map[tmp.names, "cha"]

        tmp.inds                    <- which(!is.na(hbonds[,"O-HN.1"]))
        tmp.names                   <- as.character(hbonds[tmp.inds,"O-HN.1"])
        hbonds[tmp.inds,"O-HN.1"]   <- as.numeric(tmp.map[tmp.names, "res.num"])
        hbonds[tmp.inds,"Chain2"]   <- tmp.map[tmp.names, "cha"]

        tmp.inds                    <- which(!is.na(hbonds[,"NH-O.2"]))
        tmp.names                   <- as.character(hbonds[tmp.inds,"NH-O.2"])
        hbonds[tmp.inds,"NH-O.2"]   <- as.numeric(tmp.map[tmp.names, "res.num"])
        hbonds[tmp.inds,"Chain3"]   <- tmp.map[tmp.names, "cha"]

        tmp.inds                    <- which(!is.na(hbonds[,"O-HN.2"]))
        tmp.names                   <- as.character(hbonds[tmp.inds,"O-HN.2"])
        hbonds[tmp.inds,"O-HN.2"]   <- as.numeric(tmp.map[tmp.names, "res.num"])
        hbonds[tmp.inds,"Chain4"]   <- tmp.map[tmp.names, "cha"]

        ## Set row names to "RESNUM-CHAINID"
        row.names(hbonds)           <- apply(tmp.map, 1, paste, collapse="-")
      }
    }
    else {
      hbonds <- NULL
    }

    # column numbers of phi and psi are different between
    # the old and new versions of DSSP
    phi <- as.numeric(substring(raw.lines, 104, 109))
    psi <- as.numeric(substring(raw.lines, 110, 115))
    acc <- as.numeric(substring(raw.lines, 35, 38))

    h.res <- bounds(res.num[which(sse == "H")], pre.sort=FALSE)
    g.res <- bounds(res.num[which(sse == "G")], pre.sort=FALSE)
    e.res <- bounds(res.num[which(sse == "E")], pre.sort=FALSE)
    t.res <- bounds(res.num[which(sse == "T")], pre.sort=FALSE)

    h.ind <- h.res;    g.ind <- g.res
    e.ind <- e.res;    t.ind <- t.res

    if(length(h.res) > 0) {
       res.ind <- which(sse == "H")
       h.ind[, "end"] <- res.ind[cumsum(h.res[, "length"])]
       h.ind[, "start"] <- h.ind[, "end"] - h.res[, "length"] + 1
    }

    if(length(g.res) > 0) {
       res.ind <- which(sse == "G")
       g.ind[, "end"] <- res.ind[cumsum(g.res[, "length"])]
       g.ind[, "start"] <- g.ind[, "end"] - g.res[, "length"] + 1
    }

    if(length(e.res) > 0) {
       res.ind <- which(sse == "E")
       e.ind[, "end"] <- res.ind[cumsum(e.res[, "length"])]
       e.ind[, "start"] <- e.ind[, "end"] - e.res[, "length"] + 1
    }

    if(length(t.res) > 0) {
       res.ind <- which(sse == "T")
       t.ind[, "end"] <- res.ind[cumsum(t.res[, "length"])]
       t.ind[, "start"] <- t.ind[, "end"] - t.res[, "length"] + 1
    }

    if(!resno) {
       h.res <- h.ind;    g.res <- g.ind
       e.res <- e.ind;    t.res <- t.ind
    }

    sheet = list(start=NULL, end=NULL, length=NULL, chain=NULL)
    helix = list(start=NULL, end=NULL, length=NULL, chain=NULL, type=NULL)
    turn = sheet

    ## ToDo: Add "type" for turns and strands too...

    if(length(h.res)>1) {
#      if(is.null(nrow(h.res)))
#        h.s <- as.matrix(t(h.res))
      helix$start  = c(helix$start,h.res[, "start"])
      helix$end    = c(helix$end, h.res[, "end"])
      helix$length = c(helix$length, h.res[, "length"])
      helix$chain  = c(helix$chain, cha[h.ind[, "start"]])
      helix$type   = c(helix$type, sse[h.ind[, "start"]])
    }
    if(length(g.res)>1) {
#      if(is.null(nrow(g.res)))
#        g.s <- as.matrix(t(g.res))
      helix$start  = c(helix$start,g.res[, "start"])
      helix$end    = c(helix$end, g.res[, "end"])
      helix$length = c(helix$length, g.res[, "length"])
      helix$chain  = c(helix$chain, cha[g.ind[, "start"]])
      helix$type   = c(helix$type, sse[g.ind[, "start"]])
    }
    if(length(helix$start) > 0)
       helix <- lapply(helix, function(x) {names(x) <- 1:length(helix$start); return(x)})
    if(length(e.res)>1) {
#      if(is.null(nrow(e.res)))
#        e.s <- as.matrix(t(e.res))
      sheet$start  = c(sheet$start,e.res[, "start"])
      sheet$end    = c(sheet$end, e.res[, "end"])
      sheet$length = c(sheet$length, e.res[, "length"])
      sheet$chain  = c(sheet$chain, cha[e.ind[, "start"]])
    }
    if(length(sheet$start) > 0)
       sheet <- lapply(sheet, function(x) {names(x) <- 1:length(sheet$start); return(x)})
    if(length(t.res)>1) {
#      if(is.null(nrow(t.res)))
#        t.s <- as.matrix(t(t.res))
      turn$start  = c(turn$start,t.res[, "start"])
      turn$end    = c(turn$end, t.res[, "end"])
      turn$length = c(turn$length, t.res[, "length"])
      turn$chain  = c(turn$chain, cha[t.ind[, "start"]])
    }
    if(length(turn$start) > 0)
       turn <- lapply(turn, function(x) {names(x) <- 1:length(turn$start); return(x)})

    out <- list(helix = helix, sheet = sheet, hbonds = hbonds,
                turn = turn, phi = phi, psi = psi, acc = acc,
                sse = sse, call=cl)

    class(out) <- "sse"
    return(out)
}
