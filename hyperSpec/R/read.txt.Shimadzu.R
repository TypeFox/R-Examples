  ##' Reads Shimadzu GCxGC-qMS - Spectra Files (.txt) as exported by Shimadzu Chrome Solution (v. 2.72) 
  ##' Mass Spectrometer: Shimadzu GCMS-QP 2010 Ultra (www.shimadzu.com)
  ##'
  ##' @note This is a first rough import function and the functions may change without notice.
  ##' @param filename file name and path of the .txt file
  ##' @param encoding encoding of the txt file (used by \code{\link[base]{readLines}})
  ##' @param quiet suppress printing of progress
  ##' @return list of spectra tables
  ##' @author Bjoern Egert 
  ##' @export  
  read.txt.Shimadzu <- function(filename, encoding = "", quiet = TRUE)
  {    
    
    # A file consists of several sections ([Headers])
    # Each Section consists of:
    # [Header]
    #   [MS Spectrum]
    #   [MC Peak Table]
    #   [MS Similarity Search Results for Spectrum Process Table]      
    
    impLines <- readLines(con = filename, n = -1L, ok = TRUE, warn = TRUE, encoding = encoding)
    length(impLines)
    
    # total numbers of pos1 and pos2 and pos3 are equal
    pos1 <- which(impLines == "[Header]")  														                        # row positions of Headers
    pos2 <- which(impLines == "[MC Peak Table]")  											                      # row positions of peak info tables
    pos3 <- which(impLines == "[MS Similarity Search Results for Spectrum Process Table]")    # row positions of peak annotations
    pos4 <- which(impLines == "[MS Spectrum]")													                      # row positions of peak spectra
    
    headers <- length(pos1) # number of header sections  
    
    # link spectra to headers
    pos4Li = list()
    for (i in 1:(length(pos1)))
    {
      header 			<- pos1[i] 
      headerNext	<- pos1[i+1] 
      tmp	<- (pos4>header & pos4<headerNext)	
      pos4Li[[i]] <- pos4[tmp] 
    }
    
    # treat last Header section separately ...
    headerLast <- pos1[length(pos1)]
    tmp	<- pos4>headerLast
    pos4Li[[length(pos1)]] = pos4[tmp]
    
    for (i in 1:(length(pos4Li)-1))
    {
      tmp <- length(pos4Li[[i]])
      vec	<- pos4Li[[i]]
      pos4Li[[i]] <- c(vec,pos1[i+1])
    }
    #End position
    pos4Li[[headers]] <- c(pos4Li[[headers]], length(impLines))
        
    # Check
    stopifnot(impLines[1] == "[Header]") 
    stopifnot(length(pos1) == length(pos2))
    stopifnot(length(pos2) == length(pos3))
    stopifnot(length(pos3) == length(pos4Li))
    
    
    # ----------------- 1. Import: gather [Header] informations
    
    # gather in lists
    res2Li <- list()	  # Peak Info
    res3Li <- list()    # Peak Similarity
    res4Li <- list()    # Mass Spectra
    
    for(header in 1:headers)
    {
      
      if(!quiet) cat("header:", header, "\n")
      
      # ----------------- 1a. Import: "[MC Peak Table]"
      
      start <- pos2[header]+3
      stop  <- pos3[header]-2
      
      peakMat <- read.table(file = filename, skip = start-1, nrows = stop - start, 
                            header = TRUE, sep = ";", dec =".", comment.char = "", 
                            stringsAsFactors = FALSE, quote = "\"'")
      
      # maybe faster than above ...
      #Peak#;Ret.Time;Proc.From;Proc.To;Mass;Area;Height;A/H;Conc.;Mark;Name;Ret. Index"      
      #colnames <- strsplit(impLines[start], split = ";")[[1]]
      #strsplit(impLines[(start+1):stop], split = ";")[[1]]
      #peakMat <- matrix(impLines[(start+1):stop], ncol = length(colnames), byrow = TRUE)      
      
      res2Li[[header]] <- peakMat
      
      # ----------------- 1b. Import: "[MS Similarity Search Results for Spectrum Process Table]"
      
      start <- pos3[header]+2
      stop  <- pos4Li[[header]][1]-1
      if(stop-start!=0)	# no annotation hits
      {        
        simMat <- read.table(file = filename, skip = start - 1, , nrows = stop - start, 
                            header = TRUE, sep = ";", dec =".", comment.char = "", 
                            stringsAsFactors = FALSE, quote = "")  # quote = "\"'"
      }else simMat <- NA
      res3Li[[header]] <- simMat
      
      # ----------------- 1c. Import: "[MS Spectrum]"
      
      specLi <- list()	# list of all spectra in current header section 
      for(i in 1:(length(pos4Li[[header]])-1))	
      {  		

        # extract spectra
        start <- pos4Li[[header]][i]+5  		# data starts 5 rows below
        stop  <- pos4Li[[header]][i+1]-1	  # last data row before new Spectrum 
        
        # debug purposes
        if(! quiet) cat("header:", header, "spec:", i, "start:", start, "stop:", stop, "\n")
        
        # -------------------
        # Catch expeption, when peak is reported in peakMat, but [MS Spectrum] is not available: like example: (keyrow: # of Peaks; 0)
        # [MS Spectrum]
        # # of Peaks;0
        # Raw Spectrum;38.343 (scan : 68686);Base Peak;m/z 0.00 (Inten : 0)
        # Background;38.342 <-> 38.348 (scan : 68684 <-> 68697)
        # m/z;Absolute Intensity;Relative Intensity
        # [MS Spectrum]
        
        checkRow <- pos4Li[[header]][i]+1  # check for: "# of Peaks;0"
        isEmptySpec <- impLines[checkRow] == "# of Peaks;0"
        if(isEmptySpec)
        {
          # special case: create dummy spectrum with zero intensities
          spec <- cbind("m/z" = 1:100, "Absolute Intensity" = 0, "Relative Intensity" = 0)
          
        }
        # -------------------
        
        if(!isEmptySpec)
        {
          spec  <- scan(file = filename, sep = ";", skip = start - 1, nlines = (stop - start) + 1, 
                        dec = ".", quiet = TRUE)
          spec	<- matrix(spec,ncol = 3, byrow = T)
          colnames(spec) <- c("m/z", "Absolute Intensity", "Relative Intensity")          
        }
        
        specLi[[i]] <- spec	
        
      }
      
      res4Li[[header]] <- specLi     
      
    } # for(headers)
    
    
    # ----------------- 2. combine all headers sections
    
    # res2Li --> m2
    m2 <- as.data.frame(res2Li[[1]])
    m2 <- cbind(header = 1, m2)
    for(header in 2:headers) 
    {
      tmp <- as.data.frame(res2Li[[header]])
      tmp <- cbind(header=header,tmp)
      m2  <- rbind(m2,tmp)
    }
    
    # res3Li --> m3    
    # In a header section there may be not annotation tables
    m3 <- do.call("rbind", res3Li) 
    # add header Nr. 
    tmpMat <- lapply(X = res3Li, FUN = nrow)  
    tmpMat <- as.matrix(tmpMat)
    tmp <- vector(length = 0)
    for(i in 1:nrow(tmpMat))
    {
      if(tmpMat[i,1] == "NULL" ) tmp <- c(tmp,i)
      if(tmpMat[i,1] != "NULL" ) tmp <- c(tmp,rep(x = i, times = tmpMat[i,1]))
    }    
    m3  <- cbind(header = tmp, m3) 
    m3  <- m3 [, c("header", "Spectrum.", "Hit..", "SI", "CAS..", "Name", "Mol.Weight", "Mol.Form", 
                   "Retention.Index")]  # select most important columns
    tmp <- complete.cases(m3)
    m3  <- m3[tmp,]
    
    # res4Li --> m4
    tmp <- colnames(res4Li[[1]][[1]])
    m4  <- matrix(NA, nrow = 1, ncol = length(tmp))
    colnames(m4) <- tmp
    m4  <- cbind(header = NA, spectra = NA, m4)         # header number and spectra number as first columns
    for(header in 1:headers)
    {  
      for(spectra in 1:length(res4Li[[header]]))         		    
      {
        tmp <- as.matrix(res4Li[[header]][[spectra]])	
        tmp <- cbind(header, spectra, tmp)
        m4	<- rbind(m4, tmp)
      }# spectras
    }# header
    mode(m4) <- "numeric"    
    m4 <- m4[-1,]
    
    return(list(peakInfo = m2, peakAnnotate = m3, peakMasses = m4))  
    
  } 
  
