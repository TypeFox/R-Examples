# ==================================================================
# read multialignments, allowing for different flavors of formatting
# ==================================================================

read.align <- function(file, flavor) {

  # ================================================
  # help function to make columnnames for alignments
  # ================================================
  
  makenames <- function(alignment, word) {
    
    freq <- function(col) {
      f <- sort(table(col), decreasing=T)
      if (length(f) > 1) {
        n <- paste(names(f[1:2]), collapse = "/")
      } else {
        n <- names(f)[1]
      }
      return(n)
    }
    
    names <- paste0(apply(alignment, 2, freq)
                    , " ("
                    , word
                    , ":"
                    , 1:ncol(alignment)
                    , ")"
    )
    
    return(names)
  }
  	
  # ==============
  # ==== BDPA ====
  # ==============
  
	if (flavor == "BDPA") {
		
		# METADATA
		
		meta <- scan(file
					, what = "character"
					, sep = "\n"
					, nlines = 2
					, quiet = TRUE
					)
		# extract words
		word <- sub("^.*[\"\\*](.+)[\"\\*].*$", "\\1", meta[2], perl = TRUE)
		meta <- c(meta, word)
		names(meta) <- c("dataset", "wordinfo", "word")
		meta <- as.list(meta)
		
		# ALIGNMENTS
		
		align <- read.table(file
						, sep = "\t"
						, quote = ""
						, skip = 2
						, stringsAsFactors = FALSE
						)
		
		# NAMES OF LANGUAGES
		
		doculects <- align[,1]
		doculects <- gsub("\\.*$", "", doculects)
		align <- align[, -1]
		
		# ANNOTATION OF COLUMNS
		
		selection <- grep("^[A-Z]+$", doculects)
		annotation <- align[selection, , drop = FALSE]
		if (nrow(annotation) > 0) {
			rownames(annotation) <- doculects[selection]
			align <- align[-selection, ]
			doculects <- doculects[-selection]
		} else {
			annotation = NULL
		}
	
		# MAKE COLUMNNAMES FOR ALIGN
		
		colnames(align) <- makenames(align, word)
		
		# RETURN LIST	
					
		return(list(
				meta = meta, 
				align = align, 
				doculects = doculects, 
				annotation = annotation
				))
	}
		
  # =============
  # ==== PAD ====
  # =============
  
  if (flavor == "PAD") {
    
    # METADATA
    
    meta <- scan(file
                 , what = "character"
                 , sep = "\n"
                 , quiet = TRUE
    )
    meta <- grep("^# .+:", meta, value = TRUE)
    meta <- gsub("^# (.+)", "\\1", meta)
    meta <- sapply(meta, strsplit, split = " +: +")
    names(meta) <- NULL
    meta <- unlist(meta)
    dim(meta) <- c(2,length(meta)/2)
    colnames(meta) <- meta[1,]
    meta <- as.list(meta[2,])
    
    # extract words
    word <- gsub("_\\d+\\.msa$","",file)
    word <- gsub("^.*/([^/]+)$", "\\1", word)
    meta <- c(list(word = word), meta)

    # ALIGNMENTS
    
    align <- read.table(file
                    , header = FALSE
                    , quote = ""
                    , comment.char = "#"
                    , sep = "\t"
                    )
    
    # remove dots
    
    align[,2] <- gsub(".", "", align[,2], fixed = TRUE)
 
    # ANNOTATION
    
    getAnno <- align[,1] == ":ANN"
    annotation <- align[getAnno, , drop = FALSE]
    rownames(annotation) <- annotation[,2]
    annotation <- annotation[, -c(1,2)]
    
    # DOCULECTS
    
    align <- align[!getAnno,]
    doculects <- align[,2]
    rownames(align) <- align[,1]
    align <- align[, -c(1,2)]
    
    # COLNAMES
    
    colnames(align) <- makenames(align, word)
    
    # RETURN LIST	
    
    return(list(
      meta = meta, 
      align = align, 
      doculects = doculects, 
      annotation = annotation
    ))
  }
}


