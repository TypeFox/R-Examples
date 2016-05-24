### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### textmatrix
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### dependencies: library("RStem")

### HISTORY
###
### 2014-03-21
###    * fixed snowball dependancy to snowballC
###
### 2012-07-23
###    * added html entity support and special chars support
###      for polish
###
### 2012-02-28
###    * added support for vietnamese transliterations
### 
### 2009-08-19
###    * exchanged Rstem with Snowball 
###
### 2007-11-28
###    * bugfix textvector: stemming before (!) 
###      filtering with controlled vocabulary
###    * bugfix textvector: special chars better now
###
### 2008-08-31
###    * iconv routines to solve character encoding problems
###
### 2009-01-30
###    * removed german umlauts from non-alnum summary
### 

textvector <- function (file, stemming=FALSE, language="english", minWordLength=2, maxWordLength=FALSE, minDocFreq=1, maxDocFreq=FALSE, stopwords=NULL, vocabulary=NULL, phrases=NULL, removeXML=FALSE, removeNumbers=FALSE ) {
   
   #txt = scan(file, what = "character", quiet = TRUE, encoding="UTF-8", fileEncoding="UTF-8")
   #txt = iconv(txt, to="UTF-8")
   txt = readLines(file, warn = FALSE, encoding = "UTF-8")
    
   res = try(tolower(txt), TRUE)
	if (class(res) == "try-error") {
	   stop(paste("[lsa] - could not open file ",file," due to encoding problems of the file.", sep=""))
	} else {
	   txt = res
	   res = NULL
	   gc()
	}
   
   ## Current version of R have the following bug:
   ##     R> txt <- "ae "
   ##     R> Encoding(txt)
   ##     [1] "UTF-8"
   ##     R> Encoding(gsub( "[^[:alnum:]]", " ", txt))
   ##     [1] "unknown"
   ## (The space matters.)
   ## Hence, let's save the encoding and tag it back on ...
   
	#encoding = Encoding(txt)
	
   if (removeXML) {
      txt = gsub("<[^>]*>"," ", paste(txt,collapse=" "), perl=TRUE)
      txt = gsub("<[^>]*>"," ", paste(txt,collapse=" "), perl=TRUE)
      txt = gsub("&gt;",">", txt, perl=FALSE, fixed=TRUE)
      txt = gsub("&lt;","<", txt, perl=FALSE, fixed=TRUE)
      txt = gsub("&quot;","\"", txt, perl=FALSE, fixed=TRUE)
      if (language=="german" || language=="polish") { # && l10n_info()$MBCS
         data(specialchars, envir = environment())
         specialchars = get("specialchars", envir  = environment())
         for (sc in 1:length(specialchars$entities)) {
	         txt = gsub(specialchars$entities[sc],specialchars$replacement[sc], txt, perl=FALSE, fixed=TRUE)
         }
      }
   }
   
   if (language=="arabic") {
      ## support for Buckwalter transliterations
      txt = gsub( "[^[:alnum:]'\\_\\~$\\|><&{}*`\\-]", " ", txt)
   } else if (language=="vietnamese") {
      ## support for transliterations with _ to connect words, i.e.
      ## replace everything that is non alphanumeric and not _
      txt = gsub( "[^[:alnum:]\\_]", " ", txt)
   } else if (!is.null(phrases)) {
		
		# collapse the list into a single character string
		txt = paste(txt, collapse=" ")
		
      # identify phrases in the text
      for (p in phrases) {
         # convert phrase to "word1_word2_word3"
         repl = gsub("[[:space:]]+", " ", as.character(p))
         repl = gsub("[[:space:]]+", "_", repl)
		   # replace the phrase in txt (slow but works)
         txt = gsub(p, repl, txt)
      }
		
      # filter the rest
      data(alnumx, envir = environment())
      alnumx = get("alnumx", envir  = environment())
      txt = gsub(alnumx, " ", txt)
		
      # split again by whitespaces
      txt = unlist(strsplit(txt, " "))
   } else {
      data(alnumx, envir = environment())
      alnumx = get("alnumx", envir  = environment())
      txt = gsub( alnumx, " ", txt)
   }
   
   txt = gsub("[[:space:]]+", " ", txt)
   
   # Encoding(txt) <- encoding
   
   txt = unlist(strsplit(txt, " ", fixed=TRUE))
   
   # stopword filtering?
   if (!is.null(stopwords)) txt = txt[!txt %in% stopwords]
   
   # tabulate
   tab = c(sort(table(txt), decreasing = TRUE))
   
   # stemming?
   #if (stemming) names(tab) = wordStem(names(tab), language=language)
   if (stemming) names(tab) = wordStem(names(tab), language)
	
   # vocabulary filtering?
   #if (!is.null(vocabulary)) txt = txt[txt %in% vocabulary]
	if (!is.null(vocabulary)) tab = tab[names(tab) %in% vocabulary]
   
   # bandwith for document frequency?
   tab = tab[tab >= minDocFreq]
	if (is.numeric(maxDocFreq)) tab = tab[tab <= maxDocFreq]
   
   # word-length filtering?
   tab = tab[nchar(names(tab), type="chars") >= minWordLength]
   if (is.numeric(maxWordLength)) tab = tab[nchar(names(tab), type="chars") <= maxWordLength]
   
   if (removeNumbers) {
      tab = tab[-grep("(^[0-9]+$)", names(tab), perl=TRUE)]
   }
	
   if (length(names(tab))==0) warning(paste("[textvector] - the file ", file, " contains no terms after filtering.", sep=""))
	
   return( data.frame( docs=basename(file), terms = names(tab), Freq = tab, row.names = NULL) )
   
}


textmatrix <- function( mydir, stemming=FALSE, language="english", minWordLength=2, maxWordLength=FALSE, minDocFreq=1, maxDocFreq=FALSE, minGlobFreq=FALSE, maxGlobFreq=FALSE, stopwords=NULL, vocabulary=NULL, phrases=NULL, removeXML=FALSE, removeNumbers=FALSE) {
    
	# if directory, then list its files recursively, else check
	# whether file exists and eventually append it to list.
	
	myfiles = NULL
	if ( length(mydir) > 1 ) {
	
		for (i in 1:length(mydir)) {
		
			if (file.info(normalizePath(mydir[i]))$isdir==TRUE) {
				myfiles = append(myfiles, dir(mydir[i], full.names=TRUE, recursive=TRUE))
			} else if (file.exists(normalizePath(mydir[i]))) {				
				myfiles = append(myfiles, normalizePath(mydir[i]))
			} else {
				warning( paste("[textmatrix] - WARNING: file ",mydir[i], " does not exist.", sep=""))
			}
			
		}
		
	} else if ( file.info(normalizePath(mydir))$isdir==TRUE ) {
		myfiles = dir(mydir, full.names=TRUE, recursive=TRUE)
	} else if ( file.exists(normalizePath(mydir)) ==TRUE ) {
		myfiles = normalizePath(mydir)
	} else {
		stop("[textmatrix] - ERROR: specified input file or directory does not exist.")
	}
	
    dummy = lapply( myfiles, textvector, stemming, language, minWordLength, maxWordLength, minDocFreq, maxDocFreq, stopwords, vocabulary, phrases, removeXML, removeNumbers)
    
	if (!is.null(vocabulary)) {
        
		dtm = t(xtabs(Freq ~ ., data = do.call("rbind", dummy)))
        
		if (is.numeric(minGlobFreq)){
			dtm = dtm[rowSums(lw_bintf(dtm))>=minGlobFreq,]
		}
		if (is.numeric(maxGlobFreq)) {
			dtm = dtm[rowSums(lw_bintf(dtm))<=maxGlobFreq,]
		}
		gc()
		
		result = matrix(0, nrow=length(vocabulary), ncol=ncol(dtm))
        rownames(result) = vocabulary
        result[rownames(dtm),] = dtm[rownames(dtm),]
        colnames(result) = colnames(dtm)
        dtm = result
        gc()
		
    } else {
        
		dtm = t(xtabs(Freq ~ ., data = do.call("rbind", dummy)))
		
		if (is.numeric(minGlobFreq)){
			dtm = dtm[rowSums(lw_bintf(dtm))>=minGlobFreq,]
		}
		if (is.numeric(maxGlobFreq)) {
			dtm = dtm[rowSums(lw_bintf(dtm))<=maxGlobFreq,]
		}
		gc()

    }
    
    environment(dtm) = new.env()
    class(dtm) = "textmatrix"
	
    return ( dtm )
    
}

print.textmatrix <- function ( x, bag_lines = 12, bag_cols = 10, ... ) {
    
    nc = ncol(x);
    nr = nrow(x);    
    
	# subscript out of bound bugfix 
	# by Jeff Verhulst, J&J Pharma R&D IM, 2006
    
    if ( (nc <= (3*bag_cols)) && (nr <= (3*bag_lines)) ) {
        
        y = x;
        attr(y,"class") = NULL;
        attr(y,"call") = NULL;
        environment(y) = NULL;
		ret = y
        
    } else if ( nc <= 3*bag_cols ) {
		
		redx = matrix(ncol = nc, nrow = (3*bag_lines));
		mid = round(nrow(x)/2)
			
		# top
		redx[1:bag_lines, 1:nc] = x[1:bag_lines, 1:nc]
		
		# mid
		redx[(bag_lines+1):(bag_lines*2), 1:nc] = x[mid:(mid+bag_lines-1), 1:nc]
		
		# bottom
		redx[(bag_lines*2+1):(bag_lines*3), 1:nc] = x[(nrow(x)-bag_lines+1):nrow(x), 1:nc]
				
		# dixnaxes
		rownames(redx) = c( paste(1:bag_lines,rownames(x)[1:bag_lines],sep=". "), paste(mid:(mid+bag_lines-1),rownames(x)[(mid):(mid+bag_lines-1)],sep=". "), paste((nrow(x)-bag_lines+1):nrow(x), rownames(x)[(nrow(x)-bag_lines+1):nrow(x)], sep=". "))
		colnames(redx) = paste("D", c( 1:nc ), sep="")
		docnames = paste( colnames(redx), c( colnames(x)[1:nc]), sep=" = ")	
		
		ret = NULL
		ret$matrix = round(redx,2);
		ret$legend = docnames;
		
	} else if ( nr <= 3*bag_lines ) {
		
		redx = matrix(ncol = (3*bag_cols), nrow = nr);
		midc = round(ncol(x)/2)
		
		# top = all
		redx[1:nr, 1:bag_cols] = x[1:nr, 1:bag_cols]
		redx[1:nr, (bag_cols+1):(bag_cols+bag_cols)] = x[1:nr, midc:(midc+bag_cols-1)]
		redx[1:nr, (2*bag_cols+1):(3*bag_cols)] = x[1:nr, (ncol(x)-bag_cols+1):ncol(x)]
		
		# dixnaxes
		rownames(redx) = c( paste(1:nr,rownames(x)[1:nr],sep=". "))
		colnames(redx) = paste("D", c( 1:bag_cols, midc:(midc+bag_cols-1), (ncol(x)-bag_cols+1):ncol(x) ), sep="")
		docnames = paste( colnames(redx), c( colnames(x)[1:bag_cols], colnames(x)[midc:(midc+bag_cols-1)], colnames(x)[(ncol(x)-bag_cols+1):ncol(x)] ), sep=" = ")	
		
		ret = NULL
		ret$matrix = round(redx,2);
		ret$legend = docnames;
		
	} else {
        
        redx = matrix(ncol = (3*bag_cols), nrow = (3*bag_lines));
        mid = round(nrow(x)/2)
        midc = round(ncol(x)/2)
        
        # top
        redx[1:bag_lines, 1:bag_cols] = x[1:bag_lines, 1:bag_cols]
        redx[1:bag_lines, (bag_cols+1):(bag_cols+bag_cols)] = x[1:bag_lines, midc:(midc+bag_cols-1)]
        redx[1:bag_lines, (2*bag_cols+1):(3*bag_cols)] = x[1:bag_lines, (ncol(x)-bag_cols+1):ncol(x)]
        
        # mid
        redx[(bag_lines+1):(bag_lines*2), 1:bag_cols] = x[mid:(mid+bag_lines-1), 1:bag_cols]
        redx[(bag_lines+1):(bag_lines*2), (bag_cols+1):(bag_cols+bag_cols)] = x[mid:(mid+bag_lines-1), midc:(midc+bag_cols-1)]
        redx[(bag_lines+1):(bag_lines*2), (2*bag_cols+1):(3*bag_cols)] = x[mid:(mid+bag_lines-1), (ncol(x)-bag_cols+1):ncol(x)]
        
        # bottom
        redx[(bag_lines*2+1):(bag_lines*3), 1:bag_cols] = x[(nrow(x)-bag_lines+1):nrow(x), 1:bag_cols]
        redx[(bag_lines*2+1):(bag_lines*3), (bag_cols+1):(bag_cols+bag_cols)] = x[(nrow(x)-bag_lines+1):nrow(x), midc:(midc+bag_cols-1)]
        redx[(bag_lines*2+1):(bag_lines*3), (2*bag_cols+1):(3*bag_cols)] = x[(nrow(x)-bag_lines+1):nrow(x), (ncol(x)-bag_cols+1):ncol(x)]
                
        # dixnaxes
        rownames(redx) = c( paste(1:bag_lines,rownames(x)[1:bag_lines],sep=". "), paste(mid:(mid+bag_lines-1),rownames(x)[(mid):(mid+bag_lines-1)],sep=". "), paste((nrow(x)-bag_lines+1):nrow(x), rownames(x)[(nrow(x)-bag_lines+1):nrow(x)], sep=". "))
        colnames(redx) = paste("D", c( 1:bag_cols, midc:(midc+bag_cols-1), (ncol(x)-bag_cols+1):ncol(x) ), sep="")
        docnames = paste( colnames(redx), c( colnames(x)[1:bag_cols], colnames(x)[midc:(midc+bag_cols-1)], colnames(x)[(ncol(x)-bag_cols+1):ncol(x)] ), sep=" = ")
		
		ret = NULL
		ret$matrix = round(redx,2);
		ret$legend = docnames;
		
    }
	
	print.default(ret);
	invisible(x);
	
}

summary.textmatrix <- function ( object, ... ) {
    
    s = vector(mode="numeric", length=5);
    n = vector(mode="character", length=5);
    n[1] = "vocabulary";
    s[1] = length(rownames(object));
    n[2] = "documents";
    s[2] = length(colnames(object));
    n[3] = "freqs not '0'";
    s[3] = length(which(object>0));
    n[4] = "max term length";
    s[4] = max(nchar(rownames(object),type="chars"));
    n[5] = "non-alphanumerics in terms";
    s[5] = length(which(gsub("[[:alnum:]]", "", rownames(object)) != "")); 
    names(s) = n;
    class(s) = "summary.textmatrix";
    s
    
}
