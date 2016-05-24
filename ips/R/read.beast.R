read.beast <- function(file, digits = NULL) 
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    LEFT <- grep("\\[", X)
    
   	# table of node values:
   	# ---------------------
   	tab <- extractBEASTstats(file)
   	if (!is.null(digits))
   		tab <- round(tab, digits = digits)
   	interior <- which(!is.na(tab$posterior))
 
 	# get tree as string: 'tree'
   	# ---------------------
   	RIGHT <- grep("\\]", X)
    if (length(LEFT)) {
        w <- LEFT == RIGHT
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)

	# translate:
    end <- semico[semico > i2][1]
    x <- X[(i2 + 1):end]
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    TRANS <- matrix(x, ncol = 2, byrow = TRUE)
    TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    n <- dim(TRANS)[1]
    
    start <- semico[semico > i2][1] + 1			
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    tree <- gsub("^.*= *", "", tree)
    
    # brl: vector of branch lengths
    # -----------------------------
    brl <- unlist(strsplit(tree, ":"))[-1]
	brl <- gsub("[( | ) | ;]", "", brl)
	brl <- strsplit(brl, ",")
	foo <- function(x) x <- head(x, 1)
	brl <- unlist(lapply(brl, foo))
	brl <- paste("", brl, sep = ":")
	brl <- c(brl, ";")
	
	# create list with node statistics
	# --------------------------------
	nodestats <- vector(mode = "list", length = dim(tab)[2])
	for (i in seq(along = nodestats)){
		newtree <- tree
		val <- tab[, i]
		ggg <- paste(val, brl, sep = "")
		ggg[length(ggg)] <- paste(tail(val, 1), ";", sep = "")
		for (j in interior)
			newtree <- gsub(brl[j], ggg[j], newtree)
		dt <- read.tree(text = newtree)
		# hack to suppress warning message:
		z <- dt$node.label
		z[z == "NA"] <- 9999
		z <- as.numeric(z)
		z[z == 9999] <- NA 
		nodestats[[i]] <- z
		names(nodestats)[i] <- colnames(tab)[i]	
	}
  	
  	# read tree ...
   	tr <- read.nexus(file)
   	# ... and append node statistics to phylo object
   	tr <- c(tr[1:length(tr)], nodestats[1:length(nodestats)])
   	class(tr) <- ("phylo")
   	attr(tr, "origin") <- file
   	tr   
}