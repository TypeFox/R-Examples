read.starbeast <- function(file) 
{
    # X: the scanned BEAST output
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    # phy: phylogeny
    phy <- read.nexus(file)
  
    # isolate NEWICK string
    # ------------------------------
    X <- X[grep("tree TREE1[[:space:]]+=", X)]
    X <- gsub("^.*tree TREE1[[:space:]]+= \\[&R\\] ", "", X)
    X <- gsub("[!]rotate=true,*|[!]rotate=false,*", "", X)
    X <- gsub(";$", "", X)
    
    # store stats per node in a list
    # ------------------------------
    foo <- function(x){
    	x <- gsub(",([[:lower:]])", "xxx\\1", x)
    	x <- unlist(strsplit(x, "xxx"))
    	names(x) <- sapply(x, function(x){unlist(strsplit(x, 
    		split = "="))[1]})
    	x <- sapply(x, function(x){unlist(strsplit(x, 
    		split = "="))[2]})
    	x <- gsub("[{}]", "", x)
    	x <- strsplit(x, ",")
    	lapply(x, as.numeric)
    }
    vals <- unlist(strsplit(X, "\\][[:alnum:]:.)]*\\[*&*"))
    vals <- gsub("^.*&", "", vals)
    vals <- lapply(vals, foo)
      
    #
    edges <- gsub("\\[(&[[:alnum:]_=%!.,{}-]+)\\]", "", X)
    edges <- unlist(strsplit(gsub("[()]*", "", edges), ":"))
    tips <- gsub("^[0-9]+.[0-9]+,*", "", edges)
    edges <- gsub(",$|,[[:alpha:]].+$", "", edges)
    tab <- cbind(edges, nodes = c("root", head(tips, -1)), 
        id = rep(NA, length(edges)))
    tab[tab[, 2] == "", 2] <- "internal"
    internal <- grep("internal|root", tab[, 2])
    tab[-internal, 3] <- seq(along = phy$tip.label)
    
    intnodes <- match(tab[internal, 1], 
        phy$edge.length[phy$edge[, 2] > 50], 
   	    nomatch = 0) + 51
   	tab[internal, 3] <- intnodes
   	
   	names(vals)[internal] <- paste("NODE", tab[internal, 3], sep = ":")
   	names(vals)[-internal] <- paste("TIP", tab[-internal, 3], sep = ":")
   	
   	#statsnames
   	sn <- lapply(vals, names)
   	sn <- sort(unique(unlist(sn)))
    # obj: a list for node values
   	obj <- as.list(sn)
   	names(obj) <- sn
    # loop over statsnames (sn) to create list element
    # for each node statistic
   	for (i in seq_along(sn)){
   		fo2 <- function(x, param){
   			unlist(x[param == names(x)])	
   		}
   		x <- sapply(vals, fo2, param = sn[i])
   		L <- max(sapply(x, length))
   		if (L == 1) {
   		    obj[[i]] <- unlist(x)
   		}
      if (L == 2) {
        obj[[i]] <- matrix(unlist(x), nrow = L, byrow = FALSE)
      }
      if (L > 2){
        obj[[i]] <- x
      }
   	}  	
   	tips <- vals[-internal]
   	nodes <- vals[internal]
    
    phy <- c(phy[1:length(phy)], obj[1:length(obj)])
    class(phy) <- ("phylo")
    attr(phy, "origin") <- file 
    return(phy)
}

# (&[[:alnum:]_=%!.,{}-]+) finds information between square brackets
