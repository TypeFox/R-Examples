summary.cnapath <- function(object, ..., pdb = NULL, label = NULL, col = NULL, 
                            plot = FALSE, concise = FALSE, cutoff = 0.1, normalize = TRUE) {
   
   pa <- list(object, ...)
   if(!all(sapply(pa, inherits, "cnapath")))
      stop("Input pa is not a 'cnapath' object")
   
   if(is.null(label)) label = 1:length(pa)
   if(is.null(col)) col = 1:length(pa)

   out <- list()
   
   # read node numbers on paths
   y <- lapply(pa, function(x) unlist(x$path))
   
   # store node degeneracy 
   node.deg <- lapply(y, table)
   if(normalize) {
      node.deg <- lapply(node.deg, function(x) x/max(x))
   }

   # find on-path node by the cutoff
   yy <- lapply(node.deg, function(x) x[x >= cutoff])
   onpath.node <- unique(names(unlist(yy)))
   i <- as.numeric(onpath.node)
   onpath.node <- onpath.node[order(i)]

   # generate the node degeneracy table
   o <- lapply(node.deg, function(x) {
      x <- x[match(onpath.node, names(x))]
      x[is.na(x)] <- 0
      names(x) <- onpath.node 
      x 
   } )

   # replace node id with pdb resid and resno
   if(!is.null(pdb)) {
      ca.inds <- atom.select(pdb, elety="CA", verbose = FALSE)
      resno <- pdb$atom[ca.inds$atom, "resno"]
      resid <- pdb$atom[ca.inds$atom, "resid"]
      chain <- pdb$atom[ca.inds$atom, "chain"]
      lig.inds <- atom.select(pdb, "ligand", verbose = FALSE)
      islig <- paste(chain, resno, sep="_") %in% 
           paste(pdb$atom[lig.inds$atom, "chain"], 
                 pdb$atom[lig.inds$atom, "resno"], sep="_")
      resid[!islig] <- aa321(resid[!islig])

      o <- lapply(o, function(x) {
         node <- as.numeric(names(x))

         if(length(unique(pdb$atom[, "chain"])) > 1)
            n <- paste(chain[node], paste(resid[node], resno[node], sep=""), sep="_")
         else
            n <- paste(resid[node], resno[node], sep="")

         names(x) <- n
         x 
      } )
   }

   names(o) <- label
   out$network <- label
   out$num.paths <- sapply(pa, function(x) length(x$path))
   out$hist <- lapply(pa, function(x) table(cut(x$dist, breaks=5, include.lowest = TRUE)))
   if(length(out$hist)==1) out$hist = out$hist[[1]]
   out$degeneracy <- do.call(rbind, o)
   if(normalize) out$degeneracy <- round(out$degeneracy, digits=2)
   
   if(plot) {
     
      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))
 
      layout(matrix(1:2, nrow=1), respect = TRUE)

      rgbcolors <- sapply(col, col2rgb) / 255
      rgbcolors <- rbind(rgbcolors, alpha = 0.6)
      
      ##- for path length distribution
      y1 <- lapply(pa, function(x) 
         hist(x$dist, breaks = 20, plot = FALSE) )

      par(mar=c(4, 4, 1, 1)) 
      plot(y1[[1]], freq = FALSE, col = do.call(rgb, as.list(rgbcolors[,1])), 
        border = col[1], main = "Path Length Distribution", 
        xlim = range(unlist(lapply(y1, "[[", "breaks"))),
        ylim = c(0, max(unlist(lapply(y1, "[[", "density")))), 
        xlab = "Path length", ylab = "Probability density")

      if(length(y1) > 1) 
         for(i in 2:length(y1))  {
            plot(y1[[i]], freq = FALSE, col = do.call(rgb, as.list(rgbcolors[,i])), 
            border = col[i], add = TRUE)
      }  
      legend("topleft", legend = label, bty = "n", text.col = col)

      ##- for node degeneracy 
      y2 <- lapply(pa, function(x) unlist(x$path))
      if(!is.null(pdb)) y2 <- lapply(y2, function(x) resno[x])
      if(concise) { 
         # re-number node to get more concise plot
         ii <- sort(unique(unlist(y2)))
         y2 <- lapply(y2, match, ii)
      }
      y2 <- lapply(y2, function(x)
         hist(x, breaks = c(seq(min(x), max(x), 1) - 0.5, max(x) + 0.5),
            plot = FALSE) )

      par(mar=c(4, 4, 1, 1))
      plot(y2[[1]], freq = TRUE, col = do.call(rgb, as.list(rgbcolors[,1])), 
        lty = 0, main = "Node Degeneracy", 
        xlim = range(unlist(lapply(y2, "[[", "breaks"))),
        ylim = c(0, max(unlist(lapply(y2, "[[", "counts")))), 
        xlab = "Node no", ylab = "Number of paths")

      if(length(y2) > 1) 
         for(i in 2:length(y2)) 
            plot(y2[[i]], freq = TRUE, col = do.call(rgb, as.list(rgbcolors[,i])), 
            lty = 0, add = TRUE)
        
   } 

   return(out)
}

print.cnapath <- function(x, ...) {
   dots = list(...)

   if(is.list(x) && all(sapply(x, inherits, "cnapath"))) {
      if(!"label" %in% names(dots) || is.null(dots$label)) dots$label = names(x)
      names(x) <- NULL
      args = c(x, dots)
      o <- do.call(summary, args)
   } else {
      o <- summary(x, ...)
   }

   if("plot" %in% names(dots)) plot = dots$plot
   else plot = FALSE

   if(!plot) {
      if("normalize" %in% names(dots)) normalize = dots$normalize
      else normalize = TRUE
   
      if(length(o$network) > 1) {   
         cat("Number of networks: ", length(o$network), "(", 
             paste(o$network, collapse=", "), ")\n")
      }
   
      cat("Number of paths in network(s):\n")
      if(length(o$network) > 1) {
          cat(paste("   ", o$network, ": ", o$num.paths, sep="", collapse="\n"), sep="\n")
          cat("\n")
      } else {
          cat("    ", o$num.paths, "\n\n")
      }
   
      cat("Path length distribution: \n")
      if(length(o$network) > 1) {   
         for(i in 1:length(o$network)) {
             cat("   --- ", o$network[i], " ---")
             print(o$hist[[i]])
             cat("\n")
         }
      } else {
         print(o$hist)
         cat("\n")
      }
   
      cat("Node degeneracy table: \n\n")
      if(length(o$network) == 1) rownames(o$degeneracy) = ""
      if(normalize)
         print(format(o$degeneracy, nsmall=2), quote=FALSE)
      else 
         print(o$degeneracy)
   }
}
