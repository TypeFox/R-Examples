`draw.subgraph` <-
function(mesa, handle) {
# draws a subgraph containing
# npos positive
   nitem <- ncol(mesa) - 1
   npatt <- nrow(mesa)
   edges <- matrix("n", nrow = npatt, ncol = npatt)
   rownames(edges) <- rownames(mesa)
# set up edges
# edges is set to
#   n if no edge existed
#   a if row is an ancestor of column
#   p if row is descendant of column and is pending output
#   y if row is descendant of column and has been output
   if (npatt > 1) {
   for (i in 1:(npatt - 1)) {
      for (j in (i+1):npatt) {
         differs <- which(mesa[i,1:nitem] != mesa[j,1:nitem])
         if(length(differs) == 2) {
            if ((differs[2] - differs[1]) == 1) {
               edges[i,j] <- "p"
               edges[j,i] <- "a"
            }
         }
      }
   }
   } # end if
# set up orphan
# an orphan has neither ancestors nor descendants
   orphan <- logical(npatt)
   orphan <- FALSE
   if (npatt > 1) orphan <- apply(edges, 1, function(x) sum(x != "n") == 0)
# follow routes from each node in turn
# output the links between nodes
   for (i in 1:npatt) {
      if(any(edges[i,] == "p")) {
         path <- paste("node", rownames(mesa)[i], sep = "", collapse = "")
         j <- i
         while(any(edges[j,] == "p")) {
            k <- which(edges[j,] == "p")[1]
            path <- paste(path, "->", "node", rownames(mesa)[k],
               sep = "", collapse = "")
            edges[j, k] <- "y"
            j <- k
         }
         cat(file = handle, path, "\n") # output that route
      }
   }
# now output nodes with frequencies
   for (i in 1:npatt) {   
      node <- paste("node", rownames(mesa)[i], ' [label = "',
         rownames(mesa)[i], "\\n", mesa[i, nitem + 1], '"]',
         sep = "", collapse = "")
      cat(file = handle, node, "\n")
   }
   if(any(orphan)) {
      which.orphan <- which(orphan)
      if(length(which.orphan) > 1) { 
#       need to link orphans invisibly 
#       if more than one
#       to keep them in one column and so save space
         path <- paste("node", 
            rownames(mesa)[which.orphan[1]], sep = "", collapse = "")
         for (i in 2:length(which.orphan)) {
            path <- paste(path, "->", "node", rownames(mesa)[which.orphan[i]],
               sep = "", collapse = "")
         }
         cat(file = handle, "edge [style = invis]\n", path, "\n")
      } # end if more than one orphan
   } # end if orphan
} # end of draw.subgraph

