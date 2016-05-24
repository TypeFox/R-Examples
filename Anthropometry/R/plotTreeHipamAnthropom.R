plotTreeHipamAnthropom <- function(x,main,...){
 tree.obj <- x
 n.levels <-tree.obj$n.levels

 dmat <- tree.obj$development
 max.x.clust <- max(apply(dmat, 2, function(x) sum(!is.na(unique(x)))))

 box.tuning <- 0.35

 if(ncol(dmat) > 1){
  box.size <- max(box.tuning / n.levels, box.tuning / max.x.clust)
 }else{
  box.size <- 0.2
  }

 plot(0, 0, xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "", axes = FALSE, type = "n",...)
 centre.list <- NULL

 #Draw boxes and cluster numbers; calculate box centres (needed for drawing the arrows):
 for (i in 1:n.levels){
  which.clust <- unique(dmat[,n.levels + 1 - i])
  which.clust <- which.clust[!is.na(which.clust)]

    n.clust <- length(which.clust)
     for (j in 1:n.clust){
      centre <- c((j - .5) / (n.clust),(i-.5)/(n.levels + 1))
      centre.list <- rbind(centre.list, c(which.clust[j], centre))
      make.circle.discovery(centre, box.size)
      text(centre[1], centre[2], labels = which.clust[j], cex = 2.5 / max(n.clust,n.levels)^.5, pos = 1, 
           offset = -0.45)
     }
 }

 centre.list <- centre.list[order(centre.list[,1]),]
 #Draw root box:
 root.centre <- c(0.5,(n.levels + .5) / (n.levels + 1))
 make.circle.discovery(root.centre, box.size)
 text(root.centre[1], root.centre[2], labels = c("R"), cex = 3 / max(n.clust, n.levels)^.5, pos = 1, 
      offset = -0.45)
 #Draw the arrows:
 if(ncol(dmat) > 1){
 
  for (i in 1:nrow(dmat)){
   for (j in 2:ncol(dmat)){
    if (!is.na(dmat[i,j])){
     make.arrow.circle(centre.list[dmat[i,j-1],2:3],centre.list[dmat[i,j],2:3],box.size) 
    }
   }
  }
 }

 #Draw root arrows:
 root.shoots <- unique(dmat[,1])
 for (i in root.shoots){
  make.arrow.circle(root.centre,centre.list[i,2:3],box.size)
 }
 title(main=main)
    
}