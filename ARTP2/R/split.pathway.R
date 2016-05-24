
split.pathway <- function(pathway, allele.info, group.gap){
  
  if(is.null(group.gap) || group.gap == 0 || all(is.na(allele.info$Pos))){
    return(pathway)
  }
  
  pathway <- merge(pathway, allele.info[, c('SNP', 'Pos')], by = 'SNP')
  
  if(any(is.na(pathway$Pos)) || any(pathway$Pos <= 0)){
    msg <- 'Missing or invalid SNP positions detected when trying split the pathway into independent group. Please check the fourth column in bim files. '
    stop(msg)
  }
  
  pathway <- pathway[order(pathway$Chr, pathway$Pos, pathway$Gene), ]
  chr <- sort(unique(pathway$Chr))
  Group <- NULL
  ngrp <- 0
  for(c in chr){
    path <- pathway[pathway$Chr == c, ]
    d <- diff(path$Pos)
    cut <- which(d > group.gap)
    cut <- c(cut, nrow(path))
    ncut <- length(cut)
    grp <- rep(ncut, nrow(path))
    for(j in ncut:1){
      grp[1:cut[j]] <- j
    }
    
    Group <- c(Group, grp + ngrp)
    ngrp <- ngrp + ncut
    
  }
  
  pathway <- data.frame(pathway[, c('SNP', 'Gene')], Chr = Group, stringsAsFactors = FALSE)
  rownames(pathway) <- NULL
  
  pathway
  
}
