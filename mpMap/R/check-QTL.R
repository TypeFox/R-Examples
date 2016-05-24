check_qtl <- 
function(QTL, map, n.founders)
{
  # should be formatted as:
  #qtl <- as.data.frame(matrix(data=c(
#			1, 20, .5, 0, 0, 0,
#			2, 20, 0, .5, 0, .5,
#			3, 20, 0, .5, .5, .5), nrow=3, ncol=6, byrow=T))

#names(qtl) <- c("chrm","pos", "fdr1", "fdr2", "fdr3", "fdr4")

  if (!((!is.matrix(QTL) && length(QTL) == (2+n.founders)) || 
	(is.matrix(QTL) && ncol(QTL) == (2+n.founders)))) 
          stop(paste("Model must be a matrix with ", 2+n.founders, " columns (chr, pos and effects).", sep=""))

  if (!is.matrix(QTL)) 
     QTL <- rbind(QTL)
  n.qtl <- nrow(QTL)
  
  QTL[, 2] <- QTL[, 2] + 1e-14

  map2 <- map
#  if (is.null(names(map))) names(map2) <- paste("Chr", c(1:length(map)), sep="")
  for (i in 1:length(map2))
	if (is.null(names(map2[[i]]))) names(map2[[i]]) <- 
		paste("C",i,"M", c(1:length(map2[[i]])), sep="")

  for (i in 1:n.qtl)
  {
    qch <- QTL[i,1]
    qps <- QTL[i,2]
    names(qps) <- paste("QTL", i, sep="")
    if (map2[[qch]][1]>qps) 
	map2[[qch]] <- c(qps, map2[[qch]]+qps)
    if (map2[[qch]][1]<=qps & map2[[qch]][length(map2[[qch]])] >=qps)
    {
	m1 <- which.max(map2[[qch]][map2[[qch]]<=qps])
	m2 <- which(map2[[qch]]==min(map2[[qch]][map2[[qch]]>=qps]))
	map2[[qch]] <- c(map2[[qch]][1:m1], qps, 
		map2[[qch]][m2:length(map2[[qch]])])
	tmp <- c(map2[[qch]][m1], qps)
	names(tmp) <- NULL
	if (isTRUE(all.equal(tmp[1], tmp[2]))) map2[[qch]] <- map2[[qch]][-m1]
    }
    if (map2[[qch]][length(map2[[qch]])] < qps) {
	map2[[qch]] <- c(map2[[qch]], qps)
	tmp <- map2[[qch]][(length(map2[[qch]])-1):length(map2[[qch]])]
	names(tmp) <- NULL
	if (isTRUE(all.equal(tmp[1], tmp[2]))) map2[[qch]] <- map2[[qch]][-(length(map2[[qch]])-1)]
    }
  }

  return(map2)
}
