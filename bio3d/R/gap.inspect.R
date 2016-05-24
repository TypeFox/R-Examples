"gap.inspect" <-
function(x) {
  
  # Report the number of gaps, ("-",".",NA), per
  # row (i.e. seq) and col (i.e. position) in a
  # given "alignment" 'x'

  if(is.vector(x)) {
    gaps <- is.gap(x)
    inds <- which(gaps)
    f.inds <- which(!gaps)
    gap.pos <- as.numeric(gaps)
    gap.col <- gap.pos
    gap.row <- sum(gaps)
  } else {
    if(is.list(x)) {
      if(inherits(x, "pdbs")) {
        x <- x$xyz; warning("Taking $xyz component (NOT $ali for which you should use 'gap.inspect(x$ali)')")
      } else {
        x <- x$ali
      }
    }
    
    gap.pos1 <-( as.numeric(x=="-") + as.numeric(x==".") )
    gap.pos2 <- as.numeric(is.na(gap.pos1))
    
    gap.pos<-matrix( colSums(rbind(gap.pos1,gap.pos2), na.rm=TRUE),
                    ncol=ncol(x))
  
    gap.col <- colSums(gap.pos)
    gap.row <- rowSums(gap.pos)
    inds <- which(gap.col!=0)
    ##f.inds=(1:ncol(x))[-inds]
    f.inds <- which(gap.col == 0)
  }
  output=list(t.inds=inds, f.inds=f.inds,
    row=gap.row, col=gap.col, bin=gap.pos)
}

