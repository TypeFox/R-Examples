## make a collection of cell.specs for tableplot
#
#  each argument is repeated up to the required length
#  # of patterns = n (if specified), or the maximum length of any argument

#  list of lists to data.frame and back:
#    Given pats as a list of lists: 
#      pats.df <- do.call(rbind, lapply(pats, data.frame))
#      pats.df <- do.call(rbind.data.frame, pats)    ## works best
#      pats.df <- do.call(rbind, lapply(pats, data.frame))  ## works same -- real data frame
#      pats.df <- data.frame(do.call(rbind,lapply(pats,function(x) t(as.matrix(x,ncol=10)))))
#    Given pats.df: pats <- lapply(seq(along = rownames(pats.df)),function(i) as.list(pats.df[i, ]))

make.specs0 <- function(
  n=NULL,
  as.data.frame=FALSE,
  shape=0,
  shape.col="black",
  shape.lty=1,
	shape.neg=0,
	shape.col.neg="red", 
	shape.lty.neg=1, 
  cell.fill="white",
  back.fill="white",
  label=0,
  label.size=0.7,
	label.col="black",
	ref.lines=FALSE,
  ref.col="gray80",
  scale.max=1,
	shape.lwd=1,
	frame.col="black",
	frame.lwd="0.5"
  ) {


  if(is.null(n)) {
    len <- max( unlist(lapply(list(
          shape, shape.col, shape.lty, 
          shape.neg, shape.col.neg, shape.lty.neg, 
          cell.fill, back.fill, label, label.size, label.col,
          ref.lines, ref.col,  
          scale.max,
          shape.lwd, frame.col, frame.lwd
          ), FUN=length)))
  } else {
  	len <- n
  }

  # replicate the elements in each, to required length
  shape          <- rep(shape, length.out=len)
  shape.col      <- rep(shape.col, length.out=len)
  shape.lty      <- rep(shape.lty, length.out=len)
  shape.neg      <- rep(shape.neg, length.out=len)
  shape.col.neg  <- rep(shape.col.neg, length.out=len)
  shape.lty.neg  <- rep(shape.lty.neg, length.out=len)
  cell.fill      <- rep(cell.fill, length.out=len)
  back.fill      <- rep(back.fill, length.out=len)
  label          <- rep(label, length.out=len)
  label.size     <- rep(label.size, length.out=len)
  label.col      <- rep(label.col, length.out=len)
  ref.lines      <- rep(ref.lines, length.out=len)
  ref.col        <- rep(ref.col, length.out=len)
#  ref.grid       <- rep(ref.grid, length.out=len)
  scale.max      <- rep(scale.max, length.out=len)
  shape.lwd      <- rep(shape.lwd, length.out=len)
  frame.col      <- rep(frame.col, length.out=len)
  frame.lwd      <- rep(frame.lwd, length.out=len)

#browser()
  
# easier to work with as a data.frame
	if(as.data.frame) {
  result <- data.frame(cbind(
          shape, shape.col, shape.lty, 
          shape.neg, shape.col.neg, shape.lty.neg, 
          cell.fill, back.fill, label, label.size,  label.col,
          ref.lines, ref.col,  
          scale.max,
          shape.lwd, frame.col, frame.lwd
        )
        , stringsAsFactors=FALSE)
  result <- result[1:len, ]
  }
  else { 
    result <- as.list(vector(length=len))
    for(i in 1:len) {
      item <- list(
      	shape=shape[i], shape.col=shape.col[i],  shape.lty=shape.lty[i], 
        shape.neg=shape.neg[i], shape.col.neg=shape.col.neg[i], shape.lty.neg=shape.lty.neg[i], 
        cell.fill=cell.fill[i], back.fill=back.fill[i], label=label[i],
        label.size=label.size[i], label.col=label.col[i],
        ref.lines=ref.lines[i], ref.col=ref.col[i],  
        scale.max=scale.max[i],
        shape.lwd=shape.lwd[i], frame.col=frame.col[i], frame.lwd=frame.lwd[i]
        )
      result[[i]] <- item    
      }
	class(result) <- c("specs", "list")
    }
  result  
}

### print patterns in a human readable way
#print.specs <- function(pats) {
#	result <- pats
#	if (is.list(pats))
#		result <- do.call(rbind.data.frame, pats)
#	result
#}
