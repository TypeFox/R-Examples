`vec2resno` <-
function(vec, resno) {
  ## replicate vec based on concetive
  ## similar resno entries

  if(is.pdb(resno))
    resno <- resno$atom[,"resno"]

  res.len <- rle(resno)$lengths
  if(length(vec) != length(res.len))
    stop("Length miss-match of 'vec' and concetive 'resno'")

  if( sum(res.len) != length(resno) )
    stop("Replicated length Miss-match")
  
  return( rep(vec,  times=res.len))
}

