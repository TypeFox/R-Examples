`seqaln.pair` <-
function(aln, ...) {
  cl <- match.call()
  dots <- list(...)
  dots$extra.args = paste("-matrix",
                system.file("matrices/custom.mat", package="bio3d"),
                "-gapopen -3.0 ",
                "-gapextend -0.5",
                "-center 0.0", dots$extra.args)
  args <- c(list(aln=aln), dots) 
  l <- do.call(seqaln, args)

  if(!all((seqidentity(l))==1)) {
    warning("Sequences are not identical, use seqaln()")
  }
  l$call <- cl
  return(l)
}

