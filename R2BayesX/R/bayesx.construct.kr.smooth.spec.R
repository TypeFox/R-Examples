bayesx.construct.kr.smooth.spec <- bayesx.construct.kriging.smooth.spec <- function(object, dir, prg, data)
{
  termo <- object$term
  if(length(termo) < 2L)
    stop("kriging method needs two terms!")
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
	nrknots <- object$bs.dim
  xt <- object$xt
  if(is.null(xt$full))
    term <- paste(termo[1L], "*", termo[2L], "(kriging,nrknots=", nrknots, sep = "")
  else {
    term <- paste(termo[1L], "*", termo[2L], "(kriging,full", sep = "")    
    object$xt$full <- NULL
  }
  if(!is.null(xt$knotdata)) {
    if(missing(dir))
      dir <- tempdir()
    if(!missing(dir)) {
      ok <- TRUE
	    files <- list.files(dir)
      knot.name <- paste("krknots")
      counter <- NULL
	    while(ok) {
		    knotfile <- paste(knot.name, counter, ".raw", sep = "")
		    if(any(grepl(knotfile, files))) {
          if(is.null(counter))
            counter <- 0L
          counter <- counter + 1L
        } else ok <- FALSE
      write.table(xt$knotdata, paste(dir, "/", knotfile, sep = ""), quote = FALSE, row.names = FALSE)
		  }
      knot.name <- paste(knot.name, counter, sep = "")
      term <- paste(term, ",knotdata=", knot.name, sep = "")
      object$xt$knotdata <- NULL
      if(missing(prg)) {
        prg <- paste("bayesx", counter, ".prg", sep = "")
        cat("", file = paste(dir, "/", prg, sep = ""))	  
      }
      prgfile <- paste(dir, "/", prg, sep = "")
      cat("dataset", knot.name, "\n", append = TRUE, file = prgfile)
      cat(knot.name, ".infile using ", dir, "/", knotfile, "\n", sep = "", append = TRUE, file = prgfile)
    }
  }
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

