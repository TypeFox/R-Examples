### This file contains functions for "ms".

ms <- function(nsam = NULL, nreps = 1, opts = NULL, temp.file = NULL,
               tbs.matrix = NULL){
  if(! is.null(opts) && ! is.null(nsam)){
    if(length(grep("<|>|\\|", opts))){
      stop("stdin, stdout, and pipe are all disable within opts.")
    }

    if(is.null(temp.file)){
      temp.file.ms <- tempfile("ms.")
    } else{
      temp.file.ms <- temp.file
    }

    if(nsam >= 2){
      nsam <- as.character(nsam)
      nreps <- as.character(nreps)
      new.opts <- unlist(strsplit(opts, " "))

      ### Check if "tbs" option is used in opts.
      id.tbs <- new.opts == "tbs"
      n.tbs <- sum(id.tbs)
      if(n.tbs > 0){
        if(is.null(tbs.matrix) || !is.matrix(tbs.matrix) ||
           ncol(tbs.matrix) != n.tbs){
          stop("# of tbs columns is not equal to # of tbs given in opts.")
        }

        for(i.tbs in 1:nrow(tbs.matrix)){
          tmp.opts <- new.opts
          tmp.opts[id.tbs] <- as.character(tbs.matrix[i.tbs,])
          argv <- c("ms", nsam, nreps, tmp.opts)
          .Call("R_ms_main", argv, temp.file.ms, PACKAGE = "phyclust")
        }

      } else{
        ### No "tbs" option is used in opts.
        argv <- c("ms", nsam, nreps, new.opts)
        .Call("R_ms_main", argv, temp.file.ms, PACKAGE = "phyclust")
      }

      ### Finish the calls.
      if(is.null(temp.file)){
        # ret <- scan(file = temp.file.ms,
        #             what = "character", sep = "\n", quiet = TRUE)
        # class(ret) <- "ms"
        # unlink(temp.file.ms)
        # return(ret)
        ret <- readLines(con = temp.file.ms, warn = FALSE)
        ret <- ret[ret != ""]   # Drop the empty lines.
        class(ret) <- "ms"
        unlink(temp.file.ms)
        return(ret)
      }
    }
  } else{
    temp.file.ms <- tempfile("ms.")
    argv <- c("ms", "-h")
    try(.Call("R_ms_main", argv, temp.file.ms, PACKAGE = "phyclust"),
        silent = TRUE)
    unlink(temp.file.ms)
  }

  invisible()
} # End of ms().

print.ms <- function(x, ...){
  ms <- x
  cat(ms, sep = "\n")
} # End of print.ms().
