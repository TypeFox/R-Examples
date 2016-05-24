### This file contains functions for "seq-gen".

seqgen <- function(opts = NULL, rooted.tree = NULL, newick.tree = NULL,
    input = NULL, temp.file = NULL){
  argv <- "seq-gen"

  if(! is.null(opts)){
    if(is.null(temp.file)){
      temp.file.seqgen <- tempfile("seqgen.")
    } else{
      temp.file.seqgen <- temp.file
    }
    temp.file.ms <- tempfile("ms.")

    if((! is.null(rooted.tree)) && (class(rooted.tree) == "phylo")){
      newick.tree <- write.tree(rooted.tree, digits = 12)
    }
    if((!is.null(newick.tree)) && (!is.null(input))){
      stop("rooted.tree/newick.tree and input can not work at the same time.")
    }
    if(! is.null(newick.tree)){
      write(newick.tree, file = temp.file.ms, sep = "")
    } else if(! is.null(input)){
      write(input, file = temp.file.ms, sep = "\n")
    } else{
      stop("A newick or rooted/phylo tree is required.")
    }

    argv <- c(argv, unlist(strsplit(opts, " ")), temp.file.ms)
    .Call("R_seq_gen_main", argv, temp.file.seqgen, PACKAGE = "phyclust")

    # ret <- scan(file = temp.file.seqgen,
    #             what = "character", sep = "\n", quiet = TRUE)
    # class(ret) <- "seqgen"
    # unlink(temp.file.seqgen)
    unlink(temp.file.ms)

    # return(ret)

    if(is.null(temp.file)){
      ret <- readLines(con = temp.file.seqgen, warn = FALSE)
      ret <- ret[ret != ""]   # Drop the empty lines.
      class(ret) <- "seqgen"
      unlink(temp.file.seqgen)
      return(ret)
    }
  } else{
    temp.file.seqgen <- tempfile("seqgen.")
    argv <- c(argv, "-h")
    try(.Call("R_seq_gen_main", argv, temp.file.seqgen, PACKAGE = "phyclust"),
        silent = TRUE)
    unlink(temp.file.seqgen)
  }

  invisible()
} # End of seqgen().

print.seqgen <- function(x, ...){
  seqgen <- x
  cat(seqgen, sep = "\n")
} # Enf of print.seqgen().

