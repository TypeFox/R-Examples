## ========================================
## Count number of time a subsequence appear in a sequence
## ========================================

seqeapplysub<-function(subseq,method=NULL,
                       constraint=NULL,rules=FALSE)
  {
    ## message("Event sequence analysis module is still experimental")
    if (!is.subseqelist(subseq)) {
      stop("subseq should be a subseqelist. See help on seqefsub.")
    }
    if (is.null(constraint)) {
      constraint <- subseq$constraint
    }
    if (!is.null(method)) {
      if (method=="count") {
        constraint$countMethod <- 2
      } else if (method=="presence") {
        constraint$countMethod <- 1
      } else if(method=="age") {
        constraint$countMethod <- 6
      } else if(method=="COBJ") {
        constraint$countMethod <- 1
      } else if(method=="CDIST_O") {
        constraint$countMethod <- 2
      } else if(method=="CWIN") {
        constraint$countMethod <- 3
      } else if(method=="CMINWIN") {
        constraint$countMethod <- 4
      } else if(method=="CDIST") {
        constraint$countMethod <- 5
      } else if (method%in%1:5) {
        constraint$countMethod <- method
      }
    }
    if (!inherits(constraint, "seqeconstraint")) {
      constraint = seqeconstraint()
      warning("[!] The constraint argument should be set using the seqeconstraint function. The provided constraint argument is ignored.")
    }
    if(!rules)
      {
        return(.Call("tmrmatrixsubseqinseq",
                     unlist(list(subseq$subseq)),
                     unlist(list(subseq$seqe)),
                     as.double(c(constraint$maxGap)),
                     as.double(c(constraint$windowSize)),
                     as.double(c(constraint$ageMin)),
                     as.double(c(constraint$ageMax)),
                     as.double(c(constraint$ageMaxEnd)),
                     as.double(c(constraint$countMethod)),
                     PACKAGE="TraMineR"))
      }
    else {
      return(.Call("tmrmatrixsubseqinseq",
                   unlist(list(subseq$subseq)),
                   unlist(list(subseq$subseq)),
                   as.double(c(constraint$maxGap)),
                   as.double(c(constraint$windowSize)),
                   as.double(c(constraint$ageMin)),
                   as.double(c(constraint$ageMax)),
                   as.double(c(constraint$ageMaxEnd)),
                   as.double(c(constraint$countMethod)),
                   PACKAGE="TraMineR"))
    }
  }
