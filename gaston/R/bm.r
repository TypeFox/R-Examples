setClassUnion("data.frameOrNULL",members=c("data.frame", "NULL"))
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))


pednames <- c("famid", "id", "father", "mother", "sex", "pheno")
snpnames <- c("chr", "id", "dist", "pos", "A1", "A2")

snpstatnames0 <- c("N0", "N1", "N2", "NAs", "callrate", "maf", "hz")
snpstatnames <- c("N0", "N1", "N2", "NAs", "callrate", "maf", "hz", "hwe")
pedstatnames <- c("N0", "N1", "N2", "NAs", "callrate", "hz")

is.null.df <- function(x) is.data.frame(x) & nrow(x) == 0 & ncol(x) == 0

## Class bed.matrix
setClass("bed.matrix", representation( 
                   ped = 'data.frame',
                   snps = 'data.frame', 
                   bed = 'externalptr',
                   p = 'numericOrNULL',
                   mu = 'numericOrNULL',
                   sigma = 'numericOrNULL',
                   standardize_p = 'logical',
                   standardize_mu_sigma = 'logical' ))

setValidity('bed.matrix',
           function(object) {
             errors <- character()
             if ( object@standardize_p & object@standardize_mu_sigma ) 
                errors <- c(errors, "Only one center scale parameter can be TRUE.")
             if ( object@standardize_p & is.null(object@p) ) 
                errors <- c(errors, "If 'standardize_p' is TRUE, 'p' must be defined.")
             if ( object@standardize_mu_sigma & ( is.null(object@mu) | is.null(object@sigma) ) ) 
                errors <- c(errors, "If 'standardize_mu_sigma' is TRUE, 'mu' and 'sigma' must be defined.")
             if ( !is.null(object@p) & length(object@p) != ncol(object) ) 
                errors <- c(errors, "The length of 'p' must be equal to the number of markers.")
             if ( !is.null(object@mu) & length(object@mu) != ncol(object) ) 
                errors <- c(errors, "The length of 'mu' must be equal to the number of markers.")
             if ( !is.null(object@sigma) & length(object@sigma) != ncol(object) ) 
                errors <- c(errors, "The length of 'sigma' must be equal to the number of markers.")
             if(.Call("isnullptr",  PACKAGE = "gaston", object@bed))
                errors <- c(errors, 'The externalptr is broken')
             if ( length(errors)==0 ) return(TRUE) else return(errors)
           } );


setAs("bed.matrix", "matrix",
  function(from) {
    validObject(from)
    to <- if(from@standardize_p) 
      .Call('gg_m4_as_scaled_matrix_p', PACKAGE = 'gaston', from@bed, from@p)
    else if(from@standardize_mu_sigma)
      .Call('gg_m4_as_scaled_matrix_mu_sigma', PACKAGE = 'gaston', from@bed, from@mu, from@sigma)
    else
      .Call('gg_m4_as012', PACKAGE = 'gaston', from@bed)
    colnames(to) <- from@snps$id
    rownames(to) <- if(any(duplicated(from@ped$id))) paste(from@ped$fam, from@ped$id, sep="_")
                    else from@ped$id
    to
  } );

setGeneric('as.matrix')
setMethod("as.matrix", signature="bed.matrix", definition = function(x) as(x,"matrix") )

setAs("matrix", "bed.matrix", 
  function(from) {
    bed <- .Call('gg_as_matrix4', PACKAGE = 'gaston', from)

    ped <- if(is.null(rownames(from))) 
             structure(list(), row.names = c(NA, -nrow(from)), class = "data.frame") # empty data frame with right number of lines
           else 
             data.frame(famid = rownames(from), id = rownames(from), father = 0, mother = 0, sex = 0, pheno = 0, stringsAsFactors = FALSE)

    snp <- if(is.null(colnames(from)))
             structure(list(), row.names = c(NA, -ncol(from)), class = "data.frame") #idem
           else 
             data.frame(chr = NA, id = colnames(from), dist = NA, pos = NA, A1 = NA, A2 = NA, stringsAsFactors = FALSE)

    x <- new("bed.matrix", bed = bed, snps = snp, ped = ped, p = NULL, mu = NULL,
             sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
    if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
    x
  } );

setGeneric('dim')
setMethod("dim", signature = "bed.matrix", 
    function(x) c(.Call('gg_ninds', PACKAGE = 'gaston', x@bed), .Call('gg_nsnps', PACKAGE = 'gaston', x@bed)))

setGeneric('head')
setMethod( 'head', signature(x='bed.matrix'), function(x, nrow=10, ncol=10) print( as.matrix( x[1:min( nrow, nrow(x) ),1:min( ncol, ncol(x) )] ) ) )

setMethod(show, signature("bed.matrix"), 
  function(object) {
    if(.Call("isnullptr",  PACKAGE = "gaston", object@bed))
      cat("A bed.matrix with a broken externalptr!\nHint: don't save/load bed.matrices with other functions than write.bed.matrix and read.bed.matrix\n")
    else {
      cat('A bed.matrix with ', nrow(object), ' individuals and ', ncol(object), ' markers.\n', sep='')
      if(anyDuplicated(object@snps$id)) cat("There are some duplicated SNP id's\n")
      if(all(snpstatnames0 %in% names(object@snps))) {
        cat("snps stats are set\n");
        a <- sum(object@snps$NAs == nrow(object))
        if(a > 1)  cat("  There are ", a, " SNPs with a callrate equal to zero\n");
        if(a == 1) cat("  There is one SNP with a callrate equal to zero\n");
        a <- sum(object@snps$maf == 0, na.rm = TRUE)
        if(a > 1)  cat("  There are ", a, " monomorphic SNPs\n");
        if(a == 1) cat("  There is one monomorphic SNP\n");
      } else 
        cat("snps stats are not set (or incomplete)\n")
      # ped stats et autres
      if(anyDuplicated(object@ped$id)) cat("There are some duplicated individual id's\n")
      if(all(pedstatnames %in% names(object@ped))) {
        cat("ped stats are set\n");
        a <- sum(object@ped$NAs == ncol(object))
        if(a > 1)  cat("  There are ", a, " individuals with a callrate equal to zero\n");
        if(a == 1) cat("  There is one individual with a callrate equal to zero\n");
      } else cat("ped stats are not set (or incomplete)\n")
    }
  } 
)
