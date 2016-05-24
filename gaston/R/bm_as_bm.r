as.bed.matrix <- function(x, fam, bim) {
  if ( !is.null(fam) & !is.data.frame(fam) ) stop('fam must be a data.frame or NULL.')
  if ( !is.null(bim) & !is.data.frame(bim) ) stop('bim must be a data.frame or NULL.')

  bed <- .Call('gg_as_matrix4', PACKAGE = 'gaston', x)

  ped <- if(!missing(fam)) {
            if(all(pednames %in% names(fam))) 
              fam
            else 
              stop('"fam" must contain "famid", "id", "father", "mother", "sex" and "pheno" variables')
         } else {
           if(is.null(rownames(x)))
             structure(list(), row.names = c(NA, -nrow(x)), class = "data.frame") # empty data frame with right number of lines
           else
             data.frame(famid = rownames(x), id = rownames(x), father = 0, mother = 0, sex = 0, pheno = 0, stringsAsFactors = FALSE)
         }  

  snps <- if(!missing(bim)) {
            if(all(snpnames %in% names(bim))) 
               bim
            else 
              stop('"bim" must contain "chr", "id", "dist", "pos", "A1" and "A2" variables')
          } else {
            if(is.null(colnames(x)))
             structure(list(), row.names = c(NA, -ncol(x)), class = "data.frame") #idem
           else
             data.frame(chr = NA, id = colnames(x), dist = NA, pos = NA, A1 = NA, A2 = NA, stringsAsFactors = FALSE)
          }

  x <- new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = FALSE)
  x
}

