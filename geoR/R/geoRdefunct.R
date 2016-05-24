"geoRdefunct" <-
  function()
  {
    cat("\n")
    cat("The following functions are no longer used in geoR:")
    cat("---------------------------------------------------")
    cat("\nolsfit: use variofit() instead")
    cat("\nwlsfit: use variofit() instead")
    cat("\nlikfit.old: use likfit() instead")
    cat("\n")
  }

"olsfit" <- function(...)
  stop("this function is now obsolete.\nuse variofit() instead.")

"wlsfit" <- function(...)
  stop("this function is now obsolete.\nuse variofit() instead.")


"distdiag" <-
  function(coords)
  {
    ## returns the lower triangle of the matrix with euclidean distances
    ## between pairs of points, including the diagonal.
    ##
    coords <- as.matrix(coords)
    dimc <- dim(coords)
    if(dimc[2] == 1 & dimc[1] == 2)
      return(0)
    else{
      if(dimc[2] != 2)
        stop("coords must have two columns")
      nc <- dimc[1]
      .C("distdiag",
         as.double(coords[,1]),
         as.double(coords[,2]),
         as.integer(nc),
         out = as.double(rep(0, (nc * (nc+1)/2))),
         PACKAGE = "geoR")$out
    }
  }

#"cite.geoR" <- function()
#{
#    cat("\n")
#    cat("To cite geoR in publications, use\n\n")
#    msg <- "RIBEIRO Jr., P.J. & DIGGLE, P.J. (2001) geoR: A package for geostatistical analysis. R-NEWS, Vol 1, No 2, 15-18. ISSN 1609-3631."
#    writeLines(strwrap(msg, prefix = "  "))
#    cat("\n")
#    msg <- paste("Please cite geoR when using it for data analysis!")
#    writeLines(strwrap(msg))
#    cat("\nA BibTeX entry for LaTeX users is\n\n")
#    cat("  @Article{,\n")
#    cat("     title	   = {{geoR}: a package for geostatistical analysis},\n")
#    cat("     author        = {Ribeiro Jr., P.J. and Diggle, P.J.},\n")
#    cat("     journal       = {R-NEWS},\n")
#    cat("     year	   = {2001},\n")
#    cat("     volume	   = {1},\n")
#    cat("     number	   = {2},\n")
#    cat("     pages	   = {15--18},\n")
#    cat("     issn          = {1609-3631},\n")
#    cat("     url           = {http://cran.R-project.org/doc/Rnews}\n")
#    cat("   }\n\n")
#}

#geoR.options <- function(messages = TRUE, ...)
#{
#  res <- list(...)
#  res$messages <- messages
#  .geoR.options <<- res
#  return(invisible())
#}



