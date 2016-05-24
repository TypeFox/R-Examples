# +++++++++++++++++++++++++++++++++++++++++
# Interface between R and the software CaliFloPP
# AB: 22-09-2009
# +++++++++++++++++++++++++++++++++++++++++
decodpoly <- function(X,  ficp)  {
  # FUNCTION
  # write a pair of polygons, whose analysis is required,
  # on the parameter file
  # INPUT
  # X: vector of the 2 requested polygons
  # ficp: name of the file where the parameters are written
  #       or connection to it
  # RETURN
  # none
  # CALLED BY
  # decodparam,
  # (via "apply" on the matrix of the pairs of the requested polygons)
  # -------------------------------------------
cat(X[1], X[2], "\n", file=ficp, append=TRUE)
invisible()
} # end decodpoly


# -----------------------------------------------
decodparam <- function(X, nom, ficp, nofunc)  {
  # FUNCTION
  # write a parameter and its value on the parameter file
  # INPUT
  # X: value of the  parameter
  # nom: name of parameter
  # ficp: name of  the file where the parameters are written
  #       or connection to it
  # nofunc: numbers of the requested functions
    # RETURN
  # none
  # CALLED BY
  # - itself, via "mapply" on the list of the parameters
  # which are specific to  cub, or on the list of the parameters
  # which are specific to grid
  # - califlopp, via "mapply" on the list of the parameters
  # -------------------------------------------

switch(nom,
       nfunc  = {cat("nfunc", X, "\n", file=ficp, append=TRUE)},
       func  = {cat(X, "\n", file=ficp, append=TRUE)},
         input = {cat("input", X, "\n", file=ficp, append=TRUE)},
         warn.poly = {cat("warnpoly ",  file=ficp, append=TRUE)
                      if (X==TRUE) cat("1",  file=ficp, append=TRUE) else
                    cat("0", file=ficp, append=TRUE)
                    cat("\n", file=ficp, append=TRUE)},
         warn.conv = {cat("warnconv ",  file=ficp, append=TRUE)
                      if (X==TRUE) cat("1",  file=ficp, append=TRUE) else
                    cat("0", file=ficp, append=TRUE)
                    cat("\n", file=ficp, append=TRUE)},
         send.and.receive = {cat("sendreceive ",  file=ficp, append=TRUE)
                      if (X==TRUE) cat("1",  file=ficp, append=TRUE) else
                    cat("0", file=ficp, append=TRUE)
                    cat("\n", file=ficp, append=TRUE)},
       output= {cat("output", X, "\n", file=ficp, append=TRUE)},
         verbose = {cat("verbose ", file=ficp, append=TRUE)
                    if (X==TRUE) cat("1",  file=ficp, append=TRUE) else
                    cat("0", file=ficp, append=TRUE)
                    cat("\n", file=ficp, append=TRUE)},
       delim= {cat("delim ",  file=ficp, append=TRUE);
#               cat(gsub("\"", "", deparse(X)),  file=ficp, append=TRUE);
               cat( deparse(X),  file=ficp, append=TRUE);
               cat("\n", file=ficp, append=TRUE)
               },
       poly = {
    # Decode the pairs of requested polygons
  if (is.list(X)) 
    X <- matrix(unlist(X), ncol=2, byrow=TRUE)
  if (is.vector(X)) {
    if (length(X) !=2)
      stop("Vector of polygons should be 2")
    X <- matrix(X, ncol=2)
  }
    cat("ncouples ", nrow(X), "\n", file=ficp, append=TRUE)
  apply(X, 1, decodpoly, ficp)
}, # end poly
       method = { 
                    cat("method ", file=ficp, append=TRUE)
                   if (X=="grid") cat("1", file=ficp, append=TRUE) else cat("0",file= ficp, append=TRUE)
                    cat("\n", file=ficp, append=TRUE)},
       dz  = {  cat("dz ", X, "\n", file=ficp, append=TRUE) },
       dp  = {  cat("dp ", X, "\n", file=ficp, append=TRUE) },
       tz  = {  cat("tz ", X, "\n", file=ficp, append=TRUE) },
# cubature parameters       
         cub = { if (!is.list(X))
                 stop("Component cub of param should be a list")
                 mapply( decodparam, X, names(X),
                        MoreArgs=list(ficp, nofunc))},
         maxpts = {  if (length(X) != length(nofunc))
                    stop("Component maxpts of param should have the same length than dispf")
                     for (i in 1:length(X)) {
                       if ( (X[i] <= 0) || (X[i] >800000000))
                         stop("Component maxpts  of param should be in ]0, 800000000]")
                       
           if (!is.na(X[i])) { cat("maxpts ",nofunc[i], "\n", file=ficp, append=TRUE)
                           cat(format(X[i], scientific=F),
                               "\n", file=ficp, append=TRUE) }
             } # end for
                  },
       reler  = { if (length(X) != length(nofunc))
                    stop("Component reler of param should have the same length than dispf")
                    for (i in 1:length(X)) {
           if (!is.na(X[i])) { cat("reler ",nofunc[i], "\n", file=ficp, append=TRUE)
                           cat(X[i], "\n", file=ficp, append=TRUE) }
             } # end for
                  },
       abser  = {  if (length(X) != length(nofunc)) {
                    stop("Component abser of param should have the same length than dispf")
       }
                   for (i in 1:length(X)) {
           if (!is.na(X[i])) { cat("abser ",nofunc[i], "\n", file=ficp, append=TRUE)
                           cat(X[i], "\n", file=ficp, append=TRUE) }
             } # end for
                  },
# grid parameters       
         grid = { if (!is.list(X))
                 stop("Component grid of param should be a list")
                 mapply( decodparam, X, names(X),
                        MoreArgs=list(ficp, nofunc))},
       seed =  {cat("seed", X, "\n", file=ficp, append=TRUE)},
       step =  {cat("stepx", X[1], "\nstepy", X[2],  "\n",
                    file=ficp, append=TRUE)},
       nr =  {cat("nr", X, "\n", file=ficp, append=TRUE)},

         stop(paste("Invalid name for a parameter:", nom))
         ) # end switch
  invisible()
} # end decodparam

# -----------------------------------------------        
califlopp <- function(file, dispf=c(1,2),
                      param=NULL,
                      resfile=NULL)
{
  # FUNCTION
  # call califlopp_sd (mode non-interactif)
  # (via a wrapper programme in C )
  # INPUT
  # file: pathname of the polygons file
  # dispf: the R dispersion functions
  # integer vector for the compiled functions, or
  # names of R functions
  # param: - if character, pathname of the parameter file
  #        - if list, list of parameters
  #        - if NULL: default values of all parameters
  #        Among these parameters, when param is a list
  #      poly: list of vectors of length 2 or 2-columns matrix:
  #       identifiers of the required pairs of polygons
  #       as they are noted on the polygons file.
  #       Default: all pairs of polygons
  # resfile: pathname of the results file
  #     if NULL, the function returns NULL
  #       else the function returns the main results
  # RETURN
  #   NULL
  # -------------------------------------------------------
# Verif: param must be NULL, or a character string, or a list
  if (!is.character(param) && !is.null(param)) {
    if (!is.list(param))
      stop("param should be a filename or a list of parameters")
  }

  # The dispersal functions
  dispf <- c(dispf) # put into a vector
  nfun <- length(dispf)
  if (nfun >5)
    stop("number of dispersion functions should be <=5")


  if (is.numeric(dispf)) {
    dispfc <- dispf # dispfc: indices of the functions in C
    nofunc <- dispf
    dispf <- 0 # modif 6/6/2012 NULL
  }  else    {
    dispfc <- 0
    nofunc <- seq(1, length(dispf))
  }
                  
    
  # The file which stores the parameters: nomficp 
  nomficp <- NULL
  ficp <- NULL
  if (is.character(param))
    nomficp <- param # the parameter file is already given
  else {
    # Write the parameters on a temporary file in the tmp of the user
    nomficp <- tempfile(pattern = "califlopp")

    ficp <- file(nomficp, open="w")
    
    if (!is.null(param))
      mapply( decodparam, param, names(param),
             MoreArgs=list(ficp,nofunc))
    # When the dispersal functions are not in R,
    # store the indices of thoses requested
    if (all(dispfc>0)) {
      decodparam(nfun, "nfunc", ficp, nofunc)
      decodparam(dispfc, "func", ficp, nofunc)
    }
    close(ficp)
  } # end else
  
  if (is.null(resfile))
    resfile <- ""
  else
    if (!is.character(resfile))
      stop("resfile should be the pathname of a file")

  # When all the parameters have their default value,
  # no parameter file
    if (is.null(nomficp))
      nomficp<- ""
  # Call CaliFloPP:
  .C("CALLcaliflopp",  as.integer(nfun),
     as.character(file),
   as.character(nomficp),
   as.character(resfile),
     as.integer(dispfc),
     dispf,  new.env(), PACKAGE = "RCALI")
#     dispf,  parent.frame())

  # Delete the temporary parameters file when it is not a file
  # provided by the user
  if (!is.null(ficp) && !is.character(ficp)) {
       unlink(nomficp)
     }
invisible()
}  # end califlopp
