########## R-function: add.model.flag ##########

# Determines whether or not model is a `bona fide'
# additive model

# Last changed: 02 FEB 2005

add.model.flag <- function(spm.info)
{     
   flag <- FALSE

   if ((!is.null(spm.info$lin))&(!is.null(spm.info$pen)))
      flag <- TRUE
  
   if ((!is.null(spm.info$lin))&(!is.null(spm.info$krige)))
      flag <- TRUE

   if ((!is.null(spm.info$pen))&(!is.null(spm.info$krige)))
      flag <- TRUE

   if (!is.null(spm.info$pen$x))
      if ((!is.null(spm.info$pen))&(ncol(spm.info$pen$x)>1))
         flag <- TRUE

   if (!is.null(spm.info$lin$x))
      if ((!is.null(spm.info$lin))&(ncol(spm.info$lin$x)>1))
         flag <- TRUE

   return(flag)
}

######### End of add.model.flag ########
