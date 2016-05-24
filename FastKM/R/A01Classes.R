if(!isClass("matrix or NULL")){
  setClassUnion("matrix or NULL", 
                members = c("matrix","NULL"))
}

#----------------------------------------------------------------------#
# Class for genotype data.                                             #
# Contains design matrix, kernel type to be used, and weights          #
#----------------------------------------------------------------------#
setClass(Class = "geno",
         slots = c(mat = 'matrix or NULL',
                   kernel = 'character',
                   weights = 'matrix'),
         prototype = list(mat = NULL,
                          kernel = 'linear',
                          weights = matrix(0,0,0)) )

#----------------------------------------------------------------------#
# Class for non-genotype data.                                         #
# Contains design matrix, kernel type to be used, and weights          #
#----------------------------------------------------------------------#
setClass(Class = "nongeno",
         slots = c(mat = 'matrix or NULL',
                   kernel = 'character',
                   weights = 'matrix'),
         prototype = list(mat = NULL,
                          kernel = 'linear',
                          weights = matrix(0,0,0)) )


