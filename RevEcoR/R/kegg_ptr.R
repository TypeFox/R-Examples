#'Metabolic profiles of KEGG organism Pan troglodytes (chimpanzee) (KEGG
#'organism code: ptr)
#'
#'kegg organism ptr metabolic information, which consists of enzymatic reactions
#'and metabolites.
#'
#'@details
#'
#'ptr metatolic information:
#'
#'\itemize{ \item .attrs.name: Enzymatic reactions that organism involved \item
#'substrate.name: Substrates of the corresponding reaction. \item product.name:
#'Products of the corresponding reaction. }
#'
#'@format A data frame with 1858 observations on three variables.
#'  
#'  [,1] .attrs.name, character (reaction: R)
#'  
#'  [,2] substrate.name, list (substrates: cpd)
#'  
#'  [,3] product.name, list (products: cpd)
#'  
#'  
#'  
#'  
#'@name kegg_ptr
NULL