##############
## DATASETS ##
##############

#' Autoimmune enriched regions as mapped on ImmunoChip
#'
#' Dataset. A consortium of 12 autoimmune diseases (Type 1 diabetes, Celiac
#' disease, Multiple Sclerosis, Crohns Disease, Primary Billiary Cirrhosis, 
#' Psoriasis, Rheumatoid Arthritis, Systemic Lupus Erytematosus, 
#' Ulcerative Colitis, Ankylosing Spondylitis, Autoimmune Thyroid Disease,
#' Juvenile Idiopathic Arthritis) created the ImmunoChip custom Illumina
#' iSelect microarray in order to investigate known regions from GWAS
#' associating with a p value < 5*10-8 with any of these diseases, using
#' dense mapping of locii. This object specifies the boundaries of these
#' regions, defined roughly as 0.1 centimorgan recombination distance either 
#' side of the top marker in each. The data is in the original build 36
#' coordinates as a GRanges object, but using functions in the humarray
#' package can easily be converted to build 37, 38 or RangedData/data.frame.
#' 
#' @name iChipRegionsB36
#' @docType data
#' @format An object of class GRanges
#' @keywords datasets
#' @references Cortes and Brown, Promise and pitfalls of the Immunochip (2011).
#'  Arthritis Research and Therapy 2011, 13:101
#' @seealso \code{\link{conv.36.37}} \code{\link{get.t1d.subset}}
#'  \code{\link{get.t1d.regions}} \code{\link{get.immunobase.snps}}
#' @examples
#' data(iChipRegionsB36)
#' prv(iChipRegionsB36)
#' \donttest{iChipRegionsB37 <- conv.36.37(iChipRegionsB36) }
NULL # need some sort of code to trigger devtools:document to pick up a dataset


#' ImmunoChip annotation object (built-in)
#'
#' Dataset. This object contains chromosome, position, chip SNP labels,
#' SNP rs-ids, for immunochip, plus allele codes based on a UVA/Sanger 
#' T1D dataset, which should be updated for use with a different dataset
#' with different allele coding (however the position and rs-id information
#' should still be applicable.
#' The data is in build 37 coordinates as a ChipInfo object, but using functions
#' convTo36 or convTo38 from the humarray package can easily convert this to
#' build 36 or 38.
#' 
#' @name ImmunoChipB37
#' @docType data
#' @format An object of class GRanges
#' @keywords datasets
#' @seealso \code{\link{convTo38}} \code{\link{convTo36}}
#'  \code{\link{A1}}
#' @examples
#' data(ImmunoChipB37)
#' \donttest{ImmunoChipB36 <- convTo36(ImmunoChipB37) }
NULL # need some sort of code to trigger devtools:document to pick up a dataset


