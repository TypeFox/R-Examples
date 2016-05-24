.set.attr <- function (.data, .attr, .value) {
  attr(.data, .attr) <- .value
  .data
}


#' Tables with genetic code.
#' 
#' @docType data
#' 
#' @aliases AA_TABLE AA_TABLE_REVERSED
#' 
#' @description
#' Tables with genetic code.
#' 
#' @format
#' AA_TABLE:
#' 
#' \code{Class 'table'  Named chr [1:65] "K" "N" "K" "N" ...
#' ..- attr(*, "names")= chr [1:65] "AAA" "AAC" "AAG" "AAT" ...}
#' 
#' AA_TABLE_REVERSED:
#' 
#' \code{List of 22
#' $ *: chr [1:3] "TAA" "TAG" "TGA"
#' $ A: chr [1:4] "GCA" "GCC" "GCG" "GCT"
#' $ C: chr [1:2] "TGC" "TGT"
#' $ D: chr [1:2] "GAC" "GAT"
#' ...
#' }
#' 
#' @examples
#' \dontrun{
#' AA_TABLE['ATG']  #  => "M"
#' AA_TABLE_REVERSED['K']  #  => list(K = c("AAA", "AAG"))
#' }
AA_TABLE <- table(c('TCA', 'TCG', 'TCC', 'TCT', 'TTT', 'TTC', 'TTA', 'TTG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT',
                    'TGC', 'TGA', 'TGG', 'CTA', 'CTG', 'CTC', 'CTT', 'CCA', 'CCG', 'CCC', 'CCT', 'CAT', 'CAC',
                    'CAA', 'CAG', 'CGA', 'CGG', 'CGC', 'CGT', 'ATT', 'ATC', 'ATA', 'ATG', 'ACA', 'ACG', 'ACC',
                    'ACT', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTA', 'GTG', 'GTC', 'GTT',
                    'GCA', 'GCG', 'GCC', 'GCT', 'GAT', 'GAC', 'GAA', 'GAG', 'GGA', 'GGG', 'GGC', 'GGT', 'NNN'))
AA_TABLE[c('TCA', 'TCG', 'TCC', 'TCT', 'TTT', 'TTC', 'TTA', 'TTG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT',
           'TGC', 'TGA', 'TGG', 'CTA', 'CTG', 'CTC', 'CTT', 'CCA', 'CCG', 'CCC', 'CCT', 'CAT', 'CAC',
           'CAA', 'CAG', 'CGA', 'CGG', 'CGC', 'CGT', 'ATT', 'ATC', 'ATA', 'ATG', 'ACA', 'ACG', 'ACC',
           'ACT', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTA', 'GTG', 'GTC', 'GTT',
           'GCA', 'GCG', 'GCC', 'GCT', 'GAT', 'GAC', 'GAA', 'GAG', 'GGA', 'GGG', 'GGC', 'GGT', 'NNN')] <- c('S', 'S', 'S', 'S', 'F', 'F', 'L', 'L', 'Y', 'Y', '*', '*', 'C', 'C', '*',
                                                                                                            'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R',
                                                                                                            'R', 'R', 'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S',
                                                                                                            'S', 'R', 'R', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E',
                                                                                                            'G', 'G', 'G', 'G', '~')
AA_TABLE_REVERSED <- sapply(unique(AA_TABLE), function (aa) { names(AA_TABLE)[AA_TABLE == aa] })
AA_TABLE_REVERSED <- AA_TABLE_REVERSED[order(names(AA_TABLE_REVERSED))]


#' Alphabets of TCR and Ig gene segments.
#' 
#' @docType data
#' 
#' @name segments.alphabets
#' 
#' @aliases genealphabets HUMAN_TRAV HUMAN_TRAJ HUMAN_TRBV HUMAN_TRBD HUMAN_TRBJ HUMAN_TRBV_MITCR HUMAN_TRGV HUMAN_TRGJ HUMAN_TRDV HUMAN_TRDD HUMAN_TRDJ HUMAN_IGHV HUMAN_IGHD HUMAN_IGHJ HUMAN_IGKV HUMAN_IGKJ HUMAN_IGLJ HUMAN_IGLV HUMAN_TRBV_ALS HUMAN_TRBV_FAM HUMAN_TRBV_GEN MACMUL_TRBJ MACMUL_TRBV MOUSE_TRBJ MOUSE_TRBV MOUSE_TRAV MOUSE_TRAJ MOUSE_IGKV MOUSE_IGKJ MOUSE_IGHV MOUSE_IGHD MOUSE_IGHJ MOUSE_TRDD MOUSE_TRDV MOUSE_TRDJ MOUSE_TRGV MOUSE_TRGJ MOUSE_IGLJ MOUSE_IGLV
#' 
#' @usage
#' HUMAN_TRAV
#' 
#' HUMAN_TRAJ
#' 
#' HUMAN_TRBV
#' 
#' HUMAN_TRBD
#' 
#' HUMAN_TRBJ
#' 
#' HUMAN_TRBV_MITCR
#' 
#' HUMAN_TRBV_ALS
#' 
#' HUMAN_TRGV
#' 
#' HUMAN_TRGJ
#' 
#' HUMAN_TRDV
#' 
#' HUMAN_TRDD
#' 
#' HUMAN_TRDJ
#' 
#' MOUSE_TRBV
#' 
#' MOUSE_TRBJ
#' 
#' MOUSE_TRAV
#' 
#' MOUSE_TRAJ
#' 
#' MOUSE_IGKV
#' 
#' MOUSE_IGKJ
#' 
#' MOUSE_IGHV
#' 
#' MOUSE_IGHD
#' 
#' MOUSE_IGHJ
#' 
#' MACMUL_TRBV
#' 
#' MACMUL_TRBJ
#' 
#' HUMAN_IGHV
#' 
#' HUMAN_IGHD
#' 
#' HUMAN_IGHJ
#' 
#' HUMAN_IGLV
#' 
#' HUMAN_IGLJ
#' 
#' MOUSE_IGLJ
#' 
#' MOUSE_IGLV
#' 
#' @description
#' Character vector with names for segments. With \code{tcR} we provided alphabets for all alpha, beta,
#' gamma and delta chains gene segments.
#' 
#' @format
#' Each \code{<SPECIES>_<GENES>} is a character vector. <SPECIES> is an identifier of species, <GENES> is concatenated three
#' identifiers of cell type ("TR**" for TCR, "IG**" for Ig), chain (e.g., "**A*" for alpha chains) and gene segment ("***V" for V(ariable) gene segment, 
#' "***J" for J(oining) gene segment, "***D" for D(iversity) gene segment).
#' 
#' @examples
#' \dontrun{
#' HUMAN_TRBV[1]  #  => "TRBV10-1"
#' }
HUMAN_TRAV <- c('TRAV1-1', 'TRAV1-2', 'TRAV10', 'TRAV11', 'TRAV12-1', 'TRAV12-2', 'TRAV12-3', 'TRAV13-1',
                'TRAV13-2', 'TRAV14/DV4', 'TRAV16', 'TRAV17', 'TRAV18', 'TRAV19', 'TRAV2', 'TRAV20',
                'TRAV21', 'TRAV22', 'TRAV23/DV6', 'TRAV24', 'TRAV25', 'TRAV26-1', 'TRAV26-2', 'TRAV27',
                'TRAV29/DV5', 'TRAV3', 'TRAV30', 'TRAV34', 'TRAV35', 'TRAV36/DV7', 'TRAV38-1', 'TRAV38-2/DV8',
                'TRAV39', 'TRAV4', 'TRAV40', 'TRAV41', 'TRAV5', 'TRAV6', 'TRAV7', 'TRAV8-1',
                'TRAV8-2', 'TRAV8-3', 'TRAV8-4', 'TRAV8-6', 'TRAV8-7', 'TRAV9-1', 'TRAV9-2')
HUMAN_TRAJ <- c('TRAJ10', 'TRAJ11', 'TRAJ12', 'TRAJ13', 'TRAJ14', 'TRAJ15', 'TRAJ16', 'TRAJ17',
                'TRAJ18', 'TRAJ20', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ24', 'TRAJ26', 'TRAJ27',
                'TRAJ28', 'TRAJ29', 'TRAJ3', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', 'TRAJ34',
                'TRAJ36', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ4', 'TRAJ40', 'TRAJ41', 'TRAJ42',
                'TRAJ43', 'TRAJ44', 'TRAJ45', 'TRAJ46', 'TRAJ47', 'TRAJ48', 'TRAJ49', 'TRAJ5',
                'TRAJ50', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ56', 'TRAJ57', 'TRAJ6', 'TRAJ7',
                'TRAJ8', 'TRAJ9')

HUMAN_TRBV <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-3',
                'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2',
                'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3', 'TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6',
                'TRBV6-7', 'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8',
                'TRBV7-9', 'TRBV9')
HUMAN_TRBV <- .set.attr(HUMAN_TRBV, "column", "V.gene")
HUMAN_TRBD <- c('TRBD1', 'TRBD2')
HUMAN_TRBJ <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7')

HUMAN_TRBV_MITCR <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4, TRBV12-3',
                      'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2',
                      'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                      'TRBV3-1', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                      'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3, TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-7',
                      'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV7-9', 'TRBV9')

HUMAN_TRBV_ALS <- c('TRBV10-1', 'TRBV10-2', 'TRBV10-3', 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-4', 'TRBV12-3',
                    'TRBV12-5', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 'TRBV2*01',
                    'TRBV20-1', 'TRBV21-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1',
                    'TRBV3-1*01', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 'TRBV4-3', 'TRBV5-1', 'TRBV5-4', 'TRBV5-5',
                    'TRBV5-6', 'TRBV5-8', 'TRBV6-1', 'TRBV6-3', 'TRBV6-2', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6',
                    'TRBV6-7', 'TRBV7-1', 'TRBV7-2', 'TRBV7-3', 'TRBV7-4', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8',
                    'TRBV7-9', 'TRBV9')
HUMAN_TRBV_ALS <- .set.attr(HUMAN_TRBV_ALS, "column", "V.allele")

# HUMAN_TRBV_FAM <- c()
# HUMAN_TRBV_GEN <- c()
# HUMAN_TRBV_ALS <- c()

HUMAN_TRDD = c('TRDD1', 'TRDD2', 'TRDD3')

HUMAN_TRDJ = c('TRDJ1', 'TRDJ3', 'TRDJ4', 'TRDJ2')

HUMAN_TRDV = c('TRDV3', 'TRAV38-2/DV8', 'TRDV2', 'TRAV23/DV6', 'TRAV29/DV5', 'TRAV36/DV7', 
               'TRAV14/DV4', 'TRDV1')

HUMAN_TRGJ = c('TRGJ1', 'TRGJP1', 'TRGJ2', 'TRGJP', 'TRGJP2')

HUMAN_TRGV = c('TRGV5P', 'TRGV10', 'TRGV5', 'TRGV3', 'TRGVA', 'TRGV8', 
               'TRGV4', 'TRGV11', 'TRGV2', 'TRGV9', 'TRGV1')

HUMAN_IGHD = c('IGHD4-4', 'IGHD4-23', 'IGHD2-2', 'IGHD1-7', 'IGHD3-3', 'IGHD5-18', 
               'IGHD4/OR15-4a', 'IGHD2-21', 'IGHD4/OR15-4b', 'IGHD1/OR15-1a', 'IGHD1-14', 
               'IGHD6-13', 'IGHD5-24', 'IGHD2/OR15-2a', 'IGHD5-12', 'IGHD6-19', 
               'IGHD4-11', 'IGHD5-5', 'IGHD1/OR15-1b', 'IGHD3/OR15-3a', 'IGHD4-17', 
               'IGHD6-25', 'IGHD1-1', 'IGHD3-16', 'IGHD2-8', 'IGHD3-10', 
               'IGHD3/OR15-3b', 'IGHD5/OR15-5a', 'IGHD3-9', 'IGHD7-27', 'IGHD1-20', 
               'IGHD1-26', 'IGHD3-22', 'IGHD6-6', 'IGHD2-15', 'IGHD2/OR15-2b', 
               'IGHD5/OR15-5b')

HUMAN_IGHJ = c('IGHJ4', 'IGHJ1', 'IGHJ6', 'IGHJ3', 'IGHJ2', 'IGHJ5')


HUMAN_IGHV = c('IGHV4-4', 'IGHV7-40', 'IGHV2-5', 'IGHV4-28', 'IGHV6-1', 'IGHV3/OR16-6', 
               'IGHV1-69D', 'IGHV1-2', 'IGHV5-51', 'IGHV2/OR16-5', 'IGHV3-30-52', 
               'IGHV3-53', 'IGHV3/OR16-13', 'IGHV3-74', 'IGHV1-69-2', 'IGHV3-22', 
               'IGHV4-34', 'IGHV1-45', 'IGHV3-35', 'IGHV3-64', 'IGHV4-61', 
               'IGHV3-47', 'IGHV7-34-1', 'IGHV4-59', 'IGHV3-19', 'IGHV3-15', 
               'IGHV4-39', 'IGHV2-70D', 'IGHV3-20', 'IGHV1-18', 'IGHV3-30-22', 
               'IGHV4-38-2', 'IGHV3-54', 'IGHV5-10-1', 'IGHV3-33', 'IGHV3-21', 
               'IGHV3-7', 'IGHV1/OR15-2', 'IGHV4-31', 'IGHV3/OR16-8', 'IGHV3-23D', 
               'IGHV3-30-42', 'IGHV3-62', 'IGHV3-23', 'IGHV1-8', 'IGHV3/OR16-15', 
               'IGHV1/OR15-1', 'IGHV3-43D', 'IGHV4-30-2', 'IGHV3-NL1', 'IGHV3-64D', 
               'IGHV3-13', 'IGHV3-52', 'IGHV3-11', 'IGHV3-30-33', 'IGHV1/OR15-4', 
               'IGHV5-78', 'IGHV4-55', 'IGHV3-16', 'IGHV3-33-2', 'IGHV7-4-1', 
               'IGHV2-26', 'IGHV4-30-4', 'IGHV1/OR21-1', 'IGHV1/OR15-9', 'IGHV3/OR15-7', 
               'IGHV3-9', 'IGHV3-30', 'IGHV3-29', 'IGHV3-38', 'IGHV2-70', 
               'IGHV1-NL1', 'IGHV3-30-5', 'IGHV1/OR15-5', 'IGHV1-3', 'IGHV7-81', 
               'IGHV3/OR16-14', 'IGHV1-68', 'IGHV3-63', 'IGHV3/OR16-10', 'IGHV4/OR15-8', 
               'IGHV1-58', 'IGHV3-25', 'IGHV3-69-1', 'IGHV1-38-4', 'IGHV3/OR16-16', 
               'IGHV1-69', 'IGHV3-32', 'IGHV1-46', 'IGHV1-24', 'IGHV3/OR16-9', 
               'IGHV3-38-3', 'IGHV3-30-2', 'IGHV3-48', 'IGHV1/OR15-3', 'IGHV3-30-3', 
               'IGHV3-73', 'IGHV3-49', 'IGHV3-66', 'IGHV3-43', 'IGHV2-10', 
               'IGHV3/OR16-12', 'IGHV3-71', 'IGHV3-72')

HUMAN_IGKJ = c('IGKJ5', 'IGKJ4', 'IGKJ3', 'IGKJ2', 'IGKJ1')

HUMAN_IGKV = c('IGKV2-30', 'IGKV1-17', 'IGKV1-5', 'IGKV1/OR1-1', 'IGKV1/OR2-11', 'IGKV2D-30', 
               'IGKV1-9', 'IGKV2-40', 'IGKV6D-41', 'IGKV1D-43', 'IGKV3-7', 
               'IGKV6D-21', 'IGKV2-29', 'IGKV1/OR2-1', 'IGKV1-13', 'IGKV3D-20', 
               'IGKV1-NL1', 'IGKV1D-12', 'IGKV3D-7', 'IGKV1/OR10-1', 'IGKV1/OR2-108', 
               'IGKV1D-16', 'IGKV1/OR2-118', 'IGKV1-39', 'IGKV1/OR2-9', 'IGKV1D-17', 
               'IGKV1/ORY-1', 'IGKV1/OR-4', 'IGKV1D-8', 'IGKV1D-42', 'IGKV2D-24', 
               'IGKV1-27', 'IGKV1/OR15-118', 'IGKV1/OR22-5', 'IGKV1/OR9-1', 'IGKV2/OR2-7D', 
               'IGKV2D-18', 'IGKV1/OR-3', 'IGKV1/OR9-2', 'IGKV1D-37', 'IGKV1D-39', 
               'IGKV2D-29', 'IGKV5-2', 'IGKV2D-28', 'IGKV3/OR2-268', 'IGKV2D-26', 
               'IGKV3-20', 'IGKV1-8', 'IGKV2D-40', 'IGKV6-21', 'IGKV3-11', 
               'IGKV3-15', 'IGKV1/OR2-3', 'IGKV1-12', 'IGKV2/OR22-4', 'IGKV1/OR2-2', 
               'IGKV1/OR-2', 'IGKV7-3', 'IGKV1/OR2-0', 'IGKV3D-11', 'IGKV1-16', 
               'IGKV1D-33', 'IGKV1-33', 'IGKV3D-15', 'IGKV1D-13', 'IGKV2-18', 
               'IGKV1-6', 'IGKV2-28', 'IGKV4-1', 'IGKV1-37', 'IGKV2-4', 
               'IGKV2-24')

HUMAN_IGLJ = c('IGLJ6', 'IGLJ5', 'IGLJ3', 'IGLJ4', 'IGLJ1', 'IGLJ7', 
               'IGLJ2')

HUMAN_IGLV = c('IGLV2-11', 'IGLV7-46', 'IGLV1-44', 'IGLV5-52', 'IGLV4-3', 'IGLV3-19', 
               'IGLV3-27', 'IGLV2-18', 'IGLV2-8', 'IGLV3-21', 'IGLV2-33', 
               'IGLV3-16', 'IGLV3-13', 'IGLV2-14', 'IGLV1-62', 'IGLV3-31', 
               'IGLV3-9', 'IGLV4-60', 'IGLV1-51', 'IGLV7-43', 'IGLV5-45', 
               'IGLV8-61', 'IGLV3-1', 'IGLV1-40', 'IGLV3-25', 'IGLV4-69', 
               'IGLV9-49', 'IGLV5-39', 'IGLV1-47', 'IGLV2-NL1', 'IGLV2-34', 
               'IGLV11-55', 'IGLV2-5', 'IGLV8/OR8-1', 'IGLV10-54', 'IGLV3-12', 
               'IGLV3-22', 'IGLV1-41', 'IGLV5-37', 'IGLV5-48', 'IGLV1-36', 
               'IGLV3-10', 'IGLV3-32', 'IGLV6-57', 'IGLV2-23', 'IGLV1-50')

MOUSE_TRAV <- c('TRAV1', 'TRAV10', 'TRAV10D', 'TRAV10N', 'TRAV11', 'TRAV11D', 'TRAV11N', 'TRAV12-1', 'TRAV12-2', 'TRAV12-3', 'TRAV12D-1', 'TRAV12D-2', 
                'TRAV12D-3', 'TRAV12N-1', 'TRAV12N-2', 'TRAV12N-3', 'TRAV13-1', 'TRAV13-2', 'TRAV13-3', 'TRAV13-4/DV7', 'TRAV13-5', 'TRAV13D-1', 
                'TRAV13D-2', 'TRAV13D-3', 'TRAV13D-4', 'TRAV13N-1', 'TRAV13N-2', 'TRAV13N-3', 'TRAV13N-4', 'TRAV14-1', 'TRAV14-2', 'TRAV14-3', 'TRAV14D-1', 
                'TRAV14D-2', 'TRAV14D-3/DV8', 'TRAV14N-1', 'TRAV14N-2', 'TRAV14N-3', 'TRAV15-1/DV6-1', 'TRAV15-2/DV6-2', 'TRAV15D-1/DV6D-1', 
                'TRAV15D-2/DV6D-2', 'TRAV15N-1', 'TRAV15N-2', 'TRAV16', 'TRAV16D/DV11', 'TRAV16N', 'TRAV17', 'TRAV18', 'TRAV19', 'TRAV2', 'TRAV20', 
                'TRAV21/DV12', 'TRAV3-1', 'TRAV3-3', 'TRAV3-4', 'TRAV3D-3', 'TRAV3N-3', 'TRAV4-2', 'TRAV4-3', 'TRAV4-4/DV10', 'TRAV4D-2', 'TRAV4D-3', 
                'TRAV4D-4', 'TRAV4N-3', 'TRAV4N-4', 'TRAV5-1', 'TRAV5-2', 'TRAV5-4', 'TRAV5D-4', 'TRAV5N-4', 'TRAV6-1', 'TRAV6-2', 'TRAV6-3', 'TRAV6-4', 
                'TRAV6-5', 'TRAV6-6', 'TRAV6-7/DV9', 'TRAV6D-3', 'TRAV6D-4', 'TRAV6D-5', 'TRAV6D-6', 'TRAV6D-7', 'TRAV6N-5', 'TRAV6N-6', 'TRAV6N-7', 
                'TRAV7-1', 'TRAV7-2', 'TRAV7-3', 'TRAV7-4', 'TRAV7-5', 'TRAV7-6', 'TRAV7D-2', 'TRAV7D-3', 'TRAV7D-4', 'TRAV7D-5', 'TRAV7D-6', 'TRAV7N-4',
                'TRAV7N-5', 'TRAV7N-6', 'TRAV8-1', 'TRAV8-2', 'TRAV8D-1', 'TRAV8D-2', 'TRAV8N-2', 'TRAV9-1', 'TRAV9-2', 'TRAV9-3', 'TRAV9-4', 'TRAV9D-1',
                'TRAV9D-2', 'TRAV9D-3', 'TRAV9D-4', 'TRAV9N-2', 'TRAV9N-3', 'TRAV9N-4')
MOUSE_TRAJ <- c('TRAJ11', 'TRAJ12', 'TRAJ13', 'TRAJ15', 'TRAJ16', 'TRAJ17', 'TRAJ18', 'TRAJ2', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ24', 'TRAJ26', 'TRAJ27',
                'TRAJ28', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', 'TRAJ34', 'TRAJ35', 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ4', 'TRAJ40', 'TRAJ42', 'TRAJ43', 'TRAJ45', 'TRAJ46', 
                'TRAJ48', 'TRAJ49', 'TRAJ5', 'TRAJ50', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ56', 'TRAJ57', 'TRAJ58', 'TRAJ59', 'TRAJ6', 'TRAJ7', 'TRAJ9')

MOUSE_TRBV <- c('TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 
                'TRBV14', 'TRBV15', 'TRBV16', 'TRBV17', 'TRBV19', 'TRBV2', 'TRBV20',
                'TRBV23', 'TRBV24', 'TRBV26', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 
                'TRBV4', 'TRBV5')
MOUSE_TRBJ <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-7')

MOUSE_IGKV <- c('IGKV4-54','IGKV4-59','IGKV4-60','IGKV4-61','IGKV8-16','IGKV9-123',
                'IGKV9-128','IGKV1-131','IGKV14-130','IGKV12-89','IGKV1-132',
                'IGKV3-4','IGKV4-91','IGKV8-23-1','IGKV4-77','IGKV6-15',
                'IGKV6-25','IGKV6-14','IGKV9-119','IGKV4-58','IGKV4-74',
                'IGKV3-10','IGKV1-35','IGKV9-124','IGKV3-7','IGKV8-19',
                'IGKV1-88','IGKV8-34','IGKV7-33','IGKV13-76','IGKV4-63',
                'IGKV2-a','IGKV3-3','IGKV4-62','IGKV18-36','IGKV1-133',
                'IGKV3-9','IGKV6-23','IGKV12-98','IGKV3-8','IGKV1-110',
                'IGKV4-57-1','IGKV12-47','IGKV12-44','IGKV4-70','IGKV4-90',
                'IGKV4-73','IGKV4-79','IGKV6-32','IGKV4-56','IGKV2-137',
                'IGKV14-126','IGKV5-48','IGKV4-69','IGKV1-115','IGKV4-57',
                'IGKV6-13','IGKV3-12','IGKV4-72','IGKV11-125','IGKV3-2',
                'IGKV12-38','IGKV6-29','IGKV6-b','IGKV20-101-2','IGKV6-20',
                'IGKV14-111','IGKV10-95','IGKV4-81','IGKV4-92','IGKV15-103',
                'IGKV19-93','IGKV4-52','IGKV2-109','IGKV5-45','IGKV5-39',
                'IGKV8-26','IGKV1-122','IGKV1-135','IGKV6-c','IGKV12-46',
                'IGKV17-121','IGKV6-17','IGKV13-84','IGKV4-50','IGKV1-117',
                'IGKV14-100','IGKV4-55','IGKV4-80','IGKV12-e','IGKV8-18',
                'IGKV4-68','IGKV5-43','IGKV6-d','IGKV3-1','IGKV10-94',
                'IGKV9-120','IGKV13-82','IGKV8-28','IGKV8-21','IGKV4-51',
                'IGKV12-40','IGKV8-30','IGKV9-129','IGKV15-101','IGKV4-86',
                'IGKV2-112','IGKV4-75','IGKV8-24','IGKV12-41','IGKV13-85',
                'IGKV17/OR19-2','IGKV4-83','IGKV4-71','IGKV8-27','IGKV15-102',
                'IGKV5-37','IGKV4-78','IGKV17-127','IGKV17/OR16-3','IGKV1-99',
                'IGKV3-5','IGKV4-53','IGKV10-96','IGKV16-104')

MOUSE_IGKJ <- c('IGKJ5','IGKJ2','IGKJ1','IGKJ4','IGKJ3')

MOUSE_IGHV = c('IGHV1-69','IGHV1S96','IGHV5-17','IGHV3-1','IGHV5-6-2','IGHV1-63',
               'IGHV8-13','IGHV5-6','IGHV1-28','IGHV1-47','IGHV1-18',
               'IGHV1S14','IGHV1S29','IGHV1S70','IGHV1S135','IGHV14S4',
               'IGHV5-12-4','IGHV1-27','IGHV3-4','IGHV5-6-3','IGHV12-2',
               'IGHV8-5','IGHV1-84','IGHV1S10','IGHV1S136','IGHV8-4',
               'IGHV1-85','IGHV8-2','IGHV1-14','IGHV1S32','IGHV1-31',
               'IGHV1-5','IGHV1-81','IGHV1-78','IGHV13-2','IGHV7-4',
               'IGHV6-4','IGHV1-79','IGHV5-9-1','IGHV3S1','IGHV1-17-1',
               'IGHV2-6-1','IGHV5-6-5','IGHV1S118','IGHV2-6','IGHV1-56',
               'IGHV9-3-1','IGHV1S130','IGHV1-43','IGHV1-46','IGHV2-4-1',
               'IGHV5-9-5','IGHV5-9-2','IGHV1S112','IGHV1-19','IGHV10-1',
               'IGHV3-2','IGHV5S4','IGHV8-7','IGHV5-6-1','IGHV5S21',
               'IGHV8-12','IGHV8-10','IGHV1S124','IGHV1-72','IGHV1S82',
               'IGHV1S81','IGHV8-8','IGHV12-1-1','IGHV1-26','IGHV1S22',
               'IGHV3-6','IGHV1-71','IGHV12-2-1','IGHV1S65','IGHV1-58',
               'IGHV1-23','IGHV1-7','IGHV7-2','IGHV1-53','IGHV1S56',
               'IGHV1-37','IGHV1-24','IGHV1S31','IGHV1S83','IGHV1S20',
               'IGHV4-1','IGHV3-8','IGHV13-1','IGHV8S6','IGHV1S34',
               'IGHV1S35','IGHV1-59','IGHV8S9','IGHV1S47','IGHV1-62-1',
               'IGHV10-3','IGHV1S111','IGHV1S21','IGHV1S87','IGHV1-48',
               'IGHV5-6-6','IGHV8-8-1','IGHV1S16','IGHV1-49','IGHV1S92',
               'IGHV1-20','IGHV1-70','IGHV1-42','IGHV1S100','IGHV2-6-6',
               'IGHV5-6-4','IGHV10S4','IGHV1S134','IGHV1S95','IGHV5-16',
               'IGHV3-7','IGHV1S26','IGHV1-54','IGHV1-67','IGHV2-6-5',
               'IGHV1S52','IGHV1-51','IGHV1-66','IGHV1-12','IGHV1S61',
               'IGHV1S37','IGHV14-3','IGHV1S113','IGHV2-3-1','IGHV1-83',
               'IGHV1-16','IGHV1S127','IGHV6-6','IGHV1-55','IGHV1S12',
               'IGHV1-25','IGHV1S11','IGHV1S41','IGHV1S53','IGHV1S74',
               'IGHV5-15','IGHV1-80','IGHV1-64','IGHV1S107','IGHV1S68',
               'IGHV5S12','IGHV2-5-1','IGHV2S3','IGHV5-1','IGHV4-2',
               'IGHV2-6-3','IGHV1S72','IGHV2-6-8','IGHV1S84','IGHV2-7',
               'IGHV1S46','IGHV3S7','IGHV12-1','IGHV6S2','IGHV5-12-1',
               'IGHV9-2','IGHV9-3','IGHV11-1','IGHV1-62-2','IGHV7-3',
               'IGHV1-13','IGHV1-35','IGHV12-1-2','IGHV14-4','IGHV1S19',
               'IGHV1-11','IGHV1S122','IGHV1S5','IGHV1S120','IGHV2-9',
               'IGHV1-22','IGHV1-19-1','IGHV1-60','IGHV1-61','IGHV1S44',
               'IGHV5S24','IGHV3-5','IGHV1S126','IGHV1S101','IGHV3-3',
               'IGHV9-1','IGHV1-77','IGHV2-3','IGHV5-21','IGHV1S50',
               'IGHV1S33','IGHV2-2','IGHV5-12','IGHV6-7','IGHV8-6',
               'IGHV12-3','IGHV14-2','IGHV14-1','IGHV5-12-2','IGHV6-5',
               'IGHV1S121','IGHV1S108','IGHV8-14','IGHV9-2-1','IGHV1-4',
               'IGHV10S3','IGHV8-11','IGHV5-9','IGHV1S9','IGHV16-1',
               'IGHV1-87','IGHV2-9-2','IGHV1-21-1','IGHV1-36','IGHV1-39',
               'IGHV1S78','IGHV1-62-3','IGHV1-76','IGHV2-6-2','IGHV2-6-4',
               'IGHV9S7','IGHV1-9','IGHV9-4','IGHV1S49','IGHV1-50',
               'IGHV5-9-4','IGHV1S73','IGHV6-3','IGHV2-5','IGHV1-15',
               'IGHV11-2','IGHV15-2','IGHV5-9-3','IGHV2-2-1','IGHV1S45',
               'IGHV2-4','IGHV1S103','IGHV1-82','IGHV5-2','IGHV5-4',
               'IGHV1S36','IGHV9S8','IGHV1S17','IGHV1S28','IGHV1-8',
               'IGHV7-1','IGHV1S30','IGHV1-52','IGHV1S137','IGHV2-9-1',
               'IGHV1-74','IGHV8-9','IGHV2-6-7','IGHV1S110','IGHV1S132',
               'IGHV1-32','IGHV1S18','IGHV5S9','IGHV1S15','IGHV1-21',
               'IGHV1-34','IGHV6S4','IGHV1S55','IGHV1S75','IGHV6S3',
               'IGHV1S40','IGHV2-2-2','IGHV1S51','IGHV1S67','IGHV1-75',
               'IGHV8S2')

MOUSE_IGHJ = c('IGHJ1','IGHJ3','IGHJ2','IGHJ4')

MOUSE_IGHD = c('IGHD2-10','IGHD5-4','IGHD2-13','IGHD2-4','IGHD1-3','IGHD2-2',
               'IGHD2-6','IGHD3-2','IGHD3-3','IGHD4-1','IGHD5-6',
               'IGHD2-11','IGHD2-7','IGHD2-5','IGHD5-5','IGHD6-4',
               'IGHD5-2','IGHD5-3','IGHD2-14','IGHD6-1','IGHD1-1',
               'IGHD6-3','IGHD1-2','IGHD2-9','IGHD2-1','IGHD2-8',
               'IGHD3-1','IGHD5-1','IGHD6-2','IGHD2-3','IGHD2-12')

MOUSE_IGLV = c('IGLV5','IGLV2','IGLV6','IGLV4','IGLV3','IGLV7',
               'IGLV1','IGLV8')

MOUSE_IGLJ = c('IGLJ1','IGLJ3P','IGLJ5','IGLJ3','IGLJ2','IGLJ4')

MOUSE_TRDD = c('TRDD2', 'TRDD1')

MOUSE_TRDJ = c('TRDJ1', 'TRDJ2')

MOUSE_TRDV = c('TRAV15D-1/DV6D-1', 'TRAV16D/DV11', 'TRAV4-4/DV10', 'TRAV21/DV12', 'TRDV1', 'TRAV15D-2/DV6D-2', 
               'TRAV13-4/DV7', 'TRAV6-7/DV9', 'TRDV2-1', 'TRAV15-2/DV6-2', 'TRAV14D-3/DV8', 
               'TRDV5', 'TRDV2-2', 'TRDV4', 'TRAV15-1/DV6-1')

MOUSE_TRGJ = c('TRGJ4', 'TRGJ2', 'TRGJ1', 'TRGJ3')

MOUSE_TRGV = c('TRGV4', 'TRGV6', 'TRGV5', 'TRGV3', 'TRGV7', 'TRGV1', 
               'TRGV2')

MACMUL_TRBV <- c('TRBD1', 'TRBD2', 'TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 
                 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-2P', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5',
                 'TRBJ2-6', 'TRBJ2-7', 'TRBV1-1', 'TRBV10-1', 'TRBV10-2', 'TRBV10-3', 
                 'TRBV11-1', 'TRBV11-2', 'TRBV11-3', 'TRBV12-1', 'TRBV12-2', 'TRBV12-3', 
                 'TRBV12-4', 'TRBV13', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV18', 'TRBV19', 
                 'TRBV2-1', 'TRBV2-2', 'TRBV2-3', 'TRBV20-1', 'TRBV21-1', 'TRBV22-1', 
                 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBV27', 'TRBV28', 'TRBV29-1', 
                 'TRBV3-1', 'TRBV3-2', 'TRBV3-3', 'TRBV3-4', 'TRBV30', 'TRBV4-1', 'TRBV4-2', 
                 'TRBV4-3', 'TRBV5-1', 'TRBV5-10', 'TRBV5-2', 'TRBV5-3', 'TRBV5-4', 
                 'TRBV5-5', 'TRBV5-6', 'TRBV5-7', 'TRBV5-8', 'TRBV5-9', 'TRBV6-1', 'TRBV6-2',
                 'TRBV6-3', 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-7', 'TRBV7-10', 'TRBV7-2',
                 'TRBV7-3', 'TRBV7-4', 'TRBV7-5', 'TRBV7-6', 'TRBV7-7', 'TRBV7-9', 'TRBV9')

MACMUL_TRBJ <- c('TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ2-1', 'TRBJ2-2',
                 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7')


#' Segment data.
#' 
#' @docType data
#' 
#' @aliases genesegments
#' 
#' @name segments.list
#' 
#' @description
#' \code{segments} is a list with 5 data frames with data of human alpha-beta chain segments.
#' Elements names as "TRAV", "TRAJ", "TRBV", "TRVJ", "TRVD". Each data frame consists of 5 columns:
#' 
#' - V.allelles / J.allelles / D.allelles - character column with names of V/D/J-segments.
#' 
#' - CDR3.position - position in the full nucleotide segment sequence where CDR3 starts.
#' 
#' - Full.nucleotide.sequence - character column with segment CDR1-2-3 sequence.
#' 
#' - Nucleotide.sequence - character column with segment CDR3 sequences.
#' 
#' - Nucleotide.sequence.P - character column with segment CDR3 sequences with P-insertions.
#' 
#' @format
#' \code{genesegments} is a list with data frames.
#' 
#' @examples
#' \dontrun{
#' data(genesegments)
#' genesegments$Nucleotide.sequence[segments$TRBV[,1] == "TRBV10-1"]
#' }
NULL


#' List with assembling probabilities of beta chain TCRs.
#' 
#' @docType data
#' 
#' @name beta.prob
#' 
#' @aliases beta.prob
#' 
#' @description
#' \code{beta.prob} is a list with probabilities of TCR assembling taken from
#' \code{Murugan et al. Statistical inference of the generation probability
#' of T-cell receptors from sequence repertoires}. It's a list with the following elements:
#' 
#' - P.V - matrix with 1 column and row names stands for V-beta segments. Each element is
#' a probability of choosing corresponding V-beta segment.
#' 
#' - P.del.D1 - matrix 17x17 with probabilities of choosing D5-D3 deletions for TRBD1.
#' 
#' - P.del.D1 - matrix 17x17 with probabilities of choosing D5-D3 deletions for TRBD2.
#' 
#' - P.ins.len - matrix with first columns stands for number of insertions, second and third columns filled
#' with probability values of choosing corresponding number of insertions in VD- and DJ-junctions
#' correspondingly.
#' 
#' - P.ins.nucl - data frame with probability of choosing a nucleotide in the insertion on junctions with known
#' previous nucleotide. First column with names of nucleotides, 2-5 columns are probabilities of choosing
#' adjacent nucleotide in VD-junction, 6-9 columns are probabilities of choosing adjacent nucleotide in DJ-junction.
#' 
#' - P.del.J - matrix with the first column "P.del.V" with number of deletions, and other columns with
#' names for V-segments and with probabilities of corresponding deletions.
#' 
#' - P.del.J - matrix with the first column "P.del.J" with number of deletions, and other columns with
#' names for J-segments and with probabilities of corresponding deletions.
#' 
#' - P.J.D - matrix with two columns ("TRBD1" and "TRBD2") and 13 rows with row names stands for 
#' J-beta segments. Each element is a mutual probability of choosing J-D segments.
#' 
#' @format
#' \code{beta.prob} is a list of matrices and data frames.
#' 
#' @examples
#' \dontrun{
#' # Generate 10 kmers with adjacent nucleotide probability.
#' generate.kmers.prob(rep.int(10, 10), .probs=beta.prob$P.ins.nucl[,c(1, 2:5)])
#' }
NULL


#' Twins alpha-beta chain data
#' 
#' @docType data
#' 
#' @name twinsdata
#' 
#' @aliases twa twb
#' 
#' @description
#' \code{twa.rda}, \code{twb.rda} - data frames with downsampled to the 10000 most 
#' abundant clonesets and 4 samples data of twins data (alpha and beta chains).
#' Link: http://labcfg.ibch.ru/tcr.html
#' 
#' @format
#' \code{twa} and \code{twb} are lists of 4 data frames with 10000 row in each.
#' 
#' @examples
#' \dontrun{
#' data(twa)
#' data(twb)
#' }
NULL