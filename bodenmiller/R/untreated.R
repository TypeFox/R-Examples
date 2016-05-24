#' Untreated Samples with 4 Different Treatments (Phenotypic)
#'
#' a \code{matrix} containing the signal intensity for phenotypic markers measured using CyTOF
#' on different cell types, after stimulation with BCR/FcR-XL, PMA/Ionomycin
#' and vanadate. See
#' \href{http://www.nature.com/nbt/journal/v30/n9/full/nbt.2317.html}{Bodenmiller et al 2012}
#' for details.
#'
#' @format a matrix with 15792 rows and 9 variables:
#' \itemize{
#'    \item{ CD20 }
#'    \item{ IgM }
#'    \item{ CD4 }
#'    \item{ CD33 }
#'    \item{ HLA-DR }
#'    \item{ CD14 }
#'    \item{ CD7 }
#'    \item{ CD3 }
#'    \item{ CD123 }
#' }
#' @source \url{http://reports.cytobank.org/105/v2}
"untreatedPhenoMat"

#' Untreated Samples with 4 Different Treatments (Functional)
#'
#' a \code{matrix} containing the signal intensity for functional markers measured using CyTOF
#' on different cell types, after stimulation with BCR/FcR-XL, PMA/Ionomycin
#' and vanadate. See
#' \href{http://www.nature.com/nbt/journal/v30/n9/full/nbt.2317.html}{Bodenmiller et al 2012}
#' for details.
#'
#' @format a matrix with 15792 rows and 14 variables:
#' \itemize{
#'    \item{ pStat1 }
#'    \item{ pSlp76 }
#'    \item{ pBtk }
#'    \item{ pPlcg2 }
#'    \item{ pErk }
#'    \item{ pLat }
#'    \item{ pS6 }
#'    \item{ pNFkB }
#'    \item{ pp38 }
#'    \item{ pStat5 }
#'    \item{ pAkt }
#'    \item{ pSHP2 }
#'    \item{ pZap70 }
#'    \item{ pStat3 }
#' }
#' @source \url{http://reports.cytobank.org/105/v2}
"untreatedFuncMat"

#' Untreated Samples with 4 Different Treatments (Annotation)
#'
#' A \code{data.frame} containing the source file, cell type and stimulation for every cell in
#' \code{\link{refPhenoMat}} and \code{\link{refFuncMat}}.
#' Stimulations include BCR/FcR-XL, PMA/Ionomycin and vanadate.See
#' \href{http://www.nature.com/nbt/journal/v30/n9/full/nbt.2317.html}{Bodenmiller et al 2012}
#' for details.
#'
#' @format a data.frame with 15792 rows and 2 variables:
#' \itemize{
#'    \item{Source}{the name of the source (fcs) file}
#'    \item{Treatment}{the treatment for the corresponding row in \code{\link{refPhenoMat}}
#'      or \code{\link{refFuncMat}}}
#'    \item{Cells}{the cell type for the corresponding row in \code{\link{refPhenoMat}}
#'      or \code{\link{refFuncMat}}}
#' }
#' @source \url{http://reports.cytobank.org/105/v2}
"untreatedAnnots"
