#' @importFrom R6 R6Class
#' @importClassesFrom rehh haplohh
#' @importFrom rehh calc_ehh calc_ehhs ihh2ihs
#' @importFrom methods new checkAtAssignment
stat_ihh_class <- R6Class("stat_ihh", inherit = sumstat_class,
  private = list(
    req_segsites = TRUE,
    population = NULL,
    max_snps = Inf,
    use_ihs = FALSE,
    empty_matrix = data.frame(CHR = numeric(),
                              POSITION = numeric(),
                              FREQ_a = numeric(),
                              IHHa = numeric(),
                              HHd = numeric(),
                              IES = numeric())
  ),
  public = list(
    initialize = function(name, population, max_snps,
                          calc_ihs, transformation) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      assert_that(is.numeric(max_snps))
      assert_that(length(max_snps) == 1)
      assert_that(is.logical(calc_ihs))
      assert_that(length(calc_ihs) == 1)
      private$population <- population
      private$max_snps <- max_snps
      private$use_ihs <- calc_ihs
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      assert_that(is.list(seg_sites))
      assert_that(is.model(model))
      ind <- get_population_individuals(model, private$population)
      ihh <- do.call(rbind, lapply(seq(along = seg_sites), function(locus) {
        assert_that(is_segsites(seg_sites[[locus]]))
        if (ncol(seg_sites[[locus]]) == 0) return(private$empty_matrix)

        rehh_data <- self$create_rehh_data(seg_sites[[locus]],
                                           ind, model, locus)

        if (rehh_data@nsnp == 0) return(private$empty_matrix)

        data.frame(rehh::scan_hh(rehh_data))
      }))

      if (private$use_ihs) {
        n_snps <- nrow(ihh)
        if (n_snps < 5) {
          # Standardization of iHS requires a few SNPs
          return(list(ihh = ihh,
                      iHS = data.frame(ihh[ , 1:2], iHS = rep(NA, n_snps))))
        } else {
          if (n_snps < 50) freqbin <- 0.90
          else if (n_snps < 100) freqbin <- 0.45
          else if (n_snps < 200) freqbin <- 0.225
          else if (n_snps < 400) freqbin <- 0.1
          else freqbin <- 0.05
          ihs <- suppressWarnings(
            data.frame(ihh2ihs(ihh, freqbin)$res.ihs[ , -4])
          )
          return(list(ihh = ihh, iHS = ihs))
        }
      }

      ihh
    },
    create_rehh_data = function(seg_sites, ind, model, chr_name = 1) {
      assert_that(is_segsites(seg_sites))
      assert_that(is.model(model))
      seg_sites <- seg_sites[ind, ]
      pos <- get_snp_positions(list(seg_sites), model, relative = FALSE)
      snp_mask <- self$create_snp_mask(seg_sites)

      rehh_data <- new("haplohh")
      rehh_data@haplo <- as.matrix(seg_sites[ind, snp_mask]) + 1
      rehh_data@position <- pos[[1]][snp_mask]
      rehh_data@snp.name <- as.character(seq(along = rehh_data@position))
      rehh_data@chr.name <- chr_name
      rehh_data@nhap <- length(ind)
      rehh_data@nsnp <- length(rehh_data@position)
      rehh_data
    },
    create_snp_mask = function(seg_sites) {
      if (ncol(seg_sites) < private$max_snps) return(rep(TRUE, ncol(seg_sites)))
      sort.int(sample.int(ncol(seg_sites), private$max_snps, replace = FALSE))
    }
  )
)


#' Summary Statistic: Integrated Extended Haplotype Homozygosity
#'
#' This summary statistic calculates a number of values based on
#' extended haplotype homozygosity (EHH), including iHH, iES
#' and optionally iHS.
#' Coala relies on \code{\link[rehh]{scan_hh}} from package \pkg{rehh} to
#' calculate this statistic. Please refer
#' to their documentation for detailed information on the implementation.
#' Please cite the corresponding publication (see below) if you use the
#' statistic for a publication.
#'
#' @section References:
#' \itemize{
#'  \item{Mathieu Gautier and Renaud Vitalis, rehh: an R package to detect
#'       footprints of selection in genome-wide SNP data from
#'       haplotype structure. Bioinformatics (2012) 28 (8): 1176-1177
#'       first published online March 7, 2012 doi:10.1093/bioinformatics/bts115}
#'  \item{Voight et al., A map of recent positive selection in the human
#'                        genome. PLoS Biol, 4(3):e72, Mar 2006.}
#'  }
#'
#' @inheritParams sumstat_four_gamete
#' @param max_snps The maximal number of SNPs per locus that are used for the
#'   calculation. If a locus has more SNPs, only a random subset of them will
#'   be used to increase performance. Set to \code{Inf} to use all SNPs.
#' @param calc_ihs If set to \code{TRUE}, additionally standardized iHS is
#'   calculated.
#' @return If \code{calc_ihs = FALSE}, a data.frame with values for
#'   iHH and iES is returned. Otherwise, a list of two data frames are
#'   returned, one for IHH and IES values and the other one for IHS values.
#'
#'   In all `data.frames` rows are SNPs and the colums present the following
#'   values for each SNP:
#'   \itemize{
#'    \item{CHR: The SNP's locus}
#'    \item{Positions: The SNP's absolute position on its locus}
#'    \item{FREQ_a: The SNP's absolute position on its locus}
#'    \item{IHHa: integrated EHH for the ancestral allele}
#'    \item{IHHd: integrated EHH for the derived allele}
#'    \item{IES: integrated EHHS}
#'    \item{iHS: iHS, normalized over all loci.}
#'   }
#' @export
#' @template summary_statistics
#' @examples
#'   model <- coal_model(20, 1, 1000) +
#'     feat_mutation(1000) +
#'     sumstat_ihh()
#' \dontrun{
#'     stat <- simulate(model)
#'     print(stat$ihh)}
#' @author Paul Staab
sumstat_ihh <- function(name = "ihh", population = 1,
                        max_snps = 1000, calc_ihs = FALSE,
                        transformation = identity) {
  stat_ihh_class$new(name, population, max_snps, calc_ihs, transformation)
}
