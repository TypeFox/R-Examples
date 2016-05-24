#' @docType package
#'
#' @name seroincidence
#'
#' @title
#' Estimating Infection Rates from Serological Data
#'
#' @description
#' Translates antibody levels measured in a (cross-sectional) population sample into an
#' estimate of the frequency with which seroconversions (infections) occur in the sampled population.
#'
#' @details
#' For detailed documentation type the following in the R console:\cr
#'
#' \code{vignette("installation", package = "seroincidence")}\cr
#' \code{vignette("tutorial", package = "seroincidence")}\cr
#' \code{vignette("methodology", package = "seroincidence")}
#'
#' @author
#' Author: Peter Teunis \email{<peter.teunis@@rivm.nl>}\cr
#' Contributor: Daniel Lewandowski \email{<daniel@@nextpagesoft.net>}\cr
#' Maintainer: Chantal Quinten \email{<seroincidence@@ecdc.europa.eu>}
#'
#' @references
#'
#' \strong{\emph{Methods for estimating seroincidence}}
#'
#' \itemize{
#' \item Teunis, P. F., van Eijkeren, J. C., Ang, C. W., van Duynhoven, Y. T., Simonsen, J. B., Strid, M. A., van Pelt, W.\cr
#' "Biomarker dynamics: estimating infection rates from serological data"\cr
#' Statistics in Medicine 31, no. 20 (September 9, 2012): 2240--48. doi:10.1002/sim.5322.
#'
#' \item Simonsen, J., Molbak, K., Falkenhorst, G., Krogfelt, K. A., Linneberg, A., Teunis, P. F.\cr
#' "Estimation of incidences of infectious diseases based on antibody measurements"\cr
#' Statistics in Medicine 28, no. 14 (June 30, 2009): 1882--95. doi:10.1002/sim.3592.
#' }
#'
#' \strong{\emph{Applications}}
#'
#' \itemize{
#' \item Kretzschmar, M., Teunis, P. F., Pebody, R. G.\cr
#' "Incidence and reproduction numbers of pertussis: estimates from serological and social contact data in five European countries"\cr
#' PLoS Medicine 7, no. 6 (June 1, 2010):e1000291. doi:10.1371/journal.pmed.1000291.
#'
#' \item Simonsen, J., Strid, M. A., Molbak, K., Krogfelt, K. A., Linneberg, A., Teunis, P.\cr
#' "Sero-epidemiology as a tool to study the incidence of Salmonella infections in humans"\cr
#' Epidemiology and Infection 136, no. 7 (July 1, 2008): 895--902. doi:10.1017/S0950268807009314
#'
#' \item Simonsen, J., Teunis, P. F., van Pelt, W., van Duynhoven, Y., Krogfelt, K. A., Sadkowska-Todys, M., Molbak K.\cr
#' "Usefulness of seroconversion rates for comparing infection pressures between countries"\cr
#' Epidemiology and Infection, April 12, 2010, 1--8. doi:10.1017/S0950268810000750.
#'
#' \item Falkenhorst, G., Simonsen, J., Ceper, T. H., van Pelt, W., de Valk, H., Sadkowska-Todys, M., Zota, L., Kuusi, M., Jernberg, C., Rota, M. C., van Duynhoven, Y. T., Teunis, P. F., Krogfelt, K. A., Molbak, K.\cr
#' "Serological cross-sectional studies on salmonella incidence in eight European countries: no correlation with incidence of reported cases"\cr
#' BMC Public Health 12, no. 1 (July 15, 2012): 523--23. doi:10.1186/1471-2458-12-523.
#'
#' \item Teunis, P. F., Falkenhorst, G., Ang, C. W., Strid, M. A., De Valk, H., Sadkowska-Todys, M., Zota, L., Kuusi, M., Rota, M. C., Simonsen, J. B., Molbak, K., Van Duynhoven, Y. T., van Pelt, W.\cr
#' "Campylobacter seroconversion rates in selected countries in the European Union"\cr
#' Epidemiology and Infection 141, no. 10 (2013): 2051--57. doi:10.1017/S0950268812002774.
#'
#' \item de Melker, H. E., Versteegh, F. G., Schellekens, J. F., Teunis, P. F., Kretzschmar, M.\cr
#' "The incidence of Bordetella pertussis infections estimated in the population from a combination of serological surveys"\cr
#' The Journal of Infection 53, no. 2 (August 1, 2006): 106--13. doi:10.1016/j.jinf.2005.10.020
#' }
#'
#' \strong{\emph{Quantification of seroresponse}}
#'
#' \itemize{
#' \item de Graaf, W. F., Kretzschmar, M. E., Teunis, P. F., Diekmann, O.\cr
#' "A two-phase within-host model for immune response and its application to serological profiles of pertussis"\cr
#' Epidemics 9 (2014):1--7. doi:10.1016/j.epidem.2014.08.002.
#'
#' \item Berbers, G. A., van de Wetering, M. S., van Gageldonk, P. G., Schellekens, J. F., Versteegh, F. G., Teunis, P.F.\cr
#' "A novel method for evaluating natural and vaccine induced serological responses to Bordetella pertussis antigens"\cr
#' Vaccine 31, no. 36 (August 12, 2013): 3732--38. doi:10.1016/j.vaccine.2013.05.073.
#'
#' \item Versteegh, F. G., Mertens, P. L., de Melker, H. E., Roord, J. J., Schellekens, J. F., Teunis, P. F.\cr
#' "Age-specific long-term course of IgG antibodies to pertussis toxin after symptomatic infection with Bordetella pertussis"\cr
#' Epidemiology and Infection 133, no. 4 (August 1, 2005): 737--48.
#'
#' \item Teunis, P. F., van der Heijden, O. G., de Melker, H. E., Schellekens, J. F., Versteegh, F. G., Kretzschmar, M. E.
#' "Kinetics of the IgG antibody response to pertussis toxin after infection with B. pertussis"\cr
#' Epidemiology and Infection 129, no. 3 (December 10, 2002):479. doi:10.1017/S0950268802007896.
#' }
#'
NULL
