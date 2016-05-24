#' MOdelling Tools for Reproduction and Survival Data in Ecotoxicology
#'
#' Provides tools for the analysis of survival/reproduction
#' bioassay data in quantitative environmental risk assessment. It can be
#' used to explore/visualize experimental data, and to perform an estimation
#' of \eqn{LC_{x}} (Lethal Concentration for x\% of individuals) or
#' \eqn{EC_{x}} (x\% Effect Concentration) values by fitting exposure-response
#' curves. The \eqn{LC_{x}}/\eqn{EC_{x}} and parameters of the curve are
#' provided along with an indication of the uncertainty of the estimation.
#'
#' Estimation procedures in MORSE can be used without a deep knowledge of
#' their underlying probabilistic model or inference methods. Rather, they
#' were designed to behave as well as possible without requiring a user to
#' provide values for some obscure parameters. That said, MORSE models can also
#' be used as a first step to tailor new models for more specific situations.
#'
#' The package currently handles survival and reproduction data. Functions
#' dedicated to survival (resp. reproduction) analysis start with a
#' \code{surv} (resp. \code{repro}) prefix. MORSE provides a similar
#' workflow in both cases:
#' \enumerate{
#' \item create and validate a dataset
#' \item explore a dataset
#' \item plot a dataset
#' \item fit a model on a dataset and output the expected estimates
#' }
#' Those steps are presented in more details in the "Tutorial" vignette, while
#' a more formal description of the estimation procedures are provided in the
#' vignette called "Models in MORSE package". Please refer to these documents
#' for further introduction to the use of MORSE.
#'
#' This reference manual is a detailed description of the functions exposed in
#' the package.
#'
#' \strong{Getting started} The package uses the \code{rjags} package
#' (Plummer, 2013), an R interface to the JAGS library for Bayesian model
#' estimation. Note that the \code{rjags} package does not include a copy
#' of the JAGS library: you need to install it separately. For instructions
#' on downloading JAGS, see the home page at
#' \url{http://mcmc-jags.sourceforge.net}. Once done, simply follow the steps
#' described in the tutorial vignette.
#'
#' \tabular{ll}{ Package: \tab morse\cr Type: \tab Package\cr Version: \tab
#' 2.1.1\cr Date: \tab 2015-12-21\cr License: \tab GPL (>=2)\cr }
#'
#' @name morse-package
#' @aliases morse-package morse
#' @docType package
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>, Sandrine Charles
#' <sandrine.charles@@univ-lyon1.fr>, Wandrille Duchemin
#' <wandrille.duchemin@@insa-lyon.fr>, Christelle Lopes
#' <christelle.lopes@@univ-lyon1.fr>, Guillaume Kon-Kam-king
#' <guillaume.kon-kam-king@@univ-lyon1.fr>, Philippe Veber
#' <philippe.veber@@univ-lyon1.fr>
#'
#' Maintainer: Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>
#' @seealso \code{\link[rjags]{rjags}},
#' \code{\link[ggplot2]{ggplot}}
#' @references Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S.
#' (2014) \emph{Statistical handling of reproduction data for exposure-response
#' modelling}.
#' \url{http://pubs.acs.org/doi/abs/10.1021/es502009r?journalCode=esthag}.
#'
#' Plummer, M. (2013) \emph{JAGS Version 4.0.0 user manual}.
#' \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download}.
#' @keywords package
#'
NULL

#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to
#' cadmium during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to five concentrations
#' of one metal contaminant (cadmium) during 21 days. Five concentrations were
#' tested, with four replicates per concentration. Each replicate contained 10
#' organisms. Reproduction and survival were monitored at 10 time points.
#'
#'
#' @name cadmium1
#' @docType data
#' @usage data(cadmium1)
#' @format A data frame with 200 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{integer} with
#' the replicate code (\code{1} to \code{4}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Billoir, E., Delhaye, H., Forfait, C., Clement, B.,
#' Triffault-Bouchet, G., Charles, S. and Delignette-Muller, M.L. (2012)
#' Comparison of bioassays with different exposure time patterns: The added
#' value of dynamic modelling in predictive ecotoxicology, \emph{Ecotoxicology
#' and Environmental Safety}, 75, 80-86.
#' @keywords datasets
NULL





#' Reproduction and survival datasets for snails exposed to cadmium during 56
#' days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' snails exposed to six concentrations of one metal contaminant (cadmium)
#' during 56 days. Six concentrations were tested, with six replicates per
#' concentration. Each replicate contained five organisms. Reproduction and
#' survival were monitored at 17 time points.
#'
#'
#' @name cadmium2
#' @docType data
#' @usage data(cadmium2)
#' @format A data frame with 612 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{A} to \code{F}).} \item{\code{conc}}{A vector of
#' class \code{integer} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Ducrot, V., Askem, C., Azam, D., Brettschneider, D., Brown,
#' R., Charles, S., Coke, M., Collinet, M., Delignette-Muller, M.L.,
#' Forfait-Dubuc, C., Holbech, H., Hutchinson, T., Jach, A., Kinnberg, K.L.,
#' Lacoste, C., Le Page, G., Matthiessen, P., Oehlmann, J., Rice, L.,
#' Roberts, E., Ruppert, K., Davis, J.E., Veauvy, C., Weltje, L., Wortham, R.
#' and Lagadic, L. (2014)
#' Development and validation of an OECD reproductive toxicity test guideline with
#' the pond snail Lymnaea stagnalis (Mollusca, Gastropoda),
#' \emph{Regulatory Toxicology and Pharmacology}, 70(3), 605-14.
#' @keywords datasets
NULL





#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to
#' chlordan during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to six concentrations
#' of one organochlorine insecticide during 21 days. Six concentrations were
#' tested, with 10 replicates per concentration. Each replicate contained one
#' organism. Reproduction and survival were monitored at 22 time points.
#'
#'
#' @name chlordan
#' @docType data
#' @usage data(chlordan)
#' @format A data frame with 1320 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{integer} with
#' the replicate code (\code{1} to \code{10}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with the chlordan concentrations in \eqn{\mu
#' g.L^{-1}}.} \item{\code{time}}{A vector of class \code{integer} with the
#' time points (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Manar, R., Bessi, H. and Vasseur, P. (2009) Reproductive effects
#' and bioaccumulation of chlordan in Daphnia magna, \emph{Environmental
#' Toxicology and Chemistry}, 28, 2150-2159.
#' @keywords datasets
NULL





#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to
#' copper during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to five concentrations
#' of one metal contaminant (copper) during 21 days. Five concentrations were
#' tested, with three replicates per concentration. Each replicate contained 20
#' organisms. Reproduction and survival were monitored at 16 time points.
#'
#'
#' @name copper
#' @docType data
#' @usage data(copper)
#' @format A data frame with 240 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{A} to \code{C}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with the copper concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Billoir, E., Delignette-Muller, M.L., Pery, A.R.R. and
#' Charles, S. (2008) A Bayesian Approach to Analyzing Ecotoxicological Data,
#' \emph{Environmental Science & Technology}, 42 (23), 8978-8984.
#' @keywords datasets
NULL









#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to zinc
#' during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to four concentrations
#' of one metal contaminant (zinc) during 21 days. Four concentrations were
#' tested with three replicates per concentration. Each replicate contained 20
#' organisms. Reproduction and survival were monitored at 15 time points.
#'
#'
#' @name zinc
#' @docType data
#' @usage data(zinc)
#' @format A data frame with 180 observations on the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{A} to \code{C}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with zinc concentrations in \eqn{mg.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Billoir, E.,Delignette-Muller, M.L., Pery, A.R.R. and
#' Charles S. (2008) A Bayesian Approach to Analyzing Ecotoxicological Data,
#' \emph{Environmental Science & Technology}, 42 (23), 8978-8984.
#' @keywords datasets
NULL
