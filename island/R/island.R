#' island: Stochastic Island Biogeography Theory Made Easy
#'
#' Tools to develop stochastic models based on the Theory of Island Biogeography
#' (TIB) of MacArthur and Wilson (1967) and extensions. The package allows the
#' calculation of colonization and extinction rates (including environmental
#' variables) given presence-absence data, the simulation of community assembly
#' and model selection.
#'
#' @details In the simplest stochastic model of Island Biogeography, there is a
#'   pool of species that potentially can colonize a system of islands. When we
#'   sample an island in time, we obtain a time-series of presence-absence
#'   vectors for the different species of the pool, which allows us to estimate
#'   colonization (\eqn{c}) and extinction (\eqn{e}) rates under perfect
#'   detectability. These are actual rates (in time^{-1} units). \cr The
#'   simplest stochastic model of island biogeography assumes a single
#'   colonization-extinction pair for the whole community. This model implictly
#'   assumes: first, neutrality of the species in the community, that is, all
#'   species in the community share the same values for those rates; and second,
#'   all species colonize and become extinct indepently from each other. The
#'   "species neutrality assumption" can be relaxed easily, for example,
#'   calculating different rates for different groups or on a per-species basis.
#'   In addition, we can make these rates depend on environmental variables
#'   measured at the same time that we took our samples. For more information of
#'   the basic model, please see the references.
#'
#' @section Data entry: The data should be organized in dataframes with
#'   consecutive presence-absence data of each sample ordered cronologically,
#'   being the data associated with a single species in a row. Additional
#'   columns can contain the filiations of every species to a group, i. e. a
#'   phylogenetic group or a guild.
#'
#'
#'
#' @references Alonso, D., Pinyol-Gallemi, A., Alcoverro T. and Arthur, R..
#'   (2015) Fish community reassembly after a coral mass mortality: higher
#'   trophic groups are subject to increased rates of extinction. \emph{Ecology
#'   Letters}, \bold{18}, 451--461. \cr \cr Simberloff, D. S., and Wilson, E.
#'   O.. (1969). Experimental Zoogeography of Islands: The Colonization of Empty
#'   Islands. \emph{Ecology}, \bold{50(2)}, 278--296.
#'   \url{http://doi.org/10.2307/1934856} \cr \cr Simberloff, D. S.. (1969).
#'   Experimental Zoogeography of Islands: A Model for Insular Colonization.
#'   \emph{Ecology}, \bold{50(2)}, 296--314.
#'   \url{http://doi.org/10.2307/1934857}
#'
#' @docType package
#' @name island
NULL
