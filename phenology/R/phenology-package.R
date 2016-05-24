#' Fit a parametric function that describes phenology
#'
#' \tabular{ll}{
#'  Package: \tab phenology\cr
#'  Type: \tab Package\cr
#'  Version: \tab 5.1 build 354\cr
#'  Date: \tab 2016-05-02\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title Tools to Manage a Parametric Function that Describes Phenology
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType package
#' @name phenology-package
#' @description Functions used to fit and test the phenology of species based on counts.\cr
#' Note that only the most significant changes are reported in the NEWS.\cr
#' To do:\cr
#' * There are problems with SE for fitRMU().\cr
#' * Auto-scaling for optim during fitRMU search.\cr
#' * I must adapt TCF (total clutch frequency) fit from OCF-ECF (observed clutch frequency-estimated cluth frequency) table based on:\cr
#'  Briane, J.-P., Rivalan, P., Girondot, M., 2007. The inverse problem applied to the Observed Clutch Frequency of Leatherbacks from Yalimapo beach, French Guiana. Chelonian Conservation and Biology 6, 63-69.\cr
#' Until now it is an Excel spreadsheet.\cr
#' * Fit tag-loss rate based on:\cr
#' Rivalan, P., Godfrey, M.H., Prévot-Julliard, A.-C., Girondot, M., 2005. Maximum likelihood estimates of tag loss in leatherback sea turtles. Journal of Wildlife Management 69, 540-548.\cr
#' Until now it is a RealBasic software.\cr
#' The lastest version of this package can always been installed using:\cr
#' install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/phenology.tar.gz", repos=NULL, type="source")
#' @references Girondot, M. 2010. Estimating density of animals during 
#'             migratory waves: application to marine turtles at nesting site. 
#'             Endangered Species Research, 12, 85-105.
#' @references Girondot M. and Rizzo A. 2015. Bayesian framework to integrate 
#'             traditional ecological knowledge into ecological modeling: A case 
#'             study. Journal of Ethnobiology, 35, 339-355.
#' @references Girondot, M. 2010. Editorial: The zero counts. Marine 
#'             Turtle Newsletter, 129, 5-6.
#' @keywords Seasonality Phenology Ecology
#' @seealso Girondot, M., Rivalan, P., Wongsopawiro, R., Briane, J.-P., Hulin, V.,
#'          Caut, S., Guirlet, E. & Godfrey, M. H. 2006. Phenology of marine turtle 
#'          nesting revealed by a statistical model of the nesting season. BMC Ecology, 
#'          6, 11.
#' @seealso Delcroix, E., Bédel, S., Santelli, G., Girondot, M., 2013. Monitoring 
#'          design for quantification of marine turtle nesting with limited human 
#'          effort: a test case in the Guadeloupe Archipelago. Oryx 48, 95-105.
#' @seealso Weber, S.B., Weber, N., Ellick, J., Avery, A., Frauenstein, R., 
#'          Godley, B.J., Sim, J., Williams, N., Broderick, A.C., 2014. Recovery 
#'          of the South Atlantic’s largest green turtle nesting population. 
#'          Biodiversity and Conservation 23, 3005-3018.
#' @examples
#' \dontrun{
#' library(phenology)
#' # Get the lastest version at:
#' # install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/phenology.tar.gz", 
#'      repos=NULL, type="source")
#' # Read a file with data
#' data(Gratiot)
#' # Generate a formatted list nammed data_Gratiot 
#' data_Gratiot<-add_phenology(Gratiot, name="Complete", 
#' 		reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg<-par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot<-fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' data(result_Gratiot)
#' # Plot the phenology and get some stats
#' output<-plot(result_Gratiot)
#' }

NULL

