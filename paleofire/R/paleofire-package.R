#' paleofire: A package for the Global Charcoal Database
#' 
#' The paleofire package provides tools to extract and analyse charcoal
#' sedimentary data stored in the Global Charcoal Database. Main
#' functionalities includes data extraction and sites selection, transformation
#' and homogenization of the charcoal records as well as regional to global
#' compositing.
#' 
#' \tabular{ll}{ Package: \tab paleofire\cr Type: \tab Package\cr Version: \tab
#' 1.1.8\cr Date: \tab 2015-12-21\cr License: \tab GPL (>=2)\cr }
#' 
#' @name paleofire-package
#' @aliases paleofire-package paleofire
#' @docType package
#' @section Author(s): Global Paleofire Working Group <paleofire at univ-fcomte.fr>
#' @section Maintainer: Olivier Blarquez <blarquez at gmail.com>
#' @seealso \url{http://gpwg.paleofire.org}
#' @references Blarquez, O., Vannière, B., Marlon, J. R., Daniau, A. L., Power,
#' M. J., Brewer, S., & Bartlein, P. J. (2014). paleofire: an R package to
#' analyse sedimentary charcoal records from the Global Charcoal Database to
#' reconstruct past biomass burning. Computers & Geosciences, 72, 255-261.
#' @keywords fire, global, charcoal, sediments, paleo
#' @examples
#' 
#' \dontrun{
#' ## Interactive sites selection:
#' # ID=pfInteractive()
#' 
#' ## Site selection using criterions
#' # Boreal Eastern North American sites with at least one
#' # dating point each 2500 year
#' 
#' ID=pfSiteSel(lat>50, lat<70, long>-90, long<(-50), date_int<=2500, l12==1)
#' plot(ID,zoom="world")
#' 
#' ## Modify plot
#' plot(ID,zoom="sites")
#' 
#' ## Simple test for transforming data
#' # Select site 1 (Cygnet Lake)
#' 
#' ID1=pfSiteSel(id_site==1)
#' plot(ID1)
#' 
#' # Transformation of data
#' TR=pfTransform(ID1,method=c("MinMax", "Box-Cox", "Z-Score"))
#' 
#' # Plot Transformed and raw data
#' # First retrieve raw data for Cygnet using pfExtract
#' 
#' RAW=pfExtract(ID=1)
#' 
#' dev.off()
#' par(mfrow=(c(2,1)))
#' 
#' plot(RAW[,3],RAW[,4],type="l")
#' plot(TR$Age,TR$TransData,type="l")
#' 
#' ## Transforming and Compositing
#' ## Example 1: Usage as in Power et al. 2008
#' ## Data transformation
#' 
#' TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))
#' 
#' ## Diagnostic pdf file with transformed series:
#' # pfDiagnostic(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000), 
#' # FileName = "Diagnostic.pdf")
#' 
#' ## Compositing: basic binning procedure
#' COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,8000,500))
#' plot(COMP)
#' 
#' ## The result matrix can be saved
#' # write.csv(COMP$Result,file="temp.csv")
#' 
#' ## Compositing: Using the locfit package equivalent procedure to Daniau et al. 2012
#' 
#' COMP2=pfCompositeLF(TR1, tarAge=seq(-50,8000,20), binhw=20, hw=500,nboot=100)
#' plot(COMP2)
#' 
#' ## And save
#' write.csv(COMP2$Result,file="temp2.csv")
#' }
#' 
NULL



#' Internal paleofire functions
#' 
#' Internal paleofire functions and functions waiting for man.
#' 
#' 
#' @aliases bootCI haverdist dem loessGCV bestLoess countries influx
#' @name paleofire-internal
#' 
NULL