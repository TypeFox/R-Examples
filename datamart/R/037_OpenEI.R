#' Interface to the OpenEI platform
#'
#' OpenEI is growing into a global leader in the energy data realm - specifically 
#' analyses on renewable energy and energy efficiency. The platform is a wiki, 
#' similar to Wikipedia's Wiki, and offers an SPARQL endpoint.
#'
#' The \code{openei} function provides an datamart with predefined queries on
#' this SPARQL endpoint.
#'
#' @seealso \code{\link{openei}}, \code{\link{xsparql}}
#'
#' @references
#' \href{http://chrisdavis.weblog.tudelft.nl/2010/06/12/data-mining-the-us-department-of-energy/}{Blogpost by Chris Davis},
#' \href{http://en.openei.org/wiki/OpenEI:About}{OpenEI website}
#' @export
openei <- function() datamart(
    openei.powerplants(resource="powerplants")
)

# internal function
openei.powerplants <- function(resource="powerplants") xsparql(
    resource=resource,
    url="http://en.openei.org/sparql",
    nspace=c(
        "article", "<http://enipedia.tudelft.nl/wiki/>",
        "a", "<http://enipedia.tudelft.nl/wiki/>",
        "property", "<http://enipedia.tudelft.nl/wiki/Property:>",
        "prop", "<http://enipedia.tudelft.nl/wiki/Property:>",
        "category", "<http://enipedia.tudelft.nl/wiki/Category:>",
        "cat", "<http://enipedia.tudelft.nl/wiki/Category:>",
        "rdfs", "<http://www.w3.org/2000/01/rdf-schema#>",
        "rdf", "<http://www.w3.org/1999/02/22-rdf-syntax-ns#>",
        "fn", "<http://www.w3.org/2005/xpath-functions#>",
        "afn", "<http://jena.hpl.hp.com/ARQ/function#>"
    ),
    country=function(x=NA) if(!is.na(x)) setNames(enipedia.countrynames[x], NULL) else NA,
    statement=paste(
        "select ?Name ?Generation_capacity ?Lat ?Lng ?Fuel ?City ?Plz ?Year_built ?Efficiency ?Plant_type ",
        "where {",
        "?powerPlant prop:Country a:$(country) .",
        "?powerPlant rdfs:label ?Name .",
        "?powerPlant prop:Generation_capacity_electrical_MW ?Generation_capacity . ",
        "OPTIONAL { ?powerPlant prop:Latitude ?Lat . ?powerPlant prop:Longitude ?Lng . }",
        "OPTIONAL { ?powerPlant prop:Fuel_type ?FuelURI . ?FuelURI rdfs:label ?Fuel .}",
        "OPTIONAL { ?powerPlant prop:City ?CityURI . ?CityURI rdfs:label ?City . }",
        "OPTIONAL { ?powerPlant prop:Zipcode ?Plz . }",
        "OPTIONAL { ?powerPlant prop:Year_built ?Year_built . }",
        
        "OPTIONAL { ?powerPlant prop:Operating_efficiency ?Efficiency . }",
        "OPTIONAL { ?powerPlant prop:Power_plant_type ?Plant_type . }",
        "}", sep="")
)


