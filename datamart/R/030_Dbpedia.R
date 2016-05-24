#' A class for querying Dbpedia.org
#' 
#' This class defines some resources at dbpedia.
#' See \code{queries(dbpedia())} for a list of resources.
#'
#' @seealso \code{\link{dbpedia}}, \code{\link{xsparql}}
#'
#' @examples
#' \dontrun{
#'   dbp <- dbpedia()
#'   queries(dbpedia)
#'   query(dbp, "Nuts1")
#' }
#' 
#' @name Dbpedia-class
#' @rdname Dbpedia-class
#' @exportClass Dbpedia
setClass(Class="Dbpedia", representation=representation(), contains="Xsparql")


#' Constructor for Dbpedia objects
#'
#' @param lang   two-character language code for the dbpedia, default ''
#'
#' @return a Dbpedia object, inherited from Xsparql
#' @export
#' @rdname Dbpedia-class
dbpedia <- function(lang="") xsparql(
    url=paste("http://", lang, if(lang!="") ".", "dbpedia.org/sparql", sep=""),
    nspace=c("dbo", "<http://dbpedia.org/ontology/>",
        "rdf", "<http://xmlns.com/foaf/0.1>",
        "rdfs", "<http://www.w3.org/2000/01/rdf-schema#>",
        "owl", "<http://www.w3.org/2002/07/owl#>",
        "xsd", "<http://www.w3.org/2001/XMLSchema#>",
        "dc", "<http://purl.org/dc/elements/1.1/>",
        "foaf", "<http://xmlns.com/foaf/0.1/>",
        #"", "<http://dbpedia.org/resource/>", 
        "dbpedia2", "<http://dbpedia.org/property/>", 
        "dbpedia", "<http://dbpedia.org/>",
        "skos", "<http://www.w3.org/2004/02/skos/core#>"             
    ),
    clss="Dbpedia",
    Nuts1=paste(
        "SELECT ?name, ?nuts, ?popDate, ?pop, ?area, ?gdp, ?popMetro WHERE {",
        "  ?s a <http://dbpedia.org/class/yago/StatesOfGermany>;",
        "     <http://dbpedia.org/property/nuts> ?nuts;",
        "     rdfs:label ?name . ",
        "  OPTIONAL { ?s <http://dbpedia.org/ontology/populationAsOf> ?popDate }",
        "  OPTIONAL { ?s <http://dbpedia.org/property/population> ?pop }",
        "  OPTIONAL { ?s <http://dbpedia.org/property/popMetro> ?popMetro }",
        "  OPTIONAL { ?s <http://dbpedia.org/property/gdp> ?gdp }",
        "  OPTIONAL { ?s <http://dbpedia.org/ontology/areaTotal> ?area } . ",
        "  FILTER (LANG(?name)='de') } ",
        sep=""
    ), 
    PlzAgs=paste(
      "SELECT ?name, ?plz, ?ags, ?tel, ?kfz WHERE {",
      "  ?s <http://de.dbpedia.org/property/wikiPageUsesTemplate> <http://de.dbpedia.org/resource/Vorlage:Infobox_Gemeinde_in_Deutschland> .",
      "  ?s <http://de.dbpedia.org/property/name> ?name .",
      "  ?s <http://dbpedia.org/ontology/postalCode> ?plz .",
      "  ?s <http://dbpedia.org/ontology/municipalityCode> ?ags .",
      "  OPTIONAL {?s <http://dbpedia.org/ontology/areaCode> ?tel }",
      "  OPTIONAL {?s <http://dbpedia.org/ontology/vehicleCode> ?kfz }",
      "  FILTER (LANG(?name)='de') } LIMIT 20",
      sep="\n"
    )
)


#gegen die dt Dbpedia laufen lassen!
# SELECT ?name, ?plz, ?ags WHERE {
     # ?s <http://de.dbpedia.org/property/wikiPageUsesTemplate> <http://de.dbpedia.org/resource/Vorlage:Infobox_Gemeinde_in_Deutschland> .
     # ?s <http://de.dbpedia.org/property/name> ?name .
     # ?s <http://dbpedia.org/ontology/postalCode> ?plz .
     # ?s <http://dbpedia.org/ontology/municipalityCode> ?ags .

# }


# PREFIX dbo: <http://dbpedia.org/ontology/>

# SELECT ?name, ?plz, ?ags, ?geo WHERE {
     # ?s dbpedia2:wikiPageUsesTemplate <http://dbpedia.org/resource/Template:Infobox_German_location> .
     # ?s rdfs:label ?name .
     # ?s dbo:postalCode ?plz .
     # ?s dbo:municipalityCode ?ags .
     # ?s <http://www.georss.org/georss/point> ?geo .
# FILTER(langMatches(lang(?name), "de" ))
# } LIMIT 20


# gemeinde
# PREFIX dbo: <http://dbpedia.org/ontology/>

# SELECT ?name, ?plz, ?ags, ?geo WHERE {
     # ?s rdf:type dbo:Settlement . 
     # ?s rdfs:label ?name .
     # ?s dbo:postalCode ?plz .
     # ?s dbo:municipalityCode ?ags .
     # ?s <http://www.georss.org/georss/point> ?geo .
# FILTER(langMatches(lang(?name), "de" ))
# } LIMIT 20

# all cities with more than 2 Mio inhabitants
# SELECT ?s ?population WHERE {
# ?s rdf:type <http://dbpedia.org/ontology/City>.
# ?s <http://dbpedia.org/ontology/populationUrban> ?population.
# }
# ORDER BY DESC(xsd:integer(?population))
# LIMIT 20

#fails
#SELECT * WHERE {  ?page dbpedia2:wikiPageUsesTemplate <http://dbpedia.org/resource/Template:Infobox_Gemeinde_in_Deutschland> .  ?page dbpedia2:name ?name .}

