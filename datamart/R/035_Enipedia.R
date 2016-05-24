#' Interface to the Enipedia Project
#'
#' Enipedia is an active exploration into the applications of wikis and the semantic 
#' web for energy and industry issues. Through this we seek to create a collaborative 
#' environment for discussion, while also providing the tools that allow for data from 
#' different sources to be connected, queried, and visualized from different perspectives.
#'
#' This function creates an datamart with selected queries to the SPARQL endpoint of
#' the enipedia project.
#'
#' @seealso \code{\link{enipedia}}, \code{\link{xsparql}}
#'
#' @return a Mashup object
#' @export
enipedia <- function() datamart(
    enipedia.powerplants(resource="powerplants")
)

# internal function
enipedia.powerplants <- function(resource="powerplants") xsparql(
    resource=resource,
    url="http://enipedia.tudelft.nl/sparql",
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


# internal function
enipedia.countrynames <- setNames(c("Afghanistan", "Albania", "Algeria", "American_Samoa", 
        "Andorra", "Angola", "Antarctica", "Antigua_and_Barbuda", "Argentina", 
        "Armenia", "Aruba", "Australia", "Austria", "Azerbaijan", "Bahamas", 
        "Bahrain", "Bangladesh", "Barbados", "Belarus", "Belgium", "Belize", 
        "Benin", "Bermuda", "Bhutan", "Bolivia", "Bosnia_and_Herzegovina", 
        "Botswana", "Brazil", "Brunei", "Bulgaria", "Burkina_Faso", "Burundi", 
        "Cambodia", "Cameroon", "Canada", "Cape_Verde", "Cayman_Islands", 
        "Central_African_Republic", "Chad", "Chile", "China", "Colombia", 
        "Comoros", "Congo", "Congo_Republic", "Costa_Rica", "Cote_d'Ivoire", 
        "Croatia", "Cuba", "Cyprus", "Czech_Republic", "Denmark", "Djibouti", 
        "Dominica", "Dominican_Republic", "East_Timor", "Ecuador", "Egypt", 
        "El_Salvador", "Equatorial_Guinea", "Eritrea", "Estonia", "Ethiopia", 
        "Faroe_Islands", "Fiji", "Finland", "France", "French_Guiana", 
        "French_Polynesia", "Gabon", "Gambia", "Georgia", "Germany", 
        "Ghana", "Gibraltar", "Greece", "Greenland", "Grenada", "Guam", 
        "Guatemala", "Guinea", "Guinea-Bissau", "Guyana", "Haiti", "Honduras", 
        "Hong_Kong_(China)", "Hungary", "Iceland", "India", "Indonesia", 
        "Iran", "Iraq", "Ireland", "Isle_Of_Man", "Israel", "Italy", 
        "Jamaica", "Japan", "Jordan", "Kazakhstan", "Kenya", "Kiribati", 
        "Kuwait", "Kyrgyzstan", "Laos", "Latvia", "Lebanon", "Lesotho", 
        "Liberia", "Libya", "Liechtenstein", "Lithuania", "Luxembourg", 
        "Macedonia", "Madagascar", "Malawi", "Malaysia", "Maldives", 
        "Mali", "Malta", "Marshall_Islands", "Mauritania", "Mauritius", 
        "Mayotte", "Mexico", "Micronesia", "Moldova", "Monaco", "Mongolia", 
        "Montenegro", "Morocco", "Mozambique", "Myanmar", "Namibia", 
        "Nauru", "Nepal", "Netherlands", "New_Caledonia", "New_Zealand", 
        "Nicaragua", "Niger", "Nigeria", "North_Korea", "Norway", "Oman", 
        "Pakistan", "Palau", "Palestine", "Panama", "Papua_New_Guinea", 
        "Paraguay", "Peru", "Philippines", "Poland", "Portugal", "Puerto_Rico", 
        "Qatar", "Reunion", "Romania", "Russia", "Rwanda", "Saint_Helena", 
        "Saint_Kitts_and_Nevis", "Saint_Lucia", "Samoa", "Sao_Tome_and_Principe", 
        "Saudi_Arabia", "Senegal", "Serbia", "Seychelles", "Sierra_Leone", 
        "Singapore", "Slovakia", "Slovenia", "Solomon_Islands", "Somalia", 
        "South_Africa", "South_Korea", "Spain", "Sri_Lanka", "St_Kitts_and_Nevis", 
        "St_Lucia", "St_Vincent_and_Grenadines", "Sudan", "Suriname", 
        "Swaziland", "Sweden", "Switzerland", "Syria", "Taiwan", "Tajikistan", 
        "Tanzania", "Thailand", "Timor-Leste", "Togo", "Tonga", "Trinidad_and_Tobago", 
        "Tunisia", "Turkey", "Turkmenistan", "Tuvalu", "Uganda", "Ukraine", 
        "United_Arab_Emirates", "United_Kingdom", "United_States", "Uruguay", 
        "Uzbekistan", "Vanuatu", "Venezuela", "Vietnam", "Western_Sahara", 
        "Yemen", "Zambia", "Zimbabwe"), c("AF", "AL", "DZ", 
        "AS", "AD", "AO", "AQ", "AG", "AR", "AM", "AW", "AU", "AT", "AZ", 
        "BS", "BH", "BD", "BB", "BY", "BE", "BZ", "BJ", "BM", "BT", "BO", 
        "BA", "BW", "BR", "BN", "BG", "BF", "BI", "KH", "CM", "CA", "CV", 
        "KY", "CF", "TD", "CL", "CN", "CO", "KM", "CG", "CD", "CR", "CI", 
        "HR", "CU", "CY", "CZ", "DK", "DJ", "DM", "DO", "TL", "EC", "EG", 
        "SV", "GQ", "ER", "EE", "ET", "FO", "FJ", "FI", "FR", "GF", "PF", 
        "GA", "GM", "GE", "DE", "GH", "GI", "GR", "GL", "GD", "GU", "GT", 
        "GN", "GW", "GY", "HT", "HN", "HK", "HU", "IS", "IN", "ID", "IR", 
        "IQ", "IE", "IM", "IL", "IT", "JM", "JP", "JO", "KZ", "KE", "KI", 
        "KW", "KG", "LA", "LV", "LB", "LS", "LR", "LY", "LI", "LT", "LU", 
        "MK", "MG", "MW", "MY", "MV", "ML", "MT", "MH", "MR", "MU", "YT", 
        "MX", "FM", "MD", "MC", "MN", "ME", "MA", "MZ", "MM", "NA", "NR", 
        "NP", "NL", "NC", "NZ", "NI", "NE", "NG", "KP", "NO", "OM", "PK", 
        "PW", "PS", "PA", "PG", "PY", "PE", "PH", "PL", "PT", "PR", "QA", 
        "RE", "RO", "RU", "RW", "SH", "KN", "LC", "WS", "ST", "SA", "SN", 
        "RS", "SC", "SL", "SG", "SK", "SI", "SB", "SO", "ZA", "KR", "ES", 
        "LK", "KN", "LC", "VC", "SD", "SR", "SZ", "SE", "CH", "SY", "TW", 
        "TJ", "TZ", "TH", "TL", "TG", "TO", "TT", "TN", "TR", "TM", "TV", 
        "UG", "UA", "AE", "GB", "US", "UY", "UZ", "VU", "VE", "VN", "EH", 
        "YE", "ZM", "ZW")
    )
    
# select ?Name ?Point ?Generation_capacity where {
# ?powerPlant prop:Country a:Netherlands .
# ?powerPlant rdfs:label ?Name .
# ?powerPlant prop:Point ?Point .
# ?powerPlant prop:Generation_capacity_electrical_MW ?Generation_capacity . 
# }

# "select ?plant ?Country ?GenerationMWh (?cap as ?CapacityMW) (?GenerationMWh / (365*24) / ?cap * 100 as ?LoadFactorPct) ?Fuel where {",
    # "?Country a cat:Europe .",
    # "?plant prop:Country ?Country;",
    # "prop:Annual_Energyoutput_MWh ?GenerationMWh;",
    # "prop:Generation_capacity_electrical_MW ?cap;",
    # "prop:Primary_fuel_type ?Fuel .",
    # "FILTER(?cap > 0)",
# "} order by desc(?GenerationMWh ) ",
# "limit 300",

# select ?Name ?Point ?Generation_capacity where {
# ?powerPlant prop:Country a:Netherlands .
# ?powerPlant rdfs:label ?Name .
# ?powerPlant prop:Point ?Point .
# ?powerPlant prop:Generation_capacity_electrical_MW ?Generation_capacity . 
# }


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

