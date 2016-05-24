# Stores all functions dealing directly with the API

#' .avisApiBusAvanzada
#' 
#' @param args list of arguments must have the exact names expected by the underlying API
#'
#' @return dataframe

.avisApiBusAvanzada <- function(args = list())
{
  args['formato_consulta'] <- '';
  args['tipo_consulta'] <- 'tabla';
  args['control'] <- 1;
  args['excel'] <- 1;

	queryString <- ''

	for (argName in names(args)) {
		queryString <- paste(queryString, argName, "=", args[argName], '&', sep = "")
	}
	queryString <- substr(queryString,0,nchar(queryString)-1)

	url <- paste("http://proyectoavis.com/cgi-bin/bus_avanzada.cgi", queryString, sep = "?")

	.avisVerboseMessage(paste("INFO: querying to proyectoavis.com at: ", url))

	rawData <- .avisGetURL(url)

  data <- read.csv(textConnection(rawData), sep = ";", quote = "")
  
  return (data)
}

#' .avisApiBusOrden
#' 
#' Client for the "bus_orden" endpoint at proyectoavis. API version 1.0
#' 
#' @usage .avisApiBusOrden()
#' @note This function does not allow arguments
#' @return double List id to species name

.avisApiBusOrden <- function()
{
  rawhtml<- .avisGetURL("http://proyectoavis.com/cgi-bin/bus_orden.cgi", TRUE)
  id<- grep("id_especie=[0-9]+",rawhtml)
  ids<- str_extract_all (rawhtml, "id_especie=[0-9]+")
  id_specie<- as.numeric (substring(unlist (ids [c(id)]), 12))

  name<- grep("<i>.*?</i>", rawhtml)
  names<- str_extract_all (rawhtml, "<i>.*?</i>")
  names_specie<- substring (unlist (names[c(name)]),4)
  names_species<- .avisNormalizeSpeciesName(substr(names_specie, 1, nchar(names_specie)-4))
  names(id_specie) <- names_species

  return (id_specie)
}

#' .avisApiBusEspecie
#' 
#' Client for the "bus_especie" endpoint at proyectoavis. API version 1.0
#' 
#' @usage .avisApiBusEspecie()
#' @note This function does not allow arguments
#' @return dataframe reflecting the table structure specified by the API

.avisApiBusEspecie <- function()
{
    tables<- XML::readHTMLTable ("http://proyectoavis.com/cgi-bin/bus_especie.cgi")
    table_obs<- tables[[4]]
    observ<-  table_obs[4:dim(table_obs)[1],3:6]
    names (observ)<- c("Observations", "Individuals", "UTM.10x10", "Birdwatchers")
    spsummary<- data.frame (lapply(observ, as.numeric), stringsAsFactors=FALSE)
    row.names (spsummary)<- table_obs [4:dim(table_obs)[1],2]

    return (spsummary)
}

#' .avisApiFichaUsuario
#' 
#' Client for the "ficha_usuario" endpoint at proyectoavis. API version 1.0
#' 
#' @usage .avisApiFichaUsuario(32)
#' @param usuario_id Id of the colaborator in proyectoavis
#' @return dataframe reflecting the table structure specified by the API

.avisApiFichaUsuario <- function(usuario_id)
{
  doc<-XML::htmlParse(paste ("http://proyectoavis.com/cgi-bin/ficha_usuario.cgi?id_usuario=", usuario_id, sep=""))
  nodes <- XML::getNodeSet(doc, "//table[@class=\"observaciones\"][2]//tr")

  df<-data.frame()
  
  for(node in nodes[2:length(nodes)]){
    df<-rbind(df, .avisApiFichaUsuarioExtractRowData(node))
  }

  # preliminar header spec
  names (df)<- c("SpeciesId", "Observations", "Number", "UTM.10x10", "Birdwatchers")

  return (df)  
}

#' .avisApiUsuarios
#' 
#' Client for the "usuarios" endpoint at proyectoavis. API version 1.0
#' 
#' @usage .avisApiUsuarios()
#' @note This function does not allow arguments
#' @return dataframe reflecting the table structure specified by the API

.avisApiUsuarios <- function()
{
  doc<-XML::htmlParse("http://proyectoavis.com/cgi-bin/usuarios.cgi")
  nodes <- XML::getNodeSet(doc, "//table[@class=\"observaciones1\"]/tr[@class=\"celda1\"]")

  df<- NULL
  for (node in nodes) {
    df<- rbind(df, .avisApiUsuariosExtractRowData(node))
  }

  # preliminar header
  colnames(df)<- c("UserId", "User","Observations","Species","Provinces","UTMs","Periods")

  return (df)
}


#' .avisApiFichaUsuarioExtractRowData
#'
#' Helper function to extract data from XML node
#'
#' @param node XML node
#' @return row

.avisApiFichaUsuarioExtractRowData <- function(node)
{
  clean_row_data <- XML::xmlValue(node, encoding="utf-8")
  celdas<-as.list(strsplit(gsub("\n","#", clean_row_data), "#")[[1]])
  
  # id de especie
  id<-avisSpeciesId(celdas[2])
  
  # remove species common and scientific name
  celdas<- celdas[3:length(celdas)]
  
  obsdata<- as.integer(c(id, celdas))
  
  return (obsdata)
}

#' .avisApiUsuariosExtractRowData
#'
#' Helper function to extract data from XML node
#'
#' @param node XML node
#' @return row

.avisApiUsuariosExtractRowData <- function(node)
{
  strnode <- XML::toString.XMLNode(node)   
  usu_id<-regmatches(strnode, regexpr('id_usuario=([0-9]+)', strnode))
  id<-as.integer(regmatches(usu_id, regexpr('[0-9]+', usu_id)))
  
  clean_row_data <- XML::xmlValue(node, encoding="utf-8")
  celdas<-as.list(strsplit(gsub("\n","#", clean_row_data), "#")[[1]][1:8])

  # remove username and name
  userdata<- c(id, celdas[2], as.integer(celdas[c(4:length(celdas))]))

  return (userdata)
}
