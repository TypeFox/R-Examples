# Functions related to search from http://proyectoavis.com/cgi-bin/bus_avanzada.cgi

# translate
.ravis_translated_params_map<- list( 
	id_species = 'id_especie', # string / list -> id_especie
	family = 'familia', # familia
	order = 'orden', # orden
	age = 'edad',
	sex = 'sexo', # sexo: 'macho', 'hembra', 'indeterminado', 'pareja', 'machos y hembras'
	breeding = 'reproduccion',
	habitat = 'habitat', 
	# habitat: 'terrenos agrícolas', 'roquedos de interior', 
	# 'zonas humanizadas', 'zonas húmedas interiores', 'bosque',
	# 'pastizales', 'costas', 'matorral', 'otros'
	month = 'mes', 
	year = 'ano'
)

# query parameters accepted by proyectoavis.com
.ravis_raw_search_default_params<- list( 
	formato_consulta = 'observaciones', 
	tipo_consulta = '', 
	id_observacion = '', 
	id_periodo = '', 
	id_especie = '', 
	orden = '',
	criterio = 'id_observacion', 
	familia = '', 
	edad = '', sexo = '', 
	usu = '', 
	id_ca = '',	id_provincia = '', 
	dia_ini = '', mes_ini = '', ano_ini = '', 
	dia_fin = '', mes_fin = '', ano_fin = '', 
	mes = '', ano = '', plazo = '', 
	hora_ini = '', minuto_ini = '', 
	hora_fin = '', minuto_fin = '', 
	reproduccion = '', habitat = '', codigo_habitat = '', 
	gr = '', cf = '', utm_10 = '', utm_1 = '', 
	menu = '', tipo_grafica = 'comparadas', 
	filtro_mes = '', filtro_ano = '', 
	filtro_id_especie = '', filtro_estacion = '', cobertura = '', 
	mostrar_capa = '', capa = ''
)

#' avisQuerySpecies
#' 
#' Is a wrapper for avisQuery that allows you to perform 
#' a query for more than one species at once.'names' must 
#' be either a string or a list of species names, 'args' 
#' is a list of query parameters 
#' (see avisQuery) that adds further filters to the query.
#' 
#' @usage avisQuerySpecies(names, args = list())
#' @param names Must be either a string or a list of scientific names
#' @param args A list of normalized parameters to add filters to the query. 
#' Currently in Spanish, but this might become outdated. See avisQuery.
#' @return a dataframe with the results of your specific query to Proyecto AVIS database
#' @export 
#' @examples \dontrun{
#' avisQuerySpecies("Bubo bubo")
#' avisQuerySpecies(list("Bubo bubo", "Tyto alba"), args = list(year = 2012))
#' }

avisQuerySpecies <- function (names, args = list()) 
{
	if(is.element('id_species', names(args)))
	{
		warning("id_species argument in the argument list won't be regarded")
	}

	if(!is.list(names)){ 
		names = list(names)
	}

	# check all species exists in bd before querying
	lapply(names, function(n){ 
		if(!avisHasSpecies(n)) stop(paste("Species not found: ", n)) 
	})

	df<- NULL
	for (name in names) {
		args['id_species'] <- avisSpeciesId(name)
		df<- rbind(df, avisQuery(args = args))
	}

	return (df)
}

#' avisQueryContributor
#' 
#' Is a wrapper for avisQuery that allows you to perform 
#' a query for more than one contributor at once.
#' 
#' @usage avisQueryContributor(contributor_ids, args = list())
#' @param contributor_ids must be either an integer or a list of contributors ids (integers)
#' @param args A list of normalized parameters to add filters to the query. 
#' Currently in Spanish, but this might become outdated. See avisQuery.
#' @return a dataframe with the results of your specific query to Proyecto AVIS database
#' @export 
#' @seealso avisContributorsSummary
#' @examples \dontrun{
#' avisQueryContributor(370)
#' avisQueryContributor(list(370, 399), args = list(year = 2002))
#' }
 
avisQueryContributor <- function (contributor_ids, args = list()) 
{
	if(is.element('usu', names(args))){
		warning("usu argument in the argument list won't be regarded")
	}

	ravis_username_id_list = .avisUserNameList()

	if(!is.list(contributor_ids)){ 
		contributor_ids = list(contributor_ids)
	}

	names = lapply(contributor_ids, function(id){
		name = ravis_username_id_list[as.character(id)][[1]]

		if(is.null(name) | is.na(name)){
			stop(paste("Contributor with id ",id," not found"))
		}
		return (name)
	})

	df<- NULL
	for (name in names) {
		args['usu'] <- name
		df<- rbind(df, avisQuery(args = args))
	}

	return (df)
}

#' avisQuery
#' 
#' General function for querying the database using several filters, like order, 
#' family, species, age, sex, habitat, etc. 
#' 
#' In case you set a query parameter by its name (eg: avisQuery (species="Bubo bubo"))
#' and also you set it inside the 'args' parameter (eg: avisQuery (species="Bubo bubo", args=list(species="Tyto alba")), 
#' the value setted by its name will prevail (in the example, "Bubo bubo" will apply).
#' 
#' @usage avisQuery(id_species = "", species = "", family = "", order = "", 
#' age = "", sex = "", breeding = "", habitat = "", month = "", year = "",
#' args = list())
#' @param id_species a number setting the id of the species according to proyectoavis.com database. 
#' You may get the id of a species with \code{\link{avisSpeciesId}}
#' @param species scientific name of the species (one single species): e.g. "Passer domesticus"
#' @param family To filter the data by family: e.g. "Passeridae", "Falconidae", etc.  
#' @param order To filter the data by Order: e.g. "Passeriformes", "Falconiformes", etc.
#' @param age To filter the data by age: "pollo", "juvenil", "adulto", "indeterminado".
#' @param sex To filter the data by sex: "macho", "hembra", "indeterminado", "pareja", "machos y hembras"
#' @param breeding To filter the data by breeding-migratory status: "reproducción posible", "reproducción probable", 
#' "reproducción segura", "migración", "invernada"
#' @param habitat Filter by habitat: "bosque", "matorral", "pastizales", "terrenos agrícolas", "zonas humanizadas", 
#' "zonas húmedas interiores", "roquedos de interior", "costas", "otros"
#' @param month Filter by month: 1 to 12
#' @param year Filter by year: e.g. 2001
#' @param args List of arguments accepted by www.proyectoavis.com endpoint. You may use
#' this list to set the arguments of the function (species, sex, breeding...), or 
#' you may also set all the parameters supported by the endpoint, but not normalized for its use in this package. 
#' These arguments are: id_ca, id_provincia, dia_ini, mes_ini, ano_ini, dia_fin, mes_fin, ano_fin, usu, 
#' plazo, hora_ini, minuto_ini, hora_fin, minuto_fin, codigo_habitat, gr, cf, utm_10, utm_1 (see www.proyectoavis.com)
#' @return a dataframe with the results of your specific query to Proyecto AVIS database.
#' @export 
#' @examples \dontrun{
#' # get all the observations of the species of the Order Falconiformes
#' avisQuery (order = "Falconiformes") 
#' # get all the observations of the species of the Family Falconidae
#' avisQuery(family = "Falconidae")
#' # get the observations of immatures of Iberian Imperial Eagle
#' avisQuery (species= "Aquila adalberti", age = "juvenil")
#' }

avisQuery <- function (id_species = '', species = '', family = '', order = '', age = '', 
	sex = '', breeding = '', habitat = '', month = '', year = '', args = list())
{
	if (id_species != '') args['id_species'] <- id_species
	if (species != '') args['species'] <- species
	if (family != '') args['family'] <- family
	if (order != '') args['order'] <- order
	if (age != '') args['age'] <- age
	if (sex != '') args['sex'] <- sex
	if (breeding != '') args['breeding'] <- breeding
	if (habitat != '') args['habitat'] <- habitat
	if (month != '') args['month'] <- month
	if (year != '') args['year'] <- year
	
	if(is.element('species', names(args)) && is.element('id_species', names(args)))
	{
		warning(paste("ATENTION!: you setted 'species' (", args['species'], ") and 'id_species' (", args['id_species'], ") in your query. The parameter id_species will be discarded"))
	}

	# species id
	if(is.element('species', names(args)))
	{
		args['id_species'] <- avisSpeciesId(args['species'])
		args['species']<-NULL	
	}

	rawargs <- .avisTranslateArgsToRawArgs(args)

	return (.avisQueryRaw(rawargs))
}

#' avisTranslateArgsToRawArgs
#'
#' Tranlate args (set by user) to rawargs, which can be handled by server (in spanish)
#'
.avisTranslateArgsToRawArgs<-function(args)
{
	rawargs<-args

	trans_arg_names <- names(.ravis_translated_params_map)

	for (argname in names(args)) {
		# if argname is a translated param
		if(is.element(argname, trans_arg_names)){
			raw_param_name<-.ravis_translated_params_map[argname][[1]]
			rawargs[raw_param_name] <- args[argname]

			if(!is.element(argname, .ravis_translated_params_map)){
				rawargs[argname]<-NULL
			}
		}
	}

	return (rawargs)
}

#' .avisQueryRaw
#'
#' Performs a query on observations for given arguments
#'
#' @param args list of arguments must have the exact names expected by the underlying API
#' @return dataframe
	
.avisQueryRaw <- function (args)
{
	if(!is.list(args)){
		stop("Object of type 'list' expected for query args")
	}

	args<-.avisMergeArgumentList(args, .ravis_raw_search_default_params)

	data <- .avisApiBusAvanzada(args)

	data <- .addLatLonColumns(data)

	return(data)
}

.addLatLonColumns <- function (data)
{
	utm_latlon<-.getUTMLatlong()

	x = utm_latlon$x [match (substring(data$UTM,4), utm_latlon$utm)]

	y = utm_latlon$y [match (substring(data$UTM,4), utm_latlon$utm)]
	
	data<- data.frame(data, "x"= x, "y"= y)

	return (data)
}

# merge two argument list. first argument lists overwrite seccond (default)
.avisMergeArgumentList<-function(args, defaultArgs)
{
	for (argName in names(defaultArgs)) {
		if(!is.element(argName, names(args))){
			args[argName] <- defaultArgs[argName]
		}
	}

	return (args)
}

.getUTMLatlong<- function()
{
	# hack for avoiding NOTE on check: 'no visible binding for global variable'
  	# see: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  	ravisUTMLatLong <- NULL
  	rm(ravisUTMLatLong)

	return (ravisUTMLatLong)
}