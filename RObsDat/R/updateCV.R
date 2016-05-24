updateCV <- function(){
	stopifnot(requireNamespace("SSOAP", quietly=TRUE))
	stopifnot(requireNamespace("XML", quietly=TRUE))
	
	def <- SSOAP::processWSDL("http://his.cuahsi.org/ODMCV_1_1/ODMCV_1_1.asmx?WSDL")

	ff <- SSOAP::genSOAPClientInterface(def=def)

	ans <- ff@functions$GetUnits()
	test <-     XML::xmlToList(XML::xmlParse(ans, asText = TRUE))
	sapply(test$Records, function(x){
				if(any(names(x)!="count")){
					suppressWarnings(addUnits(Name=x$UnitsName, Type=x$UnitsType, Abbreviation=x$UnitsAbbreviation))
				}

		})


	ans <- ff@functions$GetSpatialReferences()
	test <-     XML::xmlToList(XML::xmlParse(ans, asText = TRUE))
	sapply(test$Records, function(x){
				if(any(names(x)!="count")){
					if(is.null(x$SRSID)) x$SRSID="NULL"
					if(is.null(x$Notes)) x$Notes=""
					suppressWarnings(addSpatialReferences(ID=x$SpatialReferenceID, SRSID=x$SRSID, Name=x$SRSName, IsGeographic=as.logical(x$IsGeographic), Notes=x$Notes))
				}

		})

	for(tab in CVtables()){
		command <- paste("ans <- ff@functions$Get", tab, "CV()", sep = "")
		eval(parse(text=command))
		test <- XML::xmlToList(XML::xmlParse(ans, asText = TRUE))
		sapply(test$Records, function(x){
				if(any(names(x)!="count")){
				  suppressWarnings(addCV(table=tab, term=x$Term, definition=x$Definition))
				}
			})
	}


}
