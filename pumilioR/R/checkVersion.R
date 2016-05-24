checkVersion <- function(pumilio_URL, credentials = NA, pumiliologin = NA){
	
	pumilio_URL_o <- pumilio_URL
	
	if (!is.na(credentials)){
		pumilio_URL <- gsub("http://", paste("http://", credentials, "@", sep=""), pumilio_URL)
	}

	if (!is.na(pumiliologin)){
		pumilio_XML_URL <- paste(pumilio_URL, "xml.php?login=", pumiliologin, sep = "")
	}else{
		pumilio_XML_URL <- paste(pumilio_URL, "xml.php", sep = "")
	}
	
  #check valid url
	badurl <- function(...){
		stop(paste("Could not open the XML file at the address given:\n    ", paste(pumilio_URL_o, "xml.php", sep = ""), "\n  please verify the URL is correct.", sep=""))
	}
	
	#Get XML contents
	pumilio_XML <- xmlTreeParse(getURL(pumilio_XML_URL), error = badurl)
		
	pumilio_list <- xmlToList(node = pumilio_XML, addAttributes = TRUE)
	
	pumilio_version <- pumilio_list$pumilio_version
	pumilio_xml_access <- pumilio_list$pumilio_xml_access
	
	if (length(pumilio_version) == 0){
		ret = FALSE
	}else{
		pumilio_version_exp <- unlist(strsplit(pumilio_version, "\\."))
		
		if (pumilio_version_exp[1] > 2){
			ret = TRUE
		}else{
			if (pumilio_version_exp[1] == 2){
				if (pumilio_version_exp[2] >= 6){
					ret = TRUE
				}else{
					ret = FALSE
					}
			}else{
				ret = FALSE
				}
			}
	}

	cat(paste("\n Pumilio is running version ", pumilio_version, "\n   and the XML access is set to: ",  pumilio_xml_access, "\n\n", sep=""))
		
	if (ret == TRUE){
		if (length(pumilio_xml_access) == 0){
			invisible(FALSE)
		}else{
			invisible(pumilio_xml_access)
		}
	}else{
		invisible(FALSE)
	}
}	