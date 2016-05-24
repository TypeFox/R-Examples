################################################################################
# Copyright (C) 2010 by 52 North                                               #
# Initiative for Geospatial Open Source Software GmbH                          #
#                                                                              #
# Contact: Andreas Wytzisk                                                     #
# 52 North Initiative for Geospatial Open Source Software GmbH                 #
# Martin-Luther-King-Weg 24                                                    #
# 48155 Muenster, Germany                                                      #
# info@52north.org                                                             #
#                                                                              #
# This program is free software; you can redistribute and/or modify it under   #
# the terms of the GNU General Public License version 2 as published by the    #
# Free Software Foundation.                                                    #
#                                                                              #
# This program is distributed WITHOUT ANY WARRANTY; even without the implied   #
# WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU #
# General Public License for more details.                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program (see gpl-2.0.txt). If not, write to the Free Software           #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or #
# visit the Free Software Foundation web page, http://www.fsf.org.             #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2010-06-18                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#####################################################
# built-in download.file
download.file(url="http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos?request=GetCapabilities&version=1.0.0&service=SOS", destfile="text.xml")
# works, but requires saving the file locally... 

#####################################################
# http://cran.r-project.org/web/packages/httpRequest/
library("httpRequest")

# GetCapabilities request with all sections
getCapRequest <- '<?xml version="1.0" encoding="UTF-8"?><GetCapabilities xmlns="http://www.opengis.net/sos/1.0" xmlns:ows="http://www.opengis.net/ows/1.1" xmlns:ogc="http://www.opengis.net/ogc" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.opengis.net/sos/1.0 http://schemas.opengis.net/sos/1.0.0/sosGetCapabilities.xsd" service="SOS"><ows:AcceptVersions><ows:Version>1.0.0</ows:Version></ows:AcceptVersions><ows:Sections><ows:Section>OperationsMetadata</ows:Section><ows:Section>ServiceIdentification</ows:Section><ows:Section>ServiceProvider</ows:Section><ows:Section>Filter_Capabilities</ows:Section><ows:Section>Contents</ows:Section></ows:Sections></GetCapabilities>'

# variables
port <- 8080
host <- "http://giv-sos.uni-muenster.de"
host2 <- "http://mars.uni-muenster.de"
path <- "/52nSOSv3/sos?request=GetCapabilities&version=1.0.0&service=SOS"
path2 <- "OWS5SOS"
referer <- "52north.org/sos4r"

# getToHost
getcap <- getToHost(host=host, path=path, referer=referer, port=port)
# DOES NOT WORK...
# Error in make.socket(host = host, port = port, server = FALSE) : socket not established

# getToHost2
getcap <- getToHost2(host=host, path=path, referer=referer, port=port)
# Error in socketConnection(host = host, port = port, open = "a+b", blocking = TRUE) : 
#		cannot open the connection
# In addition: Warning message:
# In socketConnection(host = host, port = port, open = "a+b", blocking = TRUE) :
#  http://giv-sos.uni-muenster.de:8080 cannot be opened


# simplePostToHost
path3 <- "/52nSOSv3/sos"
data <- "request=GetCapabilities&version=1.0.0&service=SOS"
getCap <- simplePostToHost(host=host, path=path3, referer=NULL, data = data, port=port)

# postToHost
datalist <- list("request"="GetCapabilities", "version"="1.0.0", "service"="SOS")
getCap <- postToHost(host=host, path=path3, data.to.send=datalist, referer=NULL, port=port)
# Error in make.socket(host = host, port = port, server = FALSE) : socket not established


################################################################################
# http://www.omegahat.org/RCurl/
# install.packages(packageName, repos = "http://www.omegahat.org/R")
library("RCurl")
temp <- getURL("www.google.de")

# getForm
temp <- getURL(host, port=port, referer=referer, verbose=TRUE)
temp <- getURL("http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos?request=GetCapabilities&version=1.0.0&service=SOS")
# GEHT!

temp <- getForm("http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos", request="GetCapabilities", version="1.0.0", service="SOS")
# GEHT AUCH! Ist auch sauberer.


sosUrl = "http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos"
myCurlOptions = c(verbose = TRUE)
myParams <- list("request"="GetCapabilities", "version"="1.0.0", "service"="SOS")

temp <- postForm(sosUrl, .params = myParams)
# liefer nur raw[0]
rawToChar(temp)
# ""

# postForm
temp <- postForm("http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos", request="GetCapabilities", version="1.0.0", service="SOS")
# liefert exception report... unexpected element CDATA

# need to use getCapRequest
temp <- postForm("http://giv-sos.uni-muenster.de:8080/52nSOSv3/sos", "request"=getCapRequest)
# encoding error




################################################################################
# http://www.omegahat.org/SSOAP/
library("SSOAP")
library("XML")
library("XMLSchema")

#
#
#
myParseSchemaDoc <- function (url, removeComments = TRUE, 
		namespaces = c(xs = "http://www.w3.org/2001/XMLSchema"), 
		followImports = TRUE, followIncludes = followImports) {
	
	print("START...")
	
	baseURL = dirname(url)
	doc = xmlInternalTreeParse(url)
	
	print(paste("parsed doc from url:", url))
	
	if (length(namespaces) == 0 || !("xs" %in% names(namespaces))) 
		names(namespaces)[1] = "xs"
	if (followImports) {
		imports = getNodeSet(doc, "//xs:schema/xs:import", namespaces)
		imports = lapply(imports, function(node) {
					xdoc = myImportSchema(node, baseURL)
					if (is.null(xdoc)) {
						removeNodes(node)
						return(NULL)
					}
					schema = getNodeSet(xdoc, "//xs:schema", namespaces)
					sapply(schema, function(s) replaceNodes(node, xmlRoot(s)))
				})
	}
	if (followIncludes) {
		includes = getNodeSet(doc, "//xs:schema/xs:include", 
				namespaces)
		if (length(includes)) {
			sapply(includes, function(node) {
						xdoc = myImportSchema(node, baseURL)
						schema = getNodeSet(xdoc, "//xs:schema", c(xs = "http://www.w3.org/2001/XMLSchema"))
						p = xmlParent(node)
						sapply(schema, function(s) addChildren(p, kids = xmlChildren(s)))
					})
			removeNodes(includes, TRUE)
		}
	}
	if (removeComments) {
		comments = getNodeSet(doc, "//comment()", c(xs = "http://www.w3.org/2001/XMLSchema"), 
				noMatchOkay = TRUE)
		if (length(comments)) 
			removeNodes(comments)
	}
	
	print("DONE.")
	
	doc
}

myImportSchema <- function (node, baseURL) {
	print(paste("Import schema from base url ", baseURL))
	
	u = xmlGetAttr(node, "schemaLocation", NA)
	if (is.na(u)) 
		return(NULL)
	u = getRelativeURL(u, baseURL)
	myParseSchemaDoc(u)
}

myProcessWSDL <- function (fileName, handlers = WSDLParseHandlers(fileName), 
		nameSpaces = character(), useInternalNodes = TRUE, verbose = FALSE) {
	if (!is(fileName, "XMLAbstractDocument")) {
		if (useInternalNodes) {
			# wsdl = myParseSchemaDoc(fileName)
			wsdl = myParseSchemaDoc(fileName, followImports = FALSE,
					followIncludes = FALSE)
			print("Not following imports and includes!")
		}
		else wsdl = xmlTreeParse(fileName, handlers = handlers, 
					asTree = TRUE, fullNamespaceInfo = TRUE)
	}
	root = xmlRoot(wsdl)
	port = root[["service"]][["port"]]
	
	if (sum(xmlSApply(root[["service"]], xmlName) == "port") > 
			1) 
		warning("Ignoring additional <service><port> ... elements")
	
	loc = xmlGetAttr(port[["address"]], "location")
	server = SOAPServer(loc)
	
	print("Server:")
	print(server)
	
	types = processSchemaTypes(root[["types"]], root, verbose = verbose)
	tmp = root[names(root) == "binding"]
	ops = lapply(tmp, processWSDLBindings, root, types)
	names(ops) = sapply(tmp, xmlGetAttr, "name")
	if (missing(nameSpaces)) {
		schemaURIs = sapply(.SOAPDefaultNameSpaces, function(x) x["xsd"])
		uris = sapply(xmlNamespaceDefinitions(root), function(x) x$uri)
		i = match(uris, schemaURIs)
		nameSpaces = NA
		if (!all(is.na(i))) 
			nameSpaces = names(.SOAPDefaultNameSpaces)[i[!is.na(i)]]
	}
	
	SOAPServerDescription(name = xmlGetAttr(port, "name"), server = server, 
			operations = ops, types = types, nameSpaces = as.character(nameSpaces))
}

# works fine:
# tmp = processWSDL(system.file("examples", "KEGG.wsdl", package = "SSOAP"))

sosWsdlFile = "/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/data/SOS.wsdl"

# figure out why R crashes on processWSDL... might be problem with parseSchemaDoc
sosWsdl <- myParseSchemaDoc(sosWsdlFile)
# crashes R, "segfault"

sosWsdl <- myParseSchemaDoc(sosWsdlFile, followImports = FALSE,
		followIncludes = FALSE)
# goes through

# parts of schemas work, so problem probably lies in lack of memory?
gmlFeature <- myParseSchemaDoc("http://schemas.opengis.net/gml/3.1.1/base/feature.xsd")
gmlBasictype <- myParseSchemaDoc("http://schemas.opengis.net/gml/3.1.1/base/basicTypes.xsd")


w = myProcessWSDL(
		fileName = sosWsdlFile,
		verbose = TRUE,
		useInternalNodes = TRUE)
# useInternalNodes = FALSE -> Error: /home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/data/sosCommon.xsd  does not seem to be XML, nor to identify a file name
# Probably a resolvement issue.
# uesInternaLnodes = TRUE  -> crashes as well...

w = processWSDL(
		fileName = sosWsdlFile,
		verbose = TRUE,
		useInternalNodes = FALSE)


# TODO CONTINUE HERE with SOAP implementation


iface = genSOAPClientInterface(def = sosWsdl, verbose = TRUE)
# namespaces = "1.1!
# SSOAP:::.SOAPDefaultNameSpaces
# SSOAP:::.SOAPDefaultHandlers

iface

?parseSOAP
