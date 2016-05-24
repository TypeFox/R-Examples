##
## Class for NEOS communication objects
##
setClass("NeosComm", representation(url = "character", curlopts = "list", curlhandle = "CURLHandle"))
##
## Class for returned Values from requests to NEOS 
##
setClass("NeosAns", representation(ans = "character", method = "character", call = "call", nc = "NeosComm"))
##
## Class for returned XML template from NEOS 
##
setClass("NeosXml", representation(xml = "XMLNode", method = "character", call = "call", nc = "NeosComm"))
##
## Class for assigned jobnumber and password from NEOS 
##
setClass("NeosJob", representation(jobnumber = "numeric", password = "character", method = "character", call = "call", nc = "NeosComm"))
##
## Class for NEOS offset 
##
setClass("NeosOff", representation(ans = "character", offset = "integer", jobnumber = "numeric", password = "character", method = "character", call = "call", nc = "NeosComm"))

