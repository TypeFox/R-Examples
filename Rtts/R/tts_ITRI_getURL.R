tts_ITRI_getURL <- function(voice_ID){
  
  #revision for v0.32#####################
  if(is.null(voice_ID) == TRUE){
    cat("The voice ID entered is invalid. Please check.")
    return(NULL)
  }
  #revision for v0.32#####################
  
  headerFields <- 
    c(Accept = "text/xml",
      Accept = "multipart/*",
      'Content-Type' = "text/xml; charset=utf-8",
      SOAPAction = "http://tts.itri.org.tw/TTSService/GetConvertStatus")
  
  body_head <- '<?xml version="1.0" encoding="utf-8"?>
  <soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
    <soap:Body>
      <GetConvertStatus>
      <accountID>test-for-r</accountID>
      <password>test1for1r</password>
      <convertID>'
  body_tail <- '</convertID>
      </GetConvertStatus>
    </soap:Body>
  </soap:Envelope>'
  
  body <- paste(body_head, voice_ID, body_tail, sep="")
  
  # we want the actual body of the response from the web server.
  # To get this, basicTextGatherer() is used to collect the text.
  reader =  basicTextGatherer()  
  curlPerform(url = "http://tts.itri.org.tw/TTSService/Soap_1_3.php?wsdl",
              httpheader = headerFields,
              postfields = body,
              verbose = FALSE,
              writefunction = reader$update
  )
  voice_URL_raw <- reader$value()
  voice_URL <- strsplit(strsplit(voice_URL_raw, split = "</Result>")[[1]][1] ,split = "completed&amp;")[[1]][2]
  return(voice_URL)
}







