tts_ITRI_getID <- function(content, speed, volume, speaker){
  # For version 0.32 #########################
  # sometimes the ITRI server is responding very slowly. This will cause "Timed out" error for this function
  # Here I use try() function and nested structureto figure this out.
  # For version 0.32 #########################
  inner <- function(content, speed, volume, speaker){
    headerFields <- 
      c(Accept = "text/xml",
        Accept = "multipart/*",
        'Content-Type' = "text/xml; charset=utf-8",
        SOAPAction = "http://tts.itri.org.tw/TTSService/ConvertText")
    
    body_1 <- '<?xml version="1.0" encoding="utf-8"?>
    <soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
    <soap:Body>
    <ConvertText>
    <accountID>test-for-r</accountID>
    <password>test1for1r</password>
    <TTStext>'
    body_2 <- '</TTStext>
    <TTSSpeaker>'
    body_3 <- '</TTSSpeaker>
    <volume>'
    body_4 <- '</volume>
    <speed>'
    body_5 <- '</speed>
    <outType>wav</outType>
    </ConvertText>
    </soap:Body>
    </soap:Envelope>'
    
    body <- paste(body_1, content, body_2, speaker, 
                  body_3, volume, body_4, speed, body_5, 
                  sep="")
    
    # we want the actual body of the response from the web server.
    # To get this, basicTextGatherer() is used to collect the text.
    reader =  basicTextGatherer()  
    curlPerform(url = "http://tts.itri.org.tw/TTSService/Soap_1_3.php?wsdl",
                httpheader = headerFields,
                postfields = body,
                verbose = FALSE,
                writefunction = reader$update
    )
    voice_ID_raw <- reader$value()
    voice_ID <- strsplit(strsplit(voice_ID_raw, split="</Result>")[[1]][1], split = '&amp;')[[1]][3]
    return(voice_ID)
  }
  
  # the part below is modified for vision 0.32
  try(result <- inner(content, speed, volume, speaker), silent = TRUE)
  if(("result" %in% ls()) == FALSE){
    cat("\n The server is responding slowly. Please try again later.")
    return(NULL)
  }
  
  return(result)
}







