GenerateAnimationKMLFile_Track <- 
  function (sInputFile, sid, sPointsFile, 
            sOutputFile, sTrackColour) 
  {
    
    TransmitterList <- ExtractUniqueValues(sInputFile,2)
    sInputFile2 <- ExtractData(sInputFile,sQueryTransmitterList = sid)
    
    sInputFile2 <- merge(sInputFile2,sPointsFile,by.x="RECEIVERID",by.y="LOCATION")
    sInputFile2 <- sInputFile2[order(sInputFile2$DATETIME),]
    
    outfile <- file(sOutputFile, "wt")
    
    writeLines('<?xml version="1.0" encoding="UTF-8"?>',
               outfile)
    writeLines('<kml xmlns="http://www.opengis.net/kml/2.2"',
               outfile)
    writeLines('xmlns:gx="http://www.google.com/kml/ext/2.2"',
               outfile)
    writeLines('xmlns:kml="http://www.opengis.net/kml/2.2"',
               outfile)
    writeLines('xmlns:atom="http://www.w3.org/2005/Atom">',
               outfile)             
    writeLines("<Document>", outfile)
    writeLines(paste("<name>", sOutputFile, "</name>", sep = ""), outfile)
    writeLines('<StyleMap id="msn_track">', outfile)
    writeLines('<Pair>', outfile)
    writeLines('<key>normal</key>', outfile)
    writeLines('<styleUrl>#sn_track</styleUrl>', outfile)
    writeLines('</Pair>', outfile)
    writeLines('<Pair>', outfile)
    writeLines('<key>highlight</key>', outfile)
    writeLines('<styleUrl>#sh_track</styleUrl>', outfile)
    writeLines('</Pair>', outfile)
    writeLines('</StyleMap>', outfile)
    writeLines('<Style id="sn_track">', outfile)
    writeLines('<IconStyle>', outfile)
    writeLines(paste('<color>',sTrackColour,'</color>',sep=""), outfile)
    writeLines('<Icon>', outfile)
    writeLines('<href>http://maps.google.com/mapfiles/kml/shapes/track.png</href>', outfile)
    writeLines('</Icon>', outfile)
    writeLines('</IconStyle>', outfile)
    writeLines('<LabelStyle>', outfile)
    writeLines('<scale>0.6</scale>', outfile)
    writeLines(paste('<color>',sTrackColour,'</color>',sep=""), outfile)
    writeLines('</LabelStyle>', outfile)
    writeLines('<LineStyle>', outfile)
    writeLines(paste('<color>',sTrackColour,'</color>',sep=""), outfile)
    writeLines('<width>3</width>', outfile)
    writeLines('</LineStyle>', outfile)
    writeLines('</Style>', outfile)
    writeLines('<Style id="sh_track">', outfile)
    writeLines('<IconStyle>', outfile)
    writeLines('<color>ff69deb3</color>', outfile)
    writeLines('<scale>1.33</scale>', outfile)
    writeLines('<Icon>', outfile)
    writeLines('<href>http://maps.google.com/mapfiles/kml/shapes/track.png</href>', outfile)
    writeLines('</Icon>', outfile)
    writeLines('</IconStyle>', outfile)
    writeLines('<LineStyle>', outfile)
    writeLines('<color>ff69deb3</color>', outfile)
    writeLines('<width>3</width>', outfile)
    writeLines('</LineStyle>', outfile)
    writeLines('</Style>', outfile)
    writeLines('<Placemark>', outfile)
    writeLines(paste('<name>',sid,'</name>', sep = ""), outfile)
    writeLines('<styleUrl>#msn_track</styleUrl>', outfile)
    writeLines('<gx:Track>', outfile)
    for(i in 1:nrow(sInputFile2))
    {
      writeLines(paste('<when>',substr(sInputFile2$DATETIME[i], 1, 10),"T",substr(sInputFile2$DATETIME[i], 12, 19),"Z</when>",sep=""), outfile)
    }
    for(i in 1:nrow(sInputFile2))
    {
      writeLines(paste('<gx:coord>', sInputFile2$LONGITUDE[i]," ",sInputFile2$LATITUDE[i]," 0</gx:coord>",sep=""), outfile)
    }
    writeLines('</gx:Track>', outfile)
    writeLines('</Placemark>', outfile)
    writeLines("</Document>", outfile)
    writeLines("</kml>", outfile)
    close(outfile)
  }
