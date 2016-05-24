GenerateAnimationKMLFile <-
function(sInputResidenceFile,
                                     sInputNonResidenceFile,
                                     sInputPointsFile,
                                     sOutputFile,
                                     sReceiverColour)
{
  # This function constructs an animation file from the given NonResidence, Residence and Points files.
  
  # open non-residence file
  NonResidenceFile <- read.csv(sInputNonResidenceFile)
  iNonResidenceCount <- dim(NonResidenceFile)[1]
  
  # open residence file
  ResidenceFile <- read.csv(sInputResidenceFile)
  iResidenceCount <- dim(ResidenceFile)[1]
  
  # open points file
  PointFile <- read.csv(sInputPointsFile)
  iPointCount <- dim(PointFile)[1]
  
  # init First/Last time
  sFirstTIME <- "2050-01-01 00:00:00"
  sLastTime <- "2000-01-01 00:00:00"
  
  # open the output file for writing
  # we are creating a new file
  outfile <- file(sOutputFile,"wt")
  # write KML header
  writeLines("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",outfile)
  writeLines("<kml xmlns=\"http://earth.google.com/kml/2.2\">",outfile)
  writeLines("<Document>",outfile)
  writeLines(paste(GenerateTabs(1),"<name>",sOutputFile,"</name>",sep=""),outfile)
  # write styles
  
  # create a style for receivers
  writeLines(paste(GenerateTabs(1),"<Style id=\"receiver\">",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<IconStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<color>",sReceiverColour,"</color>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<scale>0.333333</scale>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<Icon>",sep=""),outfile)
  writeLines(paste(GenerateTabs(4),"<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"</Icon>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"</IconStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<LabelStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<color>00ffffff</color>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<scale>1</scale>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"</LabelStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(1),"</Style>",sep=""),outfile)
  # create a style for waypoints
  writeLines(paste(GenerateTabs(1),"<Style id=\"waypoint\">",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<IconStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<color>00ffffff</color>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<scale>0.333333</scale>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<Icon>",sep=""),outfile)
  writeLines(paste(GenerateTabs(4),"<href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"</Icon>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"</IconStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<LabelStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<color>00ffffff</color>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<scale>1</scale>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"</LabelStyle>",sep=""),outfile)
  writeLines(paste(GenerateTabs(1),"</Style>",sep=""),outfile)
  
  # write folder header
  writeLines(paste(GenerateTabs(1),"<Folder>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<name>Transmitter Residences</name>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<open>1</open>",sep=""),outfile)
  writeLines(paste(GenerateTabs(2),"<description>Exported from V-Track Software. Copyright University of Queensland 2009-2012. Authors: Ross Dwyer, Matthew Watts, Hamish Campbell and Craig Franklin.  Each point represents a residence event by a tagged animal and lines represent a non-residence event.</description>",sep=""),outfile)
  
  # write placemarks
  writeLines(paste(GenerateTabs(2),"<Folder>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<name>Residence</name>",sep=""),outfile)
  # traverse residence file, creating a placemark for each record
  for (i in 1:iResidenceCount)
  {
    sSTARTTIME <- as.character(ResidenceFile[i,1])
    sENDTIME <- as.character(ResidenceFile[i,2])
    sTRANSMITTERID <- as.character(ResidenceFile[i,4])
    sReceiver <- as.character(ResidenceFile[i,5])
    
    # update First/Last time
    if (sSTARTTIME < sFirstTIME)
      sFirstTIME <- sSTARTTIME
    if (sENDTIME < sFirstTIME)
      sFirstTIME <- sENDTIME
    if (sSTARTTIME > sLastTime)
      sLastTime <- sSTARTTIME
    if (sENDTIME > sLastTime)
      sLastTime <- sENDTIME
    
    sCoordinates <- ReturnReceiverCoordinates(sReceiver,iPointCount,PointFile)
    
    if (sCoordinates != "")
    {
      writeLines(paste(GenerateTabs(3),"<Placemark>",sep=""),outfile)
      writeLines(paste(GenerateTabs(4),"<name>",sTRANSMITTERID,"</name>",sep=""),outfile)
      writeLines(paste(GenerateTabs(4),"<TimeSpan>",sep=""),outfile)
      writeLines(paste(GenerateTabs(5),"<begin>",ConvertTimeFormat_KML(sSTARTTIME,0,0),"</begin>",sep=""),outfile)
      writeLines(paste(GenerateTabs(5),"<end>",ConvertTimeFormat_KML(sENDTIME,0,0),"</end>",sep=""),outfile)
      writeLines(paste(GenerateTabs(4),"</TimeSpan>",sep=""),outfile)
      writeLines(paste(GenerateTabs(4),"<styleUrl>#receiver</styleUrl>",sep=""),outfile)
      writeLines(paste(GenerateTabs(4),"<Point>",sep=""),outfile)
      writeLines(paste(GenerateTabs(5),"<coordinates>",sCoordinates,",0</coordinates>",sep=""),outfile)
      writeLines(paste(GenerateTabs(4),"</Point>",sep=""),outfile)
      writeLines(paste(GenerateTabs(3),"</Placemark>",sep=""),outfile)
    }
  }
  writeLines(paste(GenerateTabs(2),"</Folder>",sep=""),outfile)
  
  # traverse non-residence file, creating a series of placemarks for each record
  writeLines(paste(GenerateTabs(2),"<Folder>",sep=""),outfile)
  writeLines(paste(GenerateTabs(3),"<name>Non-Residence</name>",sep=""),outfile)
  
  for (i in 1:iNonResidenceCount)
  {
    sSTARTTIME <- as.character(NonResidenceFile[i,1])
    sENDTIME <- as.character(NonResidenceFile[i,2])
    sTRANSMITTERID <- as.character(NonResidenceFile[i,4])
    sReceiver <- as.character(NonResidenceFile[i,5])
    sRECEIVER2 <- as.character(NonResidenceFile[i,6])
    
    # update First/Last time
    if (sSTARTTIME < sFirstTIME)
      sFirstTIME <- sSTARTTIME
    if (sENDTIME < sFirstTIME)
      sFirstTIME <- sENDTIME
    if (sSTARTTIME > sLastTime)
      sLastTime <- sSTARTTIME
    if (sENDTIME > sLastTime)
      sLastTime <- sENDTIME
    
    iCoordinateIndex1 <- ReturnReceiverIndex(sReceiver,iPointCount,PointFile)
    iCoordinateIndex2 <- ReturnReceiverIndex(sRECEIVER2,iPointCount,PointFile)
    
    if (iCoordinateIndex1 > -1)
    {
      if (iCoordinateIndex2 > -1)
      {
        if (iCoordinateIndex1 != iCoordinateIndex2)
        {
          iNumberOfSteps <- iCoordinateIndex1 - iCoordinateIndex2
          
          writeLines(paste(GenerateTabs(3),"<Placemark>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"<name>",sTRANSMITTERID,"</name>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"<TimeSpan>",sep=""),outfile)
          writeLines(paste(GenerateTabs(5),"<begin>",ConvertTimeFormat_KML(sSTARTTIME,0,0),"</begin>",sep=""),outfile)
          writeLines(paste(GenerateTabs(5),"<end>",ConvertTimeFormat_KML(sENDTIME,0,0),"</end>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"</TimeSpan>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"<styleUrl>#waypoint</styleUrl>",sep=""),outfile)
          
          writeLines(paste(GenerateTabs(4),"<LineString>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"<tessellate>1</tessellate>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"<altitudeMode>clampToGround</altitudeMode>",sep=""),outfile)
          writeLines(paste(GenerateTabs(5),"<coordinates>",sep=""),outfile)
          
          
          if (iNumberOfSteps > 0)
          {
            # moving backward through the array
            for (i in iCoordinateIndex1:iCoordinateIndex2) # Step -1
            {
              sPointFile2 <- as.character(PointFile[i,2])
              sPointFile3 <- as.character(PointFile[i,3])
              writeLines(paste(GenerateTabs(6),sPointFile3,",",sPointFile2,",0",sep=""),outfile)
            }
          }else{
            # moving forward through the array
            for (i in iCoordinateIndex1:iCoordinateIndex2)
            {
              sPointFile2 <- as.character(PointFile[i,2])
              sPointFile3 <- as.character(PointFile[i,3])
              writeLines(paste(GenerateTabs(6),sPointFile3,",",sPointFile2,",0",sep=""),outfile)
            }
          }                     
          writeLines(paste(GenerateTabs(5),"</coordinates>",sep=""),outfile)
          writeLines(paste(GenerateTabs(4),"</LineString>",sep=""),outfile)                     
          writeLines(paste(GenerateTabs(3),"</Placemark>",sep=""),outfile)
        }
      }
    }
  }
  writeLines(paste(GenerateTabs(2),"</Folder>",sep=""),outfile)
  
  # write KML footer
  writeLines(paste(GenerateTabs(1),"</Folder>",sep=""),outfile)
  writeLines("</Document>",outfile)
  writeLines("</kml>",outfile)
  close(outfile)
}
