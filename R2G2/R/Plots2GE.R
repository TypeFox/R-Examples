Plots2GE <-
function(data, center, nesting = 0, customfun, goo = "Plots2GE.kml", testrun = FALSE){
  ### compute location-specific geographical coordinates
  geo = aggregate(center, by = list(nesting), mean)

  ### apply the user-defined plot function (customfun) across all locations defined in "nesting"
  for(i in levels(as.factor(nesting))){
    
    # get location-specific data
    subdata = data[ nesting == i, ]

    if(testrun == TRUE){ #run only a test plot, wihtou kml production
      x11()
      customfun(subdata)
      title(paste("Location", i))

      } else { #produce the actual plots
      png(filename = paste(goo, i, ".png", sep = ''))
      customfun(subdata)   
      title(paste("Location", i))
      dev.off()
      }
    }


  ### produce the kml file that wraps all of that
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\"\nxmlns:atom=\"http://www.w3.org/2005/Atom\">\n<Document>\n",file=goo, append = FALSE)
  cat("\t<name>R2G2 - Custom Plots</name>\n\t<description>Produced using Plots2GE R script</description>\n",file=goo, append = TRUE)

  # a kml file is organised in two blocs: styles (bloc1) and items (bloc2)
  bloc1=NULL
  bloc2=NULL

  ## prepare bloc 1 (styles)
  for(i in levels(as.factor(nesting))){
    # preparing bloc 1 (kml converting)
    bloc1=c(bloc1, paste("\t<Style id=\"graph",i,"\">\n\t\t<IconStyle>\n\t\t\t<scale>7</scale>\n\t\t\t<Icon>\n\t\t\t\t<href>",paste(goo, i, ".png", sep = ''),"</href>\n\t\t\t</Icon>\n\t\t</IconStyle>\n\t</Style>\n", sep=''))
    }

  ## prepare bloc 1 (points)
  for(i in levels(as.factor(nesting))){
    xcoord = geo[match(i, geo[,1]) , 2]
    ycoord = geo[match(i, geo[,1]) , 3]
    # preparing bloc 1 (kml converting)
    bloc2=c(bloc2, paste("\t<Placemark>\n\t\t<name>group", i,"</name>\n\t\t\t<styleUrl>#graph", i, "</styleUrl>\n\t\t<Point>\n\t\t\t<coordinates>", xcoord, ",", ycoord,",0</coordinates>\n\t\t</Point>\n\t</Placemark>", sep=''))
    }

  ## Finalise kml writing
  cat(bloc1, bloc2, file=goo, append = TRUE)
  cat("\t</Document>\n</kml>",file=goo, append = TRUE)
  }

