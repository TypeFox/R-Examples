Axis2GE <-
function(coords, maxVal, maxAlt = 1e+05, lwd = 2, apnd = TRUE, goo = "testAxis.kml"){
  
  # initiate storages
  bloc1 = NULL
  bloc2 = NULL

  ### Compute axis coordinates
  # get ticks along axis and their altitudes
  incrs = floor(log10(maxVal))
  incrs = 10^incrs
  marks = seq(0, maxVal, by = incrs)
  marks = marks[-1]
  alt.marks = maxAlt * marks / maxVal

  # get lon / lat of ticks along axis
  tmp.ticks = rbind(coords[1:2], c(coords[1] + 0.1, coords[2]))

  # prepare input for kml format
  coords.ticks = rbind(c(coords[1:2], 0), c(coords[1:2], maxAlt))
  for(i in 1:length(alt.marks)){
    tmp = cbind(tmp.ticks, rep(alt.marks[i], 2))
    coords.ticks = rbind(coords.ticks, tmp)
    }

  #### Initiate KML file
  if(apnd == FALSE){
    cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\"\nxmlns:atom=\"http://www.w3.org/2005/Atom\">\n<Document>\n",file=goo, append=F)
    cat("\t<name>R2G2 - Y-axis</name>\n\t<description>Produced using Axis2GE R script</description>\n",file=goo, append=T)
    }

  ## Draw lines
  bloc1=c(bloc1, paste("\t<Style id=\"ZStyle\">\n\t\t<LineStyle>\n\t\t\t<width>",lwd,"</width>\n\t\t</LineStyle>\n\t\t<PolyStyle>\n\t\t\t<fill>0</fill>\n\t\t</PolyStyle>\n\t</Style>\n", sep=''))

  for(i in seq(1, nrow(coords.ticks), by = 2)){
    sgmts = coords.ticks[c(i, i+1), ]

    bloc2=c(bloc2, paste("\t<Placemark>\n\t\t<name>Y-Axis</name>\n\t\t<styleUrl>#ZStyle</styleUrl>\n\t\t<Polygon>\n\t\t\t<extrude>0</extrude>\n\t\t\t<tessellate>1</tessellate>\n\t\t\t<altitudeMode>absolute</altitudeMode>\n\t\t\t<outerBoundaryIs>\n\t\t\t\t<LinearRing>\n\t\t\t\t\t<coordinates>\n", paste(paste("\t\t\t\t\t", apply(sgmts,1,paste,collapse=','),sep=''),collapse='\n'),"\n\t\t\t\t\t</coordinates>\n\t\t\t\t</LinearRing>\n\t\t\t</outerBoundaryIs>\n\t\t</Polygon>\n\t</Placemark>\n", sep=''))
    }

  ## add tickmark labels
  bloc1=c(bloc1, "\t<Style id=\"tickmark\">\n\t\t<IconStyle>\n\t\t\t<Icon></Icon>\n\t\t</IconStyle>\n\t</Style>")
  for(i in 1:length(marks)){
    sgmts = c(coords, alt.marks[i])
    labs = marks[i]
    bloc2=c(bloc2, paste("\t<Placemark>\n\t\t<name>",labs,"</name>\n\t\t<styleUrl>#tickmark</styleUrl>\n\t\t<Point>\n\t\t\t<altitudeMode>absolute</altitudeMode>\n\t\t\t\t<coordinates>",paste(paste(sgmts,collapse=','),collapse='\n'),"</coordinates>\n\t\t</Point>\n\t</Placemark>", sep = ''))
    }

  if(apnd == FALSE){
    cat(bloc1, bloc2, file=goo,append=T)
    cat("\t</Document>\n</kml>",file=goo, append=T)
    } else {
    list(bloc1 = bloc1, bloc2 = bloc2)
    }
  }
