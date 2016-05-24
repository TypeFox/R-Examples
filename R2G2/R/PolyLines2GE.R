PolyLines2GE <-
function(coords, nesting = 0, maxAlt = 1e4, goo = "Lines2GE.kml", colors = "blue", extrude = 0, fill = FALSE, lwd = 1, closepoly = FALSE){

  ## Checking inputs
  if(length(nesting) == 1) nesting = rep(1, nrow(coords))

  if(length(maxAlt) == 1) maxAlt.new = rep(maxAlt, nlevels(as.factor(nesting)))
  if(length(maxAlt) == nlevels(as.factor(nesting))){
    link = cbind(levels(as.factor(nesting)), maxAlt)
    maxAlt.new = link[match(nesting, link[,1]), 2]
    }
  if(length(maxAlt) == nrow(coords)) maxAlt.new = maxAlt

  if(length(lwd) == 1) lwd = rep(lwd, nlevels(as.factor(nesting)))
  if(is.na(match("auto", colors)) == FALSE){
    colors = rainbow(nlevels(as.factor(nesting)))
    }

  ## prepare input data
  input = data.frame(nesting, coords, maxAlt.new)


  ## get nesting params
  nesting = as.factor(nesting)
  groups = as.numeric(levels(nesting))
  nsct = nlevels(nesting) #get number of styles to be used in shapes (i.e. nb of columns in obs)

  ## Initiate KML file
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\"\nxmlns:atom=\"http://www.w3.org/2005/Atom\">\n<Document>\n",file=goo, append = FALSE)
  cat("\t<name>R2G2 - Polygons</name>\n\t<description>Produced using Lines2GE R script</description>\n",file=goo, append = TRUE)

  # a kml file is organised in two blocs: styles (bloc1) and polygons / items (bloc2)
  bloc1=NULL
  bloc2=NULL


  ## prepare bloc 1 (styles)
  #convert colors into GE format
  colors = gsub('#','',colors)
  red = substr(colors,1,2)
  blu = substr(colors,3,4)
  gre = substr(colors,5,6)
  alp = substr(colors,7,8) #alp is the transparency param, can be adapted if needed
  magic = paste(alp,blu,gre,red,sep='')

  
  # feed bloc 1
  for(j in 1:nsct){
    #styles numbering
    if(j > 1){
      nr = j - 2
      } else {
      nr = ''
      }
    # color iteration
    col = magic[j]

    # preparing bloc 1 (kml converting)
    if(fill == TRUE){
      bloc1=c(bloc1, paste("\t<Style id=\"sn_ylw-pushpin",nr,"\">\n\t\t<LineStyle id=\"khLineStyle708\">\n\t\t\t<color>ff000000</color>\n\t\t\t<width>0.5</width>\n\t\t</LineStyle>\n\t\t<PolyStyle id=\"khPolyStyle709\">\n\t\t\t<color>",magic[j],"</color>\n\t\t</PolyStyle>\n\t</Style>\n", sep=''))
      } else {
      bloc1=c(bloc1, paste("\t<Style id=\"sn_ylw-pushpin",nr,"\">\n\t\t<LineStyle id=\"khLineStyle708\">\n\t\t\t<color>",magic[j],"</color>\n\t\t\t<width>",lwd[j],"</width>\n\t\t</LineStyle>\n\t\t<PolyStyle>\n\t\t\t<fill>0</fill>\n\t\t</PolyStyle>\n\t</Style>\n", sep=''))
      }
    }


  ## Compute pie coordinates and feed bloc2 of kml file (here: polygons, make also sure that nesting into subfolders will be handy)
  for(i in groups){ #loop over all groups found into nesting param

    # open group-specific folder into GE kml file
    bloc2 = c(bloc2, paste("<Folder>\n<name>Group_",i,"</name>\n", sep=''))

    # get data of focal polygon
    sgmts = input[input[, 1] == i, 2:4]
    if(closepoly == TRUE){
      sgmts = rbind(sgmts, sgmts[1, ])
      }

    # convert these into kml format, polygon by polygon (one polygon = one pie slice)
      if(i > 1){
	nr = i - 2
	} else {
	nr = ''
	}

    # translate them into kml
    bloc2=c(bloc2, paste("\t<Placemark>\n\t\t<name>",paste("group ", nr, sep=""),"</name>\n\t\t<styleUrl>#sn_ylw-pushpin", nr,"</styleUrl>\n\t\t<Polygon>\n\t\t\t<extrude>",extrude,"</extrude>\n\t\t\t<tessellate>1</tessellate>\n\t\t\t<altitudeMode>absolute</altitudeMode>\n\t\t\t<outerBoundaryIs>\n\t\t\t\t<LinearRing>\n\t\t\t\t\t<coordinates>\n", paste(paste("\t\t\t\t\t", apply(sgmts,1,paste,collapse=','),sep=''),collapse='\n'),"\n\t\t\t\t\t</coordinates>\n\t\t\t\t</LinearRing>\n\t\t\t</outerBoundaryIs>\n\t\t</Polygon>\n\t</Placemark>\n", sep=''))
    bloc2 = c(bloc2, "</Folder>\n")
    } #end of loop over groups

  ## Finalise kml writing
  cat(bloc1, bloc2, file=goo,append = TRUE)
  cat("\t</Document>\n</kml>",file=goo, append = TRUE)
  }

