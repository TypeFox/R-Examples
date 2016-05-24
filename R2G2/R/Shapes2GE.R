Shapes2GE <-
function(center, nesting = 0, goo = "Shapes2GE.kml", nedges = 20, orient = 0, colors = "blue", maxAlt = 1e4, radius = 5e4){

  ## Checking inputs
  if(!is.null(nrow(center))){
    if(length(nesting) == 1) nesting = rep(1, nrow(center))
    if(length(nedges) == 1) nedges = rep(nedges, nrow(center))
    if(length(orient) == 1) orient = rep(orient, nrow(center))
    if(length(maxAlt) == 1) maxAlt = rep(maxAlt, nrow(center))
    if(length(radius) == 1) radius = rep(radius, nrow(center))
    if(length(colors) == 1){
      palt = colors
      colors = rep(colors, nrow(center))
      }
    if(length(colors) > 1){
      palt = colors
      }
    if(is.na(match("auto", colors)) == FALSE){
      palt = rev(rainbow(nlevels(as.factor(nesting))))
      colors = as.factor(nesting)
      levels(colors) = palt
      colors = as.character(colors)
      }
    }

  ## prepare input data
  input = data.frame(nesting, center, nedges, radius, orient, maxAlt, colors)


  ## get nesting params
  nesting = as.factor(nesting)
  groups = levels(nesting)
  nsct = nlevels(nesting) #get number of styles to be used in shapes (i.e. nb of columns in obs)

  ## Initiate KML file
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\"\nxmlns:atom=\"http://www.w3.org/2005/Atom\">\n<Document>\n",file=goo, append = FALSE)
  cat("\t<name>R2G2 - Shapes</name>\n\t<description>Produced using Shapes2GE R script</description>\n",file=goo, append = TRUE)

  # a kml file is organised in two blocs: styles (bloc1) and polygons / items (bloc2)
  bloc1=NULL
  bloc2=NULL


  ## prepare bloc 1 (styles)
  #convert colors into GE format
  colors = gsub('#','',palt)
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
    bloc1=c(bloc1, paste("\t<Style id=\"sn_ylw-pushpin",nr,"\">\n\t\t<LineStyle id=\"khLineStyle708\">\n\t\t\t<color>ff000000</color>\n\t\t\t<width>0.5</width>\n\t\t</LineStyle>\n\t\t<PolyStyle id=\"khPolyStyle709\">\n\t\t\t<color>",magic[j],"</color>\n\t\t</PolyStyle>\n\t</Style>\n", sep=''))
    }


  cnt = 1
  ## Compute pie coordinates and feed bloc2 of kml file (here: polygons, make also sure that nesting into subfolders will be handy)
  for(i in groups){ #loop over all groups found into nesting param

    # open group-specific folder into GE kml file
    bloc2 = c(bloc2, paste("<Folder>\n<name>Group_",i,"</name>\n", sep=''))

    # get data of focal group
    input.grp = input[ input[, 1] == i,]
    
    # Compute coordinates of each shape, from that group
    for(pie in 1:nrow(input.grp)){
      
      # get data of focal pie and comput coords
      input.pie = input.grp[pie, ]
      coords = GetEdges(center = as.numeric(input.pie[2:3]), 
			nedges = as.numeric(input.pie[4]), 
			radius = as.numeric(input.pie[5]), 
			orient = as.numeric(input.pie[6])) #we use the as.numeric to enhence robustness against data idiosyncraties (makes the script more robust here)

      if(any(is.na(coords))) next()

      # convert these into kml format, polygon by polygon (one polygon = one pie slice)
      if(cnt > 1){
	nr = cnt - 2
	} else {
	nr = ''
	}

      #getting actual data of polygon
      sgmts = cbind(coords, rep(input.grp[pie, 7], nrow(coords)))
      #sgmts = sgmts[nrow(sgmts):1, ] #little trick: you want the ordering of segments to be counterclockwise for proper plotting in GE!!

      # translate them into kml
      bloc2=c(bloc2, paste("\t<Placemark>\n\t\t<name>",paste("group ", nr, sep=""),"</name>\n\t\t<styleUrl>#sn_ylw-pushpin", nr,"</styleUrl>\n\t\t<Polygon>\n\t\t\t<extrude>1</extrude>\n\t\t\t<tessellate>1</tessellate>\n\t\t\t<altitudeMode>absolute</altitudeMode>\n\t\t\t<outerBoundaryIs>\n\t\t\t\t<LinearRing>\n\t\t\t\t\t<coordinates>\n", paste(paste("\t\t\t\t\t", apply(sgmts,1,paste,collapse=','),sep=''),collapse='\n'),"\n\t\t\t\t\t</coordinates>\n\t\t\t\t</LinearRing>\n\t\t\t</outerBoundaryIs>\n\t\t</Polygon>\n\t</Placemark>\n", sep=''))

    } #end of loop over pies, within groups
    bloc2 = c(bloc2, "</Folder>\n")
    cnt = cnt + 1
  } #end of loop over groups

  ## Finalise kml writing
  cat(bloc1, bloc2, file=goo,append = TRUE)
  cat("\t</Document>\n</kml>",file=goo, append = TRUE)
  }

