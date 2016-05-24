# Purpose        : Parsing "SoilProfileCollection" objects to KML
# Maintainer     : Dylan Beaudette (debeaudette@ucdavis.edu)
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Pierre Roudier (pierre.roudier@landcare.nz);
# Status         : not tested yet
# Note           : plots either a histogram or blocks (horizons);


## TODO: finish and integrate this into kml_layer.SoilProfileCollection
.SPC_to_images <- function(obj) {
	#require(Cairo)
	
	# make container dir
	cd <- 'images'
	dir.create(cd)
	
	# pre-compute some values outside of the loop
	ids <- profile_id(obj)
	fn <- paste('soil-', ids, '.png', sep='')
	fp <- paste(cd, fn, sep='/')
	
	# make links to images
	image.links <- paste('<img src="', fp, '" width=150, height=300, border=0>', sep='')
	
	# iterate over profiles and save images
	for(i in seq_along(ids)) {
# 		CairoPNG(file=fp[i], height=300, width=150)
		png(filename=fp[i], height=300, width=150)
		par(mar=c(0, 1, 0.5, 1.5))
		plot(obj[i, ], cex.names=0.75, width=0.4)
		dev.off()
	}
	
	# clean-up
	unlink(cd, recursive=TRUE)
}


kml_layer.SoilProfileCollection <- function(
  obj,
  var.name,
  var.min = 0,  
  var.scale,
  site_names = profile_id(obj),
  method = c("soil_block", "depth_function")[1],
  block.size = 100,  # in m 
  color.name,
  z.scale = 1,
  x.min,
  max.depth = 300,
  plot.points = TRUE,
  LabelScale = get("LabelScale", envir = plotKML.opts)*.7,
  IconColor = "#ff0000ff",
  shape = paste(get("home_url", envir = plotKML.opts), "circlesquare.png", sep=""),
  outline = TRUE,
  visibility = TRUE,
  extrude = TRUE,
  tessellate = TRUE,
  altitudeMode = "relativeToGround",
  camera.distance = .01,
  tilt = 90, 
  heading = 0, 
  roll = 0,
  metadata = NULL,
  html.table = NULL,
  plot.scalebar = TRUE,
  scalebar = paste(get("home_url", envir = plotKML.opts), "soilprofile_scalebar.png", sep=""),
  ...) {
  
  # deconstruct object
  h <- horizons(obj)
  s <- site(obj)
  sp <- as(obj, 'SpatialPoints')
  
  # TH: this function at the moment works only with numeric variables:
  if(method=="depth_functions" & !is.numeric(h[, var.name])) {
  	stop('numeric variable required')
  }

  # get our invisible file connection from custom evnrionment
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  # check the projection:
  prj.check <- check_projection(sp, control = TRUE) 

  # Trying to reproject data if the check was not successful
  if (!prj.check){
  	sp <- reproject(sp)
 	}
  
  # extract coordinates... for scale estimation
  LON <- as.vector(coordinates(sp)[,1])
  LAT <- as.vector(coordinates(sp)[,2])

	# library to estimate scaling factor
  if(requireNamespace("fossil", quietly = TRUE)){
    # convert meters to decimal degrees:
    new.ll <- fossil::new.lat.long(long = mean(LON), lat = mean(LAT), bearing = 90, distance = block.size/1000)
    block.size <- new.ll[2] - mean(LON)
  } else {
    block.size <- 1/240
    warning("Failed estimating 'block.size'. Install and load package 'fossil'")
  }
  if(missing(x.min)){
  	x.min = block.size/100
 	}

  if(missing(var.scale)) {   # scaling factor in x direction (estimate automatically)
    var.range <- range(h[, var.name], na.rm = TRUE, finite = TRUE)
    var.scale <- 0.003/diff(var.range)
  	if(missing(var.min)){
  		var.min <- var.range[1] - (diff(var.range)/100)
 		}
  }

  # Parsing the call for aesthetics
  aes <- kml_aes(obj, ...)

  # Read the relevant aesthetics
  # altitudeMode <- aes[["altitudeMode"]]
  balloon <- aes[["balloon"]]

  # Folder and name of the points folder
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", paste(class(obj)), parent = pl1)

  pl2c = newXMLNode("Folder", parent=pl1)
  pl3c <- newXMLNode("name", "sites", parent = pl2c)

  # Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txt <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txt, parent=pl1)
  }
  message("Writing KML...")

  if(plot.points==TRUE){
    pl2b = newXMLNode("Folder", parent=pl1)
    pl3b <- newXMLNode("name", var.name, parent = pl2b)
  }

  # Calculate 3D coordinates
  # ==========================  

  # for each site:
  prof.na <- list(NULL)
  coords <- as.list(rep("", length(obj)))
  coords.pol <- as.list(rep("", length(obj)))
  coordsB <- as.list(rep("", length(obj)))
  points_names <- as.list(rep("", length(obj)))
  soil_color <- as.list(rep("", length(obj)))
  
  for(i.site in 1:length(obj)) {
  
  # select columns of interest / mask out NA horizons:
  prof.na[[i.site]] <- which(h[[idname(obj)]] == site_names[i.site])
  if(is.integer(prof.na[[i.site]])){
    xval <- h[prof.na[[i.site]], var.name]
    htop <- h[prof.na[[i.site]], horizonDepths(obj)[1]]
    hbot <- h[prof.na[[i.site]], horizonDepths(obj)[2]]
    if(!(all(!is.na(htop))&all(!is.na(hbot))&!all(is.na(xval)))){
      warning(paste("Site '", site_names[i.site], "' missing upper and/or lower limits for one of the horizons", sep=""), immediate.=TRUE)
    }
  
    if(plot.points==TRUE){
    points_names[[i.site]] <- signif(xval, 3)
  
    # if no color is specified use standard colors:
    if(missing(color.name)){
      pal <- colorRampPalette(c("chocolate4", "darkgoldenrod1", "cornsilk")) 
      soil_color[[i.site]] <- col2kml(pal(length(xval)))
    } 
    else { 
      soil_color[[i.site]] <- col2kml(h[prof.na[[i.site]], color.name])
    }
  
    # horizon centre:
    Z <- max.depth - (htop+(hbot-htop)/2)
    
    if(method=="soil_block"){   
    	X <- LON[i.site]
    	Y <- LAT[i.site]+((block.size/2)*sqrt(2)+x.min)
    }
      else {
    	X <- round(LON[i.site]+var.scale*(xval-var.min), 6)
    	Y <- rep(LAT[i.site], length(prof.na[[i.site]]))
    }
  
    coords[[i.site]] <- paste(X, ',', Y, ',', z.scale*Z, collapse='\n ', sep = "")
  }
  
    # horizon polygons:
    if(method=="soil_block"){
      Xp <- list(NULL)
      Yp <- list(NULL)
      Zp <- list(NULL)
      XYZp <- list(NULL)
    
    if(length(xval)>0){
    for(i in 1:length(xval)){
      Xp[[i]] <- c(LON[i.site]-block.size/2, LON[i.site]+block.size/2, LON[i.site]+block.size/2, LON[i.site]-block.size/2, LON[i.site]-block.size/2)
      Yp[[i]] <- c(LAT[i.site]+((block.size/2)*sqrt(2)+x.min), LAT[i.site]+((block.size/2)*sqrt(2)+x.min), LAT[i.site]+((block.size/2)*sqrt(2)+x.min), LAT[i.site]+((block.size/2)*sqrt(2)+x.min), LAT[i.site]+((block.size/2)*sqrt(2)+x.min))
      Zp[[i]] <- c(max.depth-htop[i], max.depth-htop[i], max.depth-hbot[i], max.depth-hbot[i], max.depth-htop[i])
      XYZp[[i]] <-  paste(Xp[[i]], ',', Yp[[i]], ',', z.scale*Zp[[i]], collapse='\n ', sep = "")
    }
    coords.pol[[i.site]] <- unlist(XYZp)

    # skyscraper:
    XB <- c(LON[i.site]-block.size/2, LON[i.site]+block.size/2, LON[i.site]+block.size/2, LON[i.site]-block.size/2, LON[i.site]-block.size/2)
    YB <- c(LAT[i.site]-(block.size/2)*sqrt(2), LAT[i.site]-(block.size/2)*sqrt(2), LAT[i.site]+(block.size/2)*sqrt(2), LAT[i.site]+(block.size/2)*sqrt(2), LAT[i.site]-(block.size/2)*sqrt(2))
    ZB <- rep(max.depth, 5)
    coordsB[[i.site]] <- paste(XB, ',', YB, ',', z.scale*ZB, collapse='\n ', sep = "") 
    }
  }
  
  else {
    Xp <- round(c(as.vector(t(matrix(rep(LON[i.site]+var.scale*(xval-var.min), 2), ncol=2))), rep(LON[i.site], 2), LON[i.site]+var.scale*(xval[1]-var.min)), 6)
    Yp <- rep(LAT[i.site], length(Xp)+1)
    Zp <- c(as.vector(t(matrix(c(max.depth-htop, max.depth-hbot), ncol=2))), max.depth-hbot[length(hbot)], max.depth-htop[1], max.depth-htop[1])
    coords.pol[[i.site]] <- paste(Xp, ',', Yp, ',', z.scale*Zp, collapse='\n ', sep = "")
  }

  } 
 }
  
  ## Parse ATTRIBUTE TABLE (for each placemark):
  if ((balloon == TRUE | class(balloon) %in% c('character','numeric')) & ("horizons" %in% slotNames(obj))){
     html.table <- .df2htmltable(h[unlist(prof.na),]) 
  }

  if(plot.points==TRUE){
  # Point styles
  # =====================    
  selp <- !sapply(prof.na, function(x){length(x)==0})
  lp <- length(unlist(prof.na))
  
  txts <- sprintf('<Style id="pnt%s"><LabelStyle><scale>%.1f</scale></LabelStyle><IconStyle><color>%s</color><scale>%s</scale><Icon><href>%s</href></Icon></IconStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', paste(1:lp), rep(LabelScale, lp), rep(IconColor, lp), rep(LabelScale, lp), rep(shape, lp))
  parseXMLAndAdd(txts, parent=pl2b)  
 
  # Coordinates of points 
  # ========================== 
  txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><description><![CDATA[%s]]></description><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></Point></Placemark>', paste(unlist(points_names)), paste(1:lp), html.table, rep(as.numeric(extrude), lp), rep(altitudeMode, lp), unlist(strsplit(unlist(coords[selp]), "\n")))   
  parseXMLAndAdd(txtc, parent=pl2b)
  }
  
  # Scale bar - 59 x 2272 pixels PNG!
  # ===========================
  if(plot.scalebar==TRUE & method=="depth_function"){
  scalebar.size = 200*59/2272
  X1 = scalebar.size/2; X2 = -scalebar.size/2; X3 = scalebar.size/2; X4 = -scalebar.size/2
  Y1 = 0; Y2 = 0; Y3 = 0; Y4 = 0
  ALT1 = max.depth; ALT2 = max.depth; ALT3 = max.depth-200; ALT4 = max.depth-200 
  coords.dae <- matrix(c(X=c(X1, X2, X3, X4), Z=c(ALT1, ALT2, ALT3, ALT4), Y=c(Y1, Y2, Y3, Y4)), ncol=3)
  
  # make a COLLADA file:
  makeCOLLADA.rectangle(filename = "scalebar.dae", coords = coords.dae, href = scalebar)  
  
  # locate the scale bar in space:
  txtv <- sprintf('<Camera><longitude>%.6f</longitude><latitude>%.6f</latitude><altitude>%f</altitude><heading>%.1f</heading><tilt>%.1f</tilt><roll>%.1f</roll></Camera>', LON[1], LAT[1]-camera.distance, max.depth, heading, tilt, roll)
  txtm <- sprintf('<Placemark><name>monolith</name><Model id="%s"><altitudeMode>relativeToGround</altitudeMode><Location><longitude>%.5f</longitude><latitude>%.5f</latitude><altitude>%.0f</altitude></Location><Orientation><heading>%.1f</heading><tilt>%.1f</tilt><roll>%.1f</roll></Orientation><Scale><x>1</x><y>1</y><z>1</z></Scale><Link><href>%s</href><refreshMode>"once"</refreshMode></Link><ResourceMap><Alias><targetHref>%s</targetHref><sourceHref>%s</sourceHref></Alias></ResourceMap></Model></Placemark>', paste("scalebar", 1:length(LON), sep="_"), LON-x.min*10, LAT, rep(max.depth+100, length(LON)), rep(heading, length(LON)), rep(tilt, length(LON)), rep(roll, length(LON)), rep("scalebar.dae", length(LON)), rep(scalebar, length(LON)), rep(scalebar, length(LON))) 
  parseXMLAndAdd(txtm, parent=pl1)
  }
  
  if(method=="soil_block"){
  # Polygon styles
  # =====================
  txts <- sprintf('<Style id="poly%s"><PolyStyle><color>%s</color><fill>1</fill><outline>%s</outline></PolyStyle></Style>', paste(1:lp), paste(unlist(soil_color[selp])), rep(as.numeric(outline), lp))
  parseXMLAndAdd(txts, parent=pl2c)
 
  # Write Polygons (for each horizon)
  # ==========================  
  txtp <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><visibility>%s</visibility><Polygon><tessellate>%s</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', paste(unlist(points_names)), paste(1:lp), rep(as.numeric(visibility), lp), rep(as.numeric(tessellate), lp), rep(altitudeMode, lp), paste(unlist(coords.pol[selp])))
  parseXMLAndAdd(txtp, parent=pl2c)
  # tessellated Block:
  txtB <- sprintf('<Placemark><name>%s</name><visibility>%s</visibility><Polygon><tessellate>1</tessellate><extrude>1</extrude><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', site_names, rep(as.numeric(visibility), length(site_names)), rep(altitudeMode, length(site_names)), paste(unlist(coordsB)))
  parseXMLAndAdd(txtB, parent=pl2c)    
  }
 
    else {
  # Write Polygons (for each site) 
  txt <- sprintf('<Placemark><name>%s</name><visibility>%s</visibility><Polygon><tessellate>%s</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', site_names, rep(as.numeric(visibility), length(site_names)), rep(as.numeric(tessellate), length(site_names)), rep(altitudeMode, length(site_names)), paste(unlist(coords.pol)))     
  parseXMLAndAdd(txt, parent=pl2c)
  }


  # Add a Camera view location:
  # ==========================    
  txtv <- sprintf('<Camera><longitude>%.6f</longitude><latitude>%.6f</latitude><altitude>%f</altitude><heading>%.1f</heading><tilt>%.1f</tilt><roll>%.1f</roll></Camera>', LON[1], LAT[1]-camera.distance, max.depth, heading, tilt, roll) 
  parseXMLAndAdd(txtv, parent=pl1)

  # save results: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)
}

    
setMethod("kml_layer", "SoilProfileCollection", kml_layer.SoilProfileCollection)

# end of script;
