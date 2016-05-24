# Purpose        : Export of (spatial) metadata
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Michael Blaschek (blaschek@geographie.uni-kiel.de); 
# Dev Status     : Pre-Alpha
# Note           : Based on the US gov sp metadata standards [http://www.fgdc.gov/metadata/csdgm/], which can be converted to "ISO 19139" XML schema;


## Generate a SLD file (using the default legend):
# [http://docs.geoserver.org/stable/en/user/styling/sld-introduction.html]

metadata2SLD.Spatial <- function(obj, ...){
  
  if(xmlValue(obj@xml[["//formcont"]]) == "SpatialPixelsDataFrame"){
    metadata2SLD.SpatialPixels(obj, ...)
  }
  # ...
  ## to be continued
  else {
    stop("Format_Information_Content field in 'obj@xml' must specify an applicable sp class.")
  }
}

metadata2SLD.SpatialPixels <- function(
  obj,  # SpatialMetadata
  Format_Information_Content = xmlValue(obj@xml[["//formcont"]]),
  obj.name = normalizeFilename(deparse(substitute(obj))),
  sld.file = set.file.extension(obj.name, ".sld"),
  Citation_title = xmlValue(obj@xml[["//title"]]),
  ColorMap_type = "intervals",
  opacity = 1,
  brw.trg = 'Greys', # color scheme, according to www.colorbrewer2.org; default to 'Greys' 
  target.var, # target variable, used to calculate the class-intervals
  ...
){
  l1 = newXMLNode("StyledLayerDescriptor", attrs=c("xsi:schemaLocation" = "http://www.opengis.net/sld StyledLayerDescriptor.xsd", version="1.0.0"), namespaceDefinitions=c("http://www.opengis.net/sld", "xsi" = "http://www.w3.org/2001/XMLSchema-instance", "ogc" = "http://www.opengis.net/ogc", "gml" = "http://www.opengis.net/gml"))
  l2 <- newXMLNode("NamedLayer", parent = l1)
  l3 <- newXMLNode("Name", paste(Citation_title, "(", Format_Information_Content, ")"), parent = l2)
  l3b <- newXMLNode("UserStyle", parent = l2)
  l4 <- newXMLNode("Title", paste(obj.name, "style", sep="_"), parent = l3b)
  l4b <- newXMLNode("FeatureTypeStyle", parent = l3b)
  l5 <- newXMLNode("Rule", parent = l4b)
  l6 <- newXMLNode("RasterSymbolizer", parent = l5)
  l7 <- newXMLNode("ColorMap", attrs=c(type=ColorMap_type), parent = l6)
  ## MB: if no target variable is provided, the original ColorMapEntries (defined by spMetadata-call) are used..
  if(missing(target.var)) { 
    txt <- sprintf('<ColorMapEntry color="%s" quantity="%.2f" label="%s" opacity="%.1f"/>', obj@palette@color, obj@palette@bounds[-1], obj@palette@names, rep(opacity, length(obj@palette@color)))
  } else {
    mm <- classIntervals(target.var, ...)
    brew.p <- brewer.pal(n = length(mm$brks) - 1, name = brw.trg)
    op <- findColours(mm, pal = brew.p, under = 'under', over = 'over', between = '-', cutlabels = F)
    txt <- sprintf('<ColorMapEntry color="%s" quantity="%.2f" label="%s" opacity="%.1f"/>', attr(op, 'palette'), mm$brks[-1], attr(attr(op, 'table'), 'dimnames')[[1]], rep(opacity, length(mm$brks[-1])))  
  }
  parseXMLAndAdd(txt, l7)
  saveXML(l1, sld.file)
}


## Write the metadata dataframe to Geonetwork MEF format as specified at [http://geonetwork-opensource.org/manuals/2.6.3/developer/mef/]
#
# metadata2MEF <- function(
#    xml,  # metadata slot

#    )
#    {

# }    

# connect all methods and classes:
setMethod("metadata2SLD", "SpatialMetadata", metadata2SLD.Spatial)

# end of script;
