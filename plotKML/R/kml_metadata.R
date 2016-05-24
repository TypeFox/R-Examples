# Purpose        : Insertion of metadata into a KML file
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : Based on the US gov sp metadata standards [http://www.fgdc.gov/metadata/csdgm/];

## summary for an object of type "SpatialMetadata":
.summary.metadata <- function(object, sel, fix.enc = TRUE, full.names = "", delim.sign){
    
    if(full.names == ""){     
      full.names = read.csv(system.file("mdnames.csv", package="plotKML"))     
    }
    
    nx <- unlist(xmlToList(object@xml))
    ## convert to a table:
    if(!missing(delim.sign)){
      nxn <- gsub("\\.", delim.sign, attr(nx, "names"))
    } else {
      nxn <- attr(nx, "names")
    }
    met <- data.frame(metadata=nxn, value=paste(nx), stringsAsFactors = FALSE)
    ## add more friendly names:
    metm <- merge(x=met, y=full.names[,c("metadata","field.names")], by="metadata", all.x=TRUE)
          
    ## selected columns:
    if(missing(sel)) {
      sel = get("metadata_sel", envir = plotKML.opts)
    }
    selm <- data.frame(metadata = sel, order.no=1:length(sel))
    md <- merge(selm, metm, by="metadata", all.x=TRUE)
    md <- md[order(md$order.no),] 
    
    ## fix encoding:
    if(fix.enc==TRUE){
      md <- data.frame(lapply(md, iconv, to = "UTF8"))
    }
    
    return(md)
}

setMethod("summary", signature(object = "SpatialMetadata"), definition=.summary.metadata)


## Write metadata to a KML file:                             
setMethod("kml_metadata", signature(obj = "SpatialMetadata"), function(obj, cwidth = 150, twidth = 500, asText = FALSE){
      
    md <- .summary.metadata(obj)    
    if(!nrow(md)==0){
    # write to html:
    l1 <- newXMLNode("table", attrs=c(width=twidth, border="0", cellspacing="5", cellpadding="10"))
    l2 <- newXMLNode("caption", "Spatial metadata summary:", parent = l1)
    txt <- sprintf('<tr><th width="%.0f" scope="col"><div align="left"><strong>%s</strong></div></th><th scope="col"><div align="left"><strong>%s</strong></div></th></tr>', rep(cwidth, nrow(md)), paste(md$field.names), paste(md$value))
    parseXMLAndAdd(txt, l1)
    
    if(asText==TRUE){ 
      md.txt <- saveXML(l1)
      return(md.txt)
    }
    else {  
      return(l1) 
    }
    } else {
      warning("Metadata contains no elements. See 'spMetadata-method' for more info.")
    }
})

# end of script;