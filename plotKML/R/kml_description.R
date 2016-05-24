# Purpose        : Insertion of table data into a KML file (HTML text)
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Pierre Roudier (pierre.roudier@landcare.nz); 
# Dev Status     : Pre-Alpha
# Note           : Based on the US gov sp metadata standards [http://www.fgdc.gov/metadata/csdgm/];


## Write metadata to a KML file:                             
kml_description <- function(
    x,   # data.frame
    iframe = NULL,
    caption = "Object summary",
    fix.enc = TRUE,
    cwidth = 150,
    twidth = 300,
    delim.sign = "_",
    asText = FALSE
    ){

   if(!is.null(iframe)){
    l1 <- newXMLNode("iframe", attrs=c(src = iframe, width=twidth*2))
    }
   else {

    if(ncol(x)>2){ warning("Only first two columns will be printed in the description tag") }

    # write to html:
    l1 <- newXMLNode("table", attrs=c(width=twidth, border="0", cellspacing="5", cellpadding="10"))
    l2 <- newXMLNode("caption", caption, parent = l1)
    txt <- sprintf('<tr><th width="%.0f" scope="col"><div align="right"><strong>%s</strong>: </div></th><th scope="col"><div align="left">%s</div></th></tr>', rep(cwidth, nrow(x)), paste(x[,1]), paste(x[,2]))
    parseXMLAndAdd(txt, l1)
    }
    
    if(asText==TRUE){ 
      md.txt <- saveXML(l1)
      return(md.txt)
    }
    else {  
      return(l1) 
    }
    
}

# end of script;