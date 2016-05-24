"colorbrewer.palette" <-
  function(nclass = 5,
           type = c("qualitative", "sequential", "diverging"),
           palette = letters[1:18]){
    type <- match.arg(type)
    if(type=="sequential" && (nclass<3 || nclass>9)){
      stop("For 'sequential' type, 'nclass' must be between 3-9")
    }
    if(type=="diverging" && (nclass<3 || nclass>11)){
      stop("For 'diverging' type, 'nclass' must be between 3-11")
    }
    if(type=="qualitative" && (nclass<3 || nclass>12)){
      stop("For 'qualitative' type, 'nclass' must be between 3-12")
    }
    cd <- colorbrewer.data()
          
    palette <- match.arg(palette)
    nclass <- nclass
    cd2 <- cd[cd$type==type & cd$nclass==nclass & cd$palette==palette,]
    rgb(cd2$red, cd2$green, cd2$blue, maxColorValue = 255)
  }
