qrage <- function(links, width = 1400, height = 1000,distance=6000,nodeValue=NULL,nodeColor=NULL,linkColor='#00f',arrowColor='#f00',cut=0.25,textSize=12,linkWidth=c(1,10),nodeSize=c(5,20),linkOpacity=c(0.1,1)) {
  if (class(links) != "data.frame") {
    stop("links must be a data frame class object.");
  }
    
  df <- links;
  colnames(df) <- c("source","target","value");
  
  nv <- NULL;
  tnv <- NULL
  isSetNodeValue = F;
  if(is.null(nodeValue)){
    
  }else{
    isSetNodeValue = T;
    nv <- nodeValue;
    colnames(nv) <- c("name","nodevalue");
    tnv <- data.frame(t(nv$nodevalue));
    colnames(tnv) <- nv$name; 
  }

  td <- NULL
  if(is.null(nodeColor)){
    
  }else{
    colnames(nodeColor) <- c("name","color");
    td <- data.frame(t(as.character(t(nodeColor[,-1]))));
    colnames(td) <- nodeColor$name;
  }
  
  df <- subset(df,df$value>=(max(df$value)*cut));  
  
  data <- list(
    df = df,
    distance = distance,
    nodeColor = td,
    r = nv,
    tr = tnv,
    width = width,
    height = height,
    isSetNodeValue = isSetNodeValue,
    textSize = textSize,
    linkWidth = linkWidth,
    nodeSize = nodeSize,
    linkColor=linkColor,
    linkOpacity=linkOpacity,
    arrowColor=arrowColor
  )  
  
  # create widget
  htmlwidgets::createWidget(
    name = 'qrage',
    data,
    width = width,
    height = height,
    package = 'qrage'
  )
}

qrageOutput <- function(outputId, width = "100%", height = "400px") {
  shinyWidgetOutput(outputId, "qrage", width, height, package = "qrage");
}

renderQrage <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, qrageOutput, env, quoted = TRUE);
}