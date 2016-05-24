
#' Data layer legend
#'
#' This function creates a legend object that can be attached to a data layer.
#' 
#' @inheritParams spLayer.SpatialPoints
#' @param style the style of legend to generate.
#' Can be one of "\code{points}", "\code{lines}", "\code{polygons}", "\code{icons}", "\code{gradient}".
#' @param labels a vector to label each element of the legend.
#' @param title a title which will appear above the legend.
#' If \code{NULL} (default), it will inherit the name of the data layer during the compilation of the map.
#' If \code{NA}, the legend will appear without title.
#' @param position a character string indicating where should the legend be displayed.
#' This must be one of "\code{topleft}", "\code{topright}", "\code{bottomleft}", "\code{bottomright}", "\code{none}"
#' @param cells.range range of gradient values.
#' @param cells.col a vector of any of the three kinds of \R color specifications giving the colors of the gradient.
#' @param cells.alpha a vector of numeric values in \eqn{[0, 1]} setting the gradient opacity.
#' 
#' @return an object \code{layerlegend} which can be passed to an \code{spLayer*} function through the \code{legend} argument.
#' 
#' @export
layerLegend <- function(style,  labels = "", title = NULL, position = "bottomright",
                        png = NULL, size = 5, png.width = NULL, png.height = NULL,
                        stroke.col = 1, stroke.lwd = 1, stroke.lty = -1, stroke.alpha = 1,
                        fill.col = 2, fill.alpha = 0.5,
                        cells.range = c(1, 10), cells.col = heat.colors(12), cells.alpha = 1){
  
  position <- match.arg(position, c("bottomright", "topleft", "topright", "bottomleft", "none"), several.ok = FALSE)
  style <- match.arg(style, c("points", "lines", "polygons", "icons", "gradient"), several.ok = FALSE)
  labels <- as.character(labels)
  tab.max <- length(labels)
  
  labels <- paste0("\"", as.character(labels), "\"")
  
  if(style %in% c("points", "lines", "polygons")){
    stroke.lty <- paste0("\"", as.character(stroke.lty), "\"")
  }
  
  if(style == "points"){
    tab <- list(labels, size,
                stroke.col, stroke.lwd, stroke.lty, stroke.alpha,
                fill.col, fill.alpha)
    tab <- lapply(tab, rep, length.out = tab.max)
    names(tab) <- c("legLabels", "legSize",
                    "legStrokeCol", "legStrokeLwd", "legStrokeLty", "legStrokeAlpha",
                    "legFillCol", "legFillAlpha")
    tab$legWidth <- (tab$legSize + tab$legStrokeLwd) * 2
    tab$legHeight <- (tab$legSize + tab$legStrokeLwd) * 2
  }
  
  if(style == "icons"){
    if(is.null(png.width)){
      png.width <- size
    }
    if(is.null(png.height)){
      png.height <- size
    }
    tab <- list(labels, png, png.width, png.height)
    tab <- lapply(tab, rep, length.out = tab.max)
    names(tab) <- c("legLabels", "legPng", "legPngWidth", "legPngHeight")
    tab$legWidth <- tab$legPngWidth
    tab$legHeight <- tab$legPngHeight
  }
  
  if(style == "lines"){
    tab <- list(labels,
                stroke.col, stroke.lwd, stroke.lty, stroke.alpha)
    tab <- lapply(tab, rep, length.out = tab.max)
    names(tab) <- c("legLabels",
                    "legStrokeCol", "legStrokeLwd", "legStrokeLty", "legStrokeAlpha")
    tab$legWidth <- rep(40, tab.max)
    tab$legHeight <- rep(stroke.lwd + 5, tab.max)
  }
  
  if(style == "polygons"){
    tab <- list(labels,
                stroke.col, stroke.lwd, stroke.lty, stroke.alpha,
                fill.col, fill.alpha)
    tab <- lapply(tab, rep, length.out = tab.max)
    names(tab) <- c("legLabels",
                    "legStrokeCol", "legStrokeLwd", "legStrokeLty", "legStrokeAlpha",
                    "legFillCol", "legFillAlpha")
    tab$legWidth <- rep(30, tab.max) + (2 * tab$legStrokeLwd)
    tab$legHeight <- rep(17, tab.max) + (2 * tab$legStrokeLwd)
  }
  
  if(style == "gradient"){
    breaks <- seq(min(cells.range), max(cells.range), length.out = length(cells.col))
    tab <- list()
    tab$gradcol <- col2rgb(cells.col, alpha = FALSE)
    tab$gradalpha <- rep(cells.alpha, length.out = length(cells.col))
    tab$gradlab <- pretty(breaks, 5)
    tab$gradlab <- tab$gradlab[c(-1, -length(tab$gradlab))]
    tab$gradoffset <- round(seq(0, 100, length.out = length(cells.col)), digits = 2)
    tab$legWidth <- 15
    tab$legHeight <- 150
    tab$gradlabpos <- rev((tab$legHeight/diff(cells.range)) * (tab$gradlab - min(cells.range)))
  }
  
  if(style %in% c("points", "polygons")){
    tab$legStrokeCol <- col2rgb(col2hexa(tab$legStrokeCol, alpha.channel = TRUE, alpha = tab$legStrokeAlpha, charstring = FALSE), alpha = TRUE)
    tab$legStrokeCol["alpha", ] <- tab$legStrokeAlpha
    tab$legFillCol <- col2rgb(col2hexa(tab$legFillCol, alpha.channel = TRUE, alpha = tab$legFillAlpha, charstring = FALSE), alpha = TRUE)
    tab$legFillCol["alpha", ] <- tab$legFillAlpha
    tab$legStrokeAlpha <- tab$legFillAlpha <- NULL
  }
  
  if(style %in% c("lines")){
    tab$legStrokeCol <- col2rgb(col2hexa(tab$legStrokeCol, alpha.channel = TRUE, alpha = tab$legStrokeAlpha, charstring = FALSE), alpha = TRUE)
    tab$legStrokeCol["alpha", ] <- tab$legStrokeAlpha
  }
  
  res <- list(style = style, title = title, position = position, tab = tab)
  class(res) <- "layerlegend"
  return(res)
  
}


processLegend <- function(x, icons.legend.dir, prefix){
  legend.name <- safeVar(paste0("legend", x$layer))
  
  res <- paste0("var ", legend.name, " = L.control({position: '", x$position, "'});\n")
  res <- paste0(res, legend.name, ".onAdd = function (map) {\n\nvar div = L.DomUtil.create('div', 'info legend');\n")
  if(!is.na(x$title)){
    res <- paste0(res, "div.innerHTML += '<h1>", x$title, "</h1><hr>'\n")
  }

  if(x$style != "gradient"){
    tab.names <-  names(x$tab)
    
    if(x$style %in% c("points", "lines", "polygons")){
      x$tab$legStrokeCol <- apply(x$tab$legStrokeCol, 2, function(x) paste0("[", paste(x, collapse = ","), "]"))
    }
    if(x$style %in% c("points", "polygons")){
      x$tab$legFillCol <- apply(x$tab$legFillCol, 2, function(x) paste0("[", paste(x, collapse = ","), "]"))
    }
    if(x$style == "icons"){
      file.copy(from = levels(as.factor(x$tab$legPng)), to = icons.legend.dir)
      x$tab$legPng <- paste0("\"",
                             prefix, "_data/",
                             prefix, "_icons/",
                             prefix, "_legend_icons/",
                             gsub("(.*\\/)([^.]+\\.[[:alnum:]]+$)","\\2", x$tab$legPng), "\"")
    }
    
    tab <- lapply(x$tab, function(x) paste(x, collapse = ","))
    tab <- unlist(tab)
    tab <- paste0(tab.names, " = [", tab, "],", collapse = "\n")
    

    res <- paste0(res, "var ", tab, "\n", "labels = [];\n\n")
    
    if(x$style %in% c("points", "icons")){
      marginLeft <- paste0("var marginLeft = [", paste0((max(x$tab$legWidth) - x$tab$legWidth)/2, collapse = ","), "];\n")
      res <- paste0(res, marginLeft)
      i.sty <- "<i style=\"margin-left:'+ marginLeft[i] +'px\">"
    } else {
      i.sty <- "<i>"
    }
    
    res <- paste0(res, "for (var i = 0; i < legLabels.length; i++) {\nvar lineHeight = legHeight[i];\ndiv.innerHTML +=\n'", i.sty)
    
    if(x$style == "points"){
      svg <- "<svg height=\"'+ legHeight[i] +'\" width=\"'+ legWidth[i] + '\"><circle cx=\"'+ legWidth[i]/2 +'\" cy=\"'+ legHeight[i]/2 +'\" r=\"'+ legSize[i] +'\" style= \"stroke:rgba('+ legStrokeCol[i] +'); stroke-dasharray:'+ legStrokeLty[i] +'; stroke-width:'+ legStrokeLwd[i] +'; fill:rgba('+ legFillCol[i] +')\" /></svg>"
    }
    if(x$style == "icons"){
      svg <- "<img src=\"'+ legPng[i] +'\" style=\"width:'+ legWidth[i] + 'px;height:'+ legHeight[i] +'px\">"
    }
    if(x$style == "lines"){
      svg <- "<svg height=\"'+ legHeight[i] +'\" width=\"'+ legWidth[i] +'\"><line x1=\"0\" y1=\"'+ legHeight[i]/2 +'\" x2=\"'+ legWidth[i] +'\" y2=\"'+ legHeight[i]/2 +'\" style=\"stroke:rgba('+ legStrokeCol[i] +');stroke-dasharray:'+ legStrokeLty[i] +';stroke-width:'+ legStrokeLwd[i] +'\" /></svg>"
    }
    if(x$style == "polygons"){
      svg <- "<svg height=\"'+ legHeight[i] +'\" width=\"'+ legWidth[i] + '\"><rect x=\"'+ legStrokeLwd[i] +'\" y=\"'+ legStrokeLwd[i] +'\" width=\"'+ (legWidth[i]-2*legStrokeLwd[i]) +'\" height=\"'+ (legHeight[i]-2*legStrokeLwd[i]) +'\" style= \"stroke:rgba('+ legStrokeCol[i] +'); stroke-dasharray:'+ legStrokeLty[i] +'; stroke-width:'+ legStrokeLwd[i] +'; fill:rgba('+ legFillCol[i] +')\" />"
    }
    res <- paste0(res, svg)
    if(x$style == "lines"){
      res <- paste0(res, "</i> <p style= \"margin-top: 5px;line-height:' + lineHeight*2 + 'px;\">' + legLabels[i]  + '</p><br>';\n}")
    } else {
      res <- paste0(res, "</i> <p style= \"margin-top: 5px;line-height:' + lineHeight + 'px;\">' + legLabels[i]  + '</p><br>';\n}")
    }
  }
  
  if(x$style == "gradient"){
    tab <- x$tab
    svg <- paste0("<svg height=\"", tab$legHeight ,"\" width=\"", tab$legWidth + 10*max(nchar(tab$gradlab)), "\">")
    svg <- paste0(svg, "<defs>")
    svg <- paste0(svg, "<linearGradient id=\"grad1\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">")
    svg.def <- paste0("<stop offset=\"", tab$gradoffset, "%\" style=\"stop-color:rgb(", tab$gradcol["red",], "," , tab$gradcol["green",], ",", tab$gradcol["blue",], ");stop-opacity:", tab$gradalpha, "\"/>", collapse = "")
    svg <- paste0(svg, svg.def)
    svg <- paste0(svg, "</linearGradient>")
    svg <- paste0(svg, "</defs>")
    svg <- paste0(svg, "<rect x=\"0\" y=\"0\" width=\"", tab$legWidth, "\" height=\"", tab$legHeight, "\" stroke-width=\"0\" fill=\"url(#grad1)\" />")
    svg.tick <- paste0("<line x1=\"", tab$legWidth, "\" y1=\"", tab$gradlabpos, "\" x2=\"", tab$legWidth + 5, "\" y2=\"", tab$gradlabpos, "\" style=\"stroke:rgb(80,80,80);stroke-width:2\"/>", collapse = "")
    svg.lab <- paste0("<text style=\"dominant-baseline: central;font-size: 12px;font-style: normal\"><tspan style=\"dominant-baseline: central;\" x=\"", tab$legWidth + 10, "\" y=\"", tab$gradlabpos, "\">", tab$gradlab, "</tspan></text>", collapse = "")
    svg <- paste0(svg, svg.tick, svg.lab, "</svg>")
    
    res <- paste0(res, "div.innerHTML +=\n'<i>", svg, "</i>'")
  }
  
  res <- paste0(res, "\n\nreturn div;\n};\n\n", legend.name,".addTo(map);\n\n")
  
  res <- paste0(res, "map.on('overlayadd', function (eventLayer) {\nif (eventLayer.name === '", x$layer.name, "') {\n", legend.name, ".addTo(this);\n}\n});\n")
  res <- paste0(res, "map.on('overlayremove', function (eventLayer) {\nif (eventLayer.name === '", x$layer.name, "') {\nthis.removeControl(", legend.name, ");\n}\n});")
  
  return(res)
}

