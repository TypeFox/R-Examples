# Currently this function could only parse svg files created by the cairo
# graphics library, typically from svg() in the grDevices package (R >= 2.14.0
# required for Windows OS), and CairoSVG() in the Cairo package.
parseSVG = function(file.name)
{
    xmlParse = XML::xmlParse;
    newXMLNamespace = XML::newXMLNamespace;
    xmlRoot = XML::xmlRoot;
    getNodeSet = XML::getNodeSet;
    xmlAttrs = XML::xmlAttrs;
    xmlParent = XML::xmlParent;
    xmlName = XML::xmlName;
    xmlChildren = XML::xmlChildren;
    xmlSApply = XML::xmlSApply;
    
    svgFile = xmlParse(file.name);
    # Don't forget the name space!
    newXMLNamespace(xmlRoot(svgFile), "http://www.w3.org/2000/svg", "svg");
    # Find the first <g> child of <svg>
    pathRoot = getNodeSet(svgFile, "/svg:svg/svg:g");
    if(!length(pathRoot)) stop(sprintf("Failed in parsing file '%s'", file.name));
    pathRoot = pathRoot[[1]];
    
    # Default style for a <path> node
    defaultStyle = c("stroke" = "none",
                     "stroke-width" = "1",
                     "stroke-linecap" = "butt",
                     "stroke-linejoin" = "miter",
                     "stroke-miterlimit" = "4",
                     "stroke-opacity" = "1",
                     "fill" = "rgb(0%,0%,0%)",
                     "fill-rule" = "nonzero",
                     "fill-opacity" = "1");
    
    # Handle <path> style in named vector
    parseStyle = function(style)
    {
        if(is.null(style)) return(NULL);
        s = unlist(strsplit(style, ";"));
        val = strsplit(s, ":");
        result = sapply(val, function(x) x[2]);
        names(result) = sapply(val, function(x) x[1]);
        return(result);
    }
    # Update the attributes in "old" style with the values in "new"
    # "old" must contain "new"
    updateStyle = function(old, new)
    {
        if(is.null(new)) return(old);
        result = old;
        result[names(new)] = new;
        return(result);
    }
    # Iteratively update the style from parent nodes
    updateStyleUpward = function(node)
    {
        style = xmlAttrs(node)["style"];
        if(is.na(style)) style = NULL;
        style = parseStyle(style);
        style = updateStyle(defaultStyle, style);
        parentNode = xmlParent(node);
        # Recursively search the parent
        while(!is.null(parentNode))
        {
            parentStyle = xmlAttrs(parentNode)["style"];
            if(is.null(parentStyle) || is.na(parentStyle)) parentStyle = NULL;
            parentStyle = parseStyle(parentStyle);
            style = updateStyle(style, parentStyle);
            parentNode = xmlParent(parentNode);
        }
        return(style);
    }
    # Parse <path> and <use> nodes into structured lists
    #
    # <path style="" d="">   =====>   style=..., d=..., x=0, y=0
    #
    # <use xlink:href="#glyph0-0" x="63.046875" y="385.921875"/>
    # =====>
    # style=..., d=..., x=63.046875, y=385.921875
    #
    parseNode = function(node)
    {
        if(xmlName(node) == "use")
        {
            attrs = xmlAttrs(node);
            refID = sub("#", "", attrs["href"]);
            refPathNode = getNodeSet(svgFile, sprintf("//*[@id='%s']/svg:path", refID))[[1]];
            style = updateStyleUpward(refPathNode);
            style = updateStyle(style, updateStyleUpward(node));
            d = xmlAttrs(refPathNode)["d"];
            x = xmlAttrs(node)["x"];
            y = xmlAttrs(node)["y"];
        } else if(xmlName(node) == "path") {
            style = updateStyleUpward(node);
            d = xmlAttrs(node)["d"];
            x = y = 0;
        } else return(NULL);
        xy = as.numeric(c(x, y));
        names(d) = NULL;
        names(xy) = NULL;
        return(list(style = style, d = d, xy = xy));
    }
    # Flatten nodes
    # <g>
    #   <use />
    #   <use />
    #   <use />
    # </g>
    #
    # =====>
    #
    # <use />
    # <use />
    # <use />
    expandNode = function(node)
    {
        children = xmlChildren(node);
        res = if(!length(children)) node else children;
        return(res);
    }
    nodes = unlist(xmlSApply(pathRoot, expandNode));
    names(nodes) = NULL;
    paths = lapply(nodes, parseNode);
    path.is.null = sapply(paths, is.null);
    paths[path.is.null] = NULL;
    if(!length(paths)) stop("Unknown child node of '/svg/g'");
    return(paths);
}


#' Convert a sequence of SVG files to SWF file
#'
#' Given the file names of a sequence of SVG files, this function could
#' convert them into a Flash file (.swf).
#'
#' This function uses the XML package in R and a subset of librsvg
#' (\url{http://librsvg.sourceforge.net/}) to parse the SVG file, and
#' uses the Ming library (\url{http://www.libming.org/}) to
#' implement the conversion. Currently this function supports SVG files
#' created by \code{\link[grDevices]{svg}()} in the \pkg{grDevices}
#' package, and \code{\link[Cairo]{CairoSVG}()} in the
#' \pkg{Cairo} package.
#' @param input the file names of the SVG files to be converted
#' @param output the name of the output SWF file
#' @param bgColor background color of the output SWF file
#' @param interval the time interval (in seconds) between animation frames
#' @return The name of the generated SWF file if successful.
#' @export
#' @author Yixuan Qiu <\email{yixuan.qiu@@cos.name}>
#' @examples \dontrun{
#' if(capabilities("cairo")) {
#'   olddir = setwd(tempdir())
#'   svg("Rplot%03d.svg", onefile = FALSE)
#'   set.seed(123)
#'   x = rnorm(5)
#'   y = rnorm(5)
#'   for(i in 1:100) {
#'       plot(x <- x + 0.1 * rnorm(5), y <- y + 0.1 * rnorm(5),
#'            xlim = c(-3, 3), ylim = c(-3, 3), col = "steelblue",
#'            pch = 16, cex = 2, xlab = "x", ylab = "y")
#'   }
#'   dev.off()
#'   output = svg2swf(sprintf("Rplot%03d.svg", 1:100), interval = 0.1)
#'   swf2html(output)
#'   setwd(olddir)
#' }
#' }
#'
svg2swf = function(input, output = "movie.swf", bgColor = "white", interval = 1)
{
    # Use XML package
    # if(!require(XML))
    #     stop("svg2swf() requires XML package");
    
    xmlParse = XML::xmlParse;
    xmlAttrs = XML::xmlAttrs;
    xmlRoot = XML::xmlRoot;
    
    if(!is.character(input))
        stop("'input' must be a character vector naming the input SVG files");
    
    bg = col2rgb(bgColor, alpha = FALSE);
    bg = as.integer(bg);
    
    if(!all(file.exists(input))) stop("one or more input files do not exist");
    
    filesData = lapply(input, parseSVG);
    firstFile = xmlParse(input[1]);
    size = xmlAttrs(xmlRoot(firstFile))["viewBox"];
    size = as.numeric(unlist(strsplit(size, " ")));
    
    outfile = normalizePath(output, mustWork = FALSE);
    .Call("svg2swf", filesData, outfile, size,
          bg, as.numeric(interval), PACKAGE = "R2SWF");
    
    message("SWF file created at ", outfile);
    invisible(output);
}
