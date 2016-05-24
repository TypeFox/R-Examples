
explicitMathVariant <- function(fontfamily, fontface) {
    currentFonts <- getSVGFonts()
    stackname <- fontStackFromFontFamily(fontfamily, currentFonts)
    switch(stackname,
           sans=switch(fontface,
                       plain="sans-serif",
                       bold="bold-sans-serif",
                       italic=,
                       oblique="sans-serif-italic",
                       bold.italic="sans-serif-bold-italic"),
           serif=switch(fontface,
                        plain="normal",
                        bold="bold",
                        italic=,
                        oblique="italic",
                        bold.italic="bold-italic"),
           mono="monospace",
           switch(fontface,
                  plain="sans-serif",
                  bold="bold-sans-serif",
                  italic=,
                  oblique="sans-serif-italic",
                  bold.italic="sans-serif-bold-italic"))
}

unicode <-
c("\u0020", "\u0021", "\u2200", "\u0023", "\u2203", 
"\u0025", "\u0026", "\u220B", "\u0028", "\u0029", "\u2217", 
"\u002B", "\u002C", "\u2212", "\u002E", "\u002F", "\u0030", 
"\u0031", "\u0032", "\u0033", "\u0034", "\u0035", "\u0036", 
"\u0037", "\u0038", "\u0039", "\u003A", "\u003B", "\u003C", 
"\u003D", "\u003E", "\u003F", "\u2245", "\u0391", "\u0392", 
"\u03A7", "\u0394", "\u0395", "\u03A6", "\u0393", "\u0397", 
"\u0399", "\u03D1", "\u039A", "\u039B", "\u039C", "\u039D", 
"\u039F", "\u03A0", "\u0398", "\u03A1", "\u03A3", "\u03A4", 
"\u03A5", "\u03C2", "\u03A9", "\u039E", "\u03A8", "\u0396", 
"\u005B", "\u2234", "\u005D", "\u22A5", "\u005F", "\uF8E5", 
"\u03B1", "\u03B2", "\u03C7", "\u03B4", "\u03B5", "\u03C6", 
"\u03B3", "\u03B7", "\u03B9", "\u03D5", "\u03BA", "\u03BB", 
"\u03BC", "\u03BD", "\u03BF", "\u03C0", "\u03B8", "\u03C1", 
"\u03C3", "\u03C4", "\u03C5", "\u03D6", "\u03C9", "\u03BE", 
"\u03C8", "\u03B6", "\u007B", "\u007C", "\u007D", "\u223C", 
"\u20AC", "\u03D2", "\u2032", "\u2264", "\u2044", "\u221E", 
"\u0192", "\u2663", "\u2666", "\u2665", "\u2660", "\u2194", 
"\u2190", "\u2191", "\u2192", "\u2193", "\u00B0", "\u00B1", 
"\u2033", "\u2265", "\u00D7", "\u221D", "\u2202", "\u2022", 
"\u00F7", "\u2260", "\u2261", "\u2248", "\u2026", "\uF8E6", 
"\uF8E7", "\u21B5", "\u2135", "\u2111", "\u211C", "\u2118", 
"\u2297", "\u2295", "\u2205", "\u2229", "\u222A", "\u2283", 
"\u2287", "\u2284", "\u2282", "\u2286", "\u2208", "\u2209", 
"\u2220", "\u2207", "\uF6DA", "\uF6D9", "\uF6DB", "\u220F", 
"\u221A", "\u22C5", "\u00AC", "\u2227", "\u2228", "\u21D4", 
"\u21D0", "\u21D1", "\u21D2", "\u21D3", "\u25CA", "\u2329", 
"\uF8E8", "\uF8E9", "\uF8EA", "\u2211", "\uF8EB", "\uF8EC", 
"\uF8ED", "\uF8EE", "\uF8EF", "\uF8F0", "\uF8F1", "\uF8F2", 
"\uF8F3", "\uF8F4", "\u232A", "\u222B", "\u2320", "\uF8F5", 
"\u2321", "\uF8F6", "\uF8F7", "\uF8F8", "\uF8F9", "\uF8FA", 
"\uF8FB", "\uF8FC", "\uF8FD", "\uF8FE")

# See ~/Research/Rstuff/SVG/PlotMath/greek.R
greek <-
structure(c("\u03B1", "\u03B2", "\u03B3", "\u03B4", "\u03B5", 
"\u03B6", "\u03B7", "\u03B8", "\u03B9", "\u03BA", "\u03BB", 
"\u03BC", "\u03BD", "\u03BE", "\u03BF", "\u03C0", "\u03C1", 
"\u03C2", "\u03C3", "\u03C4", "\u03C5", "\u03D5", "\u03C7", 
"\u03C8", "\u03C9", "\u0391", "\u0392", "\u0393", "\u0394", 
"\u0395", "\u0396", "\u0397", "\u0398", "\u0399", "\u039A", 
"\u039B", "\u039C", "\u039D", "\u039E", "\u039F", "\u03A0", 
"\u03A1", "\u03A2", "\u03A3", "\u03A4", "\u03A5", "\u03A6", 
"\u03A7", "\u03A8", "\u03A9"), .Names = c("alpha", "beta", 
"gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", 
"kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", 
"", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega", 
"Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", 
"Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", 
"Pi", "Rho", "", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", 
"Omega"))

symbolNames <-
    c("..."="\u2026",
      cdots="\u22EF",
      ldots="\u2026",
      greek,
      theta1="\u03D1",
      vartheta="\u03D1",
      phi1="\u03C6",
      sigma1="\u03C2",
      varsigma="\u03C2",
      omega1="\u03D6",
      Upsilon1="\u03D2",
      aleph="\u05D0",
      infinity="\u221E",
      partialdiff="\u2202",
      nabla="\u2207",
      degree="\u00B0",
      minute="\u2032",
      second="\u2033")

# The general idea with each of these mml*() functions is to
# create a single MathML element.
# This means that, if the output is a collection of several
# elements, we wrap the whole collection in an <mrow>
mmlJuxta <- function(e, fontfamily, fontface, svgdev) {
    mrow <- newXMLNode("mrow", parent = svgDevParent(svgdev))
    svgDevChangeParent(mrow, svgdev)

    e <- e[-1]
    lapply(e, function(x) {
        toMML(x, fontfamily, fontface, svgdev)
    })

    svgDevChangeParent(xmlParent(mrow), svgdev)
}

mmlBinOp <- function(e, fontfamily, fontface, op, svgdev) {
    mrow <- newXMLNode("mrow", parent = svgDevParent(svgdev))
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    newXMLNode("mo", parent = svgDevParent(svgdev),
               newXMLTextNode(op))
    toMML(e[[3]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mrow), svgdev)
}

mmlParen <- function(e, fontfamily, fontface, svgdev) {
    mfenced <- newXMLNode("mfenced", parent = svgDevParent(svgdev))
    mrow <- newXMLNode("mrow", parent = mfenced)
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mfenced), svgdev)
}

mmlBrace <- function(e, fontfamily, fontface, svgdev) {
    mfenced <- newXMLNode("mfenced", parent = svgDevParent(svgdev),
                          attrs = list(open = "", close = ""))
    mrow <- newXMLNode("mrow", parent = mfenced)    
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mfenced), svgdev)
}

delimiters <- c(lfloor="\u230A",
                rfloor="\u230B",
                lceil="\u2308",
                rceil="\u2309")

convertDelim <- function(delim) {
    if (delim %in% names(delimiters))
        delimiters[delim]
    else
        delim
}

mmlGroup <- function(e, fontfamily, fontface, svgdev) {
    # e[[2]] and e[[4]] are the delimiters
    if (length(e) < 4)
        stop("Invalid plotmath group()")
    delim1 <- convertDelim(as.character(e[[2]]))
    delim2 <- convertDelim(as.character(e[[4]]))

    mfenced <- newXMLNode("mfenced", parent = svgDevParent(svgdev),
                          attrs = list(open = delim1, close = delim2))
    mrow <- newXMLNode("mrow", parent = mfenced)    
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[3]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mfenced), svgdev)
}

mmlSup <- function(e, fontfamily, fontface, svgdev) {
    msup <- newXMLNode("msup", parent = svgDevParent(svgdev))
    svgDevChangeParent(msup, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    toMML(e[[3]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(msup), svgdev)
}

mmlSub <- function(e, fontfamily, fontface, svgdev) {
    msub <- newXMLNode("msub", parent = svgDevParent(svgdev))
    svgDevChangeParent(msub, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    toMML(e[[3]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(msub), svgdev)
}

mmlSqrt <- function(e, fontfamily, fontface, svgdev) {
    if (length(e) > 2) {
        mroot <- newXMLNode("mroot", parent = svgDevParent(svgdev))
        svgDevChangeParent(mroot, svgdev)

        toMML(e[[2]], fontfamily, fontface, svgdev)
        toMML(e[[3]], fontfamily, fontface, svgdev)

        svgDevChangeParent(xmlParent(mroot), svgdev)
    } else {
        msqrt <- newXMLNode("msqrt", parent = svgDevParent(svgdev))
        svgDevChangeParent(msqrt, svgdev)

        toMML(e[[2]], fontfamily, fontface, svgdev)

        svgDevChangeParent(xmlParent(msqrt), svgdev)
    }
}

mmlFont <- function(e, fontfamily, fontface, svgdev) {
    toMML(e[[2]], fontfamily, fontface, svgdev)
}

mmlStyle <- function(e, fontfamily, fontface, style, svgdev) {
    displaystyle <- switch(style,
                           display="true",
                           "false")
    scriptlevel <- switch(style,
                          display=0,
                          text=0,
                          script=1,
                          scriptscript=2)

    mstyle <- newXMLNode("mstyle", parent = svgDevParent(svgdev),
                         attrs = list(displaystyle = displaystyle,
                                      scriptlevel = scriptlevel))
    svgDevChangeParent(mstyle, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mstyle), svgdev)
}

mmlSymbol <- function(e, fontfamily, fontface, svgdev) {
    newXMLNode("mtext", parent = svgDevParent(svgdev),
               attrs = list(mathvariant =
                   explicitMathVariant(fontfamily, fontface)),
               newXMLTextNode(unicode[as.integer(charToRaw(e[[2]])) - 31]))
}

mmlCSL <- function(e, fontfamily, fontface, svgdev) {
    mfenced <- newXMLNode("mfenced", parent = svgDevParent(svgdev),
                         attrs = list(open = "", close = ""))
    svgDevChangeParent(mfenced, svgdev)

    sapply(e[-1], toMML, fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mfenced), svgdev)
}

mmlAccent <- function(e, fontfamily, fontface, accent, svgdev) {
    mover <- newXMLNode("mover", parent = svgDevParent(svgdev),
                         attrs = list(accent = "true",
                                      align = "center"))
    mrow <- newXMLNode("mrow", parent = mover)
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    newXMLNode("mo", parent = mover,
               attrs = list(stretchy = "false"),
               newXMLTextNode(accent))

    svgDevChangeParent(xmlParent(mover), svgdev)
}

mmlWideAccent <- function(e, fontfamily, fontface, accent, svgdev) {
    mover <- newXMLNode("mover", parent = svgDevParent(svgdev),
                         attrs = list(accent = "true",
                                      align = "center"))
    mrow <- newXMLNode("mrow", parent = mover)
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    newXMLNode("mo", parent = mover,
               attrs = list(stretchy = "true"),
               newXMLTextNode(accent))

    svgDevChangeParent(xmlParent(mover), svgdev)
}

mmlUnderline <- function(e, fontfamily, fontface, svgdev) {
    # NOTE: <mstack> and <msline> are currently not supported
    #       by Mozilla-based browsers (2011-11-21)
    munder <- newXMLNode("munder", parent = svgDevParent(svgdev))
    mrow <- newXMLNode("mrow", parent = munder)
    svgDevChangeParent(mrow, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    newXMLNode("mo", parent = munder,
               attrs = list(stretchy = "true"),
               newXMLTextNode("\u00AF"))

    svgDevChangeParent(xmlParent(munder), svgdev)
}

mmlSpace <- function(e, fontfamily, fontface, svgdev) {
    mrow <- newXMLNode("mrow", parent = svgDevParent(svgdev))
    svgDevChangeParent(mrow, svgdev)

    if (length(e) > 2) {
        toMML(e[[2]], fontfamily, fontface, svgdev)
        newXMLNode("mtext", parent = svgDevParent(svgdev),
                   attrs = list(mathvariant =
                                explicitMathVariant(fontfamily, fontface)),
                   newXMLTextNode("\u00A0"))
        toMML(e[[3]], fontfamily, fontface, svgdev)
    } else {
        newXMLNode("mtext", parent = svgDevParent(svgdev),
                   attrs = list(mathvariant =
                                explicitMathVariant(fontfamily, fontface)),
                   newXMLTextNode("\u00A0"))
        toMML(e[[2]], fontfamily, fontface, svgdev)
    }

    svgDevChangeParent(xmlParent(mrow), svgdev)
}

mmlPhantom <- function(e, fontfamily, fontface, svgdev) {
    mphantom <- newXMLNode("mphantom", parent = svgDevParent(svgdev))
    svgDevChangeParent(mphantom, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mphantom), svgdev)
}

mmlFrac <- function(e, fontfamily, fontface, svgdev, lwd="medium") {
    mfrac <- newXMLNode("mfrac", parent = svgDevParent(svgdev),
                        attrs = list(linethickness = lwd))
    svgDevChangeParent(mfrac, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)
    toMML(e[[3]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mfrac), svgdev)
}

mmlBigOp <- function(e, fontfamily, fontface, svgdev, op=NULL) {
    mrow <- newXMLNode("mrow", parent = svgDevParent(svgdev))

    # When checking for is.null(op),
    # Either specify op special character or format the first
    # element of the expression
    if (length(e) < 3) {
        if (is.null(op)) {
            opmrow <- newXMLNode("mrow", parent = mrow)
            svgDevChangeParent(opmrow, svgdev)
            toMML(e[[1]], fontfamily, fontface, svgdev)
            svgDevChangeParent(xmlParent(opmrow), svgdev)
            newXMLNode("mtext", parent = opmrow,
                       attrs = list(mathvariant =
                           explicitMathVariant(fontfamily, fontface)),
                       newXMLTextNode("\u00A0"))
        } else {
            newXMLNode("mo", parent = mrow,
                       newXMLTextNode(op))
        }

        svgDevChangeParent(mrow, svgdev)
        toMML(e[[2]], fontfamily, fontface, svgdev)
    } else if (length(e) < 4) {
        munder <- newXMLNode("munder", parent = mrow)

        if (is.null(op)) {
            opmrow <- newXMLNode("mrow", parent = munder)
            svgDevChangeParent(opmrow, svgdev)
            toMML(e[[1]], fontfamily, fontface, svgdev)
            svgDevChangeParent(xmlParent(opmrow), svgdev)
            newXMLNode("mtext", parent = opmrow,
                       attrs = list(mathvariant =
                           explicitMathVariant(fontfamily, fontface)),
                       newXMLTextNode("\u00A0"))
        } else {
            newXMLNode("mo", parent = munder,
                       newXMLTextNode(op))
        }

        svgDevChangeParent(munder, svgdev)
        toMML(e[[3]], fontfamily, fontface, svgdev)
        svgDevChangeParent(xmlParent(munder), svgdev)
        toMML(e[[2]], fontfamily, fontface, svgdev)
    } else {
        munderover <- newXMLNode("munderover", parent = mrow)

        if (is.null(op)) {
            opmrow <- newXMLNode("mrow", parent = munderover)
            svgDevChangeParent(opmrow, svgdev)
            toMML(e[[1]], fontfamily, fontface, svgdev)
            svgDevChangeParent(xmlParent(opmrow), svgdev)
            newXMLNode("mtext", parent = opmrow,
                       attrs = list(mathvariant =
                           explicitMathVariant(fontfamily, fontface)),
                       newXMLTextNode("\u00A0"))
        } else {
            newXMLNode("mo", parent = munderover,
                       newXMLTextNode(op))
        }

        svgDevChangeParent(munderover, svgdev)
        toMML(e[[3]], fontfamily, fontface, svgdev)
        toMML(e[[4]], fontfamily, fontface, svgdev)
        svgDevChangeParent(xmlParent(munderover), svgdev)
        toMML(e[[2]], fontfamily, fontface, svgdev)
    }

    svgDevChangeParent(xmlParent(mrow), svgdev)
}

mmlFun <- function(e, fontfamily, fontface, svgdev) {
    mrow <- newXMLNode("mrow", parent = svgDevParent(svgdev))
    mtext <- newXMLNode("mtext", parent = mrow,
                        attrs = list(mathvariant = explicitMathVariant(fontfamily, fontface)),
                        newXMLTextNode(e[[1]]))
    mfenced <- newXMLNode("mfenced", parent = mrow)
    svgDevChangeParent(mfenced, svgdev)

    toMML(e[[2]], fontfamily, fontface, svgdev)

    svgDevChangeParent(xmlParent(mrow), svgdev)
}

toMML <- function(x, fontfamily, fontface, svgdev, ...) {
    UseMethod("toMML")
}

toMML.numeric <- function(x, fontfamily, fontface, svgdev, ...) {
    newXMLNode("mn", parent = svgDevParent(svgdev),
               newXMLTextNode(as.character(x)))
}

toMML.character <- function(x, fontfamily, fontface, svgdev, ...) {
    newXMLNode("mtext", parent = svgDevParent(svgdev),
               attrs = list(mathvariant =
                   explicitMathVariant(fontfamily, fontface)),
               newXMLTextNode(x))
}

toMML.name <- function(x, fontfamily, fontface, svgdev, ...) {
    # Convert special names
    if (as.character(x) %in% names(symbolNames))
        x <- symbolNames[as.character(x)]
    # R does NOT automatically italicize symbols
    newXMLNode("mtext", parent = svgDevParent(svgdev),
               attrs = list(mathvariant =
                   explicitMathVariant(fontfamily, fontface)),
               newXMLTextNode(x))
}

# A "language" object may have class "call" or "(" or "{"
"toMML.(" <- function(x, fontfamily, fontface, svgdev, ...) {
    toMML.call(x, fontfamily, fontface, svgdev, ...)
}

"toMML.{" <- function(x, fontfamily, fontface, svgdev, ...) {
    toMML.call(x, fontfamily, fontface, svgdev, ...)
}

funCallToMML <- function(x, fontfamily, fontface, svgdev) {
    funName <- as.character(x[[1]])
    switch(funName,
           "+"=,
           "/"=mmlBinOp(x, fontfamily, fontface, funName, svgdev),
           "-"=mmlBinOp(x, fontfamily, fontface, "\u2212", svgdev),
           "*"=mmlBinOp(x, fontfamily, fontface, "\u2062", svgdev),
           "%+-%"=mmlBinOp(x, fontfamily, fontface, "\u00B1", svgdev),
           "%/%"=mmlBinOp(x, fontfamily, fontface, "\u00F7", svgdev),
           "%*%"=mmlBinOp(x, fontfamily, fontface, "\u00D7", svgdev),
           "%.%"=mmlBinOp(x, fontfamily, fontface, "\u22C5", svgdev),
           "["=mmlSub(x, fontfamily, fontface, svgdev),
           "^"=mmlSup(x, fontfamily, fontface, svgdev),
           "paste"=mmlJuxta(x, fontfamily, fontface, svgdev),
           "sqrt"=mmlSqrt(x, fontfamily, fontface, svgdev),
           "("=mmlParen(x, fontfamily, fontface, svgdev),
           "{"=mmlBrace(x, fontfamily, fontface, svgdev),
           "=="=mmlBinOp(x, fontfamily, fontface, "=", svgdev),
           "!="=mmlBinOp(x, fontfamily, fontface, "\u2260", svgdev),
           "<"=mmlBinOp(x, fontfamily, fontface, "<", svgdev),
           "<="=mmlBinOp(x, fontfamily, fontface, "\u2264", svgdev),
           ">"=mmlBinOp(x, fontfamily, fontface, ">", svgdev),
           ">="=mmlBinOp(x, fontfamily, fontface, "\u2265", svgdev),
           "%~~%"=mmlBinOp(x, fontfamily, fontface, "\u2248", svgdev),
           "%=~%"=mmlBinOp(x, fontfamily, fontface, "\u2245", svgdev),
           "%==%"=mmlBinOp(x, fontfamily, fontface, "\u2261", svgdev),
           "%prop%"=mmlBinOp(x, fontfamily, fontface, "\u221D", svgdev),
           "plain"=mmlFont(x, fontfamily, "plain", svgdev),
           "bold"=mmlFont(x, fontfamily, "bold", svgdev),
           "italic"=mmlFont(x, fontfamily, "italic", svgdev),
           "bolditalic"=mmlFont(x, fontfamily, "bold.italic", svgdev),
           "symbol"=mmlSymbol(x, fontfamily, fontface, svgdev),
           "list"=mmlCSL(x, fontfamily, fontface, svgdev),
           "%subset%"=mmlBinOp(x, fontfamily, fontface, "\u2282", svgdev),
           "%subseteq%"=mmlBinOp(x, fontfamily, fontface, "\u2286", svgdev),
           "%notsubset%"=mmlBinOp(x, fontfamily, fontface, "\u2284", svgdev),
           "%supset%"=mmlBinOp(x, fontfamily, fontface, "\u2283", svgdev),
           "%supseteq%"=mmlBinOp(x, fontfamily, fontface, "\u2287", svgdev),
           "%notsupset%"=mmlBinOp(x, fontfamily, fontface, "\u2285", svgdev),
           "%in%"=mmlBinOp(x, fontfamily, fontface, "\u2208", svgdev),
           "%notin%"=mmlBinOp(x, fontfamily, fontface, "\u2209", svgdev),
           "hat"=mmlAccent(x, fontfamily, fontface, "\u005E", svgdev),
           "tilde"=mmlAccent(x, fontfamily, fontface, "\u007E", svgdev),
           "dot"=mmlAccent(x, fontfamily, fontface, "\u02D9", svgdev),
           "ring"=mmlAccent(x, fontfamily, fontface, "\u02DA", svgdev),
           # Used "macron"
           "bar"=mmlAccent(x, fontfamily, fontface, "\u00AF", svgdev),
           # FIXME:  these are just normal accents positioned as limits
           "widehat"=mmlWideAccent(x, fontfamily, fontface,
                                   "\u005E", svgdev),
           "widetilde"=mmlWideAccent(x, fontfamily, fontface,
                                     "\u007E", svgdev),
           "%<->%"=mmlBinOp(x, fontfamily, fontface, "\u2194", svgdev),
           "%->%"=mmlBinOp(x, fontfamily, fontface, "\u2192", svgdev),
           "%<-%"=mmlBinOp(x, fontfamily, fontface, "\u2190", svgdev),
           "%up%"=mmlBinOp(x, fontfamily, fontface, "\u2191", svgdev),
           "%down%"=mmlBinOp(x, fontfamily, fontface, "\u2193", svgdev),
           "%<=>%"=mmlBinOp(x, fontfamily, fontface, "\u21D4", svgdev),
           "%=>%"=mmlBinOp(x, fontfamily, fontface, "\u21D2", svgdev),
           "%<=%"=mmlBinOp(x, fontfamily, fontface, "\u21D0", svgdev),
           "%dblup%"=mmlBinOp(x, fontfamily, fontface, "\u21D1", svgdev),
           "%dbldown%"=mmlBinOp(x, fontfamily, fontface, "\u21D3", svgdev),
           "displaystyle"=mmlStyle(x, fontfamily, fontface, "display", svgdev),
           "textstyle"=mmlStyle(x, fontfamily, fontface, "text", svgdev),
           "scriptstyle"=mmlStyle(x, fontfamily, fontface, "script", svgdev),
           "scriptscriptstyle"=mmlStyle(x, fontfamily, fontface,
                                        "scriptscript", svgdev),
           "underline"=mmlUnderline(x, fontfamily, fontface, svgdev),
           "~"=mmlSpace(x, fontfamily, fontface, svgdev),
           "phantom"=mmlPhantom(x, fontfamily, fontface, svgdev),
           "over"=,
           "frac"=mmlFrac(x, fontfamily, fontface, svgdev),
           "atop"=mmlFrac(x, fontfamily, fontface, lwd="0em", svgdev),
           "sum"=mmlBigOp(x, fontfamily, fontface, svgdev, "\u2211"),
           "prod"=mmlBigOp(x, fontfamily, fontface, svgdev, "\u220F"),
           "integral"=mmlBigOp(x, fontfamily, fontface, svgdev, "\u222B"),
           "union"=mmlBigOp(x, fontfamily, fontface, svgdev, "\u22C3"),
           "intersect"=mmlBigOp(x, fontfamily, fontface, svgdev, "\u22C2"),
           "prod"=mmlBigOp(x, fontfamily, fontface, svgdev, "\u220F"),
           "lim"=mmlBigOp(x, fontfamily, fontface, svgdev),
           "min"=mmlBigOp(x, fontfamily, fontface, svgdev),
           "inf"=mmlBigOp(x, fontfamily, fontface, svgdev),
           "sup"=mmlBigOp(x, fontfamily, fontface, svgdev),
           "group"=mmlGroup(x, fontfamily, fontface, svgdev),
           "bgroup"=mmlGroup(x, fontfamily, fontface, svgdev),
           mmlFun(x, fontfamily, fontface, svgdev))
}

# Table of Unicode math ops at
# http://www.w3.org/TR/MathML2/022.html
toMML.call <- function(x, fontfamily, fontface, svgdev, ...) {
    if (is.name(x[[1]])) {
        funCallToMML(x, fontfamily, fontface, svgdev)
    } else {
        mrow <- newXMLNode("mrow", parent = svgDevParent(svgdev))
        svgDevChangeParent(mrow, svgdev)
        toMML(x[[1]], fontfamily, fontface, svgdev)
        mfenced <- newXMLNode("mfenced", parent = mrow)
        svgDevChangeParent(mfenced, svgdev) 
        toMML(x[[2]], fontfamily, fontface, svgdev)
        svgDevChangeParent(xmlParent(mrow), svgdev) 
    }
}

# fontfamily is used to set explicit 'mathvariant' when it is not
# implicit in the formula element
expr2mml <- function(e, fontfamily, fontface, svgdev) {
    math <- newXMLNode("math", parent = svgDevParent(svgdev),
                       attrs = list(display = "inline"),
                       namespaceDefinitions = "http://www.w3.org/1998/Math/MathML")
    lapply(e, function(x) {
        mrow <- newXMLNode("mrow", parent = math)
        svgDevChangeParent(mrow, svgdev)
        toMML(x, fontfamily = fontfamily, fontface = fontface,
              svgdev = svgdev)
    })

    svgDevChangeParent(xmlParent(math), svgdev)
}

