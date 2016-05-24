theme <-
function(colors=c("blue", "gray", "rose", "green", "gold", "red",
         "dodgerblue", "purple", "sienna", "orange.black",
         "gray.black", "white"),

         color.fill.bar=NULL, trans.fill.bar=NULL,
         color.fill.pt=NULL, trans.fill.pt=NULL,
         color.stroke.bar=NULL, color.stroke.pt=NULL, 
         color.bg=NULL, color.grid=NULL, color.box=NULL,
         color.heat=NULL, ghost=NULL,

         n.cat=getOption("n.cat"), lab.size=getOption("lab.size"),
         suggest=getOption("suggest"),
         quiet=getOption("quiet"), brief=getOption("brief"),

         results=getOption("results"), explain=getOption("explain"),
         interpret=getOption("interpret"), document=getOption("document"), 
         code=getOption("code"),

         width=120, show=FALSE) {

  if (!is.null(color.fill.bar)) if (color.fill.bar == "off")
    color.fill.bar <- "transparent"
  if (!is.null(color.fill.pt)) if (color.fill.pt == "off")
    color.fill.pt <- "transparent"
  if (!is.null(color.stroke.bar)) if (color.stroke.bar == "off")
    color.stroke.bar <- "transparent"
  if (!is.null(color.bg)) if (color.bg == "off") color.bg <- "transparent"
  if (!is.null(color.grid)) if (color.grid == "off" ) color.grid <- "transparent"
  if (!is.null(color.box)) if (color.box == "off") color.box <- "transparent"

  # default transparency levels
  if (!missing(colors)) {
    colors <- match.arg(colors) 
    options(colors=colors)
    if (colors == "dodgerblue")
      options(trans.fill.bar=0.25)
    else
      options(trans.fill.bar=0.00)
    options(trans.fill.pt=0.50) 
  }

  if (!is.null(ghost)) if (ghost) {
    options(trans.fill.bar = 0.50)
    if (getOption("colors") == "blue")
      options(color.fill.bar = .maketrans("lightsteelblue3", .to256("trans.fill.bar")))
    if (getOption("colors") == "gray")
      options(color.fill.bar = .maketrans("gray30", .to256("trans.fill.bar")))
    if (getOption("colors") == "green")
      options(color.fill.bar = rgb(106,127,16, alpha=.to256("trans.fill.bar"), maxColorValue=256))
    if (getOption("colors") == "rose")
      options(color.fill.bar = rgb(245,213,210, alpha=.to256("trans.fill.bar"), maxColorValue=256))
    if (getOption("colors") == "gold")
      options(color.fill.bar = .maketrans("goldenrod2", .to256("trans.fill.bar")))
    if (getOption("colors") == "red")
      options(color.fill.bar = .maketrans("firebrick2", .to256("trans.fill.bar")))
    if (getOption("colors") == "dodgerblue")
      options(color.fill.bar = .maketrans("dodgerblue3", .to256("trans.fill.bar")))
    if (getOption("colors") == "purple")
      options(color.fill.bar = .maketrans("purple1", .to256("trans.fill.bar")))
    if (getOption("colors") == "sienna")
      options(color.fill.bar = .maketrans("sienna3", .to256("trans.fill.bar")))
    if (getOption("colors") == "orange.black")
      options(color.fill.bar = rgb(249,99,2, alpha=.to256("trans.fill.bar"), maxColorValue=256))
    if (getOption("colors") == "gray.black")
      options(color.fill.bar =  .maketrans("gray80", .to256("trans.fill.bar")))
  }
  else {
    colors <- getOption("colors")
  }

  if (!is.null(trans.fill.bar)) {
    options(trans.fill.bar=trans.fill.bar)
    options(color.fill.bar = .maketrans(getOption("color.fill.bar"), .to256("trans.fill.bar")))
  }
  if (!is.null(trans.fill.pt)) {
    options(trans.fill.pt=trans.fill.pt)
    options(color.fill.pt = .maketrans(getOption("color.fill.pt"), .to256("trans.fill.pt")))
  }

  if (!is.null(color.fill.bar))
    if (color.fill.bar == "transparent")
      options(color.fill.bar = color.fill.bar) 
    else
      options(color.fill.bar = .maketrans(color.fill.bar, .to256("trans.fill.bar")))
  if (!is.null(color.fill.pt))
    if (color.fill.pt == "transparent")
      options(color.fill.pt = color.fill.pt) 
    else
      options(color.fill.pt = .maketrans(color.fill.pt, .to256("trans.fill.pt")))

  if (!is.null(color.stroke.bar)) options(color.stroke.bar = color.stroke.bar) 
  if (!is.null(color.stroke.pt)) options(color.stroke.pt = color.stroke.pt) 

  if (!is.null(color.bg)) options(color.bg=color.bg)
  if (!is.null(color.grid)) options(color.grid=color.grid)
  if (!is.null(color.box)) options(color.box=color.box)

  options(quiet=quiet)
  options(brief=brief)
  options(n.cat=n.cat)
  options(suggest=suggest)
  options(lab.size=lab.size)
  options(width=width)

  options(results=results)
  options(explain=explain)
  options(interpret=interpret)
  options(document=document)
  options(code=code)


  if (!missing(colors)) {
  # rgb(20,97,172) is dodgerblue 3.5
    theme <- options("colors")
    if (theme == "dodgerblue") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("dodgerblue3", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("dodgerblue3", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "steelblue4")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "steelblue4")
      if (is.null(color.bg)) options(color.bg = rgb(242,244,245, maxColorValue=256))
      if (is.null(color.grid)) options(color.grid = "snow3")
      if (is.null(color.heat)) options(color.heat = "dodgerblue4")
    }
    if (theme == "blue") {
      if (is.null(color.fill.bar)) 
        options(color.fill.bar = .maketrans("lightsteelblue3", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("lightsteelblue3", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "slategray")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "darkblue")
      if (is.null(color.bg)) options(color.bg = "ghostwhite")
      if (is.null(color.grid)) options(color.grid = "gray90")
      if (is.null(color.heat)) options(color.heat = "darkblue")
    }
    if (theme == "gray") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("gray35", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("gray20", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "gray60")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "black")
      if (is.null(color.bg)) options(color.bg = "gray92")
      if (is.null(color.grid)) options(color.grid = "white")
      if (is.null(color.heat)) options(color.heat = "gray5")
    }
    if (theme == "green") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = rgb(106,127,16, alpha=.to256("trans.fill.bar"),
        maxColorValue=256))
      if (is.null(color.fill.pt))
        options(color.fill.pt = rgb(106,127,16, alpha=.to256("trans.fill.pt"),
        maxColorValue=256))
      if (is.null(color.stroke.bar))
        options(color.stroke.bar = rgb(71,67,52, maxColorValue=256))
      if (is.null(color.stroke.pt))
        options(color.stroke.pt = rgb(71,67,52,  maxColorValue=256))
      if (is.null(color.bg)) options(color.bg = rgb(230,220,143, maxColorValue=256))
      if (is.null(color.grid)) options(color.grid = rgb(96,53,29, alpha=50, maxColorValue=256))
      if (is.null(color.heat)) options(color.heat = "darkgreen")
    }
    if (theme == "rose") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = rgb(245,213,210, alpha=.to256("trans.fill.bar"),
        maxColorValue=256))
      if (is.null(color.fill.pt))
        options(color.fill.pt = rgb(245,213,210, alpha=.to256("trans.fill.pt"),
        maxColorValue=256))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "mistyrose4")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "mistyrose4")
      if (is.null(color.bg)) options(color.bg = "snow1")
      if (is.null(color.grid)) options(color.grid = "snow2")
      if (is.null(color.heat)) options(color.heat = "rosybrown4")
    }
    if (theme == "gold") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("goldenrod2", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("goldenrod2", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "goldenrod4")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "goldenrod4")
      if (is.null(color.bg)) options(color.bg = rgb(255,250,245, maxColorValue=256))
      if (is.null(color.grid)) options(color.grid = rgb(220,222,200, maxColorValue=256))
      if (is.null(color.heat)) options(color.heat = "goldenrod4")
    }
    if (theme == "red") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("firebrick2", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("firebrick2", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "firebrick4")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "firebrick4")
      if (is.null(color.bg)) options(color.bg=rgb(255,251,251, maxColorValue=256))
      if (is.null(color.grid)) if (is.null(color.grid)) options(color.grid="lavenderblush2")
      if (is.null(color.heat)) options(color.heat = "darkred")
    }
    if (theme == "purple") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("purple1", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("purple1", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "purple4")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "purple4")
      if (is.null(color.bg)) options(color.bg = "lavenderblush")
      if (is.null(color.grid)) options(color.grid = "lavenderblush3")
      if (is.null(color.heat)) options(color.heat = "purple4")
    }
    if (theme == "sienna") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("sienna3", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("sienna3", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "sienna4")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "sienna4")
      if (is.null(color.bg)) options(color.bg = "seashell1")
      if (is.null(color.grid)) options(color.grid = "seashell2")
      if (is.null(color.heat)) options(color.heat = "sienna3")
    }
    if (theme == "orange.black") {
      if (is.null(color.fill.bar)) options(color.fill.bar = rgb(249,99,2, alpha=.to256("trans.fill.bar"),
        maxColorValue=256))
      if (is.null(color.fill.pt)) options(color.fill.pt = rgb(249,99,2, alpha=.to256("trans.fill.pt"),
        maxColorValue=256))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = rgb(209,87,3, maxColorValue=256))
      if (is.null(color.stroke.pt)) options(color.stroke.pt = rgb(209,87,3, maxColorValue=256))
      if (is.null(color.bg)) options(color.bg = rgb(.015,.015,.015))
      if (is.null(color.grid)) options(color.grid = rgb(100,100,100, maxColorValue=256))
      if (is.null(color.heat)) options(color.heat = "darkorange3")
    }
    if (theme == "gray.black") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = .maketrans("gray55", .to256("trans.fill.bar")))
      if (is.null(color.fill.pt))
        options(color.fill.pt = .maketrans("gray75", .to256("trans.fill.pt")))
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "gray20")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "gray95")
      if (is.null(color.bg)) options(color.bg = rgb(.015,.015,.015))
      if (is.null(color.grid)) options(color.grid = "gray30")
      if (is.null(color.heat)) options(color.heat = "gray30")
    }
    if (theme == "white") {
      if (is.null(color.fill.bar))
        options(color.fill.bar = "white")
      if (is.null(color.fill.pt))
        options(color.fill.pt = "white")
      if (is.null(color.stroke.bar)) options(color.stroke.bar = "black")
      if (is.null(color.stroke.pt)) options(color.stroke.pt = "black")
      if (is.null(color.bg)) options(color.bg = "transparent")
      if (is.null(color.grid)) options(color.grid = "gray90")
      if (is.null(color.heat)) options(color.heat = "gray70")
    }
  }

  if (!missing(ghost)) {
    options(ghost=ghost)
    if (ghost) {
      options(color.bg = "black")
      options(color.grid = "transparent")
    }
  }

  if (show) {
    cat("\n")
    cat("Note: Colors are in rgb format\n",
        "      A 4th number specifies alpha transparency\n\n", sep="")
    cat("Note: Use \"off\" to indicate transparent, eg: color.grid=\"off\"\n\n")

    cat("colors         Color:", getOption("colors"), "\n")
    cat("color.fill.bar   Bar fill color:", col2rgb(getOption("color.fill.bar"), TRUE), "\n")
    cat("color.fill.pt    Point fill color:", col2rgb(getOption("color.fill.pt"), TRUE), "\n")
    cat("trans.fill.bar Bar transparency:", getOption("trans.fill.bar"), "\n")
    cat("trans.fill.pt  Point transparency:", getOption("trans.fill.pt"), "\n")
    cat("color.stroke.bar Bar stroke color:", col2rgb(getOption("color.stroke.bar"), TRUE), "\n")
    cat("color.stroke.pt  Point stroke color:", col2rgb(getOption("color.stroke.pt"), TRUE), "\n")
    cat("color.bg         Background color:", col2rgb(getOption("color.bg")), "\n")
    cat("color.grid       Grid color:", col2rgb(getOption("color.grid")), "\n")
    cat("color.box        Graph Border color:", col2rgb(getOption("color.box")), "\n")
    cat("color.heat       Heat map color:", col2rgb(getOption("color.heat")), "\n")
    if (is.null(ghost)) ghost <- FALSE
    cat("ghost          Ghost colors (black bck, no grid lines):", getOption("ghost"), "\n")
    cat("\n")
    cat("quiet     Suppress console output for many functions:", getOption("quiet"), "\n")
    cat("brief     Reduce console output for many functions:", getOption("brief"), "\n")
    cat("suggest   Display suggests for enhanced input:", getOption("suggest"), "\n")
    cat("n.cat     Number of categories for categorical variable:", getOption("n.cat"), "\n")
    cat("lab.size  x and y axis label size:", getOption("lab.size"), "\n")
    cat("width     Column width:", getOption("width"), "\n")
    cat("\n")
  }


}
