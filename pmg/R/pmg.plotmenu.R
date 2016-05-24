### this file synchronizes with menes.R
### this one provides the widget.list stuff

par.setup.list = list(
  title = "par()",
  help = "par",
  type = NULL,                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "par(",
    ending = ")"
    ),
  arguments = list(
    Setup= list(
      bty = list(
        type="gdroplist",
        items = c('','"o"', '"l"', '"7"','"c"', '"u"',  '"]"')
        ),
      pty = list(
        type="gdroplist",
        items=c('','"s"','"m"')
        ),
      xpd=list(
        type="gdroplist",
        items=c("","TRUE","FALSE","NA")
        )
      )
    )
  )

par.axes.list = list(
  title = "par()",
  help = "par",
  type = NULL,                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "par(",
    ending = ")"
    ),
  arguments = list(
    axes=list(
      xaxt = list(
        type="gdroplist",
        items=c("",'"s"','"n"')
        ),
      yaxt = list(
        type="gdroplist",
        items=c("",'"s"','"n"')
        ),
      xlog = emptyTRUE.list,
      ylog = emptyTRUE.list,
      las = list(
        type="gdroplist",
        items=c("",0,1,2,3)
        )
      )
    )
  )
par.colors.list = list(
  title = "par()",
  help = "par",
  type = NULL,                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "par(",
    ending = ")"
    ),
  arguments = list(
    colors=list(
      bg=default.color.list,
      fg=default.color.list,
      col.main=default.color.list,
      col.sub=default.color.list,
      col.axis=default.color.list,
      col.lab=default.color.list
      )
    )
  )

par.fonts.list = list(
  title = "par()",
  help = "par",
  type = NULL,                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "par(",
    ending = ")"
    ),
  arguments = list(
    fonts = list(
      family = list(
        type="gdroplist",
        items=c("",'"serif"', '"sans"', '"mono"','"symbol"')
        ),
      font = list(
        type="gdroplist",
        items=c("",1,2,3,4)
        )
      ),
    margins=list(
      mar = list(
        type="gedit",
        text=""
        )
      )
    )
  )

par.onfigures.list = list(
  title = "par()",
  help = "par",
  type = NULL,                      # either text or graphic
  variableType = NULL,
  assignto = NULL,
  action = list(
    beginning = "par(",
    ending = ")"
    ),
  arguments = list(
    "Number of figures"=list(
      mfrow=EMPTY.list
      )
    )
  )


##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## Univariate

barplot.list = list(
  title = "univariate barplot()",
  help = "barplot",
  type = "graphic",                      # either text or graphic
  variableType = "univariatetable",
  variableTypeExtra = list(name="xlabel",value="height"),
  action = list(
    beginning = "barplot(",
    ending = ")"
    ),
  arguments = list(
    adjustments = list(
      width=list(
        type = "gedit",
        text = "1",
        coerce.with=as.numeric
        ),
      horiz= FALSE.list
      ),
    labels = labels.list
    )
  )

piechart.list = list(
  title = "pie()",
  help = "pie",
  type = "graphic",                      # either text or graphic
  variableType = "univariatetable",
  action = list(
    beginning = "pie(",
    ending = ")"
    ),
  arguments = list(
    adjustments = list(
      labels=list(
        type = "gedit",
        text = "names(x)"
        ),
      clockwise= FALSE.list
      ),
    labels = labels.list
    )
  )


univariate.boxplot.list = list(
  title = "boxplot()",
  help = "boxplot",
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  action = list(
    beginning = "boxplot(",
    ending = ")"
    ),
  arguments = list(
    adjustments = list(
      horizontal= FALSE.list,
      add= FALSE.list,
      notch = FALSE.list,
      varwidth = FALSE.list,
      col = default.color.list
      ),
    labels = labels.list
    )
 )


## other boxplots are similar
model.boxplot.list = bivariate.boxplot.list = univariate.boxplot.list
bivariate.boxplot.list$variableType="bivariate"
model.boxplot.list$variableType="model"


bivariate.boxplot.list = list(
  title = "boxplot()",
  help = "boxplot",
  type = "graphic",                      # either text or graphic
  variableType = NULL,
  action = list(
    beginning = "boxplot(",
    ending = ")"
    ),
  arguments = list(
    data = list(
      x = list(
        type="gedit"
        ),
      "..."=list(
        type="gedit"
        )
      ),
    adjustments = list(
      horizontal= FALSE.list,
      add= FALSE.list,
      notch = FALSE.list,
      varwidth = FALSE.list,
      col = default.color.list
      ),
    labels = labels.list
    )
 )






##
hist.list = list(
  title = "hist()",
  help = "hist",
  action = list(
    beginning = "hist(",
    ending = ")"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  arguments = list(
    adjustments = list(
      breaks= list(
        type="gdroplist",
        items=c("\"Sturges\"","\"Scott\"","\"Friedman-Diaconis\"")
        ),
      probability = TRUE.list,
      include.lowest = TRUE.list,
      right = TRUE.list,
      density = NULL.list,
      angle = list(
        type="gedit",
        text="45"
        ),
      border = EMPTY.list,
      col = default.color.list,
      main = EMPTY.list      
      )
    )
  )


densityplot.list = list(
  title = "density()",
  help = "density",
  action = list(
    beginning = "plot(density(",
    ending = "))"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  arguments = list(
    adjustments = list(
      bw = list(
        type = "gdroplist",
        items = c("\"nrd0\"","\"nrd\"","\"ucv\"","\"bcv\"","\"SJ\"")
        ),
      adjust = list(
        type = "gedit",
        text = 1
        ),
      kernel = list(
        type = "gdroplist",
        items = c("\"gaussian\"", "\"epanechnikov\"", "\"rectangular\"", "\"triangular\"", "\"biweight\"", "\"cosine\"", "\"optcosine\"")
        )
      )
    )
  )

qqnorm.list = list(
  title = "qqnorm()",
  help = "qqnorm",
  action = list(
    beginning = "qqnorm(",
    ending = ")"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  arguments = list(
    adjustments = list(
      xlab = EMPTY.list,
      ylab = EMPTY.list
      )
    )
  )


##################################################
## add to graphic
"add.points.list" = list(
  title = "points()",
  help = "points",
  action = list(
    beginning = "points(",
    ending = ")"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "bivariate",
  arguments = list(
    adjustments = list(
      col = default.color.list,
      pch = list(
        type = "gedit",
        text = "1"
        )
      )
    )
  )

"add.lines.list" = list(
  title = "lines()",
  help = "lines",
  action = list(
    beginning = "lines(",
    ending = ")"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "bivariate",
  arguments = list(
    adjustments = list(
      col = default.color.list,
      lty = lty.list
      )
    )

  )
"add.density.list" = list(
  title = "add density()",
  help = "density",
  action = list(
    beginning = "lines(density(",
    ending = "))"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  arguments = list(
    adjustments = list(
      bw = list(
        type = "gdroplist",
        items = c("\"nrd0\"","\"nrd\"","\"ucv\"","\"bcv\"","\"SJ\"")
      ),
      adjust = list(
        type = "gedit",
        text = 1
        ),
      kernel = list(
        type = "gdroplist",
        items = c("\"gaussian\"", "\"epanechnikov\"", "\"rectangular\"", "\"triangular\"", "\"biweight\"", "\"cosine\"", "\"optcosine\"")
        )
      )
    )

  )


## curve adds a function
"add.curve.list" = list(
  title = "curve()",
  help = "curve",
  action = list(
    beginning = "curve(",
    ending = "))"
    ),
  type = "graphic",                      # either text or graphic
  variableType = NULL,
  arguments = list(
    arguments = list(
      expr = EMPTY.list,
      label = list(
        type= "glabel",
        text = "An expression in x or name"
        ),
      add = TRUE.list,
      label = list(
        type= "glabel",
        text = "Either add or specify limits"
        ),
      from = EMPTY.list,
      to   = EMPTY.list
      )
    )
  )

rug.list = list(
  title = "rug()",
  help = "rug",
  action = list(
    beginning = "rug(",
    ending = ")"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  arguments = list(
    adjustments = list(
      ticksize = list(
        type = "gedit",
        text = "0.03"
        ),
      side = list(
        type = "gedit",
        text = "1"
        ),
      lwd = list(
        type = "gedit",
        text = "0.5"
        ),
      col = default.color.list
      )
    )
  )


title.list =  list(
  title = "title()",
  help = "title",
  action = list(
    beginning = "title(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    labels = labels.list
    )
  )



##################################################
## bivariate

scatterplot.list =  list(
  title = "plot()",
  help = "plot",
  action = list(
    beginning = "plot(",
    ending = ")"
    ),
  variableType = "bivariate",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    type = list(
      type = "gdroplist",
      items = c('"p"','"l"','"b"','"c"','"o"','"h"','"s"','"S"')
      ),
    pch = pch.list,
    cex = EMPTY.list,
    col = default.color.list,
    labels = labels.list
    )
  )
  

## sunflower plot

sunflowerplot.list =  list(
  title = "sunflowerplot()",
  help = "sunflowerplot",
  action = list(
    beginning = "sunflowerplot(",
    ending = ")"
    ),
  variableType = "bivariate",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    "plot type" = list(                   # types in genericWidget
      pch = list(
        type = "gedit",
        text = "16"
        ),
      col = default.color.list
      ),
    "Sizes" = list(
      cex = list(
        type = "gedit",
        text = "0.8"
        ),
      size = list(
        type = "gedit",
        text = "1/8"
        ),
      seg.col = list(
        type = "gedit",
        text = "2"
        ),
      seg.lwd = list(
        type = "gedit",
        text = "1.5"
        )
      ),
    labels = labels.list
    )
  )


## qqplot
qqplot.list =  list(
  title = "qqplot()",
  help = "qqplot",
  action = list(
    beginning = "qqplot(",
    ending = ")"
    ),
  variableType = "bivariate",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    "plot type" = list(                   # types in genericWidget
      pch = list(
        type = "gedit",
        text = "16"
        ),
      col = default.color.list
      ),
    "Sizes" = list(
      cex = list(
        type = "gedit",
        text = "0.8"
        )
      ),
    labels = labels.list
    )
  )


### model plots
pairs.list =  list(
  title = "pairs()",
  help = "pairs",
  action = list(
    beginning = "pairs(",
    ending = ")"
    ),
  type = "graphic",                      # either text or graphic
  variableType = "univariate",
  arguments = list(
    labels = labels.list
    )
  )




###
stripchart.list = list(
  title = "stripchart()",
  help = "stripchart",
  action = list(
    beginning = "stripchart(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    variables = list(
      x = list(
        type="geditnamedlist",
        text = ""
        )
      ),
    arguments = list(                   # types in genericWidget
      method = list(
        type = "gdroplist",
        items = c('"stack"','"jitter"','"overplot"')
        ),
      vertical = FALSE.list,
      jitter = list(
        type="gedit",
        text = 0.1
        ),
      offset = list(
        type="gedit",
        text = "1/3"
        ),
      add = FALSE.list,
      at = EMPTY.list
      ),
    "settings" = list(
      xlim = EMPTY.list,
      ylim = EMPTY.list,
      pch = pch.list,
      cex = EMPTY.list,
      col = default.color.list
      )
    )
  )

## dotchart
dotchart.list = list(
  title = "dotchart()",
  help = "dotchart",
  action = list(
    beginning = "dotchart(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      labels = NULL.list,
      groups = NULL.list,
      gdata = NULL.list
      ),
    "plot arguments" = list(
      pch = list(
        type = "gedit",
        text  = 21
        ),
      gpch = list(
        type = "gedit",
        text  = 21
        ),
      cex = EMPTY.list,
      bg = EMPTY.list,
      color = EMPTY.list,
      gcolor = EMPTY.list,
      lcolor = EMPTY.list,
      xlim = EMPTY.list
      ),
    "labels" = labels.list
    )
  )

## cumulative distribution plot. Define function to slip in plot arguments
plotecdf = function(x,ylab="Fn(x)",verticals=FALSE,col.01line="gray70") {
  plot(ecdf(x), ylab=ylab, verticals = verticals, col.01line = col.01line)
}

ecdf.list = list(
  title = "ecdf()",
  help = "ecdf",
  action = list(
    beginning = "plotecdf(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    "plot arguments" = list(
      ylab = list(
        type = "gedit",
        text  = "'Fn(x)'"
        ),
      verticals=FALSE.list,
      col.01 = "'gray70'"
      )
    )
  )

scatterplot.model.list =  list(
  title = "plot()",
  help = "plot",
  action = list(
    beginning = "plot(",
    ending = ")"
    ),
  variableType = "model",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      type = list(
        type = "gdroplist",
        items = c('"p"','"l"','"b"','"c"','"o"','"h"','"s"','"S"')
        ),
      pch = pch.list,
      cex = EMPTY.list,
      col = default.color.list
      ),
    labels = labels.list
    )
  )
  



##################################################
## lattice stuff
xyplot.list = list(
  title = "xyplot()",
  help = "xyplot",
  action = list(
    beginning = "require(lattice);xyplot(",
    ending = ")"
    ),
  variableType = "lattice",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      panel = EMPTY.list
      ),
    jitter = list(      jitter = FALSE.list,
      factor = EMPTY.list
      ),
    labels = labels.list
    )
  )

## dotplot
dotplot.list = list(
  title = "dotplot()",
  help = "dotplot",
  action = list(
    beginning = "require(lattice);dotplot(",
    ending = ")"
    ),
  variableType = "lattice",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      panel = EMPTY.list
      ),
    jitter = list(
      jitter = FALSE.list,
      factor = EMPTY.list
      ),
    labels = labels.list
    )
  )


## barchart
barchart.list = list(
  title = "barchart()",
  help = "barchart",
  action = list(
    beginning = "require(lattice);barchart(",
    ending = ")"
    ),
  variableType = "lattice",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      panel = EMPTY.list
      ),
    jitter = list(
      jitter = FALSE.list,
      factor = EMPTY.list
      ),
    labels = labels.list
    )
  )

## stripplot
stripplot.list = list(
  title = "stripplot()",
  help = "stripplot",
  action = list(
    beginning = "require(lattice);stripplot(",
    ending = ")"
    ),
  variableType = "lattice",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      panel = EMPTY.list
      ),
    jitter = list(
      jitter = FALSE.list,
      factor = EMPTY.list
      ),
    labels = labels.list
    )
  )

## bwplot
bwplot.list = list(
  title = "bwplot()",
  help = "bwplot",
  action = list(
    beginning = "require(lattice);bwplot(",
    ending = ")"
    ),
  variableType = "lattice",           # uni/bi/model/lattice
  type = "graphic",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      panel = EMPTY.list
      ),
    jitter = list(
      jitter = FALSE.list,
      factor = EMPTY.list
      ),
    labels = labels.list
    )
  )

