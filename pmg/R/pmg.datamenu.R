
## gtkDataViewer
dataViewer.list = list(
  title = "browseDataAsCSV()",
  help = "browseDataAsCSV",
  action = list(
    beginning = "browseDataAsCSV(",
    ending = ")"
    ),
  variableType = NULL,
  type = "text",                      # either text or graphic
  arguments = list(
    arguments = list(
      dataset = list(
        type = "gedit",
        text = ""
        )
      )
    )
  )


write.csv.list = list(
  title = "write.csv()",
  help = "write.csv",
  type = "text",                      # either text or graphic
  assignto = NULL,
  variableType = NULL,           # uni/bi/model/lattice
  action = list(
    beginning = "write.csv(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      x = list(
        type = "gedit"
        ),
      file = list(
        type = "gfilebrowse",
        quote=TRUE
        ),
      row.names = list(
        type = "gdroplist",
        items = c(FALSE,TRUE)
        ),
      col.names = list(
        type = "gdroplist",
        items = c(TRUE, FALSE)
        ),
      append = list(
        type = "gdroplist",
        items = c(FALSE, TRUE)
        ),
      quote =  list(
        type = "gdroplist",
        items = c(TRUE, FALSE)
        )
      )
    )
  )


## table and xtabs
table.list = list(
  title = "table()",
  help = "table",
  action = list(
    beginning = "table(",
    ending = ")"
    ),
  variableType = "bivariate",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  arguments = list(
    arguments = list(
      exclude = list(
        type="gdroplist",
        items = "c(NA,NaN)"
        ),
      deparse.level = list(
        type = "gedit",
        text = 1
        )
      )
    )
  )

xtabs.list = list(
  title = "xtabs()",
  help = "xtabs",
  action = list(
    beginning = "xtabs(",
    ending = ")"
    ),
  variableType = "model",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  arguments = list(
    arguments = list(
      drop.unused.levels = FALSE.list
      )
    )
  )

## EDA functions
stem.list = list(
  title = "stem()",
  help = "stem",
  type="text",
  assignto = NULL,
  variableType = "univariate",
  action = list(
    beginning = "stem(",
    ending=")"
    ),
  arguments = list(
    arguments=list(
      scale = list(
        type="gedit",
        text=1
        ),
      width = list(
        type="gspinbutton",
        from = 30,
        to=100,
        value = 80,
        by=1,
        digits = 0,
        horizontal=TRUE
        )
      )
    )
  )

## summary
summary.list = list(
  title = "summary()",
  help = "summary",
  type="text",
  assignto = NULL,
  variableType = NULL,
  action = list(
    beginning = "summary(",
    ending=")"
    ),
  arguments = list(
    data = list(
      object = list(
        type="gedit",
        text = ""
        )
      )
    )
  )

## mean
mean.list = list(
  title = "mean()",
  help = "mean",
  type = "text",                      # either text or graphic
  assignto = NULL,
  variableType = "univariate",           # uni/bi/model/lattice
  action = list(
    beginning = "mean(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      trim = list(
        type = "gspinbutton",
        from = 0,
        to = 0.5,
        by = 0.05,
        digits = 2
        ),
      na.rm = TRUE.list
      )
    )
  )

## median
median.list = list(
  title = "median()",
  help = "median",
  type = "text",                      # either text or graphic
  assignto = NULL,
  variableType = "univariate",           # uni/bi/model/lattice
  action = list(
    beginning = "median(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      na.rm = TRUE.list
      )
    )
  )

## sd
sd.list = list(
  title = "sd()",
  help = "sd",
  type = "text",                      # either text or graphic
  assignto = NULL,
  variableType = "univariate",           # uni/bi/model/lattice  
  action = list(
    beginning = "sd(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      na.rm = TRUE.list
      )
    )
  )

## IQR
IQR.list = list(
  title = "IQR()",
  help = "IQR",
  type = "text",                      # either text or graphic
  assignto = NULL,
  variableType = "univariate",           # uni/bi/model/lattice
  action = list(
    beginning = "IQR(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
    na.rm = TRUE.list
      )
    )
  )

## mad
mad.list = list(
  title = "mad()",
  help = "mad",
  type = "text",                      # either text or graphic
  assignto = NULL,
  variableType = "univariate",           # uni/bi/model/lattice
  action = list(
    beginning = "mad(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      constant = list(
        type = "gedit",
        text = 1.4826
        ),
      na.rm = TRUE.list,
      low = FALSE.list,
      high = FALSE.list
      )
    )
  )

## skewness
skewnessList = list(
  title = "skewness()",
  help = "",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  variableType = "univariate",           # uni/bi/model/lattice
  action = list(
    beginning = "skewness(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      na.rm = TRUE.list
      )
    )
  )

## kurtosis
kurtosisList = list(
  title = "kurtosis()",
  help = "",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  variableType = "univariate",           # uni/bi/model/lattice
  action = list(
    beginning = "kurtosis(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      na.rm = TRUE.list
      )
    )
  )

## dpqr functions
## normal
pdf.normal.list = list(
  title = "dnorm()",
  help = "dnorm",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "dnorm(",
    ending = ")"
    ),
  arguments = list(
    arugments = list(
      x = list(
        type = "gedit",
        text = ""
        )
      ),
    parameters=list(
      mean = list(
        type = "gedit",
        text = "0"
        ),
      sd = list(
        type = "gedit",
        text = "1"
        )
      ),
    others = list(
      log= FALSE.list
      )
    )
  )
  

quantile.normal.list = list(
  title = "qnorm()",
  help = "qnorm",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "qnorm(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      p = list(
        type = "gedit",
        text = ""
        )
      ),
    parameters = list(
      mean = list(
        type = "gedit",
        text = "0"
      ),
      sd = list(
        type = "gedit",
        text = "1"
        )
      ),
    others = list(
      lower.tail= TRUE.list,
      log.p= FALSE.list
      )
    )
  )
  

random.normal.list = list(
  title = "rnorm()",
  help = "rnorm",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  action = list(
    beginning = "rnorm(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
    n = list(
      type = "gedit",
      text = ""
      )
      ),
    parameters=list(
    mean = list(
      type = "gedit",
      text = "0"
      ),
    sd = list(
        type = "gedit",
        text = "1"
        )
    )
  )
  )

## t-dist
pdf.t.list = list(
  title = "dt()",
  help = "dt",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "dt(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      x = EMPTY.list
      ),
    parameters = list(
      df = list(
        type = "gedit",
        text = "0"
        ),
      ncp = list(
        type = "gedit",
        text = "0"
        )
      )
    ),
  others = list(
    log= FALSE.list
    )
  )


quantile.t.list = list(
  title = "qt()",
  help = "qt",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "qt(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      p = EMPTY.list
      ),
    parameters = list(
      df = EMPTY.list
      ),
    others = list(
      lower.tail= TRUE.list,
      log.p= FALSE.list
      )
    )
  )

random.t.list = list(
  title = "rt()",
  help = "rt",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  action = list(
    beginning = "rt(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      n = EMPTY.list
      ),
    parameters=list(
      df = list(
        type = "gedit",
        text = "0"
        )
      )
    )
  )


## exponential
pdf.exp.list = list(
  title = "dexp()",
  help = "dexp",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "dexp(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      x = EMPTY.list
      ),
    parameters=list(
      rate = list(
        type = "gedit",
        text = "1"
        )
      ),
    others = list(
      log= FALSE.list
      )
    )
  )


quantile.exp.list = list(
  title = "qexp()",
  help = "qexp",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "qexp(",
    ending = ")"
    ),
  arguments = list(
    arguments = list(
      p = EMPTY.list
      ),
    parameters=list(
      rate = list(
        type = "gedit",
        text = "1"
        )
      ),
    others = list(
      lower.tail= TRUE.list,
      log.p= FALSE.list
      )
    )
  )
random.exp.list = list(
  title = "rexp()",
  help = "rexp",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  action = list(
    beginning = "rexp(",
    ending = ")"
    ),
  arguments = list(
    arguments=list(
      n = EMPTY.list
      ),
    parameters=list(
      rate = list(
        type = "gedit",
        text = "1"
        )
      )
    )
  )


## uniform
pdf.uniform.list = list(
  title = "dunif()",
  help = "dunif",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "dunif(",
    ending = ")"
    ),
  arguments = list(
    arguments=list(
    x = EMPTY.list
      ),
    parameters=list(
      min = list(
        type = "gedit",
        text = "0"
        ),
    max = list(
      type = "gedit",
      text = "1"
      )
      ),
    others = list(
      log= FALSE.list
      )
    )
  )

quantile.uniform.list = list(
  title = "qunif()",
  help = "qunif",
  type = "text",                      # either text or graphic
  assignto = NULL,
  action = list(
    beginning = "qunif(",
    ending = ")"
    ),
  arguments = list(
    arguments=list(
    p = EMPTY.list
      ),
    parameters=list(
      min = list(
        type = "gedit",
        text = "0"
      ),
      max = list(
        type = "gedit",
        text = "1"
      )
      ),
    others = list(
      lower.tail= TRUE.list,
      log.p= FALSE.list
      )
    )
  )

random.uniform.list = list(
  title = "runif()",
  help = "runif",
  type = "text",                      # either text or graphic
  assignto = TRUE,
  action = list(
    beginning = "runif(",
    ending = ")"
    ),
  arguments = list(
    arguments=list(
      n = EMPTY.list
      ),
    parameters=list(
      min = list(
        type = "gedit",
        text = "0"
        ),
      max = list(
        type = "gedit",
        text = "1"
        )
      )
    )
  )


#######
sample.list = list(
  title = "sample()",
  help = "sample",
  action = list(
    beginning = "sample(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      size = EMPTY.list,
      replace = FALSE.list,
      prob = EMPTY.list
      )
    )
  )



##################################################
## stack, unstack, subset
 ##?? Do  theses fit in?


stack.list = list(
  title = "stack()",
  help = "stack",
  action = list(
    beginning = "stack(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice/NULL
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    Variables = list(                   # types in genericWidget
      x = list(
        type = "geditnamedlist",
        val = ""
        ),
      y= list(
        type = "glabel",
        text = "Drag value(s) to stack into x area"
        )
      )
    )
  )

unstack.list = list(
  title = "unstack()",
  help = "unstack",
  action = list(
    beginning = "unstack(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice/NULL
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    Variables = list(                   # types in genericWidget
      x = EMPTY.list
      )
    )
  )


## subset -- use pmg.subset.dialog instead
subset.list =  list(
  title = "subset()",
  help = "subset",
  action = list(
    beginning = "subset(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      subset= EMPTY.list,
      select = EMPTY.list,
      drop = FALSE.list
      )
    )
  )


## bivariate
cor.list = list(
  title = "cor()",
  help = "cor",
  action = list(
    beginning = "cor(",
    ending = ")"
    ),
  variableType = "bivariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      method = list(
        type="gdroplist",
        items = c('"pearson"', '"kendall"', '"spearman"')
        ),
      use =  list(
        type="gdroplist",
        items = c('"all.obs"', '"complete.obs"', '"pairwise.complete.obs"')
        )
      )
    )
 )

blank.list = list(
  title = "blank",
  help = "blank",
  action = list(
    beginning = "writeme(",
    ending = ")"
    ),
  variableType = "bivariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = NULL,                      # TRUE for assignto
  arguments = list(
    arguments = list(                   # types in genericWidget
      writeme = list(
        type="gedit",
        text = "Write me"
        )
      )
    )
 )

## coerce stuff
as.numeric.list =  list(
  title = "as.numeric()",
  help = "as.numeric",
  action = list(
    beginning = "as.numeric(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE                      # TRUE for assignto
  )


as.character.list = list(
  title = "as.character()",
  help = "as.character",
  action = list(
    beginning = "as.character(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE                      # TRUE for assignto
  )

data.frame.list = list(
  title = "data.frame()",
  help = "data.frame",
  action = list(
    beginning = "data.frame(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  suppressDotDotDot = TRUE,             # TRUE to suppress, NULL to handle
  arguments = list(
    variables = list(
      "..." = list(
        type="geditnamedlist",
        val = "",
        wideBody = TRUE                      # big and fat boy
        )
      ),
    "arguments"=list(
      row.names = list(
        type="gedit",
        text=""
        ),
      blank = list(
        type="glabel",
        text=""
        ),
      check.rows = TRUE.list,
      check.names = TRUE.list
      )
    )
  )

as.data.frame.list =list(
  title = "as.data.frame()",
  help = "as.data.frame",
  action = list(
    beginning = "as.data.frame(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    argument=list(
      row.names = list(
        type="gdroplist",
        items = c("","NULL")
        ),
      optional = FALSE.list
      )
    )
  )

matrix.list = list(
  title = "matrix()",
  help = "matrix",
  action = list(
    beginning = "matrix(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    argument=list(
      data = list(
        type="gedit",
        text = "NA"
        ),
      blank = list(
        type="glabel",
        text=""
        ),
      nrow = list(
        type="gedit",
        text=1
        ),
      ncol = list(
        type="gedit",
        text=1
        ),
      byrow = FALSE.list,
      dimnames = list(
        type="gedit",
        text = "NULL"
        )
      )
    )
  )

as.matrix.list =  list(
  title = "as.matrix()",
  help = "as.matrix",
  action = list(
    beginning = "as.matrix(",
    ending = ")"
    ),
  variableType = "univariate",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE                      # TRUE for assignto
  )

groupedData.list = list(
  title = "groupedData()",
  help = "groupedData",
  action = list(
    beginning = "groupedData(",
    ending = ")"
    ),
  variableType = "lattice",           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    argument=list(
      outer = list(
        type="gedit",
        text=""
        ),
      inner = list(
        type="gedit",
        text = ""
        ),
      labels = list(
        type="gedit",
        text = ""
        ),
      units = list(
        type="gedit",
        text = ""
        )
      ),
    ordered = list(
            order.groups = FALSE.list,
      FUN = list(
        type="gedit",
        text=""
        )
      )
    )
  )

factor.list = list(
  title = "factor()",
  help = "factor",
  action = list(
    beginning = "factor(",
    ending = ")"
    ),
  variableType = NULL,           # uni/bi/model/lattice
  type = "text",                        # either text or graphic
  assignto = TRUE,                      # TRUE for assignto
  arguments = list(
    variables=list(
      x = EMPTY.list
      ),
    argument=list(
      levels = list(
        type="gdroplist",
        items = c("","sort(unique.default(x)")
        ),
      labels = list(
        type="gdroplist",
        items=c("","levels")
        ),
      exclude = list(
        type="gdroplist",
        items=c("","NA")
        ),
      ordered = list(
        type="gdroplist",
        items=c("","is.ordered(x)")
        )
      )
    )
  )
