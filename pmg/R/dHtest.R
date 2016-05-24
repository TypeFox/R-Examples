## write htest interface for showing t.test, etc. with drop

## This is kinda ugly!
## TODO
## chisq.test, fisher.test, "mantelhaen.test" "mcnemar.test"
## * ---- isn't working names is giving unique guys
## 



dTestsDialog = function() {

  tests = list(
    "1-sample test of proportion" = "gui.prop.test",
    "1-sample exact test of proportion" = "gui.binom.test",
    "2-sample test of proportion" = "gui.prop.test.2sample",
    "2-sample exact test of proportion" = "gui.binom.test.2sample",
    "-----" = "-----",
    ##
    "1-sample t-test"="gui.t.test",
    "1-sample signed rank test"="gui.wilcox.signed.rank.test",
    "2-sample t-test" =   "gui.t.test.2.sample",
    "2-sample t-test var. equal" = "gui.t.test.2.sample.var.equal",
    "2-sample t-test, formula" =   "gui.t.test.2.sample.formula",
    "2-sample t-test var. equal, formula" = "gui.t.test.2.sample.var.equal.formula",
    "2-sample rank sum test" = "gui.wilcox.rank.sum.test",
    "----" = "-----",
    ##
    "Oneway ANOVA" = "gui.oneway.test",
    "Kruska-Wallis test" = "gui.kruskal.test",
    "---" = "-----",
    ##
    "Correlation test" = "gui.cor.test",
    "--" = "-----",
    ##
    "Test of variances" = "gui.var.test",
    "Ansari test" = "gui.ansari.test",
    "Bartlett test" = "gui.bartlett.test",
    "Fligner test" = "gui.fligner.test",
    "Mood test" = "gui.mood.test",
    ##
    "------" = "-----",
    "2-sample Kolmogorov-Smirnov test" = "gui.ks.test",
    "Shapiro.test" = "gui.shapiro.test"
    
    )
  
  dialogList = list()
#  for(i in names(tests)) {
#    tmp = tests[[i]]
#    if(length(tmp) == 1)
#      dialogList[[i]] <- do.call(tmp,list())
#    else
#      dialogList[[i]] <- do.call(tmp[1],tmp[2]) # tmp a list with args
#  }
  dialogList[["FirstOne"]] = glabel("Select a test from popup")
  
  
  win = pmgWC$new("Dynamic tests",
    width=400,height=300,
    handler = function(h,...) {
    for(i in dialogList) {
      ids = tag(i,"dropHandlers")
      if(!is.null(ids))
        removehandler(i,ids)
    }})
  gp = ggroup(horizontal=FALSE, container=win, raise.on.dragmotion = TRUE)
  popupGroup = ggroup(container=gp)
  addSpring(popupGroup)
  theNames = c("",names(tests))
  testPopup = gdroplist(theNames, container=popupGroup)
  
  testWindow = ggroup(container=gp)
  add(testWindow,dialogList[["FirstOne"]], expand=TRUE)
  tag(testWindow,"dialogList") <- dialogList
  tag(testWindow,"currentTest") <- dialogList[["FirstOne"]]
  
  addhandlerchanged(testPopup, handler = function(h,...) {
    popupValue = svalue(testPopup)
    if(!is.empty(popupValue) || popupValue != "-----") {
      delete(testWindow,tag(testWindow,"currentTest"))
      dialogList = tag(testWindow, "dialogList")
##      if(is.null(dialogList[[popupValue]])) {
      if(! popupValue %in% names(dialogList)) {
        dialogList[[popupValue]] <- do.call(tests[[popupValue]],list())
        tag(testWindow,  "dialogList") <- dialogList
      }
      add(testWindow,dialogList[[popupValue]], expand=TRUE)
      tag(testWindow,"currentTest") <-  dialogList[[popupValue]]
  }
  })

}




## examples of standalone usage
gui.oneway.test = function(container = NULL) {
  gui.htest("oneway.test", type="x~f",
                   template = oneway.test(mtcars$mpg ~ as.factor(mtcars$cyl)),
                   container = container)
}

gui.kruskal.test = function(container = NULL) {
  gui.htest("kruskal.test", type="x~f",
                   template = kruskal.test(mtcars$mpg ~ as.factor(mtcars$cyl)),
                   container = container)
}

gui.bartlett.test = function(container =NULL) {
  gui.htest("bartlett.test", type="x~f",
                   template = bartlett.test(mtcars$mpg ~ as.factor(mtcars$cyl)),
                   container = container)
}

gui.fligner.test = function(container = NULL) {
  gui.htest("fligner.test", type="x~f",
                   template = fligner.test(mtcars$mpg ~ as.factor(mtcars$cyl)),
                   container = container)
}


gui.t.test = function(container = NULL) {
  gui.htest("t.test",type="univariate",
            template=t.test(rnorm(100)),
            container = container)
}
gui.t.test.2.sample = function(container= NULL) {
  gui.htest("t.test",type="bivariate",
            template=t.test(rnorm(100), rnorm(100)),
            container = container)
}
gui.t.test.2.sample.formula = function(container= NULL) {
  gui.htest("t.test",type="x~f",
            template=t.test(rnorm(100) ~ factor(sample(1:2, 100,T))),
            container = container)
}
gui.t.test.2.sample.var.equal = function(container=NULL) {
  gui.htest("t.test",type="bivariate",
            template=t.test(rnorm(100), rnorm(100), var.equal=TRUE),
            extra.args = list("var.equal"=TRUE),
            container = container)
}
gui.t.test.2.sample.var.equal.formula = function(container=NULL) {
  gui.htest("t.test",type="x~f",
            template=t.test(rnorm(100) ~ factor(sample(1:2, 100,T)), var.equal=TRUE),
            extra.args = list("var.equal"=TRUE),
            container = container)
}
gui.wilcox.signed.rank.test = function(container=NULL) {
  gui.htest("wilcox.test",type="univariate",
            template=wilcox.test(rnorm(100)),
            container = container)
}
gui.wilcox.rank.sum.test = function(container=NULL) {
  gui.htest("wilcox.test",type="bivariate",
            template=wilcox.test(rnorm(100),rnorm(100)),
            container = container)
}

gui.ks.test = function(container=NULL) {
  gui.htest("ks.test",type="bivariate",
            template=ks.test(rnorm(100),rnorm(100)),
            container = container)
}

gui.ansari.test = function(container=NULL){
  gui.htest("ansari.test",type="bivariate",
            template=ansari.test(rnorm(100),rnorm(100)),
            container = container)
}

gui.cor.test = function(container=NULL) {
  gui.htest("cor.test",type="bivariate",
                   template=cor.test(rnorm(100),rnorm(100), method="pearson"),
                   extra.args=list("method"="pearson"),
                   container = container)
}

gui.mood.test = function(container=NULL) {
  gui.htest("mood.test",type="bivariate",
                   template=mood.test(rnorm(100),rnorm(100)),
                   container = container)
}

Shapiro.test = function(x,...) shapiro.test(x)
gui.shapiro.test = function(container=NULL) {
  gui.htest("Shapiro.test",type="univariate",
                   template=shapiro.test(rnorm(100)),
                   container = container)
}

gui.var.test = function(container=NULL) {
  gui.htest("var.test",type="bivariate",
                   template=var.test(rnorm(100), rnorm(100)),
                   container = container)
}

gui.prop.test = function(container=NULL) {
  gui.htest("prop.test", type="x.and.n",
                   template = prop.test(10,20,p=.5),
                   container=container)
}

gui.prop.test.2sample = function(container=NULL) {
  gui.htest("prop.test", type="x.and.n.2",
                   template = prop.test(c(10,20),c(20,30)),
                   container=container)
}

gui.binom.test  = function(container=NULL) {
  gui.htest("binom.test", type="x.and.n",
                   template = binom.test(10,20,p=.5),
                   container=container)
}

gui.binom.test.2sample = function(container=NULL) {
  gui.htest("binom.test", type="x.and.n.2",
                   template = binom.test(c(10,20),c(20,30)),
                   container=container)
}



##################################################
### Thre functions:
## * one to draw template
## * one to initialize
## * one to update on events

## sequence:
## * define template x
## * pass to HtestTemplate
## * this returns widgets. Populate widgets with defaults
## * on actinos, update widgets


### workhorse
gui.htest = function(
  FUN = "t.test",
  type = c("univariate","bivariate","x.and.n","x.and.n.2","x~f"),
  template = do.call(FUN,list(x=rnorm(100))),
  extra.args = NULL,                    # eg list("var.equal"=TRUE)
  container = NULL,
  ...) {

  type = match.arg(type)

  ##
  x = template
  obj = dHtest(x, type)
##  tag(obj,  "dropHandlers")  <- list()

  passiveComponents = tag(obj,"passiveComponents")
  
  if(!is.null(container))
    add(container,obj, expand=TRUE)

  ## now initialize
  initialize = function(obj, ...) {
    object.list = tag(obj,"object.list")
    if(!type %in% c("x.and.n","x.and.n.2")) {
      svalue(object.list[['data.name']]) <-  "Drop variable(s) here"
    }
    svalue(object.list[["null.value"]]) <-  x$null.value
    svalue(object.list[["alternative"]]) <-  "not equal to"
    svalue(object.list[["statistic.name"]]) <- names(x$statistic)
    svalue(object.list[["parameter.name"]]) <- names(x$parameter)
    svalue(object.list[["null.value.name"]]) <- names(x$null.value)
    ## confidence interval value
    x$conf.int.name = attr(x$conf.int,"conf.level")
    svalue(object.list[["conf.int.name"]]) <-  x$conf.int.name*100
     for(i in passiveComponents)
      object.list[[i]] = glabel("*")
  }

  initialize(obj)

  object.list = tag(obj,"object.list")
  ## add handlers
  if(!type %in% c("x.and.n","x.and.n.2")) {
    adddroptarget(object.list[['data.name']],
                  handler=function(h, ...) {
                    varList = tag(h$obj,"varList")
                    if(!is.list(varList)) 
                      varList = list()
                    n = length(varList)
                    varList[[n+1]] = h$dropdata
                    theName = id(h$dropdata)

                    svalue(object.list[['data.name']]) <-  theName
                    ## do I append here?
                    tag(object.list[['data.name']], "varList") <- varList
                    updatedHtest(obj)
                    ## add handler for drop data
                    ## now bind to be dynamic *if* a treeviewcolumn
                    ## NEED to make thin in gWidgets,
                    if(is.gdataframecolumn(h$dropdata)) {
                      addhandlerchanged(h$dropdata,
                                        handler=function(h,...) updatedHtest(obj)
                                        )
                    }
                  })
    addhandlerchanged(object.list[['data.name']],
                    action = obj,
                    handler=function(h,...) {
                        varList = tag(object.list[['data.name']],"varList")
                        if(!is.list(varList)) 
                          varList = list()
                        n = length(varList)
                        string = svalue(h$obj)
                        ## split on " and "
                        tmp = unlist(strsplit(string," and "))
                        for(i in 1:length(tmp)) {
                          tmp[i] = stripWhiteSpace(tmp[i])
                          if(length(grep("dropvariable",tolower(tmp[i])))==0)
                            varList[[n+i]] = tmp[i]
                        }
                        tag(object.list[['data.name']], "varList") <- varList
                        svalue(object.list[['data.name']]) <- paste(unlist(varList),sep=" and ")
                        updatedHtest(obj)
                      })
  }
  addhandlerchanged(object.list[['null.value']],
                    action = obj,
                    handler = function(h,...) {
                      updatedHtest(h$action)
                    })

  addhandlerchanged(object.list[['alternative']],
                    action = obj,
                    handler = function(h,...) {
                      updatedHtest(h$action)
                    })
  addhandlerchanged(object.list[['conf.int.name']],
                    action = obj,
                    handler = function(h,...) {
                      updatedHtest(h$action)
                    })

                       tag(obj,"FUN") <- FUN
                       tag(obj,"extra.args") <- extra.args
                       tag(obj,"type") <- type

  return(obj)
}


## This is main function
## call function, update widgets
updatedHtest = function(object, ...) {
  obj = object                          # for s3 consistency

  FUN = tag(obj,"FUN")
  extra.args = tag(obj,"extra.args")
  TYPE = tag(obj,"type")
  object.list = tag(obj,"object.list")
  dataVarList = tag(object.list[["data.name"]], "varList")
  x = list()

  h0 = as.numeric(svalue(object.list[["null.value"]]))
  hA = svalue(object.list[["alternative"]])
  if(is.empty(hA))
    hA = "not equal to"
  hA = switch(hA,                     # translate
    "not equal to" = "two.sided",
    "less than" = "less",
    "greater than" = "greater")
  conf.int.name = as.numeric(svalue(object.list[["conf.int.name"]]))/100
  
  
  ## one sample guys
  if(TYPE == "univariate") {

    ## t.test one sample
    n = length(dataVarList)
    theName = id(dataVarList[[n]])
    theValues = svalue(dataVarList[[n]])
    if(length(theValues) <= 1) {
      return()
    }

    theArgs = list(x = theValues, mu=h0, alternative=hA,
      conf.level=conf.int.name)
    x = try(do.call(FUN,c(theArgs,extra.args)), silent=TRUE)
    if(inherits(x,"try-error")) {
      cat("Error with function call:",x,"\n")
    }
    x$data.name = theName             # override
  } else if (TYPE == "bivariate") {
    ## two sample guy
    n = length(dataVarList)
    if(n == 0) {
      cat("Drop some variables\n")
      return()
    } else if(n == 1) {
      ## need another
      cat("Need another\n")
      x$data.name = Paste(id(dataVarList[[1]])," and Drop variable here")
#      tmp = Paste(id(dataVarList[[1]])," and Drop variable here")
#      svalue(object.list[["data.name"]]) <- format(tmp)
    } else  {                           # n > 1
      theXName = id(dataVarList[[n-1]])
      theXValues = svalue(dataVarList[[n-1]])
      theYName = id(dataVarList[[n]])
      theYValues = svalue(dataVarList[[n]])
      if(length(theXValues) <= 1 || length(theYValues) <= 1) {
        return()
      }

      theArgs = list(x = theXValues, y= theYValues,
        mu=h0, alternative=hA,
        conf.level=conf.int.name)

      x = try(do.call(FUN,c(theArgs, extra.args)), silent=TRUE)
      if(inherits(x,"try-error")) {
        cat("Error with function call:",x,"\n")
      }
      x$data.name = Paste(id(dataVarList[[n-1]])," and ",
        id(dataVarList[[n]]))
    }

    } else if (TYPE == "x.and.n") {
    dataVals = object.list[["data.name"]][,,drop=FALSE]
    if(length(dataVals) >=2) {
      theArgs = list(
        x = as.numeric(dataVals$x), n= as.numeric(dataVals$n),
        p=h0, alternative=hA,
        conf.level=conf.int.name)
      x = try(do.call(FUN,c(theArgs,extra.args)), silent=TRUE)
    } else {
      return()
    }
  } else if (TYPE == "x.and.n.2") {
    dataVals = object.list[["data.name"]][,,drop=FALSE] # as matrix
    if(nrow(dataVals) !=2 || ncol(dataVals) !=2 || any(is.na(dataVals))) {
      return()
    } else {
      theArgs = list(
        x = as.numeric(dataVals[,1,drop=TRUE]),
        n = as.numeric(dataVals[,2,drop=TRUE]),
        alternative=hA,
        conf.level=conf.int.name)
      x = try(do.call(FUN,c(theArgs,extra.args)), silent=TRUE)
    }
  } else if(TYPE == "x~f") {

    ## find out which is a factor, which is numeric
    ## Need two variables
    n = length(dataVarList)
    if(n == 0) {
      cat("Drop two variables. One should be a factor.\n")
      return()
    } 

    ## check for handtyped formulas in the dataVarList
    if(n >= 1) {
      typedIn = ""
      if(n >= 2) {                    # try n-1
        if(is.character(dataVarList[[n-1]]) &&
           length(grep(pattern="~",dataVarList[[n-1]]))) {
          typedIn = dataVarList[[n-1]]
        }
      } ## try n now
      if(is.character(dataVarList[[n]]) &&
         length(grep(pattern="~",dataVarList[[n]]))) {
        typedIn = dataVarList[[n]]
      }
      
      ## if typeIn then handle separately
      if(typedIn != "") {
        form = eval(parse(text=typedIn))
        tmp = all.vars(form)
        form = svalue(tmp[1]) ~ svalue(tmp[2])
        x = try(do.call(FUN,list(form)), silent=TRUE)
        if(inherits(x,"try-error")) {
          ## What to do with an error:
          cat("Error with function call:",x,"\n")
        } 
        x$data.name = Paste(tmp[1], " ~ ", tmp[2])
      } else if(n == 1) {
        tmp = Paste(id(dataVarList[[1]])," and drop variable here")
        svalue(object.list[["data.name"]]) <-  format(tmp)
        x = NULL
      } else  {                           # n > 1
        theXName = id(dataVarList[[n-1]])
        theXValues = svalue(dataVarList[[n-1]])
        theYName = id(dataVarList[[n]])
        theYValues = svalue(dataVarList[[n]])
        if(length(theXValues) <= 1 || length(theYValues) <= 1) {
          return()
        }
        
        argList = list()
        theName = ""
        ## find out which is a factor
        if(is.factor(theXValues) && is.factor(theYValues)) {
          theName = "Two factors, need one numeric variable, and one factor"
        } else if (is.factor(theXValues) && is.numeric(theYValues)) {
          argList[['formula']] = theYValues ~ theXValues
          theName = Paste(theYName, " ~ ", theXName)
        } else if (is.numeric(theXValues) && is.factor(theYValues)) {
          argList[['formula']] = theXValues ~ theYValues
          theName = Paste(theXName, " ~ ", theYName)
        } else if (is.numeric(theXValues) && is.numeric(theYValues)) {
          cat("Coercing",theYName,"to be a factor\n")
          argList[['formula']] = theXValues ~ as.factor(theYValues)
          theName = Paste(theXName, " ~ as.factor(", theYName,")")
        } else  {
          theName = "Need one numeric variable and one factor"
        }

        ## don't do this *if* FUN=...
        ## these funs don't take alter or conf.level
        if(FUN %in% c("oneway.test","kruskal.test")) {
          argList$alternative = NULL
          argList$conf.level = NULL
        } else {
          argList$alternative = hA
          argList$conf.level=conf.int.name
        }          

        
        x = try(do.call(FUN,c(argList, extra.args)), silent=TRUE)
        if(inherits(x,"try-error")) {
          ## What to do with an error:
          cat("Error with function call:",x,"\n")
        } 

        x$data.name = theName

      }
    }

  }

  if(is.null(x) || inherits(x,"try-error"))
    return()
  
#    x$statistic.name = names(x$statistic)
#    x$parameter.name = names(x$parameter)
    x$conf.int.name = attr(x$conf.int,"conf.level")
#    x$null.value.name = names(x$null.value)

    plainVals = c(
      "statistic",   
      "parameter",   
      "conf.int",
      "p.value",
      "method",    "estimate"
      )
    for(i in plainVals) {
      svalue(object.list[[i]]) <-  format(x[[i]])
    }
  if(! TYPE %in% c("x.and.n","x.and.n.2"))
    svalue(object.list[["data.name"]]) <- x$data.name

  svalue(object.list[["conf.int"]]) <- 
    Paste("(", format(x$conf.int[1]),",  ",format(x$conf.int[2]),")")
  
  svalue(object.list[["estimate"]]) <- 
    paste(names(x$estimate), format(x$estimate), sep=" = ", collapse=", ")
}

##################################################
##
##
dHtest = function(x, type=NULL, digits = 4, container = NULL, ...) {
  ## x is passed in via ...
  ## this function sets up the widgets and returns the widdgets in an object

  ## "data.name"  
  ## "statistic"   "statistic.name"
  ##" parameter"   "paramter.name"
  ## "null.value"  "null.value.name" "alternative" "method"     
  ## "p.value"     "conf.int"   
  ## "estimate"

  activeComponents = c(
    "data.name"  ,                      # data
    "null.value",                       # H_0
    "alternative",                      # H_A
    "conf.int.name")                         # alpha

  passiveComponents = c(
    "statistic",   "statistic.name",
    "parameter",   "parameter.name",
    "conf.int",
    "null.value.name",
    "p.value",
    "method",    "estimate"
    )

  group = ggroup(horizontal=FALSE, container = container)
#  obj = list(ref=group)
#  class(obj) <- c("dHtest", "gComponent", "gWidget")

  obj = group                           # avoid S3 extension
  
  tag(obj, "activeComponents") <- activeComponents
  tag(obj, "passiveComponents") <- passiveComponents
  
  object.list = list()
  for(i in passiveComponents) {
    object.list[[i]] = glabel("*")
    font(object.list[[i]]) <- c(color="red")
  }
  for(i in activeComponents) {
    object.list[[i]] = glabel("*", editable=TRUE)
    font(object.list[[i]]) <-  c(style="bold")
  }

  ## override
#  object.list[['null.value']] = gedit("")
#  object.list[['data.name']] = gbutton("Drop variables here")

  altVals = c("not equal to","less than","greater than")
  object.list[['alternative']] = gdroplist(altVals,selected=1)

#  object.list[['conf.int.name']] = gdroplist(c("80","95","99"), selected=2, editable=TRUE)

  ## data for prop,.test
  if(type == "x.and.n" || type == "x.and.n.2") {
    aDF = data.frame(x=I(c("","")),n = I(c("","")))
    rownames(aDF) = c("sample 1","sample 2")
    if(type == "x.and.n")
      aDF=aDF[1,]
    object.list[["data.name"]] <- tmp <- gdf(aDF)
    size(tmp) <- c(150, 75)
    addhandlerchanged(object.list[["data.name"]],handler=function(h,...) {
      updatedHtest(obj)
    })
  }
    
  
  ## setup window
#  text = gtext(container=container)
  tag(obj, "object.list") <- object.list
  
#  add(text,"")
#  add(text,x$method, font.attr=c("bold","large"))

  add(group, glabel(Paste("<b><i>",x$method,"</i></b>"), markup=TRUE))

  newLine = ggroup(container=group)
  add(newLine,glabel("data: "))
  add(newLine,object.list[["data.name"]], expand=TRUE) # need expand=TRUE for some widgets
  
  if(!is.null(x$statistic)) {
    newLine = ggroup(container=group)
    add(newLine, object.list[["statistic.name"]])
    add(newLine,glabel(" = "))
    add(newLine, object.list[["statistic"]])
  }

  if(!is.null(x$parameter)) {
    newLine = ggroup(container=group)
    add(newLine,object.list[["parameter.name"]])
    add(newLine,glabel(" = "))
    add(newLine, object.list[["parameter"]])
  }
  if(!is.null(x$p.value)) {
    newLine = ggroup(container=group)
    fp <- format.pval(x$p.value, digits = digits)
    inequality = if(substr(fp,1,1) == "<") "" else " = "
    add(newLine, glabel(Paste("p-value ", inequality)))
    add(newLine, object.list[["p.value"]])
  }

  if(!is.null(x$alternative)) {
    newLine = ggroup(container=group)
    add(newLine,glabel("alternative hypothesis: "))
    newLine = ggroup(container=group)
    add(newLine,glabel("   "))          # format
    if(!is.null(x$null.value)) {
      if(length(x$null.value) == 1) {
        alt.char <-
          switch(x$alternative,
                 two.sided = " not equal to ",
                 less = " less than ",
                 greater = " greater than ")
        add(newLine, glabel("true "))
        add(newLine, object.list[["null.value.name"]])
        add(newLine, glabel(" is "))
        
        add(newLine, object.list[["alternative"]])
        add(newLine, object.list[["null.value"]])
      } else {
        add(newLine, glabel(Paste(x$alternative, "\nnull values:\n")))
        add(newLine, glabel(paste(x$null.value,collapse="\t")))
      }
    } else {
      add(newLine, object.list[["alternative"]])
    }
  }

  if(!is.null(x$conf.int)) {
    newLine = ggroup(container=group)
    add(newLine,object.list[["conf.int.name"]])
    add(newLine,glabel(" percent confidence interval:"))
    newLine = ggroup(container=group)
    add(newLine,glabel("    "))         # space
    add(newLine,object.list[["conf.int"]])
  }
  if(!is.null(x$estimate)) {
    newLine = ggroup(container=group)
    add(newLine,glabel("sample estimates:"))
    newLine = ggroup(container=group)
    add(newLine,glabel("    "))         # space
    add(newLine, object.list[["estimate"]])
  }

  
  return(obj)
  
}

getGTKwidget.dHtest = function(obj, ...) obj$ref$ref


## dHtest = function(x, digits = 4, container = NULL, ...) {
##   ## x is passed in via ...
##   ## this function sets up the widgets and returns the widdgets in an object

##   ## "data.name"  
##   ## "statistic"   "statistic.name"
##   ##" parameter"   "paramter.name"
##   ## "null.value"  "null.value.name" "alternative" "method"     
##   ## "p.value"     "conf.int"   
##   ## "estimate"

##   activeComponents = c(
##     "data.name"  ,                      # data
##     "null.value",                       # H_0
##     "alternative",                      # H_A
##     "conf.int.name")                         # alpha

##   passiveComponents = c(
##     "statistic",   "statistic.name",
##     "parameter",   "parameter.name",
##     "conf.int",
##     "null.value.name",
##     "p.value",
##     "method",    "estimate"
##     )

##   object.list = list()
##   for(i in passiveComponents)
##     object.list[[i]] = glabel("*")
##   for(i in activeComponents)
##     object.list[[i]] = gbutton("")

##   ## override
##   object.list[['null.value']] = gedit("")
##   object.list[['alternative']] = gdroplist(c("two.sided","less","greater"),
##                selected=1)
##   object.list[['conf.int.name']] = gdroplist(c("80","95","99"), selected=2, editable=TRUE)
               


  
##   ## setup window
##   text = gtext(container=container)
##   obj = list(ref=text)
##   .class(obj) <- c("dHtest", "gComponent", "gWidget")
##           tag(obj, "object.list") <- object.list
  
##   add(text,"")
##   add(text,x$method, font.attr=c("bold","large"))

##   add(text,"data: ",do.newline=FALSE)
## #  add(text,x$data.name, font.attr=c("red"))
##   add(text,object.list[["data.name"]])
##   add(text, "")
  
##   if(!is.null(x$statistic)) {
##     add(text, object.list[["statistic.name"]], do.newline=FALSE)
##     add(text," = ", do.newline=FALSE)
##     add(text, object.list[["statistic"]])
##     add(text, "")
##   }

##   if(!is.null(x$parameter)) {
##     add(text,object.list[["parameter.name"]], do.newline=FALSE)
##     add(text," = ", do.newline=FALSE)
##     add(text, object.list[["parameter"]])
##     add(text, "")
##   }
##   if(!is.null(x$p.value)) {
##     fp <- format.pval(x$p.value, digits = digits)
##     inequality = if(substr(fp,1,1) == "<") "" else " = "
##     add(text, Paste("p-value ", inequality), do.newline=FALSE)
##     add(text, object.list[["p.value"]])
##     add(text, "")
##   }

##   if(!is.null(x$alternative)) {
##     add(text,"alternative hypothesis: ")
##     if(!is.null(x$null.value)) {
##       if(length(x$null.value) == 1) {
##         alt.char <-
##           switch(x$alternative,
##                  two.sided = " not equal to ",
##                  less = " less than ",
##                  greater = " greater than ")
##         add(text, "true ", do.newline=FALSE)
##         add(text, object.list[["null.value.name"]],  do.newline=FALSE)
##         add(text, " is ",  do.newline=FALSE)
        
##         add(text, object.list[["alternative"]], do.newline=FALSE)

##         add(text, "   ", do.newline=FALSE) # breathe

##         add(text, object.list[["null.value"]])
##       } else {
##         add(text, Paste(x$alternative, "\nnull values:\n"))
##         add(text, paste(x$null.value,collapse="\t"))
##       }
##     } else {
##       add(text, object.list[["alternative"]])
##     }
##     add(text, "")
##   }

## if(!is.null(x$conf.int)) {
    
##     add(text,
##         object.list[["conf.int.name"]],
##         do.newline=FALSE)
##     add(text," percent confidence interval:")
##     add(text,object.list[["conf.int"]])
##     add(text, "")
##   }
##   if(!is.null(x$estimate)) {
##     add(text,"sample estimates:")
##     add(text, object.list[["estimate"]])
##   }

  
##   return(obj)
  
## }


## print.htest <- function(x, digits = 4, quote = TRUE, prefix = "", ...)
## {
##     cat("\n")
##     writeLines(strwrap(x$method, prefix = "\t"))
##     cat("\n")
##     cat("data: ", x$data.name, "\n")
##     out <- character()
##     if(!is.null(x$statistic))
## 	out <- c(out, paste(names(x$statistic), "=",
## 			    format(round(x$statistic, 4))))
##     if(!is.null(x$parameter))
## 	out <- c(out, paste(names(x$parameter), "=",
## 			    format(round(x$parameter, 3))))
##     if(!is.null(x$p.value)) {
## 	fp <- format.pval(x$p.value, digits = digits)
## 	out <- c(out, paste("p-value",
## 			    if(substr(fp,1,1) == "<") fp else paste("=",fp)))
##     }
##     writeLines(strwrap(paste(out, collapse = ", ")))
##     if(!is.null(x$alternative)) {
## 	cat("alternative hypothesis: ")
## 	if(!is.null(x$null.value)) {
## 	    if(length(x$null.value) == 1) {
## 		alt.char <-
## 		    switch(x$alternative,
## 			   two.sided = "not equal to",
## 			   less = "less than",
## 			   greater = "greater than")
## 		cat("true", names(x$null.value), "is", alt.char,
## 		    x$null.value, "\n")
## 	    }
## 	    else {
## 		cat(x$alternative, "\nnull values:\n")
## 		print(x$null.value, ...)
## 	    }
## 	}
## 	else cat(x$alternative, "\n")
##     }
##     if(!is.null(x$conf.int)) {
## 	cat(format(100 * attr(x$conf.int, "conf.level")),
## 	    "percent confidence interval:\n",
## 	    format(c(x$conf.int[1], x$conf.int[2])), "\n")
##     }
##     if(!is.null(x$estimate)) {
## 	cat("sample estimates:\n")
## 	print(x$estimate, ...)
##     }
##     cat("\n")
##     invisible(x)
## }
