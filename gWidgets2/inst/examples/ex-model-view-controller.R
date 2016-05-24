about <- "
Even a simple form can have complicated interactions between the components.
A standard design pattern is the model-view-controller pattern. In the following
we implement this.

In this example, we create a simple form for a t.test. The user must
select a data frame, a variable `x`, and optionally a variable `y`. In addition
they may adjust various parameters.

The selection of variables is a bit tricky. The data frame is specified as a string,
but this must refer to a data frame.

The variable x must be a valid expression
within the data frame. As such, it a) shouldn't be visible until the
data frame is valid b) it should evaluate appropriately.

In addition to the restrictions on `x`, the variable `y` must be either a numeric vector
or a factor with 2 levels. 

The choice of `mu` or `conf.level` indicates if the user wants a confidence interval or
a significance test. Though R will output both, we are trying to simplify the output
so restrict to a choice of one only. If there is a value for `mu` we disable the slider.

All this can be complicated to coordinate. We will use the `objectProperties` package
by Tengfei Yin and Michael Lawrence to store our values.

This package provides several advantages:

* it allows us to set restrictions on what values are valid. See `alternative ` and `conf.level`
  for examples, though it isn't essential here.

* it allows us to listen for changes to the model through the `connect` method. We simply pass
  the changes in the GUI back to the model (see `update_var`) and put our controller logic
  through the `connect` method.

In the `TTestModel` class, we use object properties to create a model with validation and with
some model-based methods.
"



## The "model"
## objectProperties stuff
require(objectProperties)

setSingleEnum("Alternative",
              levels=c("two.sided", "less", "greater"))

setNumericWithRange("Numeric", min=0, max=1)

## A new class
setClass("DataFrameName", contains = c("character"),
         validity = function(object) {
           if(object != "" &&
              (!exists(object, .GlobalEnv) ||
               !is(get(object, .GlobalEnv), "data.frame"))
              )
             gettext("Value is not a data frame object in the global workspace")
         })

TTestModel <- setRefClass("TTestModel",
                          fields=properties(
                            fields=list(
                              dataframe="DataFrameName",
                              x="character",
                              y="character",
                              alternative="AlternativeSingleEnum",
                              mu="numeric",
                              conf.level="NumericWithMin0Max1"
                              ),
                            prototype=list(
                              dataframe="", x="", y="",
                              alternative=new("AlternativeSingleEnum", "two.sided"),
                              mu=as.numeric(NA),
                              conf.level=new("NumericWithMin0Max1", 0.95)
                              )
                            ),
                          contains = "PropertySet",
                          methods=list(
                            can_run=function() {
                              ## can we run a t-test?
                              nchar(dataframe) && nchar(x)
                            },
                            run_t_test = function(...) {
                              if(!can_run()) return()
                              out <- sapply(names(properties()), get, envir=.self, simplify=FALSE)
                              out$data <- get(out$dataframe, .GlobalEnv); out$dataframe <- NULL
                              
                              out$x <- eval(parse(text=out$x), envir=out$data)

                              if(nchar(out$y)) {
                                out$y <- eval(parse(text=out$y), envir=out$data)
                                if(is.factor(out$y)) {
                                  out$formula <- formula(sprintf("%s ~ %s", model$x, model$y))
                                  out$x
                                  out$y <- NULL
                                }
                              } else {
                                out$y <- NULL
                              }

                              
                              if(is.na(out$mu))
                                out$mu <- NULL

                              res <- do.call("t.test", out)
                              
                              ## what do we output?
                              if(!is.null(out$mu))
                                print(sprintf("The p-value is %s", res$p.value))
                              else
                                print(sprintf("The confidence interval is (%s, %s)", res$conf.int[1], res$conf.int[2]))
                            }
                          ))
model <- TTestModel$new()


##################################################
## The "View"
library(gWidgets2)
#options(guiToolkit="Qt")

w <- gwindow("t-test example", visible=FALSE)
sb <- gstatusbar(cont=w)
g <- gvbox(cont=w)
g$set_borderwidth(10)
fl <- gformlayout(cont=g)

l <- list()
l$dataframe <- gedit("", initial="data frame name", cont=fl, label="Data frame")
l$x <- gedit("", initial="x variable", cont=fl, label="x")
l$y <- gedit("", initial="y variable", cont=fl, label="y")
sapply(l[c('x','y')], function(i) enabled(i) <- FALSE)

l$alternative <- gcombobox(model$alternative@levels, cont=fl, label="alternative")
l$mu <- gedit("", initial="mean in H_0", cont=fl, label="mu", coerce.with=as.numeric)
## tcltk only works on integers! Here we multiply and divide
l$conf.level <- gslider(from=100*model$conf.level@min, to=100*model$conf.level@max,
                        by=100*0.01, value=100*model$conf.level,
                        cont=fl, label="conf.level"
                        )
l$conf.level$coerce_with <- function(x) x/100


## add buttons
bg <- ggroup(cont=g)
addSpring(bg)
gbutton("About", cont=bg, handler=function(h,...) {
  w1 <- gwindow("About", visible=FALSE, parent=w)
  g <- gvbox(cont=w1); g$set_borderwidth(10)
  glabel(about, cont=g, expand=TRUE)
  bg <- ggroup(cont=g); addSpring(bg)
  gbutton("dismiss", cont=bg, handler=function(...) dispose(w1))
  visible(w1) <- TRUE
})
run_button <- gbutton("Run", cont=bg, handler=model$run_t_test)
enabled(run_button) <- FALSE

visible(w) <- TRUE


##################################################
## The "Controller"
## connect things up. Model should not know about view and vice versa. Controllers do
## these connections

## view -> model
update_var <- function(h, ...) {
  out <- try(assign(h$action, svalue(h$obj), model), silent=TRUE)
  ## an error set invalid
  h$obj$set_invalid(inherits(out, "try-error"), out)
}
mapply(addHandlerChanged, l, handler=list(update_var), action=names(l))

## model changes -> view

## When data frame changes we update x and y values.
## changes is called on enter (not tab?)
model$dataframeChanged$connect(function() {
  ## Change to data frame updates x and y values
  df_name <- model$dataframe
  if(exists(df_name, .GlobalEnv) && is(DF <- get(df_name, .GlobalEnv), "data.frame")) {
    ## update, clear out x and y
    svalue(sb) <- ""
    sapply(l[c("x", "y")], function(widget) {
      svalue(widget) <- ""
      widget[] <- names(DF)
      enabled(widget) <- TRUE
    })
    ## move focus
    focus(l$x) <- TRUE
  } else {
    svalue(sb) <- sprintf("Data frame %s not found", df_name)
    sapply(l[c("x", "y")], function(widget) {
      svalue(widget) <- ""
      widget[] <- character(0)
      enabled(widget) <- FALSE
    })
  }
})

## y must evaluate to a numeric variable or 2-level factors
## when specified, mu does not make sense
model$yChanged$connect(function() {
  df_name = model$dataframe
  DF <- get(df_name, .GlobalEnv)
  val <- model$y
  if(val == "")
    return()

  out <- eval(parse(text=val), envir=DF)
  if(inherits(out, "try-error") ||
     (!is.numeric(out) && (is.factor(out) && length(levels(out)) != 2))) {
    svalue(sb) <- gettext("y value is invalid")
  } else {
    svalue(l$mu) <- ""
    svalue(sb) <- ""
  }
})

## if mu is specified, then conf.level not needed
model$muChanged$connect(function() {
  enabled(l$conf.level) <- model$mu == "" || !is.numeric(model$mu)
})

## If we can run we enable run button
model$changed$connect(function(name) {
  enabled(run_button) <- model$can_run()
})
