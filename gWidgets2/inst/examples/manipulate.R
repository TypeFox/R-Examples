## An implementation of the `manipulate` package of RStudio
## http://rstudio.org/docs/advanced/manipulate

## This makes it very easy to create interactive graphics.

## helpers
picker <- function(..., initial=1, label=NULL) {
  items <- list(...)
  if(!is.null(names(items)))
    items <- data.frame(value=unlist(items), label=names(items), stringsAsFactors=FALSE)
  else
    items <- unlist(items)
  list(FUN="gcombobox", items=items, selected=initial, label=label)
}

radio <- function(..., initial=1, label=NULL) {
  items <- unlist(list(...))
  list(FUN="gcombobox", items=items, selected=initial, label=label)
}


slider <- function(min, max, initial=min, step=1, ..., label=NULL)
  list(FUN="gslider",from=min, to=max, by=step, value=initial, label=label)

checkbox <- function(initial=FALSE,  label=NULL, ...)
  list(FUN="gcheckbox", checked=initial, label=label)

button <- function(label)
  list(FUN="gbutton", text=label)


## main function
manipulate <- function(.expr, ...,
                       container=gwindow("Manipulate")
                       ) {
  
  expr <- substitute(.expr)

  l <- list(...)
  pg <- gpanedgroup(container=container)
  gg <- ggroup(container=pg, expand=TRUE, fill=TRUE)
  flyt <- gframe("Controls", container=pg, horizontal=FALSE)

  if(gtoolkit() == "tcltk") {
    require(tkrplot)
    dev <- tkrplot(getToolkitWidget(gg), function() {})
    add(gg, dev, expand=TRUE)
  } else {
    ggraphics(cont=gg, expand=TRUE)
  }
  
  update_expr <- function(...) {
    values <- sapply(widgets, svalue, simplify=FALSE)
    if(gtoolkit() == "tcltk")
      tkrplot:::.my.tkdev(dev$hscale, dev$vscale)
    
    result <- withVisible(eval(expr, envir=values))
    if (result$visible) {
      eval(print(result$value))
    }

    if(gtoolkit() == "tcltk")
      .Tcl(paste("image create Rplot", dev$image))
  }
  
  make_widget <- function(nm, lst) {
    ## label
    if(lst[1] != "gbutton") {
      if(is.null(lst$label)) lst$label <- nm
      l <- glabel(paste(lst$label, ":", sep=""), cont=flyt)
      font(l) <- list(weight="bold")
    }
    ## the widget
    lst$handler <- update_expr
    lst$container <- flyt
    do.call(lst$FUN, lst[-1])
  }

  widgets <- mapply(make_widget, names(l), l, SIMPLIFY=FALSE)
  addSpring(flyt)                          # nicer layout

  update_expr()
  invisible(widgets)
}


## w <- gwindow("Manipulate example", visible=FALSE)
## manipulate({hist(rnorm(n), main=sprintf("a Histogram of %s", n))},
##            cb = picker(one=1, two=2, three=3, initial=2),
##            n=radio(c(10, 20, 100), initial=2),
##            sl = slider(0, 100, 0, 5),
##            check = checkbox(TRUE),
##            btn = button("hi there"),
##            container=w)
## visible(w) <- TRUE
