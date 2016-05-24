
#############################################################################
##
##  Stage 1: Select type of distribution and generation method
##
#############################################################################
##
##  Synopsis:
##
##    stage1 ( main )
##
##  Arguments:
##    main ... main window
##
##  Frames:
##    'Distributions'
##    'Generation Methods'
##
##  All data are stored in 'main'.
##
##  Actions (buttons)
##
##  'ok'
##     Remove all widgets from 'main' and store information about
##     distrbution and generation method in 'main' in the respective
##     lists 'distribution' and 'method'.
##     See end of function select.type() for elements of these lists.
##     At last stage 2 is started.
## 
##  'cancel'
##     Quit and destroy window.
##
##  'restart'
##     Reset all entries.
##
##
#############################################################################

## --------------------------------------------------------------------------

stage1 <- function(main) {

  ## change window title
  ## svalue(main) <- "Runuran - Select Type of Distribution and Generation Method"
  
  ## pack everything in a group
  group <- ggroup(horizontal=FALSE, spacing=10, container=main)
  tag(main,"stage1") <- group

  glabel("Stage 1:  Select type of distribution and generation method",
         container=ggroup(horizontal=TRUE,container=group))

  type.distribution(main,group)
  type.method(main,group)
  type.buttons(main,group)
}


## --------------------------------------------------------------------------
## Select type of distribution

type.distribution <- function(main,group) {

  ## the frame
  frame <- gframe(text=" Distribution ", horizontal=TRUE, container=group)

  ## radio box: select type
  type.rbx <- gradio(DISTRIBUTIONS.TYPE, container=frame,
                     handler=type.update, action=main)
  
  addSpace(frame, 5)
  gseparator(horizontal=FALSE, container=frame)
  addSpace(frame, 5)

  ## radio box: select source
  howdef.rbx <- gradio(DISTRIBUTIONS.HOWDEF, container=frame,
                       handler = type.update, action=main)

  addSpace(frame, 5)
  gseparator(horizontal=FALSE, container=frame)
  addSpace(frame, 5)

  ## droplist: continuous distributions
  cont.cbb <- gcombobox(DISTRIBUTIONS[["continuous"]]["name",], container=frame)
  tag(cont.cbb, "added") <- TRUE

  ## droplist: discrete distributions (not added by default)
  discr.cbb <- gcombobox(DISTRIBUTIONS[["discrete"]]["name",], container=frame)
  tag(discr.cbb, "added") <- FALSE
  delete(frame, discr.cbb)

  ## store all widgets in main window
  tag(main,"distr.frame")      <- frame
  tag(main,"distr.type.rbx")   <- type.rbx
  tag(main,"distr.howdef.rbx") <- howdef.rbx
  tag(main,"distr.cont.cbb")   <- cont.cbb
  tag(main,"distr.discr.cbb")  <- discr.cbb
}

## ..........................................................................
## handler: update state

type.update.distribution <- function(main) {

  ## read data
  frame      <- tag(main,"distr.frame")
  howdef     <- svalue( tag(main,"distr.howdef.rbx"))
  type       <- svalue( tag(main,"distr.type.rbx"))
  cont.cbb   <- tag(main,"distr.cont.cbb")
  discr.cbb  <- tag(main,"distr.discr.cbb")

  ## get new state
  if (howdef == "built-in") {
    if (type == "continuous")
      new.state <- c(TRUE,FALSE)
    else
      new.state <- c(FALSE,TRUE)
  } else {
    new.state <- c(FALSE,FALSE)
  }

  ## update widgets
  toggle.group(frame, cont.cbb, new.state[1])
  toggle.group(frame, discr.cbb,new.state[2])
}   


## --------------------------------------------------------------------------
## Select generation method

type.method <- function(main,group) {

  ## the frame
  frame <- gframe(text=" Generation Method ", horizontal=FALSE, container=group)

  ## radio box: type of method
  type.rbx <- gradio(METHODS.TYPE, horizontal=TRUE, container=frame,
                     handler = type.update, action=main)
  
  ## droplist: methods for continuous distributions
  cont.cbb <- gcombobox(METHODS[["continuous"]]["name",], container=frame)
  tag(cont.cbb, "added") <- FALSE
  delete(frame, cont.cbb)

  ## droplist: methods for discrete distributions
  discr.cbb <- gcombobox(METHODS[["discrete"]]["name",], container=frame)
  tag(discr.cbb, "added") <- FALSE
  delete(frame, discr.cbb)

  ## store all widgets in main window
  tag(main,"method.frame")     <- frame
  tag(main,"method.type.rbx")  <- type.rbx
  tag(main,"method.cont.cbb")  <- cont.cbb
  tag(main,"method.discr.cbb") <- discr.cbb
}

## ..........................................................................
## handler: update state

type.update.method <- function(main) {

  ## read data
  frame       <- tag(main,"method.frame")
  method.type <- svalue( tag(main,"method.type.rbx"))
  distr.type  <- svalue( tag(main,"distr.type.rbx"))
  cont.cbb    <- tag(main,"method.cont.cbb")
  discr.cbb   <- tag(main,"method.discr.cbb")
  
  ## get new state
  if (method.type == "Select method") {
    if (distr.type == "continuous")
      new.state <- c(TRUE,FALSE)
    else
      new.state <- c(FALSE,TRUE)
  }
  else {
    new.state <- c(FALSE,FALSE)
  }
  
  ## update widgets
  toggle.group(frame, cont.cbb, new.state[1])
  toggle.group(frame, discr.cbb,new.state[2])
}   


## --------------------------------------------------------------------------
## Buttons: Select 

type.buttons <- function(main,group) {

  ## the group
  buttons.grp <- ggroup(horizontal=TRUE, spacing=15, container=group)

  gbutton(action=gaction(label="Restart", icon="new",
            handler=function(h,...){type.clearup(main); stage1(main)}),
          container=buttons.grp)

  gbutton("cancel", container=buttons.grp,
          handler=function(h,...){dispose(main)})

  gbutton("help", container=buttons.grp,
          handler=show.help, action=list("stage1"))

  gbutton("ok", container=buttons.grp,
          handler=function(h,...){type.evaluate(main)})

}


## --------------------------------------------------------------------------
## Evaluate input and proceed to stage 2

type.evaluate <- function(main) {

  ## read data
  distr.type   <- svalue( tag(main,"distr.type.rbx") )
  distr.howdef <- svalue( tag(main,"distr.howdef.rbx") )
  distr.list   <- list(continuous=svalue( tag(main,"distr.cont.cbb"),  index=TRUE ),
                       discrete  =svalue( tag(main,"distr.discr.cbb"), index=TRUE ))
  
  method.type  <- svalue( tag(main,"method.type.rbx") )
  method.list  <- list(continuous=svalue( tag(main,"method.cont.cbb"), index=TRUE ),
                       discrete  =svalue( tag(main,"method.discr.cbb"),index=TRUE ))
  
  ## check settings
  if (distr.howdef == "built-in" && distr.list[[distr.type]] == 1) {
    error.message("No built-in distribution specified!")
    return()
  }

  if (method.type == "Select method" && method.list[[distr.type]] == 1) {
    error.message("No generation method specified!")
    return()
  }
  
  ## ........................................................................
  ## get data about distribution

  ## name
  distr.name <-
    if (distr.howdef == "built-in") 
      DISTRIBUTIONS[[distr.type]]["name", distr.list[[distr.type]] ]
    else
      "User-defined Distribution"

  ## symbol
  distr.symbol <- DISTRIBUTIONS[[distr.type]]["symbol", distr.list[[distr.type]] ]

  ## constructor
  distr.constr <-
    switch(distr.howdef,
           "built-in" = paste("ud",distr.symbol,sep=""),
           "user-defined" = paste("unuran.",sub("(inuous|ete)","",distr.type,perl=TRUE),".new",sep=""),
           internal.error()) 
  
  ## ........................................................................
  ## get data about method

  ## name and symbol
  if (method.type == "Select method") {
    method.name   <- METHODS[[distr.type]]["name", method.list[[distr.type]] ]
    method.symbol <- METHODS[[distr.type]]["symbol", method.list[[distr.type]] ]
  }
  else {
    type <- paste(distr.type,".auto", sep="")
    method.type  <- svalue( tag(main,"method.type.rbx"), index=TRUE )
    method.name   <- METHODS[[type]]["name", method.type]
    method.symbol <- METHODS[[type]]["symbol", method.type]
  }
  ## constructor
  method.constr <- paste(tolower(method.symbol),"d.new",sep="")

  ## ........................................................................

  ## clearup working space
  type.clearup(main)
  
  ## store information
  tag(main,"distribution") <-
    list(type=distr.type, howdef=distr.howdef,
         name=distr.name, symbol=distr.symbol, constr=distr.constr)
  tag(main,"method") <-
    list(type=distr.type,
         name=method.name, symbol=method.symbol, constr=method.constr)
  
  ## start stage 2
  stage2(main)
}


## --------------------------------------------------------------------------
## common calls for all frames in stage 1

type.update <- function(h,...) {
  main <- h$action
  type.update.distribution(main)
  type.update.method(main)
}

type.clearup <- function(main) {

  ## delete all widgets in window
  delete(main, tag(main,"stage1"))
  tag(main,"stage1")           <- NULL

  ## remove data on distribution
  tag(main,"distr.frame")      <- NULL
  tag(main,"distr.type.rbx")   <- NULL
  tag(main,"distr.howdef.rbx") <- NULL
  tag(main,"distr.cont.cbb")   <- NULL
  tag(main,"distr.discr.cbb")  <- NULL

  ## remove data on method
  tag(main,"method.frame")     <- NULL
  tag(main,"method.type.rbx")  <- NULL
  tag(main,"method.cont.cbb")  <- NULL
  tag(main,"method.discr.cbb") <- NULL
}

## --------------------------------------------------------------------------
