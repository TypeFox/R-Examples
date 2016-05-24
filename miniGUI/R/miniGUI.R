# miniGUI.R
#   Intento de GUI
# NOTAS:
# ERRORES:


##
##  load tcltk
##
require(tcltk) || stop("tcltk support is absent")


##
##  Some inits...	
##
##  .onLoad <- function(libname and pkgname;
miniGUIEnvir <- new.env()
miniGUIEnvir$miniGUIData <- list()
miniGUIEnvir$miniGUIans <- NA        ## last result


##
##  Set and get data from miniGUI
##
setMiniGUIData <- function(NAME,PART=NULL,val) 
{
  xxx <- miniGUIEnvir$miniGUIData
  if( !is.null(PART) )
    xxx[[NAME]][[PART]] <- val
  else
    xxx[[NAME]] <- val
  assign("miniGUIData",xxx,envir=miniGUIEnvir)
}


getMiniGUIData <- function(NAME,PART=NULL) 
{
  if( !is.null(PART) )
    return( miniGUIEnvir$miniGUIData[[NAME]][[PART]] )
  else
    return( miniGUIEnvir$miniGUIData[[NAME]] )
}


setMiniGUIans <- function(val) 
{
  assign("miniGUIans",val,envir=miniGUIEnvir)
}


getMiniGUIans <- function() 
{
    return( miniGUIEnvir$miniGUIans )
}





##
##  evaluation procedures
##
## miniGUIcallEval <- function(f,p)
## ok
# f = function to evaluate.
# p = params.
# RETURN:
#   The evaluation of f on params p (string) using envir of f
## {
##   ## NOTICE, all arguments are evaluated in environment(f) !!
##   tryCatch(
##     expr=do.call( f,
##                   lapply(p,
##                        function(x) eval(parse(text=x),envir=environment(f))),
##                   envir=environment(f)
##     ),
##     error=function(e) e
##   )
## }


## miniGUIcallEvalB <- function(f,p)
## ## ok
## # f = function to evaluate.
## # p = params.
## # RETURN:
## #   The call, with function f. Doest not work for real calls
## {
##   ## NOTICE, all arguments are evaluated where the call is made !!
##   do.call( f,
##            lapply(p,function(x) eval(parse(text=x))),
##            envir=environment(f)
##   )
## }

miniGUIcallEval <- function(f,p,e=environment(f))
## ok
# f = function to evaluate.
# p = params.
# e = environment where the parameters p of f are evaluated. By default .GlobalEnv 
# RETURN:
#   The evaluation of f on params p (string) using envir of f
{
  ## NOTICE, all arguments are evaluated in .GlobalEnv if e=.GlobalEnv!!
  ## NOTICE, all arguments are evaluated in environment(f) if e=environment(f)!!
  tryCatch(
    expr=do.call( f,
                  lapply( p, function(x) eval(parse(text=x), envir=e) ),
                  envir=environment(f)
    ),
    error=function(err) err
  )
}


##  set miniGUI evaluation procedure
miniGUIeval <- miniGUIcallEval


##  set miniGUI output procedure
miniGUIoutput <- function(x,mess="\nminiGUI output: \n")
{
  cat(mess)
  print(x)
}


##
##  Building the stuff
##
miniGUIgetFormals <- function(f)
## ok
# f = R function
# RETURN:
#   Filters arguments and names looking for ellipsis.
{
  res <- formals(f)
  ## avoid ellipsis param
  res <- res[names(res)!="..."]
  ## avoid arguments with ellipsis params
  elFun <- grep("[.][.][.]",as.character(res))
  if(length(elFun)>0)
    res <- res[-elFun]    
  return( res )
}


mapFuncToWidget <- function(f,frm,bttLabel="OK",STORE="ff",
                            callSubst="mini GUI call")
## ok
# f = function to display(params are labels and entry fields).
# frm = a frame where toplay trhe display.
# bttLabel = button label "OK",
# STORE = a slot in miniGUIData where to store function param vals.  
# callSubst = a substitute for the call slot in the result
#   Makes window widget to permform function computation
{
  ## get args
  ff <- miniGUIgetFormals(f)
  ## miniGUIData[[STORE]] <<- as.list(ff)
  setMiniGUIData(STORE,val=as.list(ff))
  ## widget GUI input lock 
  ## miniGUIData[["WIDGETLOCK"]] <<- TRUE
  setMiniGUIData("WIDGETLOCK",val=TRUE)
  ## tcltk stuff starts
  argsFrame <- tkframe(frm,borderwidth=2)
  fm <- tkframe(argsFrame, relief="groove", borderwidth=2)
  for(i in names(ff))
  {
    ## get parama value
    par <- deparse(ff[[i]])
    parEval <- eval(parse(text=par))  ##TODO I think we do not need this
    ## create input method widget
    if( is.miniGUIwidget(parEval) ){
      ## miniGUIData[[STORE]][[i]] <<- tclVar( par )
      setMiniGUIData(STORE,i,val=tclVar( par ))
      imw <- parEval$widget(fm,STORE,i)  
    }else{ ## text entry is the default input method widget
      ## miniGUIData[[STORE]][[i]] <<- tclVar( par )
      setMiniGUIData(STORE,i,val=tclVar( par ))
      ## imw <- miniGUIdefaultEntry(fm, textvariable=miniGUIData[[STORE]][[i]])
      imw <- miniGUIdefaultEntry(fm, textvariable=getMiniGUIData(STORE,i))
    }
    ## add the imput widget
    tkgrid(tklabel(fm, text=i), tklabel(fm, text="="),imw)
  }
  ## release GUI widget input lock
  ## miniGUIData[["WIDGETLOCK"]] <<- FALSE
  setMiniGUIData("WIDGETLOCK",val=FALSE)
  ## add execution button
  mainJob <- function (...)
  {
    ## miniGUIans <<- miniGUIeval(f,lapply(miniGUIData[[STORE]],tclvalue))
    setMiniGUIans( miniGUIeval(f,lapply(getMiniGUIData(STORE),tclvalue)) ) 
    if( "call" %in% names(getMiniGUIans()) ) 
    {
      ##miniGUIans$call <- callSubst 
      xxx <- getMiniGUIans()
      xxx$call <- callSubst
      assign("miniGUIans",xxx,envir=miniGUIEnvir)
    }
    ## show result
    miniGUIoutput( getMiniGUIans() )
  }
  ## add and pack the frames
  tkgrid(tkbutton(fm,text=bttLabel,command=mainJob))
  tkpack(fm,fill="x",expand=TRUE)
  tkpack(argsFrame)
  return(argsFrame)
}


makeWidgetCmd <- function(frmTitle,fun,baseFrame=.TkRoot,STORE="ff",GRAB=TRUE)
## ok
# frmTitle = frame title.
# fun = function to make menu command.
# baseFrame = base frame, if not suppied, it creates a stand alone window.
# STORE = Where to store call params
# GRAB = grab input widget frame, disable input from any other frame
#   Makes menu command.
{
  ## to avoid lazy eval
  fun 
  frmTitle
  ## real stuff
  res <- function()
  {
    frm <- tktoplevel(baseFrame)
    if(GRAB)tkgrab(frm) ## enable input only in this frame...
    tkwm.title(frm,paste("mini GUI:",frmTitle))
    mapFuncToWidget(fun, frm, "OK", STORE)
    quitCmd <- function()
    {
      ## Remove function storage from miniGUIData
      ## miniGUIData[[STORE]] <<- NULL
      setMiniGUIData(STORE,val=NULL)
      ## When destroying, main frame is again enabled(ungrabbed !!)
      tkdestroy(frm)
    }
    tkpack( tkbutton(frm,text=paste("Quit",frmTitle),command=quitCmd) )
  }
  return(res)
}


addMenusCmd <- function(cmdFuns,baseFrame)
## ok
# cmdFuns = command functions to add to menu
# baseFrame = base frame
#   adds functins to a menu.
{
  if(!is.null(cmdFuns))
  {
    opsMenu <- tkmenu(tkmenu(baseFrame),tearoff=TRUE)
    for(i in names(cmdFuns))
    {
      tkadd(opsMenu, "command", label=i, command=cmdFuns[[i]])
    }
    tkpopup(opsMenu,tkwinfo("pointerx", baseFrame), 
        	    tkwinfo("pointery", baseFrame))
  }
}








##
##  mini GUI 
##
# miniGUIffff <<- NA      ## storage.
# miniGUIBase <<- NA      ## main tcltk frame storage.
miniGUI <- function(mainFrameFun=evalPlugin,opFuns=NULL,title="mini GUI",
                    init=function() {},WRAPFUN=TRUE)
## ok
#  mainFrameFun = function to display(params are labels and entry fields) on
#    the main window frame or NULL.
#  opFuns= List of functions to add in the menu Ops
#  title = main window frame title
#  init = an init function to perform things after the setup. It may assume
#     miniGUIBase existence.
#  WRAPFUN = when FALSE, makeWidgetCmd is not used to create the miniGUI
#     tcltk wrapper function. 
#    Creates the gui 
{
  ##  tcltk draw main window
  miniGUIBase <- tktoplevel() ## miniGUIBase <<- tktoplevel()
  tkwm.title(miniGUIBase,title)
  ##  Some inits...
  init()
  printGUIAns <- function(...) { miniGUIoutput( getMiniGUIans() ) }
  quit <- function(...) { tkdestroy(miniGUIBase) }
  doNothing <- function(...){ }
  showGuiData <- function(...){
    res <- NULL
    if(length(miniGUIEnvir$miniGUIData)==0)
      cat('\nNo data found.')
    else{
      for(n in names(miniGUIEnvir$miniGUIData))
      res <- rbind(res,cbind(CLASS=class(miniGUIEnvir$miniGUIData[[n]]),
                             TYPE=typeof(miniGUIEnvir$miniGUIData[[n]])))
      rownames(res) <- names(miniGUIEnvir$miniGUIData)
      cat("\nMini-GUI data:\n")
      print(res)
    }
  }
  ##  file Menu function
  fileMenuCmd <- function() 
  {
    fileMenu <- tkmenu(tkmenu(miniGUIBase),tearoff=TRUE)
    tkadd(fileMenu, "command", label="GUI data", command=showGuiData)
    tkadd(fileMenu, "command", label="GUI ans.", command=printGUIAns)
    tkadd(fileMenu, "command", label="Quit", command=quit)
    tkpopup(fileMenu, 	tkwinfo("pointerx", miniGUIBase), 
			tkwinfo("pointery", miniGUIBase))
  }
  ##  ops Menu function
  if(WRAPFUN) ##when true this does not work, guess it has to do with envirs
  {
    ## avoid empty names lists
    if( is.null(names(opFuns)) )
      names(opFuns) <- paste("F",1:length(opFuns),sep="")            
    miniGUIffff <- opFuns
    for( nf in names(opFuns) )
      miniGUIffff[[nf]] <- makeWidgetCmd(nf,opFuns[[nf]],miniGUIBase)
  }else{
    miniGUIffff <- opFuns
  }
  opsMenusCmd <- function() addMenusCmd(miniGUIffff,miniGUIBase) 
  ##  adds menus
  baseMenu <- tkmenu(tkmenu(miniGUIBase))
  tkadd(baseMenu, "command", label="Basics", command=fileMenuCmd)
  if(!is.null(opFuns))
    tkadd(baseMenu, "command", label="Ops", command=opsMenusCmd)
  tkconfigure(miniGUIBase, menu=baseMenu)

  ##    add labels and fields for mainFrameFun on frame miniGUIBase
  if( is.function(mainFrameFun) )
    mapFuncToWidget(mainFrameFun,miniGUIBase,"Eval",STORE="mp")

  ##    return
  return(NA)
}
# miniGUIwidget.R
#   some predefined miniGUI widgets, of course made using tcltk
# NOTAS:
#  - A miniGUIwidget is a function that should return a function
#  of a frame and a variable name . The frame is the parent frame 
#  of the widget, while the variable whose name is given should
#  exists in the .GlobalEnv and it will be the variable that will
#  reflect the widget value.
#  - At the preset moment, only the $widget entry is used.
#  - miniGUIData[["WIDGETLOCK"]] <<- T is used to implement a lock
#  in mapFuncToWidget(), so all widget have to check 
#  miniGUIData[["WIDGETLOCK"]], i. e.:
#       if( miniGUIData[["WIDGETLOCK"]]==TRUE ) return(x)
#  to see if is called from mapFuncToWidget()
# ERRORES:


##
##  avoid building input widget when GUI is not used
##
## miniGUIData[["WIDGETLOCK"]] <- FALSE
setMiniGUIData("WIDGETLOCK",val=FALSE)

##
##  miniGUIwidget reckon
##
is.miniGUIwidget <- function(obj) "miniGUIwidget" %in% class(obj)


##
##  miniGUI reckon
##
## is.miniGUIloaded <- function() return( "package:miniGUI" %in% search() )




##
## basic GUI data entry 
##
miniGUIdefaultEntry <- tkentry
miniGUIentry <- function(x)
##
# x = init value
{
  ## if( ! miniGUIData[["WIDGETLOCK"]] ) return(x)
  if( ! getMiniGUIData("WIDGETLOCK") ) return(x)
  res <- list(widgetType="miniGUIentry",
              widget=function(FRAME,STORE,VAR)  {
                ## miniGUIData[[STORE]][[VAR]] <<- tclVar( x )
                setMiniGUIData(STORE,VAR,val=tclVar( x ))
                ## res <- tkentry(FRAME,textvariable=miniGUIData[[STORE]][[VAR]])
                res <- tkentry(FRAME,textvariable=getMiniGUIData(STORE,VAR))
                return( res )
              },  
              x=x)
  class(res) <- c(class(res),"miniGUIwidget")
  return( res )
}


##
## scale GUI data entry
##
miniGUIscale <- function(from,to,by)
##
# from, to, by = from which value, to which value, by which increment.
{
  ## if( ! miniGUIData[["WIDGETLOCK"]] ) return(from)
  if( ! getMiniGUIData("WIDGETLOCK") ) return(from)
  res <- list(widgetType="miniGUIscale",
              widget=function(FRAME,STORE,VAR)  {
                ## miniGUIData[[STORE]][[VAR]] <<- tclVar()
                setMiniGUIData(STORE,VAR,val=tclVar( ))
                res <- tkscale(FRAME,label="",from=from,to=to,resolution=by,
                           orient="horizontal",
                           ## variable=miniGUIData[[STORE]][[VAR]])
                           variable=getMiniGUIData(STORE,VAR) )
                return( res )
              },
              from=from,to=to,resolution=by)
  class(res) <- c(class(res),"miniGUIwidget")
  return( res )
}


##
## menu selection GUI data entry
##
## will use ttkcombobox(tt,textvariable=a,values=v) when available
miniGUImenusel <- function(xx)
##
# x = vector of mode numeric or character with available values.
#     x[[1]] is taken as the default value. Logicals should be
#     used as c("T","F").
{
  ## if( ! miniGUIData[["WIDGETLOCK"]] ) return(xx[[1]])
  if( ! getMiniGUIData("WIDGETLOCK") ) return(xx[[1]])
  ## to avoid lazy
  xx
  ## normal stuff
  res <- list(widgetType="miniGUImenusel",
              widget=function(FRAME,STORE,VAR)  {
                ttk85 <- as.character(tcl("info", "tclversion")) >= "8.5"
                if(ttk85) {
                  ## miniGUIData[[STORE]][[VAR]] <<- tclVar(xx[[1]])
                  setMiniGUIData(STORE,VAR,val=tclVar( xx[[1]] ))
                  res <- ttkcombobox(parent=FRAME,
                            ## textvariable=miniGUIData[[STORE]][[VAR]],
                            textvariable=getMiniGUIData(STORE,VAR),
                            values=xx)
                }else{
                  x <- "Tcl vers. < 8.5, ttkcombobox not available."
                  ## miniGUIData[[STORE]][[VAR]] <<- tclVar(xx)
                  setMiniGUIData(STORE,VAR,val=tclVar( xx ))
                  res <- tkentry(FRAME,
                            ## textvariable=miniGUIData[[STORE]][[VAR]])
                            textvariable=getMiniGUIData(STORE,VAR))
                }
                return( res )
              },
              values=xx)
  class(res) <- c(class(res),"miniGUIwidget")
  return( res )
}


# myPlugins.R
#   Algunos plugins
# NOTAS:
# ERRORES:


doNothingPlugin <- function(a)
##
{
##   cat("\ndo nothing ",a)
}


evalPlugin <- function(ev)
##
#  As objects are evaluated before giving them to functions in
#  environment(f), that's all we need.
{
  return( ev )
}
environment(evalPlugin) <- .GlobalEnv






## miniGUIAnsAssPlugin <- function(miniGUIAnsTo)
## ##
## {
##   x <- deparse(substitute(miniGUIAnsTo))
##   if(!(x==""))
##   assign(x,miniGUIans,pos=1)
## }

## lessPlugin <- function(what=miniGUIData)
## ##
## {
##   page(what,"print")
## }

## showCallPlugin <- function(f,a)
## ##
## {
##   call(f,a)
## }

