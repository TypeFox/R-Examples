## We could call this with an empty fgui() ...
fguiWindow <- function( basicMenu=TRUE, title="fgui", text="Please choose an option from the menu." ) {
  require( tcltk )

  ## cleanup from previous window...
  clearMenu()

  ## Create the main window
  window <- tktoplevel()
  tkwm.title( window, title )

  ## Add the main menu to the window
  topMenu <- tkmenu( window )
  tkconfigure( window, menu=topMenu )


  ##fileMenu <- tkmenu( topMenu, tearoff=FALSE )
  ##tkadd( topMenu, "cascade", label="File", menu=fileMenu )

  ##tkadd( fileMenu, "command", label="Quit", command=function() tkdestroy(window) )
  ##tkadd( fileMenu, "command", label="Die", command=function() print("hello") )


  ##aMenu <- tkmenu( topMenu, tearoff=FALSE )
  ##tkadd( topMenu, "cascade", label="a", menu=aMenu )

  ##abMenu <- tkmenu( aMenu, tearoff=FALSE )
  ##tkadd( aMenu, "cascade", label="a > b", menu=abMenu )

  ##tkadd( abMenu, "command", label="a > b > c", command=function(){print("hello")} )


  guiSet( "INTERNALMENU_topMenu", topMenu )
  tkfocus( window )


  ## Now set the window as created
  guiSet( "INTERNALMENU_USED", TRUE )

  ## Add a destroy event to set the window as dead
  #tcl("wm", "protocol", window, "WM_DELETE_WINDOW", quote(cat("Im staying!\n")))
  tcl( "wm", "protocol", window, "WM_DELETE_WINDOW",
       function() { guiSet("INTERNALMENU_USED",FALSE);tkdestroy(window);} )


  ## Now add the real menu
  #fguiNewMenu( c("File","me","here") )
  #fguiNewMenu( c("File","me","there") )
  #fguiNewMenu( c("File","clear") )

  ## Maybe add a text edit box
  scr <- tkscrollbar(window, repeatinterval=5,
                        command=function(...)tkyview(txt,...))
  txt <- tktext(window,yscrollcommand=function(...)tkset(scr,...))
  tkgrid(txt,scr)
  tkgrid.configure(scr,sticky="ns")
  guiSet( "INTERNALMENU_TXT", txt )

  if( basicMenu ) {
    fguiNewMenu( c("File","Clear"), command=function(){gui_tkClearTextM(txt)} )
    fguiNewMenu( c("File","Save"),
      command=function(){
        ## Get the text from the box
        text <- tclvalue( tkget( guiGet("INTERNALMENU_TXT"), '@0,0', 'end' ) )
        ## Get the filename
        fname <- tclvalue(tkgetSaveFile(filetypes="{{Text Files} {.txt}} {{All files} *}"))

        if( fname!="" ) {
          ## Write text to fname!
          f <- file( fname )
          writeLines( text, con=f )
          close(f)
        }
      } )
    ##tkadd( guiGet("INTERNALMENU_File"), "separator" )
    fguiNewMenu( c("File","SEPARATOR") )
    fguiNewMenu( c("File","Exit"), command=function(){tkdestroy(window)} )
  }



  ## And write a little something to that text edit box
  fguiWindowPrint( text )

  return( invisible() )
}

fguiWindowPrint <- function( text, endl=TRUE ) {
  txt <- guiGetSafe( "INTERNALMENU_TXT" )
  #SUCCEEDED <- FALSE
  #if( !is.na(txt[1]) ) {
  #  try( {
  #    gui_tkInsertText <- getFromNamespace( "gui_tkInsertText", "fgui" );
  #    gui_tkInsertText(txt, text);
  #    if( endl ) gui_tkInsertText(txt, "\n");
  #    SUCCEEDED <- TRUE;
  #    } )
  #}
  #if( !SUCCEEDED ) {
  #  cat( text );
  #  if( endl ) cat( "\n" );
  #}

  winExists <- guiGetSafe( "INTERNALMENU_USED" )
  if( is.na( winExists ) || !winExists ) {
    cat( text )
    if( endl ) cat( "\n" )
  }else{
    gui_tkInsertText( txt, text );
    if( endl ) gui_tkInsertText( txt, "\n" )
  }
}

## menuText is an array of the depth
##fguiNewMenu <- function( menuText, command=function(){print(menuText[length(menuText)])} ) {
fguiNewMenu <- function( menuText, command=function(){print(paste(menuText,collapse=" > "))} ) {
  ## Make sure the window existed
  winExists <- guiGetSafe( "INTERNALMENU_USED" )
  if( is.na( winExists ) || !winExists ) {
    fguiWindow()
  }

  ## Previous menu should be the top menu
  prevMenu <- guiGetSafe( "INTERNALMENU_topMenu" )
  if( is.na( prevMenu[1] ) ) {
    fguiWindow()
    prevMenu <- guiGetSafe( "INTERNALMENU_topMenu" )
  }

  for( i in 1:length(menuText) ) {
    ##print( menuText[i] )
    ##newMenu <- NA

    curName <- paste( "INTERNALMENU_", menuText[i], sep="" )
    ##print( curName )
    menu <- guiGetSafe( curName )
    ##print( menu[1] )
    if( is.na(menu[1]) ) {
      ## how should it be added?
      if( i==length(menuText) ) {
        ## Then it's the end of the menu -- should be a _command_
        ## UNLESS, it is a SEPARATOR
        if( menuText[i]!="SEPARATOR" ) {
          tkadd( prevMenu, "command", label=menuText[i], command=command )
        }else{
          tkadd( prevMenu, "separator" )
        }
      }else{
        ## Contains items -- should be _cascade_
        ## menu needs to be created!
        menu <- tkmenu( prevMenu, tearoff=FALSE )
        tkadd( prevMenu, "cascade", label=menuText[i], menu=menu )
        guiSet( curName, menu )
        registerMenu( curName ) ## so can be cleaned up...
      }
    }

    ## and reset the previous menu
    if( !is.na(menu[1]) )
      prevMenu <- menu
  }
}

## Internal -- so we can clear the menu later
registerMenu <- function( menuName ) {
  menuNames <- guiGetSafe( "INTERNALMENU_INTERNALNAMES" )
  if( is.na( menuNames[1] ) ) {
    menuNames <- menuName
  }else{
    menuNames <- c(menuNames,menuName)
  }
  guiSet( "INTERNALMENU_INTERNALNAMES", menuNames )
  ##print( "registerMenu menuNames" )
  ##print( menuNames )
}

## Internal -- clears the menu
clearMenu <- function() {
  menuNames <- guiGetSafe( "INTERNALMENU_INTERNALNAMES" )
  ##print( "clearMenu menuNames" )
  ##print( menuNames )

  if( is.na( menuNames[1] ) )
    return()

  for( m in menuNames )
    guiSet( m, NA ) ## kill it

  guiSet( "INTERNALMENU_INTERNALNAMES", NA )
}

## Here everything is the same,
##  _except_ title is now a vector representing the menu path
mgui <- function( func,
                  argOption=NULL, argFilename=NULL, argList=NULL, argSlider=NULL,
                  argCommand=NULL, argEdit=NULL, argFilter=NULL,
                  argText=NULL, argType=NULL,
                  argGridOrder=1:length(formals(func)),
                   argGridSticky=rep("a",length(formals(func))),
                   argGridFrame=rep("f",length(formals(func))),
                  title=NULL,
                  exec="OK",
                  closeOnExec=is.null(output), cancelButton=TRUE,
                  callback=NULL,
                  output='m',
                  helps='auto', helpsFunc=NULL,
                  grid=TRUE, modal=TRUE, nameFix=TRUE, getFix=TRUE,
                  verbose=FALSE ) {
  call <- match.call(expand.dots = FALSE)
  funcName <- call[[match("func", names(call))]]
  funcName <- as.character( as.expression( funcName ) )
  if( length(title)==0 )
    title <- funcName

  fguiNewMenu( menuText=title,
               command=function() {
                 gui( func=func,
                      argOption=argOption, argFilename=argFilename, argList=argList, argSlider=argSlider,
                      argCommand=argCommand, argEdit=argEdit, argFilter=argFilter,
                      argText=argText, argType=argType,
                      argGridOrder=argGridOrder,
                       argGridSticky=argGridSticky,
                       argGridFrame=argGridFrame,
                      title=title[length(title)],
                      exec=exec,
                      closeOnExec=closeOnExec, cancelButton=cancelButton,
                      callback=callback,
                      output=output,
                      helps=helps, helpsFunc=helpsFunc,
                      grid=grid, modal=modal, nameFix=nameFix, getFix=getFix,
                      verbose=verbose )
               } )
}

## debugging routines
#fguiWindowPrint( "Should go to the console." )
#####fguiWindow()
#mgui( rnorm, title=c("Random","Normal") )
#mgui( runif, title=c("Random","Uniform") )
#fguiWindowPrint( "Should go to the main window." )