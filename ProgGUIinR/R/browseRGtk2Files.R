## run examples from RGtk2 chapters
## these are in the package with pattern ex-RGtk2.*R

##' @include misc.R
NULL

##################################################

##' Makes a GUI to browse the RGtk2 examples
##'
##' @return produces the GUI
##' @export
browseRGtk2Files <- function() {
  if(!faux_require("RGtk2"))
    stop("This requires the package RGtk2")
  if(!faux_require("svMisc"))
    stop("This requires the package svMisc")

  w <- gtkWindow(show=FALSE)
  w$setTitle("Run RGtk2 Examples")
  g <- gtkVBox(); w$add(g)
  ## toolbar
  tb <- gtkToolbar()
  g$packStart(tb, expand=FALSE)
  
  ##
  pg <- gtkHPaned()
  g$packStart(pg, expand=TRUE, fill=TRUE)
  
  ## left pane for treeview
  leftPane <- gtkFrame(label="Select an example:")
  pg$add1(leftPane)
  
  ## right pane is for text view
  textView <- gtkTextView()
  textView$setLeftMargin(10)
  sw <- gtkScrolledWindow()
  sw$setPolicy("automatic", "automatic")
  sw$add(textView)
  pg$add2(sw)
  
  ## statusbar
  sb <- gtkStatusbar()
  g$packEnd(sb, expand=FALSE, fill=FALSE)
  
  w$setSizeRequest(800, 500)
  w$Show()
  pg$setPosition(250)
  
  
  ## set up tree view
  listFiles <- function() {
    ## grab from this pacakge
    dir <- system.file("Examples","ch-RGtk2", package="ProgGUIinR")
    files <- list.files(path=dir, pattern="^ex-RGtk2.*R$", full.names=TRUE)
    nms <- list.files(path=dir, pattern="^ex-RGtk2.*R$", full.names=FALSE)
    nms <- sub("^ex-RGtk2-", "", nms)
    nms <- sub(".R$", "", nms)
    nms <- gsub("-", " ", nms)
    
    ## include demos from RGtk2 package
    dir <- system.file("demo", package="RGtk2")
    files.rgtk2 <- list.files(path=dir, pattern=".R$", full.names=TRUE)
    nms.rgtk2 <- list.files(path=dir, pattern=".R$", full.names=FALSE)
    nms.rgtk2 <- sub(".R$", "", nms.rgtk2)
    nms.rgtk2 <- paste("RGtk2 demo:", nms.rgtk2)
    return(data.frame(names=c(nms, nms.rgtk2),
                      files=c(files, files.rgtk2), stringsAsFactors=FALSE))
  }
  
  
  store <- rGtkDataFrame(listFiles())
  tv <- gtkTreeView(store)
  tv$setHeadersVisible(FALSE)
  
  ## single column
  vc <- gtkTreeViewColumn()
  tv$insertColumn(vc, 0)
  cr <- gtkCellRendererText()
  vc$packStart(cr)
  vc$addAttribute(cr, "text", 0)
  
  sw <- gtkScrolledWindow()
  sw$setPolicy("automatic", "automatic")
  sw$add(tv)
  leftPane$Add(sw)
  
##################################################
  ## Create callbacks
  
  ## loadExample
  ## we mark comments in a special font
  textView$getBuffer()$createTag(tag.name="comment", weight=PangoStyle["italic"], foreground="red")
  textView$getBuffer()$createTag(tag.name="code")
  loadExample <- function(tree.view, path, column, user.data) {
    row <- as.numeric(path$toString()) + 1
    filename <- store[row, 2]
    txt <- readLines(filename)
    
    buffer <- textView$getBuffer()
    buffer$setText("")                    # clear old
    ## add
    endIter <- buffer$getEndIter()
    for(i in txt) {
      tmp <- regexpr("#", i)
      if(tmp < 0) {
        buffer$insertWithTagsByName(endIter$iter, i, "code")
      } else {
        buffer$insertWithTagsByName(endIter$iter, substr(i,1,tmp-1), "code")
        buffer$insertWithTagsByName(endIter$iter, substr(i,tmp, nchar(i)), "comment")
      }
      buffer$insert(endIter$iter, "\n")
    }
  }
  
  
  ## run example
  runCode <- function(buffer, start, end) {
    ## run code from start iter to end iter
    content <- buffer$getText(start, end)
    cmds <-  unlist(strsplit(content, "\n"))
    eval(parse(text=content), envir=.GlobalEnv)
  }
  runExample <- function(action, data) {
    buffer <- textView$getBuffer()
    bnds <- buffer$getBounds()
    runCode(buffer, bnds$start, bnds$end)
  }
  
  ## run selection
  runSelection <- function(action, data) {
    buffer <-  textView$getBuffer()
    bnds <- buffer$getSelectionBounds()
    runCode(buffer, bnds$start, bnds$end)
  }
##################################################
  ## find word under mouse. Get iter from event position
  findWord <- function(view, e, data) {
    buffer <- view$getBuffer()
    
    iter <- view$getIterAtPosition(e$getX(), e$getY())$iter 
    citer <- iter$copy()
    
    inWordRegExp <- '[a-zA-Z0-9\\_\\-\\.\\$]'
    moveIterToWordBoundary <- function(iter, direction=c("forward","backward")) {
      curText <- function(iter) {         # no iter$iter
        niter <- iter$copy()
        niter$backwardChar()
        buffer$getText(niter, iter)
      }
      flag <- TRUE
      while(flag) {
        if(!grepl(inWordRegExp, curText(iter)))
          return(iter)
        if(direction == "forward") 
          flag <- iter$forwardChar()
        else
          flag <- iter$backwardChar()
      }
      return(iter)
    }
    
    siter <- moveIterToWordBoundary(iter, "backward")
    eiter <- moveIterToWordBoundary(citer, "forward")
    word <- buffer$getText(siter, eiter)
    word <- substr(word, 1, nchar(word)-1)
    return(word)
  }
  
  ## Find "chunk" at insert point
  findChunk <- function(view, e, dt) {
    chunkRegExp <- "^### chunk number"    # marks a chunk start, end
    
    buffer <- view$getBuffer()
    insertPoint <- buffer$getMark("insert")
    iter <- buffer$getIterAtMark(insertPoint)$iter
    ## return iter pointing to current line
    moveIterToChunkLine <- function(iter, direction=c("forward","backward")) {
      curLine <- function(iter) {
        siter <- iter$copy()
        siter$backwardVisibleLine()
        eiter <- iter$copy()
        eiter$forwardToLineEnd()
        txt <- buffer$getText(siter, eiter)
        return(txt)
      }
      flag <- TRUE
      while(flag) {
        if(grepl(chunkRegExp, curLine(iter))) {
          return(iter)
        }
        if(direction=="forward")
          flag <- iter$forwardLine()
        else
          flag <- iter$backwardLine()
      }
      return(iter)
    }
    
    iter <- buffer$getIterAtMark(insertPoint) 
    siter <- moveIterToChunkLine(iter$iter, direction="backward")
    
    iter <- buffer$getIterAtMark(insertPoint)   
    eiter <- moveIterToChunkLine(iter$iter, direction="forward")
    chunk <- buffer$getText(siter, eiter)
    return(chunk)
  }
  
##################################################
  ## Toolbar items
  ## toolbar items: quit | run selection | run buffer
  quitAction <- gtkAction(name="quit",
                          label="Quit",
                          tooltip = "quit GUI",
                          stock.id="gtk-quit")
  gSignalConnect(quitAction, "activate", function(action, data) data$Destroy(), data=w)
  b <- gtkToolButton(stock.id="gtk-quit"); quitAction$connectProxy(b)
  tb$add(b)
  
  tb$add(gtkSeparatorToolItem())
  
  runBlockAction <- gtkAction(name="run block",
                              label="Run block",
                              tooltip = "Execute block containing cursorr",
                              stock.id="gtk-execute")
  gSignalConnect(runBlockAction, "activate", function(...) {
    eval(parse(text=findChunk(textView)), envir=.GlobalEnv)
  })
  b <- gtkToolButton(stock.id="gtk-execute"); runBlockAction$connectProxy(b)
  tb$add(b)
  
  
  
  runSelectedAction <- gtkAction(name="run selected",
                                 label="Run selected",
                                 tooltip="Execute selected code in text buffer",
                                 stock.id="gtk-execute")
  gSignalConnect(runSelectedAction, "activate", function(...) runSelection())
  runSelectedAction$setSensitive(FALSE)
  b <- gtkToolButton(stock.id="gtk-execute"); runSelectedAction$connectProxy(b)
  tb$add(b)
  
  runAction <- gtkAction(name="run",
                         label="Run",
                         tooltip = "Execute code in text buffer",
                         stock.id="gtk-execute")
  gSignalConnect(runAction, "activate", function(...) runExample())
  b <- gtkToolButton(stock.id="gtk-execute"); runAction$connectProxy(b)
  tb$add(b)
  
  
##################################################
  ## load text buffer on double click
  gSignalConnect(tv, "row-activated", loadExample)
  
##################################################
  ## update icon if selection
  gSignalConnect(textView$getBuffer(), "mark-set", function(buffer, iter, mark, data) {
    runSelectedAction$setSensitive(buffer['has-selection'])
  })
  
  
  gSignalConnect(textView, "motion-notify-event", function(view, e, data) {
    word <- findWord(view, e, data)
    if(word== "")
      return(FALSE)
    comps <- svMisc::completion(word)
    if(!is.data.frame(comps))
      return(FALSE)
    ind <- which(word == comps[,1])
    if(length(ind) && comps[ind, "type"] %in% c("function","variable")) {
      FUN <- comps[ind, "completion"]
      sb$push(1, svMisc::argsTip(FUN))
    } else {
      sb$push(1,"")
    }
    return(FALSE)
  })
  
}
