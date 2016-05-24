##' @include misc.R
NULL


##' GUI to display tcltk examples
##'
##' @return NULL
##' @export
browseTclTkFiles <- function() {

  if(!faux_require("tcltk"))
    stop("This function needs the tlctk package")


  ## Helper functions
  addScrollbars <- function(parent, widget,type=c("both", "x", "y")) {
    if(any(type %in% c("both","x"))) {
      xscr <- ttkscrollbar(parent, orient="horizontal",
                           command=function(...) tkxview(widget, ...))
      tkconfigure(widget,
                  xscrollcommand=function(...) tkset(xscr,...))
    }
    
    if(any(type %in% c("both","y"))) {
      yscr <- ttkscrollbar(parent, orient="vertical",
                           command=function(...) tkxview(widget, ...))
      tkconfigure(widget,
                  yscrollcommand=function(...) tkset(yscr,...))
    }
    
    ## place in grid     
    tkgrid(widget, row=0, column=0, sticky="news")
    if(any(type %in% c("both", "x"))) {
      tkgrid(xscr, row=1, column=0, sticky="ew")
      tkgrid.columnconfigure(parent, 0, weight=1)
    }
    if(any(type %in% c("both", "y"))) {
      tkgrid(yscr,row=0,column=1, sticky="ns")
      tkgrid.rowconfigure(parent, 0, weight=1)
    }
    
    
  }
  
  
  dir <- system.file("Examples","ch-tcltk", package="ProgGUIinR")
  
  Rfiles <- list.files(path=dir, pattern="^ex(.*)R$")
  
  ## basic two panel approach
  ## with top toolbar and button bar on bottom
  
  w <- tktoplevel(width=800, height=500)
  tkwm.title(w, "View example code for tcltk")
  f <- ttkframe(w, padding=c(3,3,3,12))
  tkpack(f, expand=TRUE, fill="both")
  
  ## menu if on mac to override default
  mb <- tkmenu(w)
  tkconfigure(w, menu=mb)
  
  
  
  ## paned window
  pw <- ttkpanedwindow(f, orient="horizontal")
  tkpack(pw, expand=TRUE, fill="both")
  
  ## listbox via treeview
  tvf <- ttkframe(pw)
  tcl(pw, "add", tvf, weight=1)         # no tkadd function
  tv <- ttktreeview(tvf, columns=1,show="headings", height=25)
  addScrollbars(tvf,tv, type="y")
  tcl(tv, "heading", 1, text="Example file", anchor="center")
  tcl(tv, "column", 1, width=200, stretch=TRUE, anchor="w")

  ## main text window
  tbf <- ttkframe(pw)
  tcl(pw, "add", tbf, weight=3)
  tbff <- ttkframe(tbf)
  tkpack(tbff, expand=TRUE, fill="both")
  tb <- tktext(tbff)
  addScrollbars(tbff, tb)
  
  toolBar <- ttkframe(tbf)
  tkpack(toolBar, fill="x")
  
  
  ## sash width?
  tcl(pw, "sashpos", 0, 150)              # first sash, 150 pixels
  
  
  ## summary bar
  bbf <- ttkframe(f, padding=1)
  tkpack(bbf, fill="x")
  summaryLabel <- ttklabel(bbf, text="")
  tkpack(summaryLabel, fill="x")
  
  
  ## now get to work
  
  ## first populate the Rfiles
  shade <- c("gray","none")
  for(i in seq_along(Rfiles) )
    tcl(tv, "insert", "", "end", tag=shade[1 + i %% 2],
        values=Rfiles[i])
  tktag.configure(tv, "gray", background="gray95")
  
  ## Now when double click, add to text buffer
  showFileInTextBuffer <- function(f, tb) {
    lines <- try(readLines(f), silent=TRUE)
    if(inherits(lines, "try-error")) {
      cat("Can't read", f, "\n")
      return()
    }
    ## clean out old
    tkdelete(tb, "0.0", "end")
    
    for(i in seq_along(lines)) {
      tkinsert(tb, "end", lines[i])
      tkinsert(tb, "end", "\n")
    }
    tksee(tb, "0.0")
  }
  
  tkbind(tv, "<Button-1>", function(W, x, y) {
    row <- as.character(tcl(W, "identify", "row", x, y))
    row <- 1 + as.integer(as.character(tcl(W,"index", row))) 
    f <- paste(dir,Rfiles[row], sep="/") # or get from treeview
    showFileInTextBuffer(f, tb)
  })
  
  
  ## actions
  evalBuffer <- function(W) {
    lines <- tclvalue(tkget(tb, "0.0", "end"))
    evalText(lines)
  }
  evalRegion <- function(W) {
    ## region is between lines with "^###" regexp
    cmds <- try(tclvalue(tkget(tb, "region.first", "region.last")), silent=TRUE)
    if(!inherits(cmds, "try-error"))
      evalText(cmds)
  }
  
  evalSelection <- function(W) {
    lines <- tclvalue(tkget(tb, "sel.first", "sel.last"))
    evalText(lines)
  }
  
  evalLine <- function(W) {
    line <- tclvalue(tkget(tb, "insert linestart", "insert lineend"))
    evalText(line)
  }
  
  evalText <- function(text) {
    text <- paste(text, collapse="\n")
    eval(parse(text=text), envir=.GlobalEnv)
  }
  
  
  ## add buttons to toolbar
  evalBufferButton <- ttkbutton(toolBar, text="Eval. buffer", command=evalBuffer)
  evalRegionButton <- ttkbutton(toolBar, text="Eval. chunk", command=evalRegion)
  evalSelectionButton <- ttkbutton(toolBar, text="Eval. selection", command=evalSelection, state="disabled")
  evalLineButton <- ttkbutton(toolBar, text="Eval. line", command=evalLine)
  
  tkpack(evalBufferButton, side="left")
  tkpack(evalRegionButton, side="left")
  tkpack(evalSelectionButton, side="left")
  tkpack(evalLineButton, side="left")
  
  
  ## Code to handle selection
  tkbind(tb,"<<Selection>>", function(W) {
    ## selection changed: ## update evalSelection button
    ranges <- tclvalue(tcl(W, "tag", "ranges", "sel"))
    if(length(ranges) > 1 || ranges != "") {
      cur <- tclvalue(tkget(W, "sel.first", "sel.last"))
      if(length(cur) > 0 ) {
        tkconfigure(evalSelectionButton, state="normal")
      }
    } else {
      tkconfigure(evalSelectionButton, state="disabled")
    }
  })
  
  
##################################################
  ##
  ## code to tag chunk for eval chunk
  
  tagRegion <- function(W) {
    ## region defined by "^###"
    ## make tag "region" for current region
    ## remove current "region" tag
                                        #  tktag.configure(W, "region", background="white")
    tktag.configure(W, "region", foreground="black")
    tktag.remove(tb,"region","0.0","end")
    
    ## function to get a lline from index
    getLine <- function(index) 
      tclvalue(tkget(W, paste(index,"linestart"), paste(index, "lineend")))
    inComment <- function(line) length(grep("^###", line)) > 0
    getPosFromIndex <- function(index) {
      ind <- try(tclvalue(tcl(W, "index", index)), silent=TRUE)
      if(inherits(ind, "try-error")) {
        pos <- c(-1, -1)
      } else {
        pos <- as.numeric(unlist(strsplit(ind,"\\.")))
      }
      return(pos)
    }
    atTop <- function(ctr) {
      index <- paste("insert", "-", ctr, "lines")
      pos <- getPosFromIndex(index)
      pos[1] <= 1
    }
    atBottom <- function(ctr) {
      index <- paste("insert", "+" , ctr, "lines")
      pos <- getPosFromIndex(index)
      endPos <- getPosFromIndex("end")
      pos[1] >= endPos[1]
    }
    ## get top
    ctr <- 0
    if(atTop(ctr)) {
      topCtr <- 0
    } else {
      curLine <- getLine("insert"); 
      while(!inComment(curLine) && !atTop(ctr)) {
        ctr <- ctr + 1
        curLine <- getLine(paste("insert", "-", ctr, "lines"))
      }
      if(atTop(ctr)) {
        topCtr <- ctr
      } else {
        ## search for comment
        prevLine <- getLine(paste("insert", "-", ctr,"lines"))
        while(inComment(prevLine) && !atTop(ctr)) {
          ctr <- ctr + 1
          prevLine <- getLine(paste("insert", "-", ctr,"lines"))
        }
        if(atTop(ctr))
          topCtr <- ctr
        else
          topCtr <- ctr - 1
      }
    }
    
    ## get Bottom
    ctr <- 0
    curLine <- getLine("insert")
    ## might be in a comment, clear this first
    if(inComment(curLine)) {
      while(inComment(curLine) && !atBottom(ctr)) {
        ctr <- ctr + 1
        curLine <-  getLine(paste("insert", "+" , ctr, "lines"))
      }
    }
    if(atBottom(ctr)) {
      endCtr <- ctr
    } else {
      ctr <- ctr + 1
      nextLine <-  getLine(paste("insert", "+" , ctr, "lines"))
      while(!inComment(nextLine) && !atBottom(ctr)) {
        ctr <- ctr + 1
        nextLine <-  getLine(paste("insert", "+" , ctr, "lines"))
      }
      if(atBottom(ctr))
        endCtr <- ctr
      else
        endCtr <- ctr - 1
    }
    
    tktag.add(tb, "region", paste("insert", "-", topCtr, "lines linestart", sep=" "),
              paste("insert", "+", endCtr, "lines lineend", sep=" "))
    ## highlight
    tktag.configure(tb, "region", foreground="gray40") # inidicate through shading chunk
  }
  
  tkbind(tb, "<ButtonPress-1>", function(W,K,k) {
    tagRegion(W)
  })
  tkbind(tb, "<KeyRelease>", function(W,K,k) {
    tagRegion(W)
  })
  
##################################################
  ## Motion binding to give quick defn. of object
  ## summary of object mouse is over
  shortSummary <- function(x) UseMethod("shortSummary")
  shortSummary.default <- function(x, len=80) {
    if(inherits(x, "try-error")) {
      return("")
    } else {
      out <- paste("An object of class ", paste(class(x), collaspe=","),",", sep="")
      return(out)
    }
  }
  
  shortSummary.function <- function(x) {
    out <- format(x)[1]
    return(out)
  }
  
  tkbind(tb, "<Motion>", function(W, x, y) {
    ## update status line with summary
    cur <- tkget(W,"current  wordstart", "current wordend") # misses . part of functions
    cur <- tclvalue(cur)
    cur <- try(get(cur, envir=.GlobalEnv), silent=TRUE)
    out <- shortSummary(cur)
    tkconfigure(summaryLabel, text=out)
    
    ## highlight region
    
  })
  
}

## run it
##showTclTkExamples()

