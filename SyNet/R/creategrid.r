creategrid <- function(dnpoint){  
  require(tcltk) 
  require(tkrplot)
  if (is.null(class(dnpoint)) | class(dnpoint) != "dnpoint") {
        cat("Argument is not of class 'dnpoint' \n")
        return(invisible())
  }
  nsp <- length(dnpoint$Label)
  long <- dnpoint$Points[,2]
  lat <- dnpoint$Points[,3]
  xlims <- range(long)
  xlims <- xlims + c(-0.1*diff(xlims), 0.1*diff(xlims))
  ylims <- range(lat)
  ylims <- ylims + c(-0.1*diff(ylims), 0.1*diff(ylims))
  #tablainput <- data.frame(x = long, y = lat, idlab = rep(dnpoint$Label, tapply(dnpoint$Points[,1],dnpoint$Points[,1],length)), idnum = dnpoint$Points[,1])
  origin <- NULL
  xorvar <- tclVar(round(xlims[1], 3))
  yorvar <- tclVar(round(ylims[2], 3))
  latres <- longres <- 0
  ncols <- tclVar("10")
  nlins <- tclVar("10")
  longgr <- tclVar(round(diff(xlims)/10, 3))
  latgr <- tclVar(round(diff(ylims)/10, 3))
  currentx <- tclVar()
  currenty <- tclVar()
  maxrich <- tclVar()
  gridx <- gridy <- c()
  modpar <- function(ruta) {
  	options(warn = -1)
    if(ruta > 0) {
      n <- as.integer(tclvalue(ncols))
      if(is.na(n)) {tclvalue(ncols) <- 1; tkrreplot(img); return()}
      if(ruta == 1) tclvalue(ncols) <- as.character(aux <- n + 1) 
      if(ruta == 2) tclvalue(ncols) <- as.character(aux <- n - 1) 
      if(ruta == 3) tclvalue(ncols) <- as.character(aux <- n)
      if(aux < 1) tclvalue(ncols) <- "1"
      if(aux > 9999) tclvalue(ncols) <- "9999" 
      tkrreplot(img)
    }
    if(ruta < 0) {
      n <- as.integer(tclvalue(nlins))
      if(is.na(n)) {tclvalue(nlins) <- 1; tkrreplot(img); return()}
      if(ruta == -1) tclvalue(nlins) <- as.character(aux <- n + 1) 
      if(ruta == -2) tclvalue(nlins) <- as.character(aux <- n - 1)  
      if(ruta == -3) tclvalue(nlins) <- as.character(aux <- n) 
      if(aux < 1) tclvalue(nlins) <- "1"
      if(aux > 9999) tclvalue(nlins) <- "9999"  
      tkrreplot(img)
    }
  }
  
  curpos <- function (x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
    if(any(c(x < 0, x > width, y < 0, y > height))) {return()}
    rangeX <- diff(range(gridx))
    rangeY <- diff(range(gridy))
    auxx <- x/width 
    nx <- gridx[1] + auxx*rangeX
    auxy <- y/height
    ny <- gridy[1] - auxy*rangeY
    return(c(nx, ny))
  }
  
  drawgrid <- function() {
    gridx <<- seq(xlims[1], xlims[2], length.out = as.integer(tclvalue(ncols)) + 1)
    gridy <<- seq(ylims[2], ylims[1], length.out = as.integer(tclvalue(nlins)) + 1) 
    tclvalue(latgr) <<- as.character(round(diff(ylims)/as.integer(tclvalue(nlins)), 3))
    tclvalue(longgr) <<- as.character(round(diff(xlims)/as.integer(tclvalue(ncols)), 3))
    par(bg = "cyan")
    plot.new()
    par(plt = c(0,1,0,1))
    par(usr = c(xlims, ylims))
    points(long, lat, pch = 19, xlab = "", ylab = "", col = "red")
    abline(v = gridx, h = gridy)
    if(!is.null(origin)) abline(v = origin[1], h = origin[2], col = 3, lwd = 2)
  }
  
  movegrid <- function(x, y) {
    origin <<- c(x, y)
    tclvalue(xorvar) <<- as.character(round(x, 3))
    tclvalue(yorvar) <<- as.character(round(y, 3))
    tkrreplot(img)
  }
  
  fixnewor <- function() {
    if(is.na(as.numeric(tclvalue(xorvar)))) {
      tkmessageBox(message= "Your NW_corner (x) is ambiguous. Previous value recovered", icon="warning")
      tclvalue(xorvar) <<- as.character(round(xlims[1], 3))
    }
    if(is.na(as.numeric(tclvalue(yorvar)))) {
      tkmessageBox(message= "Your NW_corner (y) is ambiguous. Previous value recovered", icon="warning")
      tclvalue(yorvar) <<- as.character(round(ylims[2], 3))
    }
    auxx <- as.numeric(tclvalue(xorvar))
    auxy <- as.numeric(tclvalue(yorvar))
    if(auxx > max(long)) {
      tkmessageBox(message= "Your NW_corner is beyond the right limit of data. Previous value recovered", icon="warning")
      tclvalue(xorvar) <<- as.character(round(xlims[1], 3))    
    }
    if(auxy < min(lat)) {
      tkmessageBox(message= "Your NW_corner is below the bottom of data. Previous value recovered", icon="warning")
      tclvalue(yorvar) <<- as.character(round(ylims[2], 3))    
    }
    xlims[1] <<- round(as.numeric(tclvalue(xorvar)), 3)
    ylims[2] <<- round(as.numeric(tclvalue(yorvar)), 3)    
    origin <<- NULL 
    tkrreplot(img)
  }
  
  changecell <- function(ruta) {
    options(warn = -1)
    if(ruta == 1) {
      celh <- as.numeric(tclvalue(longgr))
      if(is.na(celh)) {modpar(3); return()}
      if(celh < (aux <- max(long) - xlims[1])) {
        p <- ceiling(aux/celh)
        xlims[2] <<- xlims[1] + p*celh 
        tclvalue(ncols) <<- as.character(p)
        longres <<- celh
      }
      modpar(3)
    }
    if(ruta == 2) {
      celv <- as.numeric(tclvalue(latgr))
      if(is.na(celv)) {modpar(-3); return()}
      if(celv < (aux <- ylims[2] - min(lat))) {
      p <- ceiling(aux/celv)
      ylims[1] <<- ylims[2] - p*celv 
      tclvalue(nlins) <<- as.character(p) 
      latres <<- celv
      }
      modpar(-3)
    }
  }
  
  modalDialog <- function(title, question, entryInit, entryWidth = 20,
      returnValOnCancel = "ID_CANCEL") {
      dlg <- tktoplevel()
      tkwm.deiconify(dlg)
      tkgrab.set(dlg)
      tkfocus(dlg)
      tkwm.title(dlg, title)
      textEntryVarTcl <- tclVar(paste(entryInit))
      textEntryWidget <- tkentry(dlg, width = paste(entryWidth),
          textvariable = textEntryVarTcl)
      tkgrid(tklabel(dlg, text = "       "))
      tkgrid(tklabel(dlg, text = question), textEntryWidget)
      tkgrid(tklabel(dlg, text = "       "))
      ReturnVal <- returnValOnCancel
  
      onOK <- function() {
          ReturnVal <<- tclvalue(textEntryVarTcl)
          tkgrab.release(dlg)
          tkdestroy(dlg)
          tkfocus(tt)
      }
      onCancel <- function() {
          ReturnVal <<- returnValOnCancel
          tkgrab.release(dlg)
          tkdestroy(dlg)
          tkfocus(tt)
      }
      OK.but <- tkbutton(dlg, text = "   OK   ", command = onOK)
      Cancel.but <- tkbutton(dlg, text = " Cancel ", command = onCancel)
      tkgrid(OK.but, Cancel.but)
      tkgrid(tklabel(dlg, text = "    "))
  
      tkfocus(dlg)
      tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg); tkfocus(tt)})
      tkbind(textEntryWidget, "<Return>", onOK)
      tkwait.window(dlg)
  
      return(ReturnVal)
  }
  
  outtable <- function() {
      ngrx <- as.integer(tclvalue(ncols))
      ngry <- as.integer(tclvalue(nlins))
      if(diff(xlims) > 0) gridx <- seq(xlims[1], xlims[2], length.out = ngrx + 1) else gridx <- xlims + c(-1, 1)
      if(diff(ylims) > 0) gridy <- seq(ylims[2], ylims[1], length.out = ngry + 1) else gridy <- ylims + c(1, -1)
      difx <- diff(gridx)
      dify <- diff(gridy)
      cornersi <- expand.grid(1:ngrx, 1:ngry)
      dntable <- matrix(0, nrow = nsp, ncol = nrow(cornersi))
      rownames(dntable) <- dnpoint$Label
      for(pts in 1:length(long)) {
        auxx <- long[pts] - gridx
        auxy <- lat[pts] - gridy
        celx <- which(diff((auxx > 0) + 2*(auxx < 0))!=0)
        cely <- which(diff((auxy > 0) + 2*(auxy < 0))!=0)
        celocc <- celx[1] + (cely - 1)*ngrx
        if(length(celx) == 2) celocc <- c(celocc, celx[2] + (cely - 1)*ngrx)  
        dntable[dnpoint$Points[pts, 1],celocc] <- 1
      }
      sel <- which(apply(dntable, 2, sum)>0)
      cornersi[sel,1]-> xleft
      cornersi[sel,2] + 1 -> ybottom
      cornersi[sel, 1] + 1 -> xright
      cornersi[sel,2] -> ytop
      outobj <- list(dntable = dntable[, sel, drop = FALSE], xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, ncells = c(ngrx, ngry))
      return(outobj)
  }
  
  launchout <- function() {
      ReturnVal <- modalDialog("Save data into an object", "Enter the variable name", "outgr")
      if (ReturnVal == "ID_CANCEL") return()
      out <- outtable()
      outobj <- data.frame(rbind(colid = out$xleft, rowid = out$ytop, out$dntable))
      assign(ReturnVal, outobj, envir = .GlobalEnv)
      #create distributional table
      tkmessageBox(title = "Out", message = paste("The object ", ReturnVal, "holds your grid data"))
  }
  
  drawrich <- function(occmap){
      gridx <- seq(xlims[1], xlims[2], length.out = as.integer(tclvalue(ncols)) + 1)
      gridy <- seq(ylims[2], ylims[1], length.out = as.integer(tclvalue(nlins)) + 1) 
      par(bg = "white")
      plot.new()
      par(plt = c(0,1,0,1))
      par(usr = c(xlims, ylims))
      totrich <- apply(occmap$dntable, 2, sum)
      rel <- totrich/max(totrich)
      tclvalue(maxrich) <<- max(totrich)
      rect(gridx[occmap$xleft], gridy[occmap$ybottom], gridx[occmap$xright], gridy[occmap$ytop], col = gray(1- rel))
      abline(v = gridx, h = gridy, col = 2)
  }
  
  constructpdf <- function(){
    filename <- tclvalue(tkgetSaveFile(initialfile = "Rplots.pdf",
                defaultextension = ".pdf", title = "Save graph...",
                filetypes = "{PDF {.pdf}} {{All Files} {*.*}}"))
    if (filename != "") pdf(file = filename) else { tkrreplot(richplot); return()}
    outdn <- outtable()
    plot(rep(1, 4), 1:4, type = "n", xlab = "", ylab = "", axes = FALSE, main = "MAP OF OCCUPIED CELLS")
    text(1, 1, paste("Grid of ", tclvalue(ncols), "cols and ", tclvalue(nlins), "rows"))
    text(1, 2, paste("Northwest corner at: (", tclvalue(xorvar), ", ", tclvalue(yorvar), ")"))
    text(1, 3, paste("Cell Width = ", tclvalue(longgr), "units and Cell Height = ", tclvalue(latgr), "units"))
    text(1, 4, paste("Gray tones proportional to richness by cell\nHighest richness (black cell[s]) = ", tclvalue(maxrich))) 
    drawrich(outdn)
    dev.off()
  }
  
  tt<-tktoplevel()
  tkwm.title(tt, "Grid creator")
  leftframe <- tkframe(tt, relief = "groove", borderwidth = 2)
  rightframe <- tkframe(tt, relief = "groove", borderwidth = 2)
  gridframe <- tkframe(leftframe, relief = "groove", borderwidth = 2)
  outputframe <- tkframe(leftframe, relief = "groove", borderwidth = 2)
  plotframe <- tkframe(leftframe)
  cols <- tkentry(gridframe, width="6", textvariable = ncols)
  lins <- tkentry(gridframe, width="6", textvariable = nlins)
  xorigen <- tkentry(gridframe, width = "6", textvariable = xorvar)
  yorigen <- tkentry(gridframe, width = "6", textvariable = yorvar)
  xsize <- tkentry(gridframe, width="6", textvariable = longgr) 
  ysize <- tkentry(gridframe, width="6", textvariable = latgr) 
  addcol <- tkbutton(gridframe, text = "+", command = function(...) modpar(1))
  subtractcol <- tkbutton(gridframe, text = "-", command = function(...)  modpar(2))
  addlin <- tkbutton(gridframe, text = "+", command = function(...) modpar(-1))
  subtractlin <- tkbutton(gridframe, text = "-", command = function(...) modpar(-2))
  img <- tkrplot(plotframe, fun = drawgrid, hscale=1.5,vscale=1.5)
  richplot <- tkrplot(rightframe, function() drawrich(outtable()), hscale=1.2,vscale=1.2)
  outButton <- tkbutton(outputframe, text = "Create object", command = launchout)
  refreshButton <- tkbutton(outputframe, text = "Full extent", 
                            command = function(){xlims <<- range(long)
                                                 xlims <<- xlims + c(-0.1*diff(xlims), 0.1*diff(xlims))
                                                 ylims <<- range(lat)
                                                 ylims <<- ylims + c(-0.1*diff(ylims), 0.1*diff(ylims))
                                                 tclvalue(xorvar) <<- as.character(round(xlims[1],3))
                                                 tclvalue(yorvar) <<- as.character(round(ylims[2], 3))
                                                 tkrreplot(img)
                                                 })
  richButton <- tkbutton(rightframe, text = "Refresh richness map", command = function() tkrreplot(richplot))
  pdfButton <- tkbutton(rightframe, text =  "  Save pdf report  ", command = function()  {tkrreplot(richplot); constructpdf()})
  clipButton <- tkbutton(outputframe, text = "Copy to clipboard", command = function() tkrreplot(img))                                               
  xpos <- tklabel(outputframe, textvariable = currentx)
  ypos <- tklabel(outputframe, textvariable = currenty)
  mrich <- tklabel(rightframe, textvariable = maxrich, font = 13)
  tkgrid(tklabel(gridframe, text = "GRID PARAMETERS", font = "Times 14", foreground = "blue"), columnspan = 4)
  tkgrid(tklabel(gridframe, text = " # rows ="), lins, subtractlin, addlin)
  tkgrid(tklabel(gridframe, text = " # cols ="),  cols, subtractcol, addcol)
  tkgrid.configure(subtractlin,sticky="e", ipadx = 4)
  tkgrid.configure(addlin,sticky="w", ipadx = 2)
  tkgrid.configure(subtractcol,sticky="e", ipadx = 4)
  tkgrid.configure(addcol,sticky="w", ipadx = 2)
  tkgrid(tklabel(gridframe, text = "Cell Size"), columnspan = 4)
  tkgrid(tklabel(gridframe, text = "width ="), xsize, tklabel(gridframe, text = "height ="), ysize)
  tkgrid(tklabel(gridframe, text = "NW Corner"), columnspan = 4)
  tkgrid(tklabel(gridframe, text = "    x ="), xorigen, tklabel(gridframe, text = "     y ="), yorigen)
  tkgrid(tklabel(outputframe, text = "GENERAL ACTIONS", font = "Times 14", foreground = "blue"), columnspan = 2)
  tkgrid(outButton, columnspan = 2)
  tkgrid(refreshButton, columnspan = 2)
  tkgrid(clipButton, columnspan = 2)
  tkgrid(tklabel(outputframe, text = "Current Position (x, y)", font = "Times 12"), columnspan = 2)
  tkgrid(xpos, ypos) 
  tkgrid(gridframe,outputframe)
  tkgrid(img)
  tkgrid(plotframe, columnspan = 2)                            
  tkgrid(tklabel(rightframe, text = "OUTPUT DETAILS", font = "Times 14", foreground = "blue"))
  tkgrid(richButton)
  tkgrid(pdfButton)
  tkgrid(tklabel(rightframe, text = ""))
  tkgrid(richplot)
  tkgrid(tklabel(rightframe, text = "Maximum recorded richness (black cells)", font = 13))
  tkgrid(mrich)
  tkgrid(leftframe, rightframe)
  tkbind(img, "<B1-Motion>", function(x, y) {
                             aa <- curpos(x, y)
                             if(!is.null(aa)) movegrid(aa[1], aa[2])})
  tkbind(img, "<ButtonRelease-1>", fixnewor)
  tkbind(img, "<Motion>", function(x, y) {
                          aa <- curpos(x, y)
                          tclvalue(currentx) <<- round(aa[1], 3)
                          tclvalue(currenty) <<- round(aa[2], 3)}) 
  tkbind(img, "<Leave>", function() {tclvalue(currentx) <<- ""; tclvalue(currenty) <<- ""})
  tkbind(cols, "<KeyPress-Return>", function() {modpar(3); tkfocus(lins)})
  tkbind(lins, "<KeyPress-Return>", function() {modpar(-3); tkfocus(cols)})
  tkbind(cols, "<FocusOut>", function() modpar(3))
  tkbind(lins, "<FocusOut>", function() modpar(-3))
  tkbind(xsize, "<KeyPress-Return>", function() {changecell(1); tkfocus(ysize)})
  tkbind(ysize, "<KeyPress-Return>", function() {changecell(2); tkfocus(xsize)})
  tkbind(xsize, "<FocusOut>", function() changecell(1))
  tkbind(ysize, "<FocusOut>", function() changecell(2))
  tkbind(xorigen, "<KeyPress-Return>", function() {fixnewor(); tkfocus(yorigen)})
  tkbind(yorigen, "<KeyPress-Return>", function() {fixnewor(); tkfocus(xorigen)}) 
  tkbind(xorigen, "<FocusOut>", fixnewor)
  tkbind(yorigen, "<FocusOut>", fixnewor) 
}
