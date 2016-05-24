viewsp <- function () {
      tf <- tktoplevel()
      tkwm.title(tf, "Group Selected")
      tkwm.geometry(tf, "750x600")
      tkwm.resizable(tf, 0, 0)
      ptn <- NULL
      listofsp <- NULL
      targetgr <- NULL
      part <- tclVar()
      msg <- tclVar()
      done <- tclVar(0)
      targetsp <- NULL
      frameList <- tkframe(tf, relief = "groove", borderwidth = 2)
      tlbgr <- tklistbox(frameList)
      scrgr <- tkscrollbar(frameList, repeatinterval = 5, command = function(...) tkyview(tlbgr,
          ...))
      tkconfigure(tlbgr, yscrollcommand = function(...) tkset(scrgr, ...))
      tlbsp <- tklistbox(frameList)
      scrsp <- tkscrollbar(frameList, repeatinterval = 5, command = function(...) tkyview(tlbsp,
          ...))
      tkconfigure(tlbsp, yscrollcommand = function(...) tkset(scrsp,
          ...))
      loadelements <- function () {
        tkdelete(tlbsp, 0, "end")
        listofsp <<- NULL
        numc <- tclvalue(tkcurselection(tlbgr))
        if (numc == "") return()
        targetgr <<- tclvalue(tkget(tlbgr, numc))
        apply(ptn[[2]], 1, function(x) {
                  if (x[2] == targetgr) {tkinsert(tlbsp, "end", x[1])
                  listofsp <<- c(listofsp, x[1])
                  }})
        tkrreplot(mapview) #call to drawmap()
      } #End loadelements function
      loadpartition <- function() {
        selectpt(entryPartition)
        aux <- as.character(tclvalue(part))
        if (aux %in% ls(globalenv())) ptn <<- get(aux) else return()
        if (!is(ptn, "nampartition")) return()
        tkdelete(tlbgr, 0, "end")
        tkdelete(tlbsp, 0, "end")
        listofsp <<- NULL
        targetgr <<- NULL
        for (i in unique(ptn[[2]][,2])) tkinsert(tlbgr, "end", i)
        tkrreplot(mapview) #call to drawmap()
      } #End loadpartition function
      drawmap <- function () {
        par(bg = "white")
        if(is.null(ptn)) {plot.new(); return()}
        if(ptn$kind == "none") {plot.new(); return()}
        else if(ptn$kind == "grids") {
            x <- as.numeric(ptn$coords[,1])
            y <- as.numeric(ptn$coords[,2])
            v <- seq(0.5, max(x) + 0.5, by = 1)
            h <- seq(-0.5, -1*max(y) - 0.5, by = -1)
            plot(x, - 1*y, xlab = "", ylab = "", main = "SPATIAL EXPRESSION")
            abline(v = v)
            abline(h = h)
            rect(x - 0.5, -1*y - 0.5, x + 0.5, -1*y + 0.5, col = gray(0.9))
            if (is.null(listofsp)) return()
            tclvalue(msg) <- paste(targetgr," (None Element Selected) ")
            whforgr <- unique(unlist(ptn$occupancy[listofsp]))
            rect(x[whforgr] - 0.5, -1*y[whforgr] - 0.5, x[whforgr] + 0.5, -1*y[whforgr] + 0.5, col = 1)
            numc <- tclvalue(tkcurselection(tlbsp))
            if (numc == "") return()
            targetsp <- tclvalue(tkget(tlbsp, numc))
            whforsp <- unlist(ptn$occupancy[targetsp])
            rect(x[whforsp] - 0.5, -1*y[whforsp] - 0.5, x[whforsp] + 0.5, -1*y[whforsp] + 0.5, col = 2)
            tclvalue(msg) <- paste(targetgr," (", targetsp, " )")
        }
        else if(ptn$kind == "points") {
            plot (ptn$coords, xlab = "LONGITUDE", ylab = "LATITUDE", main = "SPATIAL EXPRESSION",
                  col = gray(0.9), pch = 19, asp = 1)
            grsel <- unique(unlist(ptn$occupancy[listofsp]))
            points(ptn$coords[grsel, ,drop = FALSE], pch = 19)
            tclvalue(msg) <- paste(targetgr," (None Element Selected)")
            numc <- tclvalue(tkcurselection(tlbsp))
            if (numc == "") return()
            targetsp <- tclvalue(tkget(tlbsp, numc))
            spsel <- unlist(ptn$occupancy[targetsp])
            points(ptn$coords[spsel, ,drop = FALSE], col = 2, pch = 19)
            tclvalue(msg) <- paste(targetgr," (", targetsp, " )")
        }
      } #End drawmap function
      partButton <- tkbutton(frameList, text = "Search", command = loadpartition)
      loadButton <- tkbutton(frameList, text = "Load", command = loadelements, padx = 10)
      butFrame <- tkframe(frameList)
      cancel.but <- tkbutton(butFrame, text = "Cancel", command = function() tkdestroy(tf), padx = 10)
      submit.but <- tkbutton(butFrame, text = "View", default = "active",
                    command = function() tkrreplot(mapview), padx = 15)
      ref <- tklabel(tf, textvariable = msg, font = "Times 14", foreground = "blue")
      entryPartition <- tkentry(frameList, textvariable = part)
      tkgrid(entryPartition, partButton)
      tkgrid.configure(partButton, sticky = "w")
      tkgrid(tklabel(frameList, text =""))
      tkgrid(tklabel(frameList, text = "PARTITION", foreground = "blue",
             font = "Times 12"), columnspan = 2)
      tkgrid(tlbgr, scrgr)
      tkgrid.configure(scrgr, sticky = "wns")
      tkgrid(tklabel(frameList, text = "||"), columnspan = 2)
      tkgrid(loadButton, columnspan = 2)
      tkgrid(tklabel(frameList, text = "\\/"), columnspan = 2)
      tkgrid(tlbsp, scrsp)
      tkgrid.configure(scrsp, sticky = "wns")
      tkgrid(tklabel(frameList, text = "ELEMENTS", foreground = "blue",
             font = "Times 12"), columnspan = 2)
      tkgrid(tklabel(frameList, text =""))
      tkgrid(submit.but, cancel.but)
      tkgrid(butFrame, columnspan = 2)
      tkgrid(tklabel(frameList, text = "NOTE: Double-click can\nbe used over the lists"),
             columnspan = 2)
      mapview <- tkrplot(tf, fun = drawmap, 1.4, 1.4)
      tkgrid(frameList, mapview)
      tkgrid(tklabel(tf, text = ""), ref)
      tkbind(tlbsp, "<Double-ButtonPress-1>", function() tkrreplot(mapview))
      tkbind(tlbgr, "<Double-ButtonPress-1>", loadelements)
      tkbind(tf, "<Destroy>", function() tkdestroy(tf))
      tkbind(tf, "<KeyPress-Return>", function() tkrreplot(mapview))
      tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
} 

