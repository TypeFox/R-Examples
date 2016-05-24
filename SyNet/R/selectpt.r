  selectpt <- function (entryPartition) {
      tf <- tktoplevel()
      tkwm.title(tf, "Choose Partition")
      done <- tclVar(0)
      tlb <- tklistbox(tf)
      scr <- tkscrollbar(tf, repeatinterval = 5, command = function(...) tkyview(tlb,
          ...))
      tkconfigure(tlb, yscrollcommand = function(...) tkset(scr,
          ...))
      frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
      cancel.but <- tkbutton(frame1, text = "Cancel", command = function() tkdestroy(tf))
      submit.but <- tkbutton(frame1, text = "Select", default = "active",
          command = function() tclvalue(done) <- 1)
      tkpack(cancel.but, submit.but, side = "left")
      tkpack(frame1, side = "bottom")
      tkpack(tlb, side = "left", fill = "both", expand = TRUE)
      tkpack(scr, side = "right", fill = "y")
      obj <- ls(globalenv())
      for (i in obj) {
          xobj <- get(i, envir = globalenv())
          if (class(xobj)== "nampartition") tkinsert(tlb, "end", i)
      }
      tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done) <- 1)
      tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
      tkbind(tf, "<KeyPress-Return>", function() tclvalue(done) <- 1)
      tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
      tkwait.variable(done)
      if (tclvalue(done) == "2")
          return(0)
      numc <- tclvalue(tkcurselection(tlb))
      if (numc == "") {
          tkdestroy(tf)
          return(0)
      }
      choix <- tclvalue(tkget(tlb, numc))
      tkdelete(entryPartition, 0, "end")
      tkinsert(entryPartition, "end", choix)
      tkdestroy(tf)
  }
