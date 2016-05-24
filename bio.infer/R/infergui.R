"infergui" <-

function() {

  inferguienv <- environment(NULL)


  getfile <- function() {
    name <- tclvalue(tkgetOpenFile())
    if (name == "") return;

    inferguienv[["bcnt"]] <- read.delim(name)

    if (ncol(inferguienv[["bcnt"]]) != 3) {
      tkmessageBox(message="Please make sure your benthic count file is tab-delimited and has only the following three fields: site ID, taxon name, taxon abundance",icon="error",type="ok")
      inferguienv[["bcnt"]] <- NULL
    }

  }

  getfile.r <- function() {
    getfile()
    if (! is.null(inferguienv[["bcnt"]])) {
      #tkinsert(inferguienv[["txt"]], "end", "Biological data loaded.\n")
      cat("Biological data loaded.\n")
      statevec <- c("normal", "disabled" ,"disabled")
      tkdestroy(inferguienv[["tt"]])
      draw.toplev(statevec)
    }
  }

  select.te <- function() {
    inferguienv[["tefname"]] <- sel.coeffile()

    if (length(inferguienv[["tefname"]]) > 0) {
#      tkinsert(inferguienv[["txt"]], "end", "Taxon-environment coefficient file selected.\n")
      data(list = inferguienv[["tefname"]], envir = inferguienv)
      cat("Taxon-environment coefficient file selected.\n")

#      tkgrab.release(inferguienv[["txt"]])
#      tkfocus(inferguienv[["msgwin"]])
      tkdestroy(inferguienv[["tt"]])
      draw.toplev(c("normal", "normal", "disabled"))
    }
    else {
      tkmessageBox(message = "Select coefficient file to continue",
                   icon = "error", type = "ok")
    }
  }

  compute.inf <- function() {
    tkconfigure(inferguienv[["tt"]], cursor = "wait")
    coef <- get(inferguienv[["tefname"]])
#    data(itis.ttable)
#    tkinsert(inferguienv[["txt"]], "end", "Merging biological data with ITIS...\n")
#    tkfocus(inferguienv[["msgwin"]])
    tkfocus(inferguienv[["tt"]])
    bcnt.tax <- get.taxonomic(inferguienv[["bcnt"]])
#    tkinsert(inferguienv[["txt"]], "end", "Assigning operational taxonomic units...\n")
#    tkfocus(inferguienv[["msgwin"]])
    flush.console()
    bcnt.otu <- get.otu(bcnt.tax, coef, gui = TRUE)
#    tkinsert(inferguienv[["txt"]], "end", "Computing inferences...\n")
#    tkfocus(inferguienv[["msgwin"]])
    flush.console()
    ss <- makess(bcnt.otu)
    inferguienv[["inf.out"]] <- mlsolve(ss, coef)
#    tkinsert(inferguienv[["txt"]], "end", "Inferences computed!")
    cat("Inferences computed!\n")
#    tkfocus(inferguienv[["msgwin"]])
    tkconfigure(inferguienv[["tt"]], cursor = "arrow")
    flush.console()
    if (! is.null(inferguienv[["inf.out"]])) {
      tkdestroy(inferguienv[["tt"]])
      draw.toplev(c("normal", "normal", "normal"))
    }
  }

  export.res <- function() {
    name <- tclvalue(tkgetSaveFile())
    if (name != "") {
      write.table(inferguienv[["inf.out"]], file = name, sep = "\t",
                  row.names = FALSE)
    }
  }

  quitinf <- function() {
    tkgrab.release(inferguienv[["tt"]])
    tkdestroy(inferguienv[["tt"]])
#    tkdestroy(inferguienv[["msgwin"]])
  }

  draw.toplev <- function(statevec) {
    inferguienv[["tt"]] <- tktoplevel()

    tkwm.title(inferguienv[["tt"]], "PECBO")
    button.widget1 <- tkbutton(inferguienv[["tt"]], text = "Load biological data",
                               command = getfile.r)

    button.widget2 <- tkbutton(inferguienv[["tt"]],
                               text = "Select taxon-environment relationships",
                               command = select.te, state = statevec[1])

    button.widget3 <- tkbutton(inferguienv[["tt"]], text = "Compute inferences",
                               command = compute.inf, state = statevec[2])

    button.widget4 <- tkbutton(inferguienv[["tt"]], text = "Export inferences",
                               command = export.res, state = statevec[3])
    button.widget5 <- tkbutton(inferguienv[["tt"]], text= "Quit", command = quitinf)

    tkgrid(tklabel(inferguienv[["tt"]], text = " "))
    tkgrid(tklabel(inferguienv[["tt"]], text = "  "), button.widget1, tklabel(inferguienv[["tt"]], text = "  "))
    tkgrid(tklabel(inferguienv[["tt"]], text = " "))
    tkgrid(tklabel(inferguienv[["tt"]], text = "  "),button.widget2, tklabel(inferguienv[["tt"]], text = "  "))
    tkgrid(tklabel(inferguienv[["tt"]], text = " "))
    tkgrid(tklabel(inferguienv[["tt"]], text = "  "),button.widget3, tklabel(inferguienv[["tt"]], text = "  "))
    tkgrid(tklabel(inferguienv[["tt"]], text = " "))
    tkgrid(tklabel(inferguienv[["tt"]], text = "  "),button.widget4, button.widget5,
           tklabel(inferguienv[["tt"]], text = "  "))
    tkgrid(tklabel(inferguienv[["tt"]], text = " "))
  }

  inferguienv[["bcnt"]] <- NULL
  inferguienv[["tefname"]] <- character(0)
  statevec <- c("disabled", "disabled", "disabled")
  # Set up status window

  draw.toplev(statevec)


}


