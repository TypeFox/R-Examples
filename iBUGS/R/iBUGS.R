#' An Intelligent Interface to WinBUGS/OpenBUGS/JAGS
#'
#' Create a GUI by \pkg{gWidgetsRGtk2} and all options for \pkg{R2WinBUGS} and
#' \pkg{R2jags} can be set in the GUI. The letter \sQuote{i} can be interpreted
#' as \sQuote{\bold{i}ntelligent} or \sQuote{\bold{i}nterface} -- depends on
#' what you think.
#'
#' \pkg{iBUGS} will try to find the directories of WinBUGS/OpenBUGS/JAGS in your
#' system and use them when calling \code{\link[R2WinBUGS]{bugs}} or
#' \code{\link[R2jags]{jags}} in \pkg{R2WinBUGS} or \pkg{R2jags}. For most
#' users, this search will succeed, unless WinBUGS/OpenBUGS/JAGS were
#' not installed in the default directory.
#'
#' \pkg{iBUGS} will also try to get the data object names from the current R
#' session and guess the parameter names in a BUGS model. Click the
#' \dQuote{Preference} button and you will see the lists of names.
#'
#' More intelligence is to be added.
#'
#' @return Invisible \code{NULL}.
#' @author Yihui Xie <\url{http://yihui.name}>
#' @seealso \code{\link[R2WinBUGS]{bugs}}, \code{\link[R2jags]{jags}}, \code{\link{bugs.options}}
#' @keywords utilities
#' @export
#' @examples
#'
#' \dontrun{
#'
#' iBUGS()
#'
#' }
#'
iBUGS = function() {
  options(guiToolkit = "RGtk2")
  auto = "No"
  fn = "bugsdata.txt"
  g = ggroup(horizontal = FALSE, container = gw0 <- gwindow("iBUGS - Intelligent (Open|Win)BUGS Interface"))
  g1 = ggroup(container = g, expand = TRUE)
  if (.Platform$OS.type == "windows") {
    g7 = gframe(container = g, text = "Program")
    items = c("OpenBUGS", "WinBUGS", "JAGS")
    rb = gradio(items, horizontal = TRUE, container = g7)
    bugs.options(program = svalue(rb))
    addHandlerClicked(rb, handler = function(h, ..) {
      bugs.options(program = svalue(h$obj))
    })
  } else bugs.options(program = "JAGS")
  g2 = ggroup(container = g)
  txt = gtext("model\n{\n\t## likelihood\n\tfor (i in 1:N) {\n\t\t\n\t}\n\t## prior\n\n}", container = g1, wrap = FALSE, 
              font.attr = c(family = "monospace", size = "large", styles = "normal", weights = "light"), expand = TRUE)
  gbutton("Open", container = g2, handler = function(h, ...) {
    s = gfile("Open Model", initialfilename = basename(bugs.options("model.file")))
    if (!is.na(s)) {
      svalue(txt) = readLines(s)
      tag(txt, "src.file") = s
      tooltip(txt) = s
      bugs.options(model.file = s)
    }
    focus(txt)
  })
  gbutton("Save", container = g2, handler = function(h, ...) {
    s = tag(txt, "src.file")
    if (is.null(s)) 
      s = gfile("Save the Model", type = "save", initialfilename = basename(bugs.options("model.file")))
    if (!is.na(s)) {
      writeLines(svalue(txt), s)
      bugs.options(model.file = s)
      tag(txt, "src.file") = s
      tooltip(txt) = s
      galert(paste("BUGS model saved to", s))
    }
  })
  gbutton("Preferences", container = g2, handler = function(h, ...) {
    g = ggroup(horizontal = FALSE, container = gw <- gwindow("Options"))
    ## data & parameters
    g1 = ggroup(container = g, expand = TRUE)
    ## other options
    g2 = glayout(container = ggroup(container = g), expand = TRUE, spacing = 2)
    # JAGS can be auto-updated until converge using R2jags::autojags;
    if (bugs.options("program") == "JAGS") {
      g8 = gframe(container = g, text = "Auto-update until the model converges?")
      items.auto = c("No", "Yes")
      rb.auto = gradio(items.auto, horizontal = TRUE, container = g8)
      addHandlerClicked(rb.auto, handler = function(h, ..) {
        auto = svalue(h$obj)
      })
    }
    ## buttons
    g3 = ggroup(container = g)
    g.data = gtable(data.frame(data = unique(c(bugs.options("data"), unlist(sapply(grep("^[^(package:)]", search(), value = TRUE), 
                                                                                   ls))))), multiple = TRUE, container = g1, expand = TRUE)
    svalue(g.data) = bugs.options("data")
    addHandlerMouseMotion(g.data, handler = function(h, ...) focus(h$obj))
    con = textConnection(svalue(txt))
    model.line = gsub("^[[:space:]]+|[[:space:]]+$", "", readLines(con))
    b = which(model.line == "}")
    par.list = ""
    if (length(b) >= 2) {
      par.list = gsub("([[:space:]]+|)(~|<-).*", "", grep("~|<-", model.line[(b[length(b) - 1] + 1):(b[length(b)] - 1)], value = TRUE))
    } else galert("I cannot find any parameters in your model!", "Warning")
    par.list = unique(c(bugs.options("parameters.to.save"), par.list))
    par.list = par.list[nzchar(par.list)]
    g.parameters.to.save = gtable(data.frame(parameters.to.save = par.list, stringsAsFactors = FALSE), multiple = TRUE, container = g1, 
                                  expand = TRUE)
    svalue(g.parameters.to.save) = bugs.options("parameters.to.save")
    close(con)
    addHandlerMouseMotion(g.parameters.to.save, handler = function(h, ...) focus(h$obj))
    addhandlerdoubleclick(g.parameters.to.save, handler = function(h, ...) {
      gi = ginput("Well, you beat me! It seems I have not figured out all the parameter names in your model... Please manually input them here:", 
                  title = "Parameters to Save", icon = "question")
      if (!is.na(gi) && gi != "") 
        g.parameters.to.save[, ] = c(g.parameters.to.save[], gi)
      g.parameters.to.save[, ] = g.parameters.to.save[][!duplicated(g.parameters.to.save[, , drop = FALSE])]
      g.parameters.to.save[, ] = g.parameters.to.save[, ][nzchar(g.parameters.to.save[, ])]
    })
    opt.names = setdiff(names(bugs.options()), c("data", "parameters.to.save"))
    if (bugs.options("program") == "JAGS") 
      opt.names = setdiff(opt.names, c("n.sims", "bin", "debug", "codaPkg", "bugs.directory", "program", "clearWD", "bugs.seed", 
                                       "summary.only", "save.history", "over.relax", "progress.bar")) else opt.names = setdiff(opt.names, c("jags.seed", "refresh", "progress.bar"))
    for (i in 1:ceiling(length(opt.names)/3)) {
      for (j in 1:3) {
        if (3 * i - (3 - j) <= length(opt.names)) {
          g2[i, 2 * j - 1] = opt.names[3 * i - (3 - j)]
          g2[i, 2 * j, expand = TRUE] = (tmp <- gedit(bugs.options(opt.names[3 * i - (3 - j)]), container = g2, expand = TRUE))
          id(tmp) = opt.names[3 * i - (3 - j)]
          addhandlerblur(tmp, handler = function(h, ...) {
            v = svalue(h$obj)
            warn.level = getOption("warn")
            options(warn = -1)
            eval(parse(text = sprintf("bugs.options(%s = if (!nzchar(v)) NULL else ifelse(!is.na(as.integer(v)),  as.integer(v), ifelse(!is.na(as.logical(v)), as.logical(v), v)))", 
                                      id(h$obj))))
            options(warn = warn.level)
            ## need validation here
          })
          addHandlerMouseMotion(tmp, handler = function(h, ...) focus(h$obj))
        }
      }
    }
    gbutton("ok", container = g3, handler = function(h, ...) {
      bugs.options(data = as.character(svalue(g.data)), parameters.to.save = as.character(svalue(g.parameters.to.save)))
      dispose(gw)
    })
    gbutton("help", container = g3, handler = function(h, ...) {
      gw1 = gwindow(paste("Help on ", ifelse(bugs.options("program") == "JAGS", "jags", "bugs")), visible = FALSE)
      size(gw1) = c(500, 500)
      g4 = ggroup(horizontal = FALSE, container = gw1)
      g5 = ggroup(container = g4, expand = TRUE)
      helpWidget = ghelp(container = g5, expand = TRUE)
      visible(gw1) = TRUE
      add(helpWidget, ifelse(bugs.options("program") == "JAGS", "R2jags::jags", "R2WinBUGS:::bugs"))
      g6 = ggroup(container = g4)
      gbutton("cancel", container = g6, handler = function(h, ...) {
        dispose(gw1)
      })
    })
    gbutton("cancel", container = g3, handler = function(h, ...) {
      dispose(gw)
    })
    size(g1) = c(200, 200)
  })
  gbutton("Execute", container = g2, handler = function(h, ...) {
    writeLines(svalue(txt), bugs.options("model.file"))
    if (file.exists("bugsdata.txt")) 
      source(file = "bugsdata.txt")
    out = with(bugs.options(), {
      if (bugs.options("program") == "JAGS") 
        jags(data, if (!is.null(bugs.options("inits"))) 
          eval(parse(text = bugs.options("inits"))) else NULL, parameters.to.save, model.file, n.chains, n.iter, n.burnin, n.thin, DIC, digits, working.directory = NULL, jags.seed, 
             refresh, progress.bar = "none") else bugs(data, if (!is.null(bugs.options("inits"))) 
               eval(parse(text = bugs.options("inits"))) else NULL, parameters.to.save, model.file, n.chains, n.iter, n.burnin, n.thin, n.sims, bin, debug, DIC, digits, codaPkg, bugs.directory, 
                                                       program, working.directory, clearWD, useWINE = FALSE, WINE = NULL, newWINE = FALSE, WINEPATH = NULL, bugs.seed, summary.only, 
                                                       save.history, over.relax)
    })
    # working.directory and progress.bar are fixed for jags
    if (bugs.options("program") == "JAGS" && auto == "Yes" && max(get(bugs.options("model.name"))$BUGSoutput$summary[, "Rhat"]) >= 
          1.1) 
      out = autojags(out)
    message(sprintf("(*) Returned values saved to the R object '%s';\n    you may play with it now, e.g. press the buttons of Print and Plot", 
                    bugs.options("model.name")))
    dput(out, bugs.options()$model.name)
  })
  gbutton("Print", container = g2, handler = function(h, ...) {
    if (file.exists(bugs.options()$model.name)) 
      print(dget(bugs.options()$model.name)) else galert("Please execute a model first!", "Warning")
  })
  gbutton("Plot", container = g2, handler = function(h, ...) {
    if (file.exists(bugs.options()$model.name)) 
      plot(dget(bugs.options()$model.name)) else galert("Please execute a model first!", "Warning")
  })
  gbutton("Demo", container = g2, handler = function(h, ...) {
    if (isTRUE(gconfirm("I will overwrite the current model and show the demo. Do you want to continue?"))) {
      Ex.Y = rnorm(100, mean = rnorm(1, mean = 5))
      Ex.N = 100
      svalue(txt) = c("## I have created two objects ", "##   in your R session: Ex.Y, Ex.N", "model", "{", "  ## sampling distribution: Normal(mu, 1)", 
                      "  for(i in 1:Ex.N){", "    Ex.Y[i] ~ dnorm(mu,1)", "  }", "  ## prior: Normal(5, 1)", "  mu ~ dnorm(5, 1)", "}")
      bugs.options(data = c("Ex.Y", "Ex.N"), parameters.to.save = "mu")
      dump(c("Ex.Y", "Ex.N"), file = fn)
    }
  })
  gbutton("Help", container = g2, handler = function(h, ...) {
    if (isTRUE(gconfirm(paste(c(1:3, rep("*", 3)), ". ", c("Write the BUGS model in the textbox;", "Click the 'Preference' button and set the options;", 
                                                           "Execute! Done!\n\n", "Or you can try the demo first;", "I might crash R if certain options were not correctly specified! (to be improved...)", 
                                                           sprintf("Do you want to view the package vignette (%s) now?", system.file("doc", "iBUGS.pdf", package = "iBUGS"))), sep = "", 
                              collapse = "\n"), "Help"))) 
      system(paste(getOption("pdfviewer"), system.file("doc", "iBUGS.pdf", package = "iBUGS")))
  })
  invisible(NULL)
  gbutton("cancel", container = g2, handler = function(h, ...) {
    dispose(gw0)
  })
}
