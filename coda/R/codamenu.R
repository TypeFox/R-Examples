"codamenu" <- function () 
{
  coda.options(default=TRUE)
  file.menu <- c("Read BUGS output files", 
                 "Use an mcmc object", 
                 "Quit")
  pick <- menu(file.menu, title = "CODA startup menu")
  if (pick == 0 || pick == 3) 
    return(invisible())
  else if (pick == 1) {
    coda.dat <- read.coda.interactive()
    if (is.null(coda.dat)) {
      return(invisible())
    }
  }
  else if (pick == 2) {
    msg <- "\nEnter name of saved object (or type \"exit\" to quit)"
    repeat {
      cat(msg, "\n")
      outname <- read.and.check(what = character())
      if (outname == "exit" || outname == "\"exit\"") {
        return(invisible())
      }
      else if (!exists(outname)) 
        msg <- "Can't find this object"
      else {
        coda.dat <- eval(parse(text = outname))
        if (is.mcmc(coda.dat)) {
          coda.dat <- mcmc.list(coda.dat)
        }
        if (!is.mcmc.list(coda.dat)) {
          msg <- "Not an mcmc or mcmc.list object"
        }
        else {
          break
        }
      }
    }
  }
  else stop("Invalid option")
  if (is.null(chanames(coda.dat))) {
    chanames(coda.dat) <- chanames(coda.dat, allow.null = FALSE)
  }
  if (is.matrix(coda.dat[[1]]) && is.null(varnames(coda.dat))) {
    varnames(coda.dat) <- varnames(coda.dat, allow.null = FALSE)
  }
  
  ## Check for variables that are linear functions of the
  ## iteration number
  is.linear <- rep(FALSE, nvar(coda.dat))
  for (i in 1:nchain(coda.dat)) {
      for (j in 1:nvar(coda.dat)) {
          lm.out <- lm(as.matrix(coda.dat[[i]])[,j] ~ time(coda.dat))
          if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
              is.linear[j] <- TRUE
          }
      }
  }
  if (any(is.linear)) {
      cat("Dropping the following variables, which are linear\n")
      cat("functions of the iteration number\n")
      print(varnames(coda.dat)[is.linear])
      inset <- varnames(coda.dat)[!is.linear]
      coda.dat <- coda.dat[, inset, drop=FALSE]
  }

  ## Sample size test
  cat("Checking effective sample size ...")
  ess <- lapply(gelman.transform(coda.dat), effectiveSize)
  warn.small <- FALSE
  for (i in 1:length(ess))
    {
      if (any(ess[[i]] < 200))
        warn.small <- TRUE
    }
  if (warn.small)
    {
      cat("\n")
      cat("*******************************************\n")
      cat("WARNING !!!                              \n")
      cat("Some variables have an effective sample  \n")
      cat("size of less than 200 in at least one    \n")
      cat("chain.                                   \n")
      cat("This is too small, and may cause errors  \n")
      cat("in the diagnostic tests                  \n")
      cat("HINT:                                    \n")
      cat("Look at plots first to identify variables\n")
      cat("with slow mixing.  (Choose menu Output   \n")
      cat("Analysis then Plots)                     \n")
      cat("Re-run your chain with a larger sample   \n")
      cat("size and thinning interval. If possible, \n")
      cat("reparameterize your model to improve mixing\n")
      cat("*******************************************\n")
    }
  else
    {
      cat("OK\n")
    }
  current.menu <- "codamenu.main"
  old.opt <- options(warn=-1, show.error.messages=FALSE)
  on.exit(options(old.opt))
  ## Create working data, a subset of coda.dat
  work.dat <- coda.dat

  repeat {
    next.menu <- try(do.call(current.menu, list(work.dat, coda.dat)))
    if (!is.null(class(next.menu)) && class(next.menu) == "try-error")
      {
        if (current.menu == "codamenu.main")
          {
            cat("A crash has occurred in the main menu\nBailing out\n")
            return(invisible());
          }
        else
          {
            cat("\n\n")
            cat("**********************\n")
            cat("An error has occurred\n")
            cat("Returning to main menu\n")
            cat("**********************\n")
            current.menu <- "codamenu.main"
          }
      }
    else
      {
        if (is.list(next.menu) && !is.null(next.menu[["work.dat"]])) {
          work.dat <- next.menu$work.dat
          next.menu <- next.menu[[1]]
        }
        if (next.menu == "quit") {
          if(read.yesno("Are you sure you want to quit", FALSE))
            break
        }
        else current.menu <- next.menu
      }
  }
  invisible()
}

"codamenu.anal" <- function (work.dat, ...) 
{
  next.menu <- "codamenu.anal"
  choices <- c("Plots", "Statistics", "List/Change Options", 
               "Return to Main Menu")
  next.menu.list <- c("plots", "summary", "codamenu.options", 
                      "codamenu.main")
  cat("\n")
  pick <- menu(choices, title = "CODA Output Analysis menu")
  if (pick == 0) 
    next.menu <- "quit"
  else if (next.menu.list[pick] == "summary") {
    if (coda.options("combine.stats")) {
      print(summary(work.dat, quantiles = coda.options("quantiles"), 
                    digits = coda.options("digits")))
    }
    else for (i in 1:nchain(work.dat)) {
      cat(chanames(work.dat, allow.null = FALSE)[i], "\n")
      print(summary(work.dat[[i]], quantiles = coda.options("quantiles"), 
                    digits = coda.options("digits")))
    }
  }
  else if (next.menu.list[pick] == "plots") {
    auto.layout <- !coda.options("user.layout") 
    ask <- TRUE
    repeat {
      if (coda.options("combine.plots")) 
        plot(work.dat, trace = coda.options("trace"), 
             density = coda.options("densplot"),
             smooth = coda.options("lowess"), 
             auto.layout = auto.layout, bwf = coda.options("bandwidth"), 
             combine.chains = !coda.options("combine.plots"), 
             ask = ask)
      else for (i in 1:nchain(work.dat)) {
        plot(work.dat[[i]], trace = coda.options("trace"), 
             density = coda.options("densplot"),
             smooth = coda.options("lowess"), 
             auto.layout = auto.layout, bwf = coda.options("bandwidth"), 
             combine.chains = coda.options("combine.plots"), 
             ask = ask)
      }
      codamenu.ps()
      if (names(dev.cur()) == "postscript") 
        ask <- FALSE
      else break
    }
  }
  else next.menu <- next.menu.list[pick]
  return(next.menu)
}

"codamenu.diags" <- function (work.dat, ...) 
{
  next.menu <- "diags"
  while (next.menu == "diags") {
    choices <- c("Geweke", "Gelman and Rubin", "Raftery and Lewis", 
                 "Heidelberger and Welch", "Autocorrelations",
                 "Cross-Correlations", "List/Change Options",
                 "Return to Main Menu")
    next.menu.list <- c("codamenu.diags.geweke", "codamenu.diags.gelman", 
                        "codamenu.diags.raftery", "codamenu.diags.heidel", 
                        "codamenu.diags.autocorr", "codamenu.diags.crosscorr", 
                        "codamenu.options", "codamenu.main")
    pick <- menu(choices, title = "CODA Diagnostics Menu")
    if (pick == 0) 
      return("quit")
    else next.menu <- next.menu.list[pick]
  }
  return(next.menu)
}

"codamenu.diags.autocorr" <- function (work.dat, ...) 
{
  next.menu <- "codamenu.diags"
  codamenu.output.header("AUTOCORRELATIONS WITHIN EACH CHAIN:", work.dat)
  print(autocorr(work.dat), digits = coda.options("digits"))
  choices <- c("Plot autocorrelations", "Return to Diagnostics Menu")
  pick <- menu(choices, title = "Autocorrelation Plots Menu")
  if (pick == 0) 
    next.menu <- "quit"
  else if (pick == 1) {
    ask <- TRUE
    repeat {
      autocorr.plot(work.dat, auto.layout = !coda.options("user.layout"), 
                    ask = ask)
      codamenu.ps()
      if (names(dev.cur()) == "postscript") 
        ask <- FALSE
      else break
    }
  }
  return(next.menu)
}

"codamenu.diags.crosscorr" <- function (work.dat, ...) 
{
  next.menu <- "codamenu.diags.crosscorr"
  crosscorr.out <- if (coda.options("combine.corr")) {
    crosscorr(work.dat)
  }
  else lapply(work.dat, crosscorr)
  if (coda.options("combine.corr") & nchain(work.dat) > 1) 
    cat("Pooling over chains:", chanames(work.dat, allow.null = FALSE), 
        sep = "\n", collapse = "\n")
  print(crosscorr.out, digits = coda.options("digits"))
  cat("\n")
  choices <- c("Change options",
               "Plot Cross Correlations", 
               "Return to Diagnostics Menu")
  pick <- menu(choices, title = "Cross correlation plots menu")
  if (pick == 0) 
    next.menu <- "quit"
  else
    switch(pick,
           change.tfoption("Combine chains", "combine.corr"), 
           {
             repeat {
               if (coda.options("combine.corr")) 
                 crosscorr.plot(work.dat)
               else {
                 opar <- par(ask = TRUE)
                 lapply(work.dat, crosscorr.plot)
                 par(opar)
               }
               codamenu.ps()
               if (names(dev.cur()) != "postscript") 
                 break
             }
           },
           next.menu <- "codamenu.diags")
  return(next.menu)
}

"codamenu.diags.heidel" <- function (work.dat, ...) 
{
  this.menu <- "codamenu.diags.heidel"
  next.menu <- "codamenu.diags"
  title <- "HEIDELBERGER AND WELCH STATIONARITY AND INTERVAL HALFWIDTH TESTS"
  codamenu.output.header(title, work.dat)
  cat("Precision of halfwidth test =", coda.options("halfwidth"), "\n\n")
  heidel.out <- heidel.diag(work.dat, eps = coda.options("halfwidth"))
  print(heidel.out, digits = coda.options("digits"))
  choices <- c("Change precision", "Return to diagnostics menu")
  pick <- menu(choices)
  if (pick == 0) 
    next.menu <- "quit"
  else if (pick == 1) 
    next.menu <- codamenu.options.heidel(this.menu)
  return(next.menu)
}

"codamenu.diags.raftery" <- function (work.dat, ...) 
{
  next.menu <- this.menu <- "codamenu.diags.raftery"
  codamenu.output.header("RAFTERY AND LEWIS CONVERGENCE DIAGNOSTIC", work.dat)
  print(raftery.diag(work.dat, q = coda.options("q"), r = coda.options("r"), 
                     s = coda.options("s")), digits = coda.options("digits"))
  choices <- c("Change parameters", "Return to diagnostics menu")
  pick <- menu(choices)
  next.menu <- if (pick == 0) 
    "quit"
  else if (pick == 1) {
    codamenu.options.raftery(this.menu)
  }
  else "codamenu.diags"
  return(next.menu)
}

"codamenu.main" <- function (work.dat, ...) 
{
  choices <- c("Output Analysis", "Diagnostics", "List/Change Options", "Quit")
  next.menu.list <- c("codamenu.anal", "codamenu.diags", "codamenu.options", 
                      "quit")
  pick <- menu(choices, title = "CODA Main Menu")
  if (pick == 0) 
    next.menu <- "quit"
  else next.menu <- next.menu.list[pick]
  return(next.menu)
}

"codamenu.diags.gelman" <- function (work.dat, ...) 
{
  next.menu <- this.menu <- "codamenu.diags.gelman"
  if (nchain(work.dat) == 1) {
    cat("\nError: you need more than one chain.\n\n")
    return(next.menu = "codamenu.diags")
  }
  else if (niter(work.dat) <= 50) {
    cat("\nError: you need > 50 iterations in the working data\n")
    return(next.menu = "codamenu.diags")
  }
  z <- window(work.dat, start = niter(work.dat)/2)
  for (i in 2:nchain(z)) {
    for (j in 1:(i - 1)) {
      if (any(apply(as.matrix(z[[i]] - z[[j]]), 2, var)) < 1e-08) {
        cat("\nError: 2nd halves of",
            chanames(z, allow.null = FALSE)[c(j, i)],
            "are identical for at least one variable\n")
        return(next.menu = "codamenu.diags")
      }
    }
  }
  codamenu.output.header("GELMAN AND RUBIN DIAGNOSTIC", work.dat)
  print(gelman.diag(work.dat, transform = TRUE),
        digits = coda.options("digits"))
  choices <- c("Shrink Factor Plots", "Change bin size for shrink plot", 
               "Return to Diagnostics Menu")
  action.list <- c("ShrinkPlot", "ChangeBin", "Return")
  while (next.menu == "codamenu.diags.gelman") {
    pick <- menu(choices, title = "Gelman & Rubin menu")
    if (pick == 0) 
      next.menu <- "quit"
    else switch(action.list[pick], ShrinkPlot = {
      ask <- TRUE
      repeat {
        gelman.plot(work.dat, max.bins = coda.options("gr.max"), 
                    bin.width = coda.options("gr.bin"),
                    auto.layout = !coda.options("user.layout"), 
                    ask = ask)
        codamenu.ps()
        if (names(dev.cur()) == "postscript") 
          ask <- FALSE
        else break
      }
    }, ChangeBin = {
      codamenu.options.gelman(NULL, work.dat)
    }, Return = {
      next.menu <- "codamenu.diags"
    })
  }
  return(next.menu)
}

"codamenu.diags.geweke" <- function (work.dat, ...) 
{
  next.menu <- "codamenu.diags.geweke"
  codamenu.output.header("GEWEKE CONVERGENCE DIAGNOSTIC (Z-score)", work.dat)
  geweke.out <- geweke.diag(work.dat, frac1 = coda.options("frac1"), 
                            frac2 = coda.options("frac2"))
  print(geweke.out, digits = coda.options("digits"))
  choices <- c("Change window size", "Plot Z-scores",
               "Change number of bins for plot", 
               "Return to Diagnostics Menu")
  action.list <- c("ChangeWindow", "Plot", "ChangeBin", "Return")
  while (next.menu == "codamenu.diags.geweke") {
    pick <- menu(choices, title = "Geweke plots menu")
    if (pick == 0) 
      return("quit")
    switch(action.list[pick], ChangeWindow = {
      codamenu.options.geweke.win(NULL)
      geweke.out <- geweke.diag(work.dat,
                                frac1 = coda.options("frac1"), 
                                frac2 = coda.options("frac2"))
      print(geweke.out, digits = coda.options("digits"))
    }, Plot = {
      ask <- TRUE
      repeat {
        if(start(work.dat) >= end(work.dat)) {
          cat("Chain too short: end iteration must be at least twice\n")
          cat("the start iteration\n")
          break
        }
        geweke.plot(work.dat, frac1 = coda.options("frac1"), 
                    frac2 = coda.options("frac2"),
                    nbins = coda.options("geweke.nbin"),
                    auto.layout = !coda.options("user.layout"), 
                    ask = ask)
        codamenu.ps()
        if (names(dev.cur()) == "postscript") 
          ask <- FALSE
        else break
      }
    }, ChangeBin = {
      codamenu.options.geweke.bin(NULL)
    }, Return = {
      next.menu <- "codamenu.diags"
    })
  }
  return(next.menu)
}

"codamenu.options" <- function (work.dat, ...) 
{
  next.menu <- "codamenu.options"
  choices <- c("List current options", "Data Options", "Plot Options", 
               "Summary Statistics Options", "Diagnostics Options", 
               "Output Analysis", "Diagnostics", "Main Menu")
  action.list <- c("ListOptions", "codamenu.options.data", 
                   "codamenu.options.plot", "codamenu.options.stats",
                   "codamenu.options.diag", "codamenu.anal",
                   "codamenu.diags", "codamenu.main")
  pick <- menu(choices, title = "CODA main options menu")
  if (pick == 0) 
    return("quit")
  if (action.list[pick] == "ListOptions") {
    display.working.data(work.dat)
    display.coda.options(stats = TRUE, plots = TRUE, diags = TRUE)
    next.menu <- "codamenu.options"
  }
  else next.menu <- action.list[pick]
  return(next.menu)
}

"codamenu.options.data" <- function (work.dat, coda.dat) 
{
  next.menu <- "codamenu.options.data"
   
  work.vars <- varnames(work.dat)
  work.chains <- chanames(work.dat)
  work.start <- start(work.dat)
  work.end <- end(work.dat)
  work.thin <- thin(work.dat)
  
  choices <- c("List current data options", "Select variables for analysis", 
               "Select chains for analysis", "Select iterations for analysis", 
               "Select thinning interval", "Return to main options menu")
  action.list <- c("ListDataOptions", "SelectVars", "SelectChains", 
                   "SelectIters", "SelectThinInterval", "MainOptionsMenu")
  pick <- menu(choices, title = "CODA data options menu")
  if (pick == 0) 
    return("quit")
  switch(action.list[pick], ListDataOptions = {
    display.working.data(work.dat)
  }, SelectVars = {
    work.vars <- multi.menu(varnames(coda.dat, allow.null = FALSE), 
                            "Select variables for analysis",
                            c("VARIABLE NUMBER", "VARIABLE NAME"),
                            allow.zero = FALSE)
  }, SelectChains = {
    work.chains <- multi.menu(chanames(coda.dat, allow.null = FALSE), 
                              "Select chains for analysis:",
                              c("CHAIN NUMBER", "CHAIN NAME"),
                              allow.zero = FALSE)
  }, SelectIters = {
    cat("\nIterations available = ", start(coda.dat), ":", 
        end(coda.dat), "\n", sep = "")
    work.start <- read.and.check("Enter iteration you wish to start at", 
                                 lower = start(coda.dat),
                                 upper = end(coda.dat),
                                 default = start(work.dat))
    work.end <- read.and.check("Enter iteration you wish to end at", 
                               lower = work.start,
                               upper = end(coda.dat),
                               default = end(work.dat))
  }, SelectThinInterval = {
    cat("\nThinning interval of full data = ", thin(coda.dat), "\n", sep = "")
    work.thin <- read.and.check("Enter thinning interval:", 
                                lower = thin(coda.dat),
                                default = thin(work.dat))
  }, MainOptionsMenu = {
    next.menu <- "codamenu.options"
  })
  if (action.list[pick] != "ListDataOptions" && action.list[pick] != 
      "MainOptionsMenu") {
    cat("Recreating working data...\n")
    wd <- window(coda.dat[, work.vars, drop = FALSE], start = work.start, 
                 end = work.end, thin = work.thin)
    work.dat <- wd[work.chains, drop=FALSE]
  }
  return(list(next.menu, "work.dat"=work.dat))
}

"codamenu.options.diag" <- function (work.dat, ...) 
{
  next.menu <- this.menu <- "codamenu.options.diag"
  choices <- c("Display current diagnostic options",
               "Window sizes for Geweke's diagnostic", 
               "Bin size for plotting Geweke's diagnostic",
               "Bin size for plotting Gelman & Rubin's diagnostic", 
               "Parameters for Raftery & Lewis' diagnostic",
               "Halfwidth precision for Heidelberger & Welch's diagnostic", 
               "Combine chains to calculate correlation matrix",
               "Return to main options menu")
  pick <- menu(choices, title = "CODA diagnostics options menu")
  if (pick == 0) 
    return("quit")
  switch(pick, display.coda.options(diags = TRUE),
         next.menu <- codamenu.options.geweke.win(this.menu), 
         next.menu <- codamenu.options.geweke.bin(this.menu), 
         next.menu <- codamenu.options.gelman(this.menu, work.dat),
         next.menu <- codamenu.options.raftery(this.menu), 
         next.menu <- codamenu.options.heidel(this.menu),
         {
           change.tfoption("Do you want to combine all chains to calculate correlation matrix",  "combine.corr")
         }, next.menu <- "codamenu.options")
  return(next.menu)
}

"codamenu.options.gelman" <- function (last.menu, work.dat) 
{
  choices <- c("Default: bin width = 10; maximum number of bins = 50", 
               "User-specified bin width",
               "User-specified total number of bins")
  pick <- menu(choices, title = "Options for defining bin size to plot Gelman-Rubin-Brooks diagnostic")
  if (pick == 0) 
    return("quit")
  switch(pick, {
    coda.options(gr.max = 50)
    coda.options(gr.bin = 10)
  }, {
    coda.options(gr.max = Inf)
    default <- if (coda.options("gr.bin") == 0) 
      10
    else coda.options("gr.bin")
    msg <- "Enter required bin width:"
    coda.options(gr.bin = read.and.check(msg, lower = 1, 
                   upper = niter(work.dat) - 50, default = default))
  }, {
    coda.options(gr.bin = 0)
    default <- if (is.infinite(coda.options("gr.max"))) 
      50
    else coda.options("gr.max")
    msg <- "Enter total number of bins required:"
    coda.options(gr.max = read.and.check(msg, lower = 1, 
                   upper = niter(work.dat) - 50, default = default))
  })
  return(last.menu)
}

"codamenu.options.geweke.bin" <- function (last.menu) 
{
  msg <- "Enter number of bins for Geweke-Brooks plot"
  ans <- read.and.check(msg, what=numeric(), lower=1,
                        default=coda.options("geweke.nbin"))
  coda.options(geweke.nbin = ans)
  return(last.menu)
}

"codamenu.options.geweke.win" <-
function (last.menu) 
{
  msg1 <- "Enter fraction of chain to include in 1st window:"
  msg2 <- "Enter fraction of chain to include in 2nd window:"
  ans1 <- ans2 <- 1
  while (ans1 + ans2 >= 1) {
    ans1 <- read.and.check(msg1, lower = 0, upper = 1,
                           default = coda.options("frac1"))
    ans2 <- read.and.check(msg2, lower = 0, upper = 1,
                           default = coda.options("frac2"))
    ## Check that sum of fractions doesn't exceed 1.0
    if (ans1 + ans2 >= 1) 
      cat("Error: Sum of fractions in 1st and 2nd windows must be < 1.0\n")
  }
  coda.options(frac1 = ans1, frac2 = ans2)
  return(last.menu)
}

"codamenu.options.heidel" <- function (last.menu) 
{
  coda.options(halfwidth = read.and.check("Enter precision for halfwidth test",
                 lower = 0, default = coda.options("halfwidth")))
  return(last.menu)
}

"codamenu.options.plot" <- function (...) 
{
  next.menu <- "codamenu.options.plot"
  choices <- c("Show current plotting options",
               "Plot trace of samples", 
               "Plot kernel density estimate",
               "Add smooth line through trace plot", 
               "Combine chains",
               "Single plot per page",
               "Specify page layout for plots", 
               "Select bandwidth function for kernel smoothing",
               "Return to main options menu")
  pick <- menu(choices, title = "CODA plotting options menu")
  if (pick == 0) 
    return("quit")
  switch(pick,
         display.coda.options(plots = TRUE),
         change.tfoption(choices[2], "trace"),
         change.tfoption(choices[3], "densplot"),
         change.tfoption(choices[4], "lowess"),
         change.tfoption(choices[5], "combine.plots"), 
         {
           ans <- read.yesno(choices[6], default=TRUE)
           if(ans) {
             coda.options(user.layout = TRUE)
             par(mfrow = c(1,1))
           }
         },
         {
           change.tfoption("Do you want to specify your own page layout for the plots", "user.layout")
           if (coda.options("user.layout")) {
             mrows <- read.and.check("Enter number of rows per page",
                                     lower = 1, upper = 7)
             mcols <- read.and.check("Enter number of columns per page",
                                     lower = 1, upper = 8)
             par(mfrow = c(mrows, mcols))
           }
         }, {
           next.menu <- "codamenu.options.plot.kernel"
         }, NULL)
  if (pick == length(choices)) 
    next.menu <- "codamenu.options"
  return(next.menu)
}

"codamenu.options.plot.kernel" <- function (...) 
{
  if (!coda.options("densplot")) {
    cat("\nNo density plots requested - this option is irrelevant\n")
  }
  else {
    kernel.menu <- c("Smooth (0.25 * sample range)",
                     "Coarse (Silverman 1986 eqn. 3.28 & 3.30)", 
                     "User-defined function",
                     "Return to Plotting Options Menu")
    pick1 <- menu(kernel.menu, title = "Select kernel bandwidth function")
    if (pick1 == 0) 
      return("quit")
    switch(pick1, {
      bwf <- function(x) {
        (max(x) - min(x))/4
      }
      coda.options(bandwidth = bwf)
    }, {
      bwf <- function(x) {
        1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2
      }
      coda.options(bandwidth = bwf)
    }, {
      func.OK <- FALSE
      while (!func.OK) {
        cat("Enter bandwidth as an expression in terms of x,\n")
        cat("the vector of sampled values, e.g. \n")
        cat("(max(x) - min(x)) / 4\n")
        ans <- scan(what = character())
        if (length(ans) > 0) {
          bwf <- "function(x){"
          for (i in 1:length(ans)) {
            bwf <- paste(bwf, ans[i], sep = "")
          }
          bwf <- paste(bwf, "}", sep = "")
          bwf <- try(eval(parse(text = bwf)), silent=TRUE)
          if (inherits(bwf, "try-error")) {
            cat("Invalid expression\n")
          }
          else {
            ## Carry out simple test to check whether the
            ## function entered makes sense
            ##
            bw <- try(bwf(1:10), silent=FALSE)
            if (inherits(bw, "try-error")) {
               cat("Error calling function with input 1:10\n")
            }
            else {
               func.OK <- is.numeric(bw) && (length(bw) == 1)
               if(!func.OK) {
                 cat("This is not a suitable function: it must return a\n")
                 cat("single numeric value given a numeric vector x.\n")
               }
            }
          }
        }
      }
      coda.options(bandwidth = bwf)
    }, NULL)
  }
  return("codamenu.options.plot")
}

"codamenu.options.raftery" <- function (last.menu) 
{
  coda.options(q = read.and.check("Enter quantile to be estimated:", 
                 lower = 0, upper = 1, default = coda.options("q")))
  coda.options(r = read.and.check("Enter required precision:", 
                 upper = coda.options("q"), default = coda.options("r")))
  coda.options(s = read.and.check("Enter required probability:", 
                 lower = 0, upper = 1, default = coda.options("s")))
  return(last.menu)
}

"codamenu.options.stats" <- function (...) 
{
  next.menu <- "codamenu.options.stats"
  choices <- c("Display current statistics options", "Combine chains for summary statistics", 
               "Quantiles for summary statistics", "Number of significant digits for printing", 
               "Return to main options menu")
  pick <- menu(choices, title = "CODA options for summary statistics")
  if (pick == 0) 
    return("quit")
  switch(pick, display.coda.options(stats = TRUE), {
    mssg <- "Do you want to combine all chains when calculating summary statistics"
    change.tfoption(mssg, "combine.stats")
  }, {
    mssg <- paste("Enter quantiles required, separated by commas\n(Default =", 
                  paste(coda.options("quantiles"), collapse = ", "))
    repeat {
      cat("\n", mssg, "\n")
      if (is.R()) {
        ans <- as.numeric(scan(what = character(), sep = ",", nlines = 1,
                               quiet = TRUE))
      }
      else {
        ans <- as.numeric(scan(what = character(), sep = ",", nlines = 1))
      }
      if (length(ans) == 0) 
        ans <- coda.options("quantiles")
      if (any(is.na(ans))) 
        mssg <- "You must enter numeric values"
      else if (any(ans >= 1) || any(ans <= 0)) 
        mssg <- "You must enter values between 0 and 1"
      else break
    }
    if (length(ans) > 0) 
      coda.options(quantiles = sort(ans))
  }, {
    mssg <- "Enter number of significant digits to be printed"
    ans <- read.and.check(mssg, what = integer(), lower = 0, 
                          default = coda.options("digits"))
    coda.options(digits = ans)
  }, {
    next.menu <- "codamenu.options"
  })
  return(next.menu)
}

"display.working.data" <-  function (data)
{
  cat("WORKING DATA\n")
  cat("============\n")
  cat("Variables selected : ",
      paste(varnames(data, allow.null = FALSE), collapse=", ")
      ,"\n", sep="")
  cat("Chains selected    : ",
      paste(chanames(data, allow.null = FALSE), collapse=", ")
      , "\n", sep="")
  cat("Iterations - start : ", start(data), "\n", sep="")
  cat("               end : ", end(data), "\n", sep="")
  cat("Thinning interval  : ", thin(data), "\n", sep="")
  cat("\n")
}

"display.coda.options" <-
  function (stats = FALSE, plots = FALSE, diags = FALSE) 
{
  cat("\nCurrent option settings:")
  cat("\n=======================\n\n")
  if (stats) {
    cat("SUMMARY STATISTICS OPTIONS\n")
    cat("==========================\n\n")
    cat("Combine chains     : ", coda.options("combine.stats"), "\n", sep="")
    cat("Quantiles          : ",
        paste(coda.options("quantiles") * 100, "%", sep="", collapse = ", "),
        "\n", sep="")
    cat("Significant digits : ", coda.options("digits"), "\n", sep="")
    cat("\n")
  }
  if (plots) {
    cat("PLOTTING OPTIONS\n")
    cat("================\n\n")
    cat("Trace               : ", coda.options("trace"),         "\n", sep="")
    cat("Density             : ", coda.options("densplot"),      "\n", sep="")
    cat("Smooth lines        : ", coda.options("lowess"),        "\n", sep="")
    cat("Combine chains      : ", coda.options("combine.plots"), "\n", sep="")
    cat("User-defined layout : ", coda.options("user.layout"),   "\n", sep="")
    if(coda.options("user.layout")) {
      cat("                    : ", paste(par("mfrow"), collapse=" X "), "\n", sep="")
    }
    cat("Bandwidth function  :\n")
    print(coda.options("bandwidth"))
    cat("\n")
  }
  if (diags) {
    cat("DIAGNOSTICS OPTIONS\n")
    cat("===================\n\n")
    cat("Geweke\n")
    cat("------\n")
    cat("Window 1 fraction  : ", coda.options("frac1"),      "\n", sep="")
    cat("Window 2 fraction  : ", coda.options("frac2"),      "\n", sep="")
    cat("Number of bins     : ", coda.options("geweke.nbin"), "\n", sep="")
    cat("\n")
    
    cat("Gelman & Rubin\n")
    cat("--------------\n")
    cat("Bin width          : ", coda.options("gr.bin"), "\n", sep="")
    cat("Max number of bins : ", coda.options("gr.max"), "\n", sep="")
    cat("\n")
    
    cat("Raftery & Lewis\n")
    cat("---------------\n")
    cat("Quantile (q)       : ", coda.options("q"), "\n", sep="")
    cat("Precision (+/- r)  : ", coda.options("r"), "\n", sep="")
    cat("Probability (s)    : ", coda.options("s"), "\n", sep="")
    cat("\n")
           
    cat("Cross-correlations\n")
    cat("------------------\n")
    cat("Combine chains     : ", coda.options("combine.corr"), "\n", sep="") 
    cat("\n")
  }
  invisible()
}

"read.coda.interactive" <- function () 
{
  repeat {
    cat("Enter CODA index file name\n")
    cat("(or a blank line to exit)\n")
    if (is.R()) {
      index.file <- scan(what = character(), sep = "\n", strip.white = TRUE, 
                         nlines=1, quiet=TRUE)
    }
    else {
      index.file <- scan(what = character(), sep = "\n", strip.white = TRUE)
    }
    if (length(index.file) == 0)
      return(invisible())
    cat("Enter CODA output file names, separated by return key\n")
    cat("(leave a blank line when you have finished)\n")
    if (is.R()) {
      output.files <- scan(what = character(), sep = "\n", strip.white = TRUE, 
                           quiet = TRUE)
    }
    else {
      output.files <- scan(what = character(), sep = "\n", strip.white = TRUE)
    }
    all.files <- c(index.file, output.files)
    if (any(!file.exists(all.files))) {
      cat("The following files were not found:\n")
      cat(paste(all.files[!file.exists(all.files)], 
                collapse = "\n"), "\n\n")
    }
    else break
  }
  nfiles <- length(output.files)
  chains <- vector("list", nfiles)
  names(chains) <- output.files
  for (i in 1:nfiles) chains[[i]] <- read.coda(output.files[i], index.file)
  return(mcmc.list(chains))
}

"codamenu.ps" <- function () 
{
  if (names(dev.cur()) == "postscript") {
    dev.off()
  }
  else {
    cat("\nSave plots as a postscript file (y/N) ?\n")
    ans <- readline()
    if (length(ans) == 0) 
      ans <- "n"
    if (ans == "Y" | ans == "y") {
      repeat {
        mssg <- "Enter name you want to call this postscript file"
        ps.name <- read.and.check(mssg, what = character(), 
                                  default = "Rplots.ps")
        if (file.exists(ps.name)) {
          pick <- menu(title = "File exists", choices = c("overwrite", 
                                                "choose another file name"))
          if (pick == 1) 
            break
        }
        else {
          break
        }
      }
      postscript(file = ps.name)
    }
  }
  return(dev.cur())
}

"codamenu.output.header" <- function (title, data) 
{
  ##
  ## A short header: common to most codamenu output
  ##
  cat("\n", title, sep = "")
  cat("\n", paste(rep("=", nchar(title)), collapse = ""), "\n\n", 
      sep = "")
  cat("Iterations used = ", start(data), ":", end(data), 
      "\n", sep = "")
  cat("Thinning interval =", thin(data), "\n")
  cat("Sample size per chain =", niter(data), "\n\n")
  invisible()
}
