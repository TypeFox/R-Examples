# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.
set.max.replicates <- function (first = TRUE) {
  value <- .cur.db@Prefs@Bootgam.prefs$n
  if (first == TRUE) {
    cat("Number of (maximum) bootstrap replicates to perform \n")
    if (is.null(value)) {
      cat("The current value is NULL...\n")
    }
    else {
            cat("The current number is", value, "...\n")
          }
    cat("\nPlease type the new number: ")
  }
  ans <- as.numeric(readline())
  if (ans == "NULL" || ans == "null") {
      Recall(first = FALSE)
  } else {
    if (ans > 0) {
      .cur.db@Prefs@Bootgam.prefs$n <- ans
      c1 <- call("assign",pos = 1, ".cur.db", .cur.db)
      eval(c1)
      invisible()
      return()
    } else {
      cat("Please enter a numeric value larger than 0 ")
      Recall(first = FALSE)
    }
  }
}

bootgam.conv.crit1 <- function() {
  cat("\nType the critical value of the fluctuation ratio\n")
  cat("you want to use (0 to exit)\n")
  ans <- readline()
  if (ans == 0) {
    return()
  }
  if (ans < 1) {
    cat("The number must be greater than 1\n")
    Recall()
  }
  .cur.db@Prefs@Bootgam.prefs$crit1.conv <- ans
  c2 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c2)
  invisible()
  return()
}

bootgam.conv.crit2 <- function(skip=F) {
  if(!skip) {
    cat("\nType the lowest important relative inclusion frequency\n")
    cat("or return for the default (0.2):")
    ans <- readline()
    if(ans == "") {
      .cur.db@Prefs@Bootgam.prefs$crit2.liif <- 0.2
    } else  if(ans < 0 || ans > 1) {
      cat("The number must be greater than 0 and lower than 1\n")
      Recall()
    } else {
      .cur.db@Prefs@Bootgam.prefs$crit2.liif <- ans
    }
  }

  cat("\nType the lowest absulute joint inclusion frequency\n")
  cat("or return for the default (25):")
  ans2 <- readline()
  if(ans2 == "") {
    .cur.db@Prefs@Bootgam.prefs$crit2.ljif.conv <- 25
  } else if(ans2 < 1) {
    cat("The number must be greater than 1\n")
    Recall(skip=T)
  } else {
    .cur.db@Prefs@Bootgam.prefs$crit2.ljif.conv <- ans2
  }
  c3 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c3)
  invisible()
  return()
}

specify.start.model <- function () {
  if (any(is.null(covs <- xvardef("covariates", .cur.db)))) {
    cat("\nThe current data base has no covariates defined\n")
    invisible()
    return()
  }
  cat("\nThe following covariates are defined in the current data base:\n")
  cat(covs, fill = 60)
  cat("\nType the names of the covariates that should be included in the\n")
  cat("\nmodel and end with a blank line:\n")
  st.covs <- scan(what = character())
  if(length(st.covs) == 0) {
    st.covs <- NULL
  }
  .cur.db@Prefs@Bootgam.prefs$start.mod <- st.covs
  c4 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c4)
  invisible()
  return()
}

normalize.median <- function () {
  if (.cur.db@Prefs@Bootgam.prefs$median.norm == FALSE) {
    .cur.db@Prefs@Bootgam.prefs$median.norm <- TRUE
    cat ("\nNormalize to median is now set to ON\n")
  } else {
    .cur.db@Prefs@Bootgam.prefs$median.norm <- FALSE
    cat ("\nNormalize to median is now set to OFF\n")
  }
  c5 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c5)
  
  invisible()
  return()
}

bg.conv.crit2 <- function(skip = FALSE) {
  if(!skip) {
    cat("\nType the lowest important relative inclusion frequency\n")
    cat("or return for the default (0.2):")
    ans <- readline()
    if(ans == "") {
      .cur.db@Prefs@Bootgam.prefs$liif <- 0.2
    } else  if(ans < 0 || ans > 1) {
      cat("The number must be greater than 0 and lower than 1\n")
      Recall()
    } else {
      .cur.db@Prefs@Bootgam.prefs$liif <- ans
    }
  }
  cat("\nType the lowest absulute joint inclusion frequency\n")
  cat("or return for the default (25):")
  ans2 <- readline()
  if(ans2 == "") {
    .cur.db@Prefs@Bootgam.prefs$ljif.conv <- 25
  } else if(ans2 < 1) {
    cat("The number must be greater than 1\n")
    Recall(skip=T)
  } else {
    .cur.db@Prefs@Bootgam.prefs$ljif.conv <- ans2
  }
  c6 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c6)
  invisible()
  return()
}

bg.conv.crit1 <- function() {
  cat("\nType the critical value of the fluctuation ratio\n")
  cat("you want to use (0 to exit)\n")
  ans <- readline()
  if(ans == 0) {
    return()
  }
  if(ans < 1) {
    cat("The number must be greater than 1\n")
    Recall()
  }
  .cur.db@Prefs@Bootgam.prefs$conv.value <- ans
  c7 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c7)
  invisible()
  return()
}

change.conv.alg <- function () {
  cat("\nSpecify the algorithm for convergence calculations\n")
  choices <- c("Return to previous menu",
	       "Fluctuation ratio",
	       "Lowest absolute joint inclusion frequency"
	       )
  pick <- menu(choices)
  switch(pick,
	 return(),
	 { .cur.db@Prefs@Bootgam.prefs$algo <- "fluct.ratio"
           c8 <- call("assign",pos = 1, ".cur.db", .cur.db)
           eval(c8)
           bg.conv.crit1()
         },
	 { .cur.db@Prefs@Bootgam.prefs$algo <- "liif"
           c9 <- call("assign",pos = 1, ".cur.db", .cur.db)
           eval(c9)
  	   bg.conv.crit2()
         })
  Recall()
  return()
}

specify.start.check <- function () {
  cat("\nType the iteration number at which you want to start checking\n")
  cat("convergence (0 to exit):")
  ans <- readline()
  if(ans == 0) {
    return()
  }
  if(ans < 0) {
    cat("The number must be positive.\n")
    Recall()
  }
  .cur.db@Prefs@Bootgam.prefs$start.check <- ans
  c10 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c10)
  invisible()
  return()
}

specify.interval <- function () {
  cat("\nType the interval at which the convergence should be \n")
  cat("checked (0 to exit):")
  ans <- readline()
  if (ans == 0) {
    return()
  }
  if (ans < 0) {
    cat ("The number must be positive.\n")
    Recall()
  }
  .cur.db@Prefs@Bootgam.prefs$check.interval <- ans
  c11 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c11)
  invisible()
  return()
}

exclude.individuals <- function () {
  cat("Please type the ID number of the individuals you want to exclude\n")
  cat("from the bootgam analysis and finish with a blank line:\n")
  inds <- scan(what = character())
  if(length(inds) == 0) {
    inds <- NULL
  }
  .cur.db@Prefs@Bootgam.prefs$excluded.ids <- inds
  c12 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c12)
  invisible()
  return()
}

set.seed.number <- function () {
  cat ("Type a seed number between 1 and 1000 (return to exit):")
  ans <- readline()
  if (ans == "") {
    return()
  }
  if (as.numeric(ans) >= 1 && as.numeric(ans) <= 1000) {
    .cur.db@Prefs@Bootgam.prefs$seed <- ans
  } else {
    cat("The number must be between 1 and 1000\n")
    Recall()
  }
  c13 <- call("assign",pos = 1, ".cur.db", .cur.db)
  eval(c13)
}


bootgam.settings.menu <- function () {
  choices <- c("Return to the previous menu ->",
               "List current bootGAM settings",
               "Set maximum number of bootstrap replicates",
	       "Change convergence algorithm",
	       "Specify iteration to start check convergence at",
	       "Specify at what interval to check the convergence",
	       "Set seed number",
	       "Exclude individuals",
	       "Specify starting model"
               )
  title = "\nBootGAM SETTINGS MENU\n  \\main\\covariate model\\BootGAM\\settings for the BootGAM"
  pick <- menu(choices, title = title)
  qx <- 0
  switch(pick + 1, qx <- 2, qx <- 1,
         list.bootgam.settings(eval(as.name(".cur.db"))),
         set.max.replicates(),
         change.conv.alg(),
         specify.start.check(),
         specify.interval(),
         set.seed.number(),
         exclude.individuals(),
         specify.start.model()
         )
  if (qx == 2) {
    return(invisible(2))
  }
  else {
    if (qx == 1) {
      return(invisible(0))
    }
    else {
      Recall()
    }
  }
}

list.bootgam.settings <- function (object) {
  if (exists("object")) {
    cat(paste("\nThe current run number is ", object@Runno, ".\n\n", sep = ""))
    if (!any(is.null(object@Prefs@Xvardef$parms))) {
      cat("\nMaximum number of bootstrap replicates: ", object@Prefs@Bootgam.prefs$n)
      cat(":", object@Prefs@Bootgam.prefs$n)
      cat("\nConvergence algorithm to use: ", object@Prefs@Bootgam.prefs$algo)
      if (object@Prefs@Bootgam.prefs$algo == "fluct.ratio") {
        cat("\nConvergence criterion: ", object@Prefs@Bootgam.prefs$conv.value)
      } else {
        cat("\nLowest importan relative inclusion freq: ", object@Prefs@Bootgam.prefs$liif)
        cat("\nCritical value (ljif): ", object@Prefs@Bootgam.prefs$ljif.conv)
      }
      cat("\nStarting model: ", object@Prefs@Bootgam.prefs$n)
    }
 } else {
    cat("The current run number is", object@Runno, "but no matching database was found.\n")
  }
}

bootgam.menu <- function () {
  choices <- c("Return to previous menu ->",
               "Run a bootGAM",
               "Summarize bootGAM",
               "Plot bootGAM results ->",
               "Settings for the BootGAM ->",
               "Settings for the GAM ->")
  title = "\nBootGAM MENU\n  \\main\\covariate model\\BootGAM\n\n*** Note that the bootGAM also uses the settings from the GAM!\n    Please go the GAM settings menu to alter these.\n"
  pick <- menu(choices, title = title)
  qx <- 0
  switch(pick + 1, qx <- 2, qx <- 1,
         xp.bootgam (eval(as.name(".cur.db")), overwrite = FALSE),
         bootgam.print(),
         qx <- bootgam.plot.menu(),
         qx <- bootgam.settings.menu(),
         qx <- gam.settings.menu())
  if (qx == 2) {
    return(invisible(2))
  }
  else {
    if (qx == 1) {
      return(invisible(0))
    }
    else {
      Recall()
    }
  }
}

bootgam.plot.menu <- function () {
  choices <- c("Return to previous menu",
               "Inclusion frequencies",
               "Most common 2-covariate combinations",
               "Distribution of model size",
               "Inclusion stability - covariates",
               "Inclusion index of covariates",
               "Inclusion index of covariates/individuals",
               "Compare index of covariates/individuals"
               )
    title = "\nBootGAM plot MENU\n  \\main\\covariate model\\BootGAM\\Plot"
    pick <- menu(choices, title = title)
    qx <- 0
    switch(pick + 1, qx <- 2, qx <- 1,
           print(xp.inc.prob()),
           print(xp.inc.prob.comb.2()),
           print(xp.distr.mod.size()),
           print(xp.inc.stab.cov()),
           print(xp.incl.index.cov()),
           print(xp.incl.index.cov.ind()),
           print(xp.incl.index.cov.comp())
           )
    if (qx == 2) {
        return(invisible(2))
    }
    else {
        if (qx == 1) {
            return(invisible(0))
        }
        else {
            Recall()
        }
    }
}

bootscm.plot.menu <- function () {
  choices <- c("Return to previous menu ->",
               "Inclusion frequencies",
               "Most common 2-covariate combinations",
               "Distribution of model size",
               "Inclusion stability - covariates",
               "Inclusion index covariates",
               "Inclusion index of covariates/indidividuals",
               "Compare index of covariates/individuals",
               "Bias parameter estimates (hurricane plot) ",
               "Correlation in parameters covariate effects",
               "Difference in OFV final models (optimism plot)"
               )
  title = "\nBOOTSCM PLOT MENU\n  \\main\\covariate model\\BootSCM\\Plot menu\n\n"
  pick <- menu(choices, title = title)
  qx <- 0
  switch(pick + 1, qx <- 2, qx <- 1,
         print(xp.inc.prob()),
         print(xp.inc.prob.comb.2()),
         print(xp.distr.mod.size()),
         print(xp.inc.stab.cov()),
         print(xp.incl.index.cov()),
         print(xp.incl.index.cov.ind()),
         print(xp.incl.index.cov.comp()),
         print(xp.boot.par.est()),
         print(xp.boot.par.est.corr(ask.covs=TRUE)),
         print(xp.dofv.plot())
         )
  if (qx == 2) {
    return(invisible(2))
  }
  else {
    if (qx == 1) {
      return(invisible(0))
    }
    else {
      Recall()
    }
  }
}

bootscm.menu <- function () {
  choices <- c("Return to previous menu ->",
               "Import bootSCM data (from PsN folder)",
               "Summarize bootSCM",
               "Plot menu"
               )
  title = "\nBOOTSCM MENU\n  \\main\\covariate model\\BootSCM\n\nThe BootSCM is implemented in PsN. Xpose can import its output and\ngenerate plots similar to those for the BootGAM\n"
  pick <- menu(choices, title = title)
  qx <- 0
  switch(pick + 1, qx <- 2, qx <- 1,
         bootscm.import(),
         bootgam.print(),
         qx <- bootscm.plot.menu()
         )
  if (qx == 2) {
    return(invisible(2))
  }
  else {
    if (qx == 1) {
      return(invisible(0))
    }
    else {
      Recall()
    }
  }
}


