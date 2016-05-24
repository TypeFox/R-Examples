# Copyright (C) 2010-2015 Johannes Ranke {{{
# Contact: jranke@uni-bremen.de

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/> }}}

mkinmod <- function(..., use_of_ff = "min", speclist = NULL, quiet = FALSE, verbose = FALSE)
{
  if (is.null(speclist)) spec <- list(...)
  else spec <- speclist
  obs_vars <- names(spec)

  # Check if any of the names of the observed variables contains any other
  for (obs_var in obs_vars) {
    if (length(grep(obs_var, obs_vars)) > 1) stop("Sorry, variable names can not contain each other")
    if (grepl("_to_", obs_var)) stop("Sorry, names of observed variables can not contain _to_")
    if (obs_var == "sink") stop("Naming a compound 'sink' is not supported")
  }

  if (!use_of_ff %in% c("min", "max"))
    stop("The use of formation fractions 'use_of_ff' can only be 'min' or 'max'")

  # The returned model will be a list of character vectors, containing {{{
  # differential equations (if supported), parameter names and a mapping from
  # model variables to observed variables. If possible, a matrix representation
  # of the differential equations is included
  # Compiling the functions from the C code generated below only works if the
  # implicit assumption about differential equations specified below
  # is satisfied
  parms <- vector()
  # }}}

  # Do not return a coefficient matrix mat when FOMC, IORE, DFOP or HS is used for the parent {{{
  if(spec[[1]]$type %in% c("FOMC", "IORE", "DFOP", "HS")) {
    mat = FALSE 
  } else mat = TRUE
  #}}}

  # Establish a list of differential equations as well as a map from observed {{{
  # compartments to differential equations
  diffs <- vector()
  map <- list()
  for (varname in obs_vars)
  {
    # Check the type component of the compartment specification {{{
    if(is.null(spec[[varname]]$type)) stop(
      "Every part of the model specification must be a list containing a type component")
    if(!spec[[varname]]$type %in% c("SFO", "FOMC", "IORE", "DFOP", "HS", "SFORB")) stop(
      "Available types are SFO, FOMC, IORE, DFOP, HS and SFORB only")
    if(spec[[varname]]$type %in% c("FOMC", "DFOP", "HS") & match(varname, obs_vars) != 1) {
        stop(paste("Types FOMC, DFOP and HS are only implemented for the first compartment,", 
                   "which is assumed to be the source compartment"))
    }
    #}}}
    # New (sub)compartments (boxes) needed for the model type {{{
    new_boxes <- switch(spec[[varname]]$type,
      SFO = varname,
      FOMC = varname,
      IORE = varname,
      DFOP = varname,
      HS = varname,
      SFORB = paste(varname, c("free", "bound"), sep = "_")
    )
    map[[varname]] <- new_boxes
    names(map[[varname]]) <- rep(spec[[varname]]$type, length(new_boxes)) #}}}
    # Start a new differential equation for each new box {{{
    new_diffs <- paste("d_", new_boxes, " =", sep = "")
    names(new_diffs) <- new_boxes
    diffs <- c(diffs, new_diffs) #}}}
  } #}}}

  # Create content of differential equations and build parameter list {{{
  for (varname in obs_vars)
  {
    # Get the name of the box(es) we are working on for the decline term(s)
    box_1 = map[[varname]][[1]] # This is the only box unless type is SFORB
    # Turn on sink if this is not explicitly excluded by the user by
    # specifying sink=FALSE
    if(is.null(spec[[varname]]$sink)) spec[[varname]]$sink <- TRUE
    if(spec[[varname]]$type %in% c("SFO", "IORE", "SFORB")) { # {{{ Add decline term
      if (use_of_ff == "min") { # Minimum use of formation fractions
        if(spec[[varname]]$type == "IORE" && length(spec[[varname]]$to) > 0) {
           stop("Transformation reactions from compounds modelled with IORE\n",
                "are only supported with formation fractions (use_of_ff = 'max')")
        }
        if(spec[[varname]]$sink) {
          # If sink is required, add first-order/IORE sink term
          k_compound_sink <- paste("k", box_1, "sink", sep = "_")
          if(spec[[varname]]$type == "IORE") {
            k_compound_sink <- paste("k__iore", box_1, "sink", sep = "_")
          }
          parms <- c(parms, k_compound_sink)
          decline_term <- paste(k_compound_sink, "*", box_1)
          if(spec[[varname]]$type == "IORE") {
            N <- paste("N", box_1, sep = "_")
            parms <- c(parms, N)
            decline_term <- paste0(decline_term, "^", N)
          }
        } else { # otherwise no decline term needed here
          decline_term = "0" 
        }
      } else {
        k_compound <- paste("k", box_1, sep = "_")
        if(spec[[varname]]$type == "IORE") {
          k_compound <- paste("k__iore", box_1, sep = "_")
        }
        parms <- c(parms, k_compound)
        decline_term <- paste(k_compound, "*", box_1)
        if(spec[[varname]]$type == "IORE") {
          N <- paste("N", box_1, sep = "_")
          parms <- c(parms, N)
          decline_term <- paste0(decline_term, "^", N)
        }
      }
    } #}}}
    if(spec[[varname]]$type == "FOMC") { # {{{ Add FOMC decline term
      # From p. 53 of the FOCUS kinetics report, without the power function so it works in C
      decline_term <- paste("(alpha/beta) * 1/((time/beta) + 1) *", box_1)
      parms <- c(parms, "alpha", "beta")
    } #}}}
    if(spec[[varname]]$type == "DFOP") { # {{{ Add DFOP decline term
      # From p. 57 of the FOCUS kinetics report
      decline_term <- paste("((k1 * g * exp(-k1 * time) + k2 * (1 - g) * exp(-k2 * time)) / (g * exp(-k1 * time) + (1 - g) * exp(-k2 * time))) *", box_1)
      parms <- c(parms, "k1", "k2", "g")
    } #}}}
    HS_decline <- "ifelse(time <= tb, k1, k2)" # Used below for automatic translation to C
    if(spec[[varname]]$type == "HS") { # {{{ Add HS decline term
      # From p. 55 of the FOCUS kinetics report
      decline_term <- paste(HS_decline, "*", box_1)
      parms <- c(parms, "k1", "k2", "tb")
    } #}}}
    # Add origin decline term to box 1 (usually the only box, unless type is SFORB)#{{{
    diffs[[box_1]] <- paste(diffs[[box_1]], "-", decline_term)#}}}
    if(spec[[varname]]$type == "SFORB") { # {{{ Add SFORB reversible binding terms
      box_2 = map[[varname]][[2]]
      if (use_of_ff == "min") { # Minimum use of formation fractions
        k_free_bound <- paste("k", varname, "free", "bound", sep = "_")
        k_bound_free <- paste("k", varname, "bound", "free", sep = "_")
        parms <- c(parms, k_free_bound, k_bound_free)
        reversible_binding_term_1 <- paste("-", k_free_bound, "*", box_1, "+",
          k_bound_free, "*", box_2)
        reversible_binding_term_2 <- paste("+", k_free_bound, "*", box_1, "-",
          k_bound_free, "*", box_2)
      } else { # Use formation fractions also for the free compartment
        stop("The maximum use of formation fractions is not supported for SFORB models")
        # The problems were: Calculation of dissipation times did not work in this case
        # and the coefficient matrix is not generated correctly by the code present 
        # in this file in this case
        f_free_bound <- paste("f", varname, "free", "bound", sep = "_")
        k_bound_free <- paste("k", varname, "bound", "free", sep = "_")
        parms <- c(parms, f_free_bound, k_bound_free)
        reversible_binding_term_1 <- paste("+", k_bound_free, "*", box_2)
        reversible_binding_term_2 <- paste("+", f_free_bound, "*", k_compound, "*", box_1, "-",
          k_bound_free, "*", box_2)
      }
      diffs[[box_1]] <- paste(diffs[[box_1]], reversible_binding_term_1)
      diffs[[box_2]] <- paste(diffs[[box_2]], reversible_binding_term_2)
    } #}}}

    # Transfer between compartments#{{{
    to <- spec[[varname]]$to
    if(!is.null(to)) {
      # Name of box from which transfer takes place
      origin_box <- box_1

      # Number of targets
      n_targets = length(to)

      # Add transfer terms to listed compartments
      for (target in to) {
        if (!target %in% obs_vars) stop("You did not specify a submodel for target variable ", target)
        target_box <- switch(spec[[target]]$type,
          SFO = target,
          IORE = target,
          SFORB = paste(target, "free", sep = "_"))
        if (use_of_ff == "min" && spec[[varname]]$type %in% c("SFO", "SFORB"))
        {
          k_from_to <- paste("k", origin_box, target_box, sep = "_")
          parms <- c(parms, k_from_to)
          diffs[[origin_box]] <- paste(diffs[[origin_box]], "-", 
            k_from_to, "*", origin_box)
          diffs[[target_box]] <- paste(diffs[[target_box]], "+", 
            k_from_to, "*", origin_box)
        } else {
          # Do not introduce a formation fraction if this is the only target
          if (spec[[origin_box]]$sink == FALSE && n_targets == 1) {
            diffs[[target_box]] <- paste(diffs[[target_box]], "+",
                                         decline_term)
          } else {
            fraction_to_target = paste("f", origin_box, "to", target, sep = "_")
            parms <- c(parms, fraction_to_target)
            diffs[[target_box]] <- paste(diffs[[target_box]], "+", 
                fraction_to_target, "*", decline_term)
          }
        }
      }
    } #}}}
  } #}}}

  model <- list(diffs = diffs, parms = parms, map = map, spec = spec, use_of_ff = use_of_ff)

  # Create coefficient matrix if appropriate#{{{
  if (mat) {
    boxes <- names(diffs)
    n <- length(boxes)
    m <- matrix(nrow=n, ncol=n, dimnames=list(boxes, boxes))

    if (use_of_ff == "min") { # {{{ Minimum use of formation fractions
      for (from in boxes) {
        for (to in boxes) {
          if (from == to) { # diagonal elements
            k.candidate = paste("k", from, c(boxes, "sink"), sep = "_")
            k.candidate = sub("free.*bound", "free_bound", k.candidate)
            k.candidate = sub("bound.*free", "bound_free", k.candidate)
            k.effective = intersect(model$parms, k.candidate)
            m[from,to] = ifelse(length(k.effective) > 0,
                paste("-", k.effective, collapse = " "), "0")

          } else {          # off-diagonal elements
            k.candidate = paste("k", from, to, sep = "_")
	    if (sub("_free$", "", from) == sub("_bound$", "", to)) {
              k.candidate = paste("k", sub("_free$", "_free_bound", from), sep = "_")
	    }
	    if (sub("_bound$", "", from) == sub("_free$", "", to)) {
              k.candidate = paste("k", sub("_bound$", "_bound_free", from), sep = "_")
	    }
            k.effective = intersect(model$parms, k.candidate)
            m[to, from] = ifelse(length(k.effective) > 0,
                k.effective, "0")
          }
        }
      } # }}}
    } else { # {{{ Use formation fractions where possible
      for (from in boxes) {
        for (to in boxes) {
          if (from == to) { # diagonal elements
            k.candidate = paste("k", from, sep = "_")
            m[from,to] = ifelse(k.candidate %in% model$parms,
                paste("-", k.candidate), "0")
            if(grepl("_free", from)) { # add transfer to bound compartment for SFORB
              m[from,to] = paste(m[from,to], "-", paste("k", from, "bound", sep = "_"))
            }
            if(grepl("_bound", from)) { # add backtransfer to free compartment for SFORB
              m[from,to] = paste("- k", from, "free", sep = "_")
            }
            m[from,to] = m[from,to]
          } else {          # off-diagonal elements
            f.candidate = paste("f", from, "to", to, sep = "_")
            k.candidate = paste("k", from, to, sep = "_")
            # SFORB with maximum use of formation fractions not implemented, see above
            m[to, from] = ifelse(f.candidate %in% model$parms,
              paste(f.candidate, " * k_", from, sep = ""), 
              ifelse(k.candidate %in% model$parms, k.candidate, "0"))
            # Special case: singular pathway and no sink
            if (spec[[from]]$sink == FALSE && length(spec[[from]]$to) == 1 && to %in% spec[[from]]$to) {
              m[to, from] = paste("k", from, sep = "_")
            }
          }
        }
      }
    } # }}}
    model$coefmat <- m
  }#}}}

  # Try to create a function compiled from C code if more than one observed {{{
  # variable and gcc is available
  if (length(obs_vars) > 1) {
    if (Sys.which("gcc") != "") {

      # Translate the R code for the derivatives to C code
      diffs.C <- paste(diffs, collapse = ";\n")
      diffs.C <- paste0(diffs.C, ";")

      # HS
      diffs.C <- gsub(HS_decline, "(time <= tb ? k1 : k2)", diffs.C, fixed = TRUE)

      for (i in seq_along(diffs)) {
        state_var <- names(diffs)[i]

        # IORE
        if (state_var %in% obs_vars) {
          if (spec[[state_var]]$type == "IORE") {
            diffs.C <- gsub(paste0(state_var, "^N_", state_var),
                            paste0("pow(y[", i - 1, "], N_", state_var, ")"),
                            diffs.C, fixed = TRUE)
          }
        }

        # Replace d_... terms by f[i-1]
        # First line
        pattern <- paste0("^d_", state_var)
        replacement <- paste0("\nf[", i - 1, "]")
        diffs.C <- gsub(pattern, replacement, diffs.C)
        # Other lines
        pattern <- paste0("\\nd_", state_var)
        replacement <- paste0("\nf[", i - 1, "]")
        diffs.C <- gsub(pattern, replacement, diffs.C)

        # Replace names of observed variables by y[i],
        # making the implicit assumption that the observed variables only occur after "* "
        pattern <- paste0("\\* ", state_var)
        replacement <- paste0("* y[", i - 1, "]")
        diffs.C <- gsub(pattern, replacement, diffs.C)
      }

      derivs_sig <- signature(n = "integer", t = "numeric", y = "numeric",
                              f = "numeric", rpar = "numeric", ipar = "integer")
      
      # Declare the time variable in the body of the function if it is used
      derivs_code <- if (spec[[1]]$type %in% c("FOMC", "DFOP", "HS")) {
        paste0("double time = *t;\n", diffs.C)
      } else {
        diffs.C
      }

      # Define the function initializing the parameters
      npar <- length(parms)
      initpar_code <- paste0(
        "static double parms [", npar, "];\n",
        paste0("#define ", parms, " parms[", 0:(npar - 1), "]\n", collapse = ""),
        "\n",
        "void initpar(void (* odeparms)(int *, double *)) {\n",
        "    int N = ", npar, ";\n",
        "    odeparms(&N, parms);\n",
        "}\n\n")

      # Try to build a shared library
      cf <- try(cfunction(list(func = derivs_sig), derivs_code,
                          otherdefs = initpar_code,
                          verbose = verbose,
                          convention = ".C", language = "C"),
                silent = TRUE)

      if (!inherits(cf, "try-error")) {
        if (!quiet) message("Successfully compiled differential equation model from auto-generated C code.")
        model$cf <- cf
      }
    }
  }
  # }}}

  class(model) <- "mkinmod"
  return(model)
}

print.mkinmod <- function(x, ...) {
  cat("<mkinmod> model generated with\n")
  cat("Use of formation fractions $use_of_ff:", x$use_of_ff, "\n")
  cat("Specification $spec:\n")
  for (obs in names(x$spec)) {
    cat("$", obs, "\n", sep = "")
    spl <- x$spec[[obs]]
    cat("$type:", spl$type)
    if (!is.null(spl$to) && length(spl$to)) cat("; $to: ", paste(spl$to, collapse = ", "), sep = "")
    cat("; $sink: ", spl$sink, sep = "")
    if (!is.null(spl$full_name)) if (!is.na(spl$full_name)) cat("; $full_name:", spl$full_name)
    cat("\n")
  }
  if (is.matrix(x$coefmat)) cat("Coefficient matrix $coefmat available\n")
  if (!is.null(x$cf)) cat("Compiled model $cf available\n")  
}
# vim: set foldmethod=marker ts=2 sw=2 expandtab:
