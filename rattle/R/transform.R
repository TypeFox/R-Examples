# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2016-01-07 08:12:50 gjw>
#
# TRANSFORM TAB
#
# Copyright (c) 2009-2011 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

########################################################################
# CALLBACKS

# When a radio button on the TRANSFORM tab is selected, display the
# appropriate option.

on_impute_radiobutton_toggled <- function(button)
{
  if (button$getActive()) 
  {
    crv$TRANSFORM$setCurrentPage(crv$TRANSFORM.IMPUTE.TAB)
  }
  setStatusBar()
}

# TODO 080423 Change to RESCALE
on_normalise_radiobutton_toggled <- function(button)
{
  if (button$getActive()) 
  {
    crv$TRANSFORM$setCurrentPage(crv$TRANSFORM.NORMALISE.TAB)
  }
  setStatusBar()
}

on_normalise_interval_radiobutton_toggled <- function(button)
{
  active <- button$getActive()
  theWidget("normalise_interval_numgroups_label")$setSensitive(active)
  theWidget("normalise_interval_numgroups_spinbutton")$setSensitive(active)
}

on_remap_radiobutton_toggled <- function(button)
{
  if (button$getActive()) 
  {
    crv$TRANSFORM$setCurrentPage(crv$TRANSFORM.REMAP.TAB)
  }
  setStatusBar()
}

on_cleanup_radiobutton_toggled <- function(button)
{
  if (button$getActive()) 
  {
    crv$TRANSFORM$setCurrentPage(crv$TRANSFORM.CLEANUP.TAB)
  }
  setStatusBar()
}

on_impute_constant_radiobutton_toggled <- function(button)
{
  theWidget("impute_constant_entry")$setSensitive(button$getActive())
}

########################################################################
# CAPTURE A RECORD OF TRANSFORMS
#
# 090605 Record a list of transforms. Treat as an object. Previously
# we simply had a list of strings and recorded the parameters of the
# transform within the string. This was not scalable. Now we have a
# list of lists. Each element in the list is named with the name of
# the new variable. The sub list then records at least the source
# variable name and the type of the transform. Other transform
# parameters that completely define how the transform is calculated
# are also recorded.
#
# RRC_AGE
#	orig	Age
#	type	RRC
#	mean	38.622
#	sd	13.58475

union.transform <- function(tr, nm, lst)
{
  # 090605 Add the new transform, resulting in a variable named NM
  # with paramaters LST (which has a minimum of source and type), into
  # the list TR, optionally overwriting any previous transform of the
  # same type. The list of transforms, TR, is returned, with the new
  # transform appedned.

  # Remove any transforms of the same variable and type. Note that the
  # names already embody the variable name and the transform type.
  
  tr[which(nm == names(tr))] <- NULL

  # Add the new transform, and give it the name of the variable being
  # transformed.
  
  tr <- c(tr, list(lst))
  names(tr)[length(tr)] <- nm

  # Return the augmented list of transforms.
  
  return(tr)
}

########################################################################
# UTILITIES

modalvalue <- function(x, na.rm=FALSE)
{
  # Determine the modal value of the data x.
  
  x = unlist(x)
  if(na.rm) x = x[!is.na(x)]
  u = unique(x)
  n = length(u)
  frequencies = rep(0, n)
  for(i in seq_len(n))
  {
    if(is.na(u[i]))
    {
      frequencies[i] = sum(is.na(x))
    } else
    {
      frequencies[i] = sum(x==u[i], na.rm=TRUE)
    }
  }
  u[which.max(frequencies)]
}

# RESCALE BY GROUP
#
# 110529 Take a vector of numbers and rescale it appropriately. The
# rescaling is either done in the context of the population of
# observations when by is NULL, or grouped into subpopulations as
# defined by the by argument. Note that in the transform code there is
# not much to be gained by calling this function when there is no by
# argument, and in fact it is probably more educational to show the
# direct call to the other commands in the log tab.
#
# TYPE == irank
#
# Rescale the data to integers from 0 to itop-1, with any missing
# values mapped to the midpoint. Original idea from Tony Nolan.
#
# TYPE == recenter
#
# The usual z-score rescaling, subtracting the mean and dividing by
# the standard deviation.
#
# We plan to generalise this function to do any kind of
# rescaling by group, and then get Rattle to simply call this function
# to do its work.
#
# Examples:
#
# rescale.by.group(crs$dataset[["Age"]], crs$dataset[["Gender"]])

rescale.by.group <- function(x, by=NULL, type="irank", itop=100)
{
  # 110529 TODO Check that by is a factor. 160107 Make sure it is a
  # factor... Any consequences?

  bylevels <- levels(as.factor(by))

  # Initialise the result.

  y <- rep(0, length(x))

  # Choose the operation

  cmd <- switch(type,
                irank='round(rescaler(x[elts], "range") * (itop-1))',
                recenter='scale(x[elts])[,1]',
                range='rescaler(x[elts], "range")',
                rank='rescaler(x[elts], "rank")',
                robust='rescaler(x[elts], "robust")')
  
  if (is.null(by))
    y <- eval(parse(text=sub('\\[elts\\]', '', cmd)))
  else
    for (vl in bylevels)
    {
      elts <- sapply(by==vl, isTRUE)
      y[elts] <- eval(parse(text=cmd))
    }

  # 101007 We used to map to itop-1 but now itop/2. This case is when
  # there is only a single observation in a group (or when the
  # original data is missing), and itop/2 makes more sense than either
  # 0 or itop-1.
  
  if (type == "irank") y[is.nan(y)] <- round(itop/2)

  return(y)
}

########################################################################
# DISPATCH

executeTransformTab <- function()
{
  # We cannot do any transforms if there is no dataset.

  if (noDatasetLoaded()) return()

  # Dispatch to the appropriate option.

  # TODO 080423 Change NORMALISE to RESCALE
  
  if (theWidget("normalise_radiobutton")$getActive())
    executeTransformNormalisePerform()
  else if (theWidget("impute_radiobutton")$getActive())
    executeTransformImputePerform()
  else if (theWidget("remap_radiobutton")$getActive())
    executeTransformRemapPerform()
  else if (theWidget("cleanup_radiobutton")$getActive())
    executeTransformCleanupPerform()
}

########################################################################
# RESCALE

executeTransformNormalisePerform <- function(variables=NULL,
                                             action=NULL,
                                             vprefix=action)
{
  # TODO 080609 We should rename this in line with the interface change,
  # since it is not necessarily normalisation but is rescaling.

  # 090323 The paramaters were added to the function call so that we
  # can override what the interface says. This is primarily useful for
  # my Sweave work where I'm automatically controlling Rattle to get
  # screensshots.
  
  # First determine which normalisation option has been chosen and the
  # prefix of the new variable that will be introduced.  Default to
  # NULL in the hope of picking up an error if something has gone wrong.

  # TODO 071124 The radio buttons could be checkbuttons, and we do
  # multiple transformations for the selected variables, but for now,
  # stay with radio buttons as it is simple, without loss of
  # functionality.

  # 110226 Allow numeric and categoric for recenter, scale01,
  # medianad, and rank. ByGroup variables will be created, as
  # suggested by Tony. This will then replace the By Group optoin
  # altogether.

  bygroup <- FALSE
  byvname <- NULL

  # If the action is not passed through by the command line then
  # determine what action is required from the GUI.
  
  if (is.null(action))
  {
    if (theWidget("normalise_recenter_radiobutton")$getActive())
    {
      action <- "recenter"
      vprefix <- "RRC"
    }
    else if (theWidget("normalise_scale01_radiobutton")$getActive())
    {
      action <- "scale01"
      vprefix <- "R01"
    }
    else if (theWidget("normalise_rank_radiobutton")$getActive())
    {
      action <- "rank"
      vprefix <- "RRK"
    }
    else if (theWidget("normalise_medianad_radiobutton")$getActive())
    {
      action <- "medianad"
      vprefix <- "RMD"
    }
    else if (theWidget("normalise_interval_radiobutton")$getActive())
    {
      action <- "interval"
      vprefix <- "RIN" # 110530 Changed from "RBG" for "bygroup" to "interval"
    }
    else if (theWidget("rescale_matrix_radiobutton")$getActive())
    {
      action <- "matrix"
      vprefix <- "RMA"
    }
    else if (theWidget("rescale_log_radiobutton")$getActive())
    {
      action <- "log"
      vprefix <- "RLG"
      #remap.comment <- "Log transform."
    }
    else if (theWidget("rescale_log10_radiobutton")$getActive())
    {
      action <- "log10"
      vprefix <- "R10"
    }
  }

  # 110530 DEPRECATED The "bygroup" action is now "interval". This has
  # been renamed globally in Rattle but Sweave documents may still
  # be using "bygroup", so for now map "bygroup" to "interval" and
  # issue a warning.

  if (action == "bygroup")
  {
    warnDialog("DEPRECATED transform action of 'bygroup'. Use 'interval' instead.")
    action <- "interval"
  }
  
  # Obtain the list of selected variables from the treeview.

  if (is.null(variables))
  {
    selected <- theWidget("impute_treeview")$getSelection()
    selected$selectedForeach(function(model, path, iter, data)
    {
      variables <<- c(variables, model$get(iter, 1)[[1]])
    }, TRUE)
    if (length(variables)) Encoding(variables) <- "UTF-8"
  }

  # 110226 If no variables are selected then there is nothing to do.
  
  if (!length(variables))
  {
    warnDialog(Rtxt("No variables have been selected for rescaling.",
                    "Please select some variables and Execute again."))
    setStatusBar(Rtxt("No variables selected to be rescaling."))
    return(FALSE)
  }

  ## GROUP BY 110530 If a categoric is also selected, and the action
  ## allows it (i.e., not "matrix" or "log") the set up things for
  ## groupby to work. We can only handle a single categoric for the
  ## groupby so put up an error dialogue and return FALSE. We used to
  ## simply remove the categorics from the list of variables to be
  ## rescaled and continue with that but now that we allow groupby it
  ## makes more sense to do nothing if not sure what to do.

  # Identify the class of each variable.  080328 For any ordered
  # factors class returns two values (since the object inherits first
  # from ordered and then factor), so we remove the "ordered" from the
  # list to hopefully get back to the right length for classes (i.e.,
  # one class for each variable). NOTE objects can inherit from
  # multiple classes, and the order presented is the order in which
  # they inherit. TODO Maybe before we unlist we need to turn multiple
  # results into one, like "ordered_factor".

  # TODO Allow multiple categorics and then group across all the
  # cateogircs: MaleMarried MaleDivorced FemaleMarried etc.
  
  classes <- unlist(lapply(variables, function(x) class(crs$dataset[[x]])))
  classes <- classes[classes!="ordered"]

  # If there are multiple factors amongst the variables then for now bail out.

  nfactors <- sum("factor" == classes)

  # Currently, only support grouping by a single categoric. TODO
  # Support a group by of multiple categorics.
      
  if (nfactors > 1)
  {
    infoDialog(Rtxt("We only support By Group with a single categoric",
                    "variable for now.\n\nPlease select just one",
                    "categoric."))
    return(FALSE)
  }
  else if (nfactors == 1)
  {
    if (action %in% c("matrix", "log", "log10"))
    {
      # 110220 Choosing a categoric for one of these does not make
      # sense. So remove it. 100428 BUG When using Rtxt in sprintf,
      # and substituting a UTF-8 encoded variable we get garbage, so
      # convert to unknown then back again. Making the Rtxt UTF-8 does
      # not fix it.

      Encoding(variables) <- "unknown"
      infoDialog(sprintf(Rtxt("We cannot rescale using '%s'",
                              "on a categoric variable.",
                              "Ignoring: %s."),
                         action, paste(variables[which("factor" == classes)],
                                       collapse=", ")))
      Encoding(variables) <- "UTF-8"
      variables <- variables[-which("factor" == classes)] # Remove the factors.
      if (length(variables) == 0) return(FALSE)
    }
    else
    {
      # 110220 If a categoric is selected for one of the other
      # transforms, then turn the operation into a bygroup
      # transform.

      bygroup <- TRUE
      numnumerics <- sum(classes=="numeric" | classes=="integer")

      # The prefix for the new variable will be BYC for recenter, BG1
      # for scale01, BGK for rank, BGD for medianad, and BGN for
      # interval.
     
      vprefix <- paste("BG", substr(x=vprefix, start=3, stop=3), sep="")

      # Remove the categoric (if any) from the list of variables and
      # store its name in byvname. This allows us to continue to use
      # the loop below, having just the numeric variables in the list.
      
      byvname <- variables[which("factor" == classes)]
      variables <- variables[-which("factor" == classes)]

      # Ensure we have at least one numeric variable.
    
      if (numnumerics == 0)
      {
        infoDialog(Rtxt("We must have a numeric variable to normalize for the",
                        "By Group option.\n\nPlease select one numeric variable."))
        return(FALSE)
      }
    }
  }
  
  startLog(Rtxt("Transform variables by rescaling."))
  
  # Make sure we have the reshape library from where the rescaler
  # function comes.
  
  if (action %in% c("scale01", "rank", "medianad", "interval"))
  {
    if (! packageIsAvailable("reshape", Rtxt("normalize data"))) return()
    lib.cmd <- "library(reshape, quietly=TRUE)"
    appendLog(packageProvides("reshape", "rescaler"), lib.cmd)
    eval(parse(text=lib.cmd))
  }

  # Record the current variable roles so that we can maintain these,
  # modified appropriately by ignore'ing the imputed variables, and
  # input'ing the newly imputed variables.
  
  input <- getSelectedVariables("input")
  target <- getSelectedVariables("target")
  risk <- getSelectedVariables("risk")
  ident <- getSelectedVariables("ident")
  ignore <- getSelectedVariables("ignore")
  weight <- getSelectedVariables("weight")

  # For MATRIX obtain the matrix total first and then divide each
  # variable by this. We do this here so that it is done only once,
  # rather than through each loop below.

  if (action == "matrix")
  {
    matrix.total <- 0
    total.cmd <- sprintf(paste("matrix.total <- sum(crs$dataset[,",
                               'c("%s")],',
                               "na.rm=TRUE)"),
                         paste(variables, collapse='", "'))
    appendLog(Rtxt("Calculate the matrix total."), total.cmd)
    eval(parse(text=total.cmd))
  }

  # Loop through each of the supplied variables.
    
  for (v in variables)
  {
    norm.score.command <- NULL
    norm.score.comment <- NULL
    
    # Create the new variable name.

    if (action %in% c("interval"))
    {
      num.groups <- theWidget("normalise_interval_numgroups_spinbutton")$getValue()
      if (bygroup)
        vname <- paste(vprefix, byvname, v, num.groups, sep="_")
      else
        vname <- paste(vprefix, v, num.groups, sep="_")
    }
    else
    {
      if (bygroup)
        vname <- paste(vprefix, byvname, v, sep="_")
      else
        vname <- paste(vprefix, v, sep="_")
    }
    
    # Check variable specific preconditions, and if we fail then
    # proceed to next variable.

    if (action == "medianad")
    {
      # 080609 For audit$Deductions this returns all NaN or Inf
      # because the median is 0. We can see this with
      # rescaler(crs$dataset[["Deductions"]], "robust") So check for
      # this and do nothing!

      median.cmd <- sprintf('median(crs$dataset[["%s"]], na.rm=TRUE)', v)
      if (eval(parse(text=median.cmd)) == 0)
      {
        # 100428 BUG When using Rtxt in sprintf, and substituting a
        # UTF-8 encoded variable we get garbage, so convert to unknown
        # then back again. Making the Rtxt UTF-8 does not fix it.
        Encoding(v) <- "unknown"
        warnDialog(sprintf(Rtxt("The variable '%s' has a median of 0.",
                                "We cannot compute the Median/MAD Rescaler",
                                "for this variable."), v))
        Encoding(v) <- "UTF-8"
        next()
      }
    }
         
    # Generate the command to copy the current variable into a new
    # variable, prefixed appropraitely.

    copy.cmd <- sprintf('crs$dataset[["%s"]] <- crs$dataset[["%s"]]',
                        vname, v)
    cl <- class(crs$dataset[[v]])

    # Take a copy of the variable to be imputed.

    # 100428 BUG When using Rtxt in sprintf, and substituting a UTF-8
    # encoded variable we get garbage, so convert to unknown then back
    # again. Making the Rtxt UTF-8 does not fix it.
    Encoding(v) <- "unknown"
    appendLog(sprintf(Rtxt("Rescale %s."), v), copy.cmd)
    Encoding(v) <- "UTF-8"
    eval(parse(text=copy.cmd))
    
    # Determine what action to perform.
    
    if (action == "recenter")
    {
      if (bygroup)
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n',
                                  '   rescale.by.group(crs$dataset[["%s"]],',
                                  'crs$dataset[["%s"]],',
                                  'type="recenter")'),
                            vname, v, byvname)
      else
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n',
                                  '   scale(crs$dataset[["%s"]])[,1]'),
                            vname, v)
      
      norm.comment <- Rtxt("Recenter and rescale the data around 0.")

      # Record the transformation for inclusion in PMML.

      # 090605 New transforms data structure
      
      lst <- list(orig=v, type=vprefix, status="active",
                  mean=mean(crs$dataset[[vname]], na.rm=TRUE),
                  sd=sd(crs$dataset[[vname]], na.rm=TRUE))
      crs$transforms <- union.transform(crs$transforms, vname, lst)

      # For the log, record the command to use when scoring the data.
      
      norm.score.command <- sprintf(paste('crs$dataset[["%s"]] <-',
                                          '(crs$dataset[["%s"]] -',
                                          '%f)/%f'),
                                    vname, v,
                                    mean(crs$dataset[[vname]], na.rm=TRUE),
                                    sd(crs$dataset[[vname]], na.rm=TRUE))
    }
    else if (action == "scale01")
    {
      if (bygroup)
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n',
                                  '   rescale.by.group(crs$dataset[["%s"]],',
                                  'crs$dataset[["%s"]],',
                                  'type="range")'),
                            vname, v, byvname)
      else
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <- ',
                                  'rescaler(crs$dataset[["%s"]], "range")'),
                            vname, v)

      norm.comment <- Rtxt("Rescale to [0,1].")

      # Record the transformation for inclusion in PMML.

      # 090606 New transforms data structure
      
      lst <- list(orig=v, type=vprefix, status="active",
                  min=min(crs$dataset[[vname]], na.rm=TRUE),
                  max=max(crs$dataset[[vname]], na.rm=TRUE))
      crs$transforms <- union.transform(crs$transforms, vname, lst)

      # For the log, record the command to use when scoring the data.
      
      norm.score.command <- sprintf(paste('crs$dataset[["%s"]] <-',
                                          '(crs$dataset[["%s"]] -',
                                          '%f)/abs(%f - %f)'),
                                    vname, v,
                                    min(crs$dataset[[vname]], na.rm=TRUE),
                                    max(crs$dataset[[vname]], na.rm=TRUE),
                                    min(crs$dataset[[vname]], na.rm=TRUE))
    }
    else if (action == "rank")
    {
      if (bygroup)
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n',
                                  '   rescale.by.group(crs$dataset[["%s"]],',
                                  'crs$dataset[["%s"]],',
                                  'type="rank")'),
                            vname, v, byvname)
      else
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <- ',
                                  'rescaler(crs$dataset[["%s"]], "rank")'),
                            vname, v)
      
      norm.comment <- Rtxt("Convert values to ranks.")

      # How would we rank a new item? Thus, can we actually use a rank
      # in a transform?

      # 090606 Record the transformation for inclusion in PMML.

      lst <- list(orig=v, type=vprefix, status="active")
      crs$transforms <- union.transform(crs$transforms, vname, lst)

    }
    else if (action == "medianad")
    {
      if (bygroup)
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n',
                                  '   rescale.by.group(crs$dataset[["%s"]],',
                                  'crs$dataset[["%s"]],',
                                  'type="robust")'),
                            vname, v, byvname)
      else
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <- ',
                                  'rescaler(crs$dataset[["%s"]], "robust")'),
                            vname, v)
      
      norm.comment <- Rtxt("Rescale by subtracting median and dividing",
                           "by median abs deviation.")

      # Record the transformation for inclusion in PMML.
      
      # 090606 New transforms data structure
      
      lst <- list(orig=v, type=vprefix, status="active",
                   median=median(crs$dataset[[vname]], na.rm=TRUE),
                   mad=mad(crs$dataset[[vname]], na.rm=TRUE))
      crs$transforms <- union.transform(crs$transforms, vname, lst)

      # For the log, record the command to use when scoring the data.
      
      norm.score.command <- sprintf(paste('crs$dataset[["%s"]] <-',
                                          '(crs$dataset[["%s"]] -',
                                          '%f)/%f'),
                                    vname, v,
                                    median(crs$dataset[[vname]], na.rm=TRUE),
                                    mad(crs$dataset[[vname]], na.rm=TRUE))
    }
    else if (action == "interval")
    {
      # v <- current numeric variable name from variables
      # byvname <- categoric variable name (no longer in variables)
      # vname <-  the new variable name set up as above

      num.groups <- theWidget("normalise_interval_numgroups_spinbutton")$getValue()

      # 110529 Considerably simplified by making use of the new
      # rescale.by.group function.
      
      if (bygroup)
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n    rescale.by.group(',
                                  'crs$dataset[["%s"]], crs$dataset[["%s"]],\n',
                                  '                     type="irank", itop=%s)', sep=""),
                            vname, v, byvname, num.groups)
      else
        norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <-\n    rescale.by.group(',
                                  'crs$dataset[["%s"]], ',
                                  'type="irank", itop=%s)',
                                  sep=""),
                            vname, v, num.groups)
      
      norm.comment <- sprintf(Rtxt("Rescale to 0-%s within each group."),
                              num.groups-1)

      # 090606 Record the transformation for inclusion in PMML.
      
      lst <- list(orig=v, type=vprefix, status="active", group=byvname)
      crs$transforms <- union.transform(crs$transforms, vname, lst)

    }
    else if (action == "matrix")
    {
      norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <- ',
                                'crs$dataset[["%s"]]/matrix.total'),
                          vname, v)
      norm.comment <- Rtxt("Divide variable values by matrix total.")

      # 090606 Is this still needed with the new data structure?
      # 090117 Remove any old instances of the same transform. Could
      # this be a general test outside the loop?

      ## present <- grep(vname, crs$transforms)
      ## if (length(present) >0) crs$transforms <- crs$transforms[-present]
      
      # 090117 Record the transformation for inclusion in PMML. Note
      # that we only need matrix.total, but all other members of
      # .TRANSFORMS.NORM.CONTINUOUS have two paramters, so include the
      # 1 to keep the group consistent.

      # 090605 New transforms data structure
      
      lst <- list(orig=v, type=vprefix, status="active",
                  sum=matrix.total, vars=variables)
      crs$transforms <- union.transform(crs$transforms, vname, lst)

      # For the log, record the command to use when scoring the data.
      
      norm.score.command <- sprintf(paste('crs$dataset[["%s"]] <-',
                                          'crs$dataset[["%s"]]/%f'),
                                    vname, v, matrix.total)
    }
    else if (action == "log")
    {
      norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <- ',
                                'log(crs$dataset[["%s"]])',
                                '\n  crs$dataset[crs$dataset[["%s"]] == -Inf &',
                                '! is.na(crs$dataset[["%s"]]), "%s"] <- NA'),
                          vname, v, vname, vname, vname)
      norm.comment <- Rtxt("Take a log transform of the variable - treat -Inf as NA.")

      # Record the transformation for inclusion in PMML.

      lst <- list(orig=v, type=vprefix, status="active")
      crs$transforms <- union.transform(crs$transforms, vname, lst)

      # For the log, record the command to use when scoring the data.
      
      norm.score.command <- norm.cmd
    }
    else if (action == "log10")
    {
      norm.cmd <- sprintf(paste('crs$dataset[["%s"]] <- ',
                                'log10(crs$dataset[["%s"]])',
                                '\n  crs$dataset[crs$dataset[["%s"]] == -Inf &',
                                '! is.na(crs$dataset[["%s"]]), "%s"] <- NA'),
                          vname, v, vname, vname, vname)
      norm.comment <- Rtxt("Take a log10 transform of the variable - treat -Inf as NA.")

      # Record the transformation for inclusion in PMML.

      lst <- list(orig=v, type=vprefix, status="active")
      crs$transforms <- union.transform(crs$transforms, vname, lst)

      # For the log, record the command to use when scoring the data.
      
      norm.score.command <- norm.cmd
    }

    appendLog(norm.comment,
              "if (building)\n{\n  ",
              norm.cmd,
              "\n}")
    eval(parse(text=norm.cmd))
    if (! is.null(norm.score.command))
      appendLog(Rtxt("When scoring transform using the training data parameters."),
                "if (scoring)\n{\n  ", norm.score.command, "\n}")
    
    # Now update the variable roles.
    
    if (v %in% input)
    {
      input <- setdiff(input, v)
      input <- union(input, vname)
    }
    else if (v %in% target)
    {
      target <- setdiff(target, v)
      target <- union(target, vname)
    }
    else if (v %in% risk)
    {
      risk <- setdiff(risk, v)
      risk <- union(risk, vname)
    }
    else if (v %in% ident)
    {
      ident <- setdiff(ident, v)
      ident <- union(ident, vname)
    }
    else
    {
      # If the source variable was ignore, then leave it as such, and
      # put the new variable in as input.

      input <- union(input, vname)
    }
    ignore <- union(ignore, v)
  }
  
  if (length(variables) > 0)
  {
    
    # Reset the dataset views keeping the roles unchanged except for
    # those that have been normalised, which have just been added as
    # inputs, with the originals now ignored.

    resetDatasetViews(input, target, risk, ident, ignore, weight)
    resetTestTab()
    
    # Update the status bar

    setStatusBar(sprintf(Rtxt("Normalized variables added to the dataset",
                              "with '%s' prefix."), vprefix))
  }
  else
  {
    warnDialog(Rtxt("No variables have been selected for rescaling.",
                    "Please select some variables and Execute again."))
    setStatusBar(Rtxt("No variables selected to be rescaling."))
  }
}  

########################################################################
# IMPUTE

executeTransformImputePerform <- function()
{
  # First determine which imputation option has been chosen and the
  # prefix of the new variable that will be introduced.  Default to
  # NULL so that if the value is not changed, we may get error (it
  # should be an error if the value is not changed).

  # TODO 071124 The rdaio buttons could be checkbuttons, and we do
  # multiple imputations for the selected variables, but for now, stay
  # with radio buttons as it is simply, without loss of functionality.
  
  action <- NULL
  vprefix <- NULL
  if (theWidget("impute_zero_radiobutton")$getActive())
  {
    action <- "zero"
    vprefix <- "IZR" # May want to distinguish ZERO and MISSING
  }
  else if (theWidget("impute_mean_radiobutton")$getActive())
  {
    action <- "mean"
    vprefix <- "IMN"
  }
  else if (theWidget("impute_median_radiobutton")$getActive())
  {
    action <- "median"
    vprefix <- "IMD"
  }
  else if (theWidget("impute_mode_radiobutton")$getActive())
  {
    action <- "mode"
    vprefix <- "IMO"
  }
  else if (theWidget("impute_constant_radiobutton")$getActive())
  {
    action <- "constant"
    vprefix <- "ICN"
  }
  
  # Obtain the list of selected variables from the treeview.

  imputed <- NULL
  selected <- theWidget("impute_treeview")$getSelection()
  selected$selectedForeach(function(model, path, iter, data)
  {
    imputed <<- c(imputed, model$get(iter, 1)[[1]])
  }, TRUE)
  if (length(imputed)) Encoding(imputed) <- "UTF-8"


  if (is.null(imputed)) 
    warnDialog(Rtxt("No variables have been selected for imputation.",
                    "Please select some variables and Execute again."))

  # We check here if the action is mean or median, and we have any
  # categoric variables to be imputed. If so put up an info dialogue
  # and remove the cateorigcals from the list of variables to be
  # imputed. We can't impute a mean or median for a categoric
  # variable.

  classes <- unlist(lapply(imputed, function(x) class(crs$dataset[[x]])))

  if (action %in% c("mean", "median") && "factor" %in% classes)
  {
    # 100429 BUG When using Rtxt in sprintf, and substituting a UTF-8
    # encoded variable we get garbage, so convert to unknown then back
    # again. Making the Rtxt UTF-8 does not fix it.
    Encoding(imputed) <- "unknown"
    infoDialog(sprintf(Rtxt("We cannot impute the %s for a",
                            "categoric variable. Ignoring: %s."),
                       action, paste(imputed[which("factor" %in% classes)],
                                     collapse=", ")))
    Encoding(imputed) <- "UTF-8"
    imputed <- imputed[-which("factor" %in% classes)] # Remove the factors.
  }
  
  # Record the current variable roles so that we can maintain these,
  # modified appropriately by ignore'ing the imputed variables, and
  # input'ing the newly imputed variables.
  
  input <- getSelectedVariables("input")
  target <- getSelectedVariables("target")
  risk <- getSelectedVariables("risk")
  ident <- getSelectedVariables("ident")
  ignore <- getSelectedVariables("ignore")
  weight <- getSelectedVariables("weight")

  if (length(imputed) > 0) startLog(Rtxt("Perform missing value imputation."))

  # [TODO 071124] The following code could be tidied up quite a
  # bit. It has evolved. Bits of the code handling the categorics
  # were copied from the numeric parts and vice versa, and they do it
  # different ways. Should try to do it the same way. Works for now!
      
  startLog(Rtxt("Transform variables by imputing missing values."))
  
  for (z in imputed)
  {
    # Generate the command to copy the current variable into a new
    # variable, prefixed appropraitely.
    
    vname <- paste(vprefix, z, sep="_")
    
    copy.cmd <- sprintf('crs$dataset[["%s"]] <- crs$dataset[["%s"]]',
                            vname, z)
    cl <- class(crs$dataset[[z]])
    
    # 110313 Note that cl could be "ordered" "factor". We have not
    # considered whether handling an "ordered" is different to
    # handling a "factor", in terms of adding "Missing" as another
    # level.
    
    if ("factor" %in% cl)
    {
      # Mean and median are not supported for categorics!

      if (action == "zero")
      {

        # Take a copy of the variable to be imputed.

        # 100429 BUG When using Rtxt in sprintf, and substituting a
        # UTF-8 encoded variable we get garbage, so convert to unknown
        # then back again. Making the Rtxt UTF-8 does not fix it.
        Encoding(z) <- "unknown"
        appendLog(sprintf(Rtxt("Impute %s."), z), copy.cmd)
        Encoding(z) <- "UTF-8"
        eval(parse(text=copy.cmd))
             
        # If "Missing" is not currently a category for this variable,
        # add it in.

        if ("Missing" %notin% levels(crs$dataset[[vname]]))
        {
          levels.cmd <- sprintf(paste('levels(crs$dataset[["%s"]]) <-',
                                      'c(levels(crs$dataset[["%s"]]),',
                                      '"Missing")'),
                                vname, vname)
          appendLog(Rtxt("Add a new category 'Missing' to the variable"),
                    levels.cmd)
          eval(parse(text=levels.cmd))
        }
      
        # Change all NAs to Missing.
      
        missing.cmd <- sprintf(paste('crs$dataset[["%s"]][is.na(',
                                     'crs$dataset[["%s"]])] <- "Missing"',
                                     sep=""),
                               vname, z)
        appendLog(Rtxt("Change all NAs to 'Missing'"), missing.cmd)
        eval(parse(text=missing.cmd))

        # 090630 I seemed to have only recorded these transforms for
        # non-categorics so now record the transformation for
        # inclusion in PMML.

        lst <- list(orig=z, type=vprefix, status="active", impute="Missing")
        crs$transforms <- union.transform(crs$transforms, vname, lst)
      }
      else if (action == "mode")
      {
        # Take a copy of the variable to be imputed.
    
        appendLog(sprintf("IMPUTE %s.", z), copy.cmd)
        eval(parse(text=copy.cmd))
             
        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 ' <- modalvalue(crs$dataset[["%s"]], ',
                                 "na.rm=TRUE)", sep=""), vname, z, z)
        appendLog(Rtxt("Change all NAs to the modal value (not advisable)."),
                  imp.cmd)
        eval(parse(text=imp.cmd))
      }
      else if (action == "constant")
      {
        # Take a copy of the variable to be imputed.
    
        # 100428 BUG When using Rtxt in sprintf, and substituting a UTF-8
        # encoded variable we get garbage, so convert to unknown then back
        # again. Making the Rtxt UTF-8 does not fix it.
        Encoding(z) <- "unknown"
        appendLog(sprintf(Rtxt("Impute %s."), z), copy.cmd)
        Encoding(z) <- "UTF-8"
        eval(parse(text=copy.cmd))
             
        val <- theWidget("impute_constant_entry")$getText()

        # If val is not currently a category for this variable, add it
        # in.

        if (val %notin% levels(crs$dataset[[vname]]))
        {
          levels.cmd <- sprintf(paste('levels(crs$dataset[["%s"]]) <-',
                                      'c(levels(crs$dataset[["%s"]]),',
                                      sprintf('"%s")', val)),
                                vname, vname)
          appendLog(sprintf(Rtxt("Add a new category '%s' to the variable"), val), 
                    levels.cmd)
          eval(parse(text=levels.cmd))
        }

        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 ' <- "%s"', sep=""), vname, z, val)
        appendLog(sprintf(Rtxt("Change all NAs to the constant value: %s"), val),
                  imp.cmd)
        eval(parse(text=imp.cmd))
      }
      else
      {
        # 100429 BUG When using Rtxt in sprintf, and substituting a
        # UTF-8 encoded variable we get garbage, so convert to unknown
        # then back again. Making the Rtxt UTF-8 does not fix it.
        Encoding(z) <- "unknown"
        infoDialog(sprintf(Rtxt("The option to impute the %s for the",
                                "categoric variable (%s) is not (yet)",
                                "available."), action, z))
        Encoding(z) <- "UTF-8"
      }
    }
    else
    {
      imp.val <- Rtxt("Not yet implemented.")
      
      # Take a copy of the variable to be imputed.
    
      # 100429 BUG When using Rtxt in sprintf, and substituting a
      # UTF-8 encoded variable we get garbage, so convert to unknown
      # then back again. Making the Rtxt UTF-8 does not fix it.
      Encoding(z) <- "unknown"
      appendLog(sprintf(Rtxt("Impute %s."), z), copy.cmd)
      Encoding(z) <- "UTF-8"
      eval(parse(text=copy.cmd))
      
      # Determine what action to perform.

      if (action == "zero")
      {
        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 " <- 0", sep=""), vname, z)
        imp.comment <- Rtxt("Change all NAs to 0.")
        imp.val <- 0

        # Record the transformation for inclusion in PMML.

        lst <- list(orig=z, type=vprefix, status="active", impute=imp.val)
        crs$transforms <- union.transform(crs$transforms, vname, lst)
      }
      else if (action == "mean")
      {
        # Note that if z is an integer (e.g. audit$Age) then the
        # imputed variable will be numeric.
        
        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 ' <- mean(crs$dataset[["%s"]], ',
                                 "na.rm=TRUE)", sep=""), vname, z, z)
        imp.comment <- Rtxt("Change all NAs to the mean value (not advisable).")
        imp.val <- mean(crs$dataset[[z]], na.rm=TRUE)

        # Record the transformation for inclusion in PMML.

        lst <- list(orig=z, type=vprefix, status="active", impute=imp.val)
        crs$transforms <- union.transform(crs$transforms, vname, lst)
      }
      else if (action == "median")
      {
        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 ' <- median(crs$dataset[["%s"]], ',
                                 "na.rm=TRUE)", sep=""), vname, z, z)
        imp.comment <- Rtxt("Change all NAs to the median (not advisable).")
        imp.val <- median(crs$dataset[[z]], na.rm=TRUE)

        # Record the transformation for inclusion in PMML.

        lst <- list(orig=z, type=vprefix, status="active", impute=imp.val)
        crs$transforms <- union.transform(crs$transforms, vname, lst)
      }
      else if (action == "mode")
      {
        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 ' <- modalvalue(crs$dataset[["%s"]], ',
                                 "na.rm=TRUE)", sep=""), vname, z, z)
        imp.comment <- Rtxt("Change all NAs to the modal value (not advisable).")
        imp.val <- modalvalue(crs$dataset[[z]], na.rm=TRUE)

        # Record the transformation for inclusion in PMML.

        lst <- list(orig=z, type=vprefix, status="active", impute=imp.val)
        crs$transforms <- union.transform(crs$transforms, vname, lst)
      }
      else if (action == "constant")
      {
        val <- theWidget("impute_constant_entry")$getText()
        if (is.na(as.numeric(val)))
        {
          # 100429 BUG When using Rtxt in sprintf, and substituting a
          # UTF-8 encoded variable we get garbage, so convert to
          # unknown then back again. Making the Rtxt UTF-8 does not
          # fix it.
          Encoding(z) <- "unknown"
          errorDialog(sprintf(Rtxt("The supplied value of '%s' for the variable '%s'",
                                   "is not numeric. Please supply a numeric value."),
                              val, z))
          Encoding(z) <- "UTF-8"
          next
        }
        imp.cmd <- sprintf(paste('crs$dataset[["%s"]]',
                                 '[is.na(crs$dataset[["%s"]])]',
                                 ' <- %s ', sep=""), vname, z, val)
        imp.comment <- sprintf(Rtxt("Change all NAs to the constant: %s."), val)
        imp.val <- val

        # Record the transformation for inclusion in PMML.

        lst <- list(orig=z, type=vprefix, status="active", impute=imp.val)
        crs$transforms <- union.transform(crs$transforms, vname, lst)
      }
        
      appendLog(imp.comment, "if (building)\n{\n  ",  
                imp.cmd, "\n}")
      eval(parse(text=imp.cmd))
      appendLog(Rtxt("When scoring, transform using the training data parameters:"),
                "if (scoring)\n{\n",
                sprintf(paste('  crs$dataset[["%s"]]',
                              '[is.na(crs$dataset[["%s"]])] <- %s',
                              sep=""), vname, z, imp.val),
                "\n}")
    }

    if (z %in% input)
    {
      input <- setdiff(input, z)
      input <- union(input, vname)
    }
    else if (z %in% target)
    {
      target <- setdiff(target, z)
      target <- union(target, vname)
    }
    else if (z %in% risk)
    {
      risk <- setdiff(risk, z)
      risk <- union(risk, vname)
    }
    else if (z %in% ident)
    {
      ident <- setdiff(ident, z)
      ident <- union(ident, vname)
    }
    else
    {
      # If the source variable was ignore, then leave it as such, and
      # put the new variable in as input.
      input <- union(input, vname)
    }
    ignore <- union(ignore, z)
  }
  
  if (length(imputed) > 0)
  {
    # Reset the dataset views keeping the roles unchanged except for
    # those that have been imputed, wich have just been added as
    # inputs, with the originals now ignored.

    resetDatasetViews(input, target, risk, ident, ignore, weight)
    resetTestTab()

    # Update the status bar

    setStatusBar(sprintf(Rtxt("Imputed variables added to the dataset",
                              "with '%s_' prefix."), vprefix))
  }
  else
  {
    setStatusBar(Rtxt("No variables selected to be imputed."))
  }
}  

executeTransformRemapPerform <- function(vars=NULL,
                                         action=NULL,
                                         num.bins=4,
                                         remap.prefix=paste(action, "_",
                                           sep=""),
                                         remap.comment="Remap")
{
  # Remap variables in some way:
  # quantiles, kmeans, eqwidth, log, indicator, joincat, asfactor, asnumeric

  # Obtain the list of selected variables from the treeview.

  if (is.null(vars))
  {
    selected <- theWidget("impute_treeview")$getSelection()
    selected$selectedForeach(function(model, path, iter, data)
    {
      vars <<- c(vars, model$get(iter, 1)[[1]])
    }, TRUE)
    if (length(vars)) Encoding(vars) <- "UTF-8"
  }
  
  if (length(vars) == 0)
  {
    infoDialog(Rtxt("Please select some variables to remap first. Then Execute."))
    return()
  }

  # Record the current variable roles so that we can maintain these,
  # modified appropriately by ignore'ing the binned variables, and
  # input'ing the newly binned variables.
  
  input <- getSelectedVariables("input")
  target <- getSelectedVariables("target")
  risk <- getSelectedVariables("risk")
  ident <- getSelectedVariables("ident")
  ignore <- getSelectedVariables("ignore")
  weight <- getSelectedVariables("weight")

  # Determine the action requested.

  if (is.null(action))
  {
    if (theWidget("remap_quantiles_radiobutton")$getActive())
    {
      action <- "quantiles"
      num.bins <- theWidget("remap_bins_spinbutton")$getValue()
      remap.prefix <- sprintf("BQ%d", num.bins)
      remap.comment <- sprintf(Rtxt("Bin the variable(s) into %d bins",
                                    "using quantiles."), num.bins)
    }
    else if (theWidget("remap_kmeans_radiobutton")$getActive())
    {
      action <- "kmeans"
      num.bins <- theWidget("remap_bins_spinbutton")$getValue()
      remap.prefix <- sprintf("BK%d", num.bins)
      remap.comment <- sprintf(Rtxt("Bin the variable(s) into %d bins",
                                    "using kmeans."), num.bins)
    }
    else if (theWidget("remap_eqwidth_radiobutton")$getActive())
    {
      action <- "eqwidth"
      num.bins <- theWidget("remap_bins_spinbutton")$getValue()
      remap.prefix <- sprintf("BE%d", num.bins)
      remap.comment <- sprintf(Rtxt("Bin the variable(s) into %d bins",
                                    "using equal widths."), num.bins)
    }
    else if (theWidget("remap_indicator_radiobutton")$getActive())
    {
      action <- "indicator"
      remap.prefix <- "TIN"
      remap.comment <- Rtxt("Turn a factor into indicator variables.")
    }
    else if (theWidget("remap_joincat_radiobutton")$getActive())
    {
      action <- "joincat"
      remap.prefix <- "TJN"
      remap.comment <- Rtxt("Turn two factors into one factor.")
    }
    else if (theWidget("remap_asfactor_radiobutton")$getActive())
    {
      action <- "asfactor"
      remap.prefix <- "TFC"
      remap.comment <- Rtxt("Transform into a factor.")
    }
    else if (theWidget("remap_asnumeric_radiobutton")$getActive())
    {
      action <- "asnumeric"
      remap.prefix <- "TNM"
      remap.comment <- Rtxt("Transform into a numeric.")
    }
  }
  
  # 090603 Check if it is an indicator transform, and more than a
  # single variable selected. Cannot handle this yet - do them one at
  # a time for now.

  if (action %in% c("indicator") && length(vars) > 1)
  {
    errorDialog(Rtxt("The Indicator Variable transform can only be performed",
                     "one variable at a time. Please select just one variable",
                     "to transform."))
    return(FALSE)
  }
  
  # Check if the action is one that only works on numeric data, and we
  # have any categoric variables selected. If so put up an info
  # dialogue and remove the categorics from the list of variables to
  # be imputed.

  classes <- unlist(lapply(vars, function(x) class(crs$dataset[[x]])))
  if (action %in% c("quantiles", "kmeans", "eqwidth", "log", "asfactor")
      && "factor" %in% classes)
  {
    # 100429 BUG When using Rtxt in sprintf, and substituting a UTF-8
    # encoded variable we get garbage, so convert to unknown then back
    # again. Making the Rtxt UTF-8 does not fix it.
    Encoding(vars) <- "unknown"
    infoDialog(sprintf(Rtxt("Only numeric data is permissible for the %s transform.",
                            "\n\nIgnoring: %s."),
                       switch(action,
                              quantiles=Rtxt("Quantiles"),
                              kmeans=Rtxt("KMeans"),
                              eqwidth=Rtxt("Equal Width"),
                              log=Rtxt("Log"),
                              asfactor=Rtxt("As Categoric")),
                       paste(vars[which("factor" %in% classes)], collapse=", ")))
    Encoding(vars) <- "UTF-8"
    vars <- vars[-which("factor" %in% classes)] # Remove the factors.
  }
  if (action %in% c("indicator", "joincat", "asnumeric")
      && ("numeric" %in% classes || "integer" %in% classes))
  {
    # 100429 BUG When using Rtxt in sprintf, and substituting a UTF-8
    # encoded variable we get garbage, so convert to unknown then back
    # again. Making the Rtxt UTF-8 does not fix it.
    Encoding(vars) <- "unknown"
    infoDialog(sprintf(Rtxt("Only non numeric data is permissible for the %s",
                            " transform.\n\nIgnoring: %s."),
                       switch(action,
                              indicator=Rtxt("Indicator Variable"),
                              joincat=Rtxt("Join Categorics"),
                              asnumeric=Rtxt("As Numeric")),
                       paste(vars[which(classes == "numeric" |
                                        classes == "integer")],
                             collapse=", ")))
    Encoding(vars) <- "UTF-8"
    vars <- vars[-which(classes == "numeric" | classes == "integer")]
  }

  # If, as a result of removing variables from consideration we end up
  # with no variables left, silenty exit as we have already popped up
  # a meassage about removing the categoric variables.
  
  if (length(vars) == 0) return()

  # Now that we know which variables we are remapping, we can specify
  # the actions. 080406 We set ordered=FALSE here for now because
  # randomForest does not handle them. Andy is working on it.
  
  if (action == "quantiles")
  {
    # 111025 For now, if a weight is selected in the data tab, then
    # silently do weighted binning.
    
    if(length(weight <- getSelectedVariables("weight")))
    {
      remap.cmd <- paste(sprintf(paste('  crs$dataset[["%s_%s"]] <- binning(crs$',
                                       'dataset[["%s"]], %d, method="wtd.quantile", ',
                                       'ordered=FALSE, weights=crs$dataset[["%s"]])',
                                       sep=""),
                                 remap.prefix, vars, vars, num.bins, weight),
                         collapse="\n")
    }
    else
    {
      remap.cmd <- paste(sprintf(paste('  crs$dataset[["%s_%s"]] <- binning(crs$',
                                       'dataset[["%s"]], %d, method="quantile", ',
                                       'ordered=FALSE)',
                                       sep=""),
                                 remap.prefix, vars, vars, num.bins),
                         collapse="\n")
    }
  }
  else if (action == "kmeans")
  {
    remap.cmd <- paste(sprintf(paste('  set.seed(23456)\n',
                                     '  crs$dataset[["%s_%s"]] <- binning(crs$',
                                     'dataset[["%s"]], %d, method="kmeans", ',
                                     'ordered=FALSE)',
                                     sep=""),
                               remap.prefix, vars, vars, num.bins),
                       collapse="\n")
  }
  else if (action == "eqwidth")
  {
    remap.cmd <- paste(sprintf(paste('  crs$dataset[["%s_%s"]] <- cut(crs$',
                                     'dataset[["%s"]], %d)',
                                     sep=""),
                               remap.prefix, vars, vars, num.bins),
                       collapse="\n")
  }
  else if (action == "indicator")
  {
    remap.cmd <- paste(sprintf(paste('  crs$dataset[, make.names(paste("%s_%s_", ',
                                     'levels(',
                                     'crs$dataset[["%s"]]), sep=""))] ',
                                     '<- diag(nlevels(',
                                     'crs$dataset[["%s"]]))[crs$dataset',
                                     '[["%s"]],]',
                                     sep=""),
                               remap.prefix, vars, vars, vars, vars),
                       collapse="\n")
  }
  else if (action == "joincat")
  {
    if (length(vars) > 2)
    {
      infoDialog(Rtxt("A join of only two categoric variables at a time is allowed.",
                      "Please select just two categoric variables, then Execute."))
      return()
    }
    if (length(vars) < 2)
    {
      infoDialog(Rtxt("A join of categoric variables requires two categoric variables.",
                      "Please select two categoric variables, then Execute."))
      return()
    }
    
      
    remap.cmd <- sprintf(paste('  crs$dataset[, "%s_%s_%s"] <- ',
                               'interaction(paste(crs$dataset[["%s"]], "_",',
                               'crs$dataset[["%s"]], sep=""))\n',
                               '  crs$dataset[["%s_%s_%s"]]',
                               '[grepl("^NA_|_NA$", crs$dataset[["%s_%s_%s"]])]',
                               ' <- NA\n',
                               '  crs$dataset[["%s_%s_%s"]] <- ',
                               'as.factor(as.character(crs$dataset[["%s_%s_%s"]]))',
                               sep=""),
                         remap.prefix, vars[1], vars[2],
                         vars[1], vars[2],
                         remap.prefix, vars[1], vars[2],
                         remap.prefix, vars[1], vars[2],
                         remap.prefix, vars[1], vars[2],
                         remap.prefix, vars[1], vars[2])
  }
  else if (action == "asfactor")
  {
    remap.cmd <- paste(sprintf(paste('  crs$dataset[["%s_%s"]] <- ',
                                     'as.factor(crs$',
                                     'dataset[["%s"]])',
                                     sep=""),
                               remap.prefix, vars, vars),
                       collapse="\n")

    # 090718 Remap the levels to correspond to how the transform will
    # appear when exported to PMML.

    ol <- NULL # 090808 Keep "R check" happy.

    relevel.cmd <- paste(sprintf('\n  ol <- levels(crs$dataset[["%s_%s"]])',
                                 remap.prefix, vars),
                         "  lol <- length(ol)",
                         paste('  nl <- c(sprintf("[%s,%s]", ol[1], ol[1]),',
                               'sprintf("(%s,%s]", ol[-lol], ol[-1]))'),
                         sprintf('  levels(crs$dataset[["%s_%s"]]) <- nl',
                                 remap.prefix, vars),
                         sep="\n", collapse="\n")
    remap.cmd <- paste(remap.cmd, relevel.cmd, sep="\n")
  }
  else if (action == "asnumeric")
  {
    remap.cmd <- paste(sprintf(paste('  crs$dataset[["%s_%s"]] <- ',
                                     'as.numeric(crs$',
                                     'dataset[["%s"]])', sep=""),
                               remap.prefix, vars, vars),
                       collapse="\n")
  }
  
  # Perform the remapping.

  startLog(Rtxt("Remap variables."))
  appendLog(remap.comment,
            # 090601 build/score difference only for some remap ops.
            ifelse(action %in% c('quantiles', 'kmeans', "eqwidth"),
                   sprintf("if (building)\n{\n %s\n}\n", remap.cmd),
                   remap.cmd))
  eval(parse(text=remap.cmd))

  # Record the transformation as well as reporting to the log.

  vname <- paste(remap.prefix, vars, sep="_")

  if (action %in% c('quantiles', 'kmeans', "eqwidth"))
  {
    # 090109 At this stage the remapping has been performed. Why not
    # get the breaks from the levels of the variables instead of
    # relying on other information.

    breaks <- sapply(vname, function(x)
                      sort(unique(as.numeric(unlist(strsplit(
                      gsub(",", " ",
                      gsub("\\(|\\]|\\[", "",
                      paste(levels(crs$dataset[[x]]),
                            collapse=" "))), " "))))))

    # 090108 The binning command returns as an attribute the
    # breaks. We extract the breaks, add that to the crs$transforms.

    ## lst <- paste(remap.prefix, vars, sep="_")
    ## if (action == "eqwidth")
    ##   breaks <- sapply(vars, function(x)
    ##                   sort(unique(as.numeric(unlist(strsplit(
    ##                   gsub(",", " ",
    ##                   gsub("\\(|\\]|\\[", "",
    ##                   paste(levels(cut(crs$dataset[[x]], num.bins)),
    ##                         collapse=" "))), " "))))))
    ## else
    ##   breaks <- sapply(lst, function(x) sort(attr(crs$dataset[[x]], "breaks")))

    # 090606 New transforms data structure

    ## lst <- paste(lst, apply(breaks, 2, paste, collapse="_"), sep="_")
    ## crs$transforms <- union(crs$transforms, lst)

    lst <- list(orig=vars, type=substr(remap.prefix, 1, 2), status="active",
                breaks=as.vector(breaks))
    crs$transforms <- union.transform(crs$transforms, vname, lst)
    
    appendLog(Rtxt("When scoring, use the training data parameters to bin new data."),
              "if (scoring)\n{\n",
              # Print the transforms based on the training parameters
              paste(sprintf(paste('  crs$dataset[["%s"]] <- ',
                                  'cut(crs$dataset[["%s"]],\n    ',
                                  '%s,\n    include.lowest=TRUE)', sep=""),
                            paste(remap.prefix, vars, sep="_"),
                            vars,
                            paste("c(", apply(breaks, 2, paste,
                                              collapse=','),
                                  ")", sep="")),
                    collapse="\n"),
              # Comment the transforms based on the test parameters.
              Rtxt("\n\n# Alternatively, use the min/max from the new data.\n\n"),
              paste(sprintf(paste('#  crs$dataset[["%s"]] <- ',
                                  'cut(crs$dataset[["%s"]],\n#    ',
                                  '%s,\n#    include.lowest=TRUE)',
                                  sep=""),
                            paste(remap.prefix, vars, sep="_"),
                            vars,
                            sprintf(sub(',[^),]*)$',
                                        ', max(crs$dataset[["%s"]], na.rm=TRUE))',
                                        sub('c\\([^,]*,',
                                            'c(min(crs$dataset[["%s"]], na.rm=TRUE),',
                                            paste("c(", apply(breaks, 2, paste,
                                                              collapse=','),
                                                  ")", sep=""))), vars, vars)),
                    collapse="\n"),
              "\n}")
    
  }
  else if (action == "indicator")
  {
    # 090606 Record each of the newly added variables in the list of
    # transforms.

    sapply(levels(crs$dataset[,vars]), function(x) 
           {
             lst <- list(orig=vars, type=remap.prefix, status="active", level=x)
             crs$transforms <- union.transform(crs$transforms,
                                               make.names(paste(vname, x,
                                                                sep="_")), lst)
           })
  }
  else if (action == "joincat")
  {
    lst <- list(orig=vars, type=remap.prefix, status="active",
                levels=list(levels(crs$dataset[,vars[1]]),
                  levels(crs$dataset[,vars[2]])))
    names(lst$levels) <- vars
    crs$transforms <- union.transform(crs$transforms,
                                      paste(remap.prefix, vars[1], vars[2], sep="_"),
                                      lst)
  }
  else if (action == "asfactor")
  {
    # 090707 Add in a breaks attribute. We have already transformed
    # the variable so we can get the information from it. Also note
    # that unlike the other binning transforms, the name of the bins
    # here are the values of the variable, not ranges. We also need to
    # repeat the first entry, so that when we represent the ranges in
    # the XML we get the first range as, for example, [0,0].
    #
    # 090718 The levels have been converted to be ranges, to conform
    # with the bins, thus the breaks may now need to be determined
    # from the second value of each level. However, we used the
    # variable "ol" above to record the old levels, so we can use that
    # as the breaks. 100417 But this use of ol does npot work when
    # there is more than a single transform. So get the breaks from
    # the original transformation.

    for (v in vars)
    {
      # 090718 breaks <- as.numeric(levels(crs$dataset[[vname]]))
      breaks <- as.numeric(levels(as.factor(crs$dataset[[v]])))
      # 100417 breaks <- as.numeric(ol)
      breaks <- c(breaks[1], breaks)
      lst <- list(orig=v, type=remap.prefix, status="active", breaks=breaks)
      crs$transforms <- union.transform(crs$transforms,
                                        paste(remap.prefix, v, sep="_"),
                                        lst)
    }
  }
  else if (action == "asnumeric")
  {
    for (v in vars)
    {
      lst <- list(orig=v, type=remap.prefix, status="active",
                  levels=levels(crs$dataset[[v]]))
      crs$transforms <- union.transform(crs$transforms,
                                        paste(remap.prefix, v, sep="_"),
                                        lst)
    }
  }
  
  # Record the new variables as having an INPUT role. 090110
  # Previously implemented no other changes as the original variables
  # are probably still required for modelling, but that was
  # "suprising" so now ignore the originals, like the other
  # transforms. Particularly because the new and old variables will be
  # correlated, and so lead to singularities in regression models.

  if (action == "joincat")
  {
    input <- union(input, paste(remap.prefix, vars[1], vars[2], sep="_"))
    ignore <- union(ignore, vars)
  }
  else if (action == "indicator")
  {
    # 090603 This needs work to support multiple variables at the one
    # time. 090724 Add the make.names to catch levels that have spaces
    # in their names. 101120 By default all but the first are Input,
    # and the first is Ignore.

    new.vars <- paste(remap.prefix, vars,
                      make.names(levels(crs$dataset[[vars]])), sep="_")
    
    input <- union(input, new.vars[-1])
    ignore <- union(ignore, union(vars, new.vars[1]))
  }
  else if (action %in% c('quantiles', 'kmeans', "eqwidth", "asfactor", "asnumeric"))
  {
    # 090722 Added asfactor and asnumeric here - avoid surprising the
    # user as they are probably expecting the old variables to be
    # ignored, like most other operations. I can't see any reason for
    # keeping the originals as input, except if the user really wanted
    # to, which they can purposfully do.
    
    input <- setdiff(input, vars) # 090110 Added
    ignore <- union(ignore, vars) # 090110 Added
    input <- union(input, paste(remap.prefix, vars, sep="_"))
  }
  else
    input <- union(input, paste(remap.prefix, vars, sep="_"))

  # 090731 Remove any vars as a target/risk/ident, since they have
  # been ignored and we don't want both!
  
  target <- setdiff(target, vars)
  risk <- setdiff(risk, vars)
  ident <- setdiff(ident, vars)

  # Reset the dataset views keeping the roles unchanged except for
  # those that have been created, wich have just been added as inputs.

  resetDatasetViews(input, target, risk, ident, ignore, weight)
  resetTestTab()
  
  # Update the status bar
  
  setStatusBar(sprintf(Rtxt("Remapped variables added to the dataset",
                            "with '%s_' prefix."), remap.prefix))
}

#-----------------------------------------------------------------------

executeTransformCleanupPerform <- function()
{
  # First, record the current variable roles so that we can maintain
  # these, modified appropriately.
  
  input <- getSelectedVariables("input")
  target <- getSelectedVariables("target")
  risk <- getSelectedVariables("risk")
  ident <- getSelectedVariables("ident")
  ignore <- getSelectedVariables("ignore")
  weight <- getSelectedVariables("weight")

  startLog("CLEANUP the Dataset")

  if (theWidget("delete_ignored_radiobutton")$getActive())
  {
    if (variablesHaveChanged(Rtxt("deleting the selected ignored variables"))) return()
    to.delete <- getSelectedVariables("ignore")
  }    
  else if (theWidget("delete_selected_radiobutton")$getActive())
  {

    # Obtain the list of selected variables from the treeview.

    to.delete <- NULL
    selected <- theWidget("impute_treeview")$getSelection()
    selected$selectedForeach(function(model, path, iter, data)
    {
      to.delete <<- c(to.delete, model$get(iter, 1)[[1]])
    }, TRUE)
    if (length(to.delete)) Encoding(to.delete) <- "UTF-8"

    if (length(to.delete) == 0)
    {
      infoDialog(Rtxt("Please select some variables to delete first. Then Execute."))
      return()
    }
  }
  else if (theWidget("delete_navars_radiobutton")$getActive())
  {
    # Get a list of all variables. For now (and perhaps always),
    # ignore the role.
    
    to.delete <- names(crs$dataset)

    # Remove from the list any variables that do not have missing
    # values.
    
    for (v in to.delete)
      if (sum(is.na(crs$dataset[[v]])) == 0)
        to.delete <- setdiff(to.delete, v)
  }
  else if (theWidget("delete_naents_radiobutton")$getActive())
  {
    # Here, ignore the variables that have a role of Ignore, so we
    # only delete entities that have missing values for non-ignored
    # variables.

    if (! length(ignore))
    {
      cases <- complete.cases(crs$dataset)
      del.cmd <- "crs$dataset <- crs$dataset[complete.cases(crs$dataset),]"
    }
    else
    {
      cases <- complete.cases(crs$dataset[,-getVariableIndicies(ignore)])
      del.cmd <- sprintf(paste("crs$dataset <- crs$dataset[complete.cases(",
                               "crs$dataset[,-%s]),]", sep=""),
                         simplifyNumberList(getVariableIndicies(ignore)))
    }
    
    if (! questionDialog(sprintf(Rtxt("We are about to delete %d",
                                      "entities from the in-memory dataset.",
                                      "These have missing values for some of the",
                                      "non-Ignore variables.\n\nAre you sure you",
                                      "want to delete these entities?"),
                                 sum(!cases))))
      return()

    # Perform the deletions.
  
    appendLog(Rtxt("Remove rows with missing values."), del.cmd)
    eval(parse(text=del.cmd))

  }

  if (!theWidget("delete_naents_radiobutton")$getActive())
  {
    mapped <- FALSE
    if ("UTF-8" %in% Encoding(to.delete))
    {
      mapped <- TRUE
      Encoding(to.delete) <- "unknown"
    }    
    if (! questionDialog(sprintf(Rtxt("We are about to delete the following variables",
                                      "from the in-memory dataset.",
                                      "This will permanently remove them from",
                                      "the memory copy of the data, but will not",
                                      "affect any file system copy.\n\n",
                                      "Delete: %s",
                                      "\n\nAre you sure you want to delete these",
                                      "variables?"),
                                 paste(to.delete, collapse=", "))))
      return()
    if (mapped) Encoding(to.delete) <- "UTF-8"

    del.cmd <- paste(sprintf('crs$dataset$%s <- NULL', to.delete),
                     collapse="\n")
    del.comment <- Rtxt("Remove specific variables from the dataset.")

    # Perform the deletions.
  
    appendLog(del.comment, del.cmd)
    eval(parse(text=del.cmd))

    # Ensure any deleted variables are no longer included in the list
    # of transformed variables. 090606 Modified to work with new
    # transforms data structure. Note that we are only removing
    # deleted transformed variables. What about when we delete a
    # variable that a transform depends on! 090701 Perhaps we don't
    # delete them from the list of transforms. Is there anything to
    # lose keeping the information there, in particular if there are
    # other transforms derived from any transforms about to be
    # deleted. 090801 Restore the removal of deleted transforms since
    # otherwise they get used to identify source variables to
    # export. We may like to add some logic here to not allow any
    # variable to be deleted if it is used in a transform
    # (irrespective of whether that tranform is ignored).

    crs$transforms[names(crs$transforms) %in% to.delete] <- NULL

  }
  
  # Reset the dataset views keeping the roles unchanged except for
  # those that have been delete.

  resetDatasetViews(input, target, risk, ident, ignore, weight)
  resetTestTab()

  # Update the status bar

  setStatusBar(Rtxt("The deletions from the dataset have been completed."))
}
