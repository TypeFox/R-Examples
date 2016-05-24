## UNDER DEVELOPEMENT - HOW TO MERGE WITH pmml.transforms.R

# PMML: Predictive Modelling Markup Language
#
# Part of the Rattle package for Data Mining
#
# Time-stamp: <2012-12-04 06:11:03 Graham Williams>
#
# Copyright (c) 2009-2013 Togaware Pty Ltd
# Copyright (c) 2009 Information Builders Inc
#
# This files is part of the Rattle and RStat suite for Data Mining in R.
#
# This file will become freely available under an open source license
# three years from 17 May 2009.
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
# TRANSFORMS FRAMEWORK

# GLOBAL CONSTANTS

# 090517 List here all transforms that we know about from the rattle
# package. In using the pmml package with rattle these are used to
# identify tranformations. Information about how the transformation is
# performed (i.e., the transformation parameters) is separately
# recorded within rattle (crs$transforms) and must be used separately
# to generate the appropriate PMML code. Eventually I expect to see
# the handling of transform variables incorporated into the pmml
# package.

.TRANSFORMS.NORM.CONTINUOUS <- c("RRC", "R01", "RMD", "RMA")
.TRANSFORMS.IMPUTE <- c("IZR", "IMN", "IMD", "IMO", "ICN")
.TRANSFORMS.APPLY <- c("RLG", "R10")
.TRANSFORMS.BIN <- c("BQ", "BK", "BE", "TFC")
.TRANSFORMS.INDICATOR <- c("TIN")

.TRANSFORMS.OTHER.NUM <- c("RRK", "RIN", "BGC", "BG1", "BGD", "BGK")
.TRANSFORMS.OTHER.CAT <- c("TJN", "TNM")
.TRANSFORMS.OTHER <- c(.TRANSFORMS.OTHER.NUM, .TRANSFORMS.OTHER.CAT)

.TRANSFORMS <- c(.TRANSFORMS.NORM.CONTINUOUS,
                 .TRANSFORMS.APPLY,
                 .TRANSFORMS.IMPUTE,
                 .TRANSFORMS.BIN,
                 .TRANSFORMS.INDICATOR,
                 .TRANSFORMS.OTHER)

.TRANSFORMS.TO.CATEGORIC <- c(.TRANSFORMS.BIN, "TJN")

.supportTransformExport <- function(transforms=NULL)
{
  # 090208 Returns TRUE/FALSE. The default is that we don't have
  # Rattle code to generate the transforms in the pmml - i.e.,
  # pmml.transforms is not defined. A package which does will define a
  # pmml.transforms.
  
  return(exists("pmml.transforms") &&
         # length(getAnywhere("pmml.transforms")$objs) > 0 &&
         ! is.null(transforms))
}

.activateTransforms <- function(transforms)
{
  # 090813 For each transform ensure the stats is set to active. This
  # is needed when we have multiple transforms and we have remove from
  # the list of transforms those which are truely ignored (i.e.,
  # inactive). Those that remain with an inactive status are then
  # assumed to be used for active transforms. So make them active.

  return(lapply(transforms, function(x) {x$status <- "active"; return(x)}))
}

.activateDependTransforms <- function(transforms)
{
  # 090813 For each transform, starting from the last (assuming the
  # order of multiple transforms goes from the beginning to the end -
  # i.e., RRC_TNM_EDUCATION appears later than TNM_EDUCATION),
  # activate any inactive transforms that are needed for later active
  # transforms.

  for (i in rev(seq_along(transforms)))
  {
    if (i == length(transforms)) next
    if (names(transforms)[i] %in%
        sapply(transforms[seq(i+1, length(transforms))],
               function(x) ifelse(x$status=="active", x$orig, NA)))
      transforms[[i]]$status <- "active"
  }
  return(transforms)
}

.unifyTransforms <- function(field, transforms, keep.first=TRUE)
{
  # 090102 Unify the list of variables in FIELD based on the known
  # TRANSFORMS, so that the variables in the TRANSFORMS are removed
  # from FIELD and replaced with the original variables from which the
  # TRANFROMS are derived. 090617 The final list of fields then are those
  # variables which need to be supplied to run the model.
  #
  # The variable TRANSFORMS is a list of transform structures (lists),
  # and each list has at least named elements called type, original
  # name, and status.

  # 090111 We make sure that the FIELD class is reset appropriately,
  # in special cases. So a binning transform, which is now a
  # categoric, was originally a numeric! Change it back to numeric.

  # 090806 Obtain a list of those varaibles that are listed as the
  # original variable within an active transform. This is used to
  # determine which transforms to ignore. Some inactive transforms are
  # needed because there is an active transform derived from it.

  active.by.orig <- unique(unlist(sapply(transforms, function(x)
                                         if (x$status=="active") x$orig)))
  
  for (i in seq_along(transforms))
  {
    var <- transforms[i][[1]]
    type <- var$type
    vname <- names(transforms)[i] # The transformed variable name

    # 090801 If the status is inactive and there is no active
    # transform that uses this vname (i.e., it is not in the
    # active.by.orig vector), then skip its processing - the variable
    # is not required hence it's source variable is not required.

    status <- transforms[i][[1]]$status
    if (status == "inactive" &&
        (is.null(active.by.orig) ||
         ! (vname  %in% active.by.orig))) next
    
    # 090102 The vname should be in the list of variables in field,
    # where we replace it with the original name. 090617 If it is not
    # in the list of variables then it should be an "intermediate"
    # transform, which is used later in another transform. We need to
    # handle this.

    # 090617 if (vname %in% field$name)
    # 090617 {
    index <- which(vname == field$name)
    if (length(index))
    {
          
      # 081229 I should probably be testing if index is integer(0) -
      # 090617 Now doing so.

      # 090607 REMOVE oname <- sub("^([^_]*)_", "", field$name[index])

      # 090111 For the binning operations, be sure to change the class
      # from categoric back to numeric. 091105 For TIN transform, be
      # sure to change the class to categoric since the original
      # variable will have been a categoric variable. Does this
      # capture them all?

      if (type %in% .TRANSFORMS.BIN)
        field$class[index] <- "numeric"
      else if (type %in% c("TIN", "TNM"))
        field$class[index] <- "factor"

      # 090102 If the original variable (the untransformed variable)
      # for this transformed variable vname is already in the input
      # variables or is itself a transformed variable, then simply
      # remove the entry naming the transformed variable vname,
      # otherwise replace the entry naming the transformed var with
      # the original variable name.

      for (v in var$orig)
      {
        if (v %in% union(field$name, names(transforms)))
        {
          # 090819 Not 100% this is correct, but it captures where
          # TJN_EDUCATION_GENDER is a transform and GENDER is also the
          # target. the second time through this loop v is GENDER, but
          # we do not want to remove the "index" item since that will
          # have been set to EDUCATION the first time through the
          # loop!
          if (length(var$orig) == 1)
          {
            field$name <- field$name[-index]
            field$class <- field$class[-index]
          }
        }
        else
        {
          # 090724 For TJN with two orig variables we loop through
          # here twice, so the second time through (i.e., when v is
          # not the same as var$orig[1]) add the dependent variable to
          # the end of the list and record a class for it - assumed to
          # always be factor.
          if (v != var$orig[1]) index <- length(field$name) + 1
          field$name[index] <- v
          if (v != var$orig[1]) field$class[index] <- "factor"
        }
        
      }
    }
    else if (length(which(var$orig == field$name)) == 0) # Orig not in field
    {
      field$name <- c(field$name, var$orig)
      # 090724 Add this rep to account for TJN which has multiple orig vars.
      field$class <- c(field$class, rep(.origVarType(var$type), length(var$orig)))
    }
  }

  # 090102 Reset the field$class names to correspond to the new
  # variables.

  names(field$class) <- field$name

  # 090801 If this is called in the context of a running Rattle, then
  # sort the fields to be the same order they appear in
  # names(crs$dataset). We need to keep the first field as is if it is
  # the target field to ensure the follow on logic for predictive
  # models. Not required for descriptive models.

  ## Remove this for now - remove the dependency between pmml and
  ## rattle
  ##
  ## if (exists("crs") && ! is.null(crs$dataset))
  ## {
  ##   if (keep.first)
  ##   {
  ##     orig <- field
  ##     f1 <- field$class[1]
  ##     field$class <- field$class[-1]
  ##     field$name <- field$name[-1]
  ##     field$class <- field$class[names(crs$dataset)[names(crs$dataset) %in% field$name]]
  ##     field$class <- c(f1, field$class)
  ##     field$class <- c(field$class, orig$class[! (orig$name %in% names(field$class))])
  ##     field$name <- names(field$class)
  ##   }
  ##   else
  ##   {
  ##     orig <- field
  ##     field$class <- field$class[names(crs$dataset)[names(crs$dataset) %in% field$name]]
  ##     field$class <- c(field$class, orig$class[! (orig$name %in% names(field$class))])
  ##     field$name <- names(field$class)
  ##   }
  ## }
  
  return(field)
}

pmmlCanExport <- function(vname)
{
  # 090517 In general this support function will return TRUE if the
  # variable name supplied is either not a transform (hence it is
  # exportable in general, since it is a supplied variable) or if it
  # is a transform, it is a supported transformation, if there is such
  # a list available. This is intended for the implementation of
  # transforms into PMML. An extermal package will provide a
  # pmml.transforms function and will override the simple default here
  # to test whether the transform is one that it can handle.
  return(! .isTransformVar(vname))
}

.isTransformVar <- function(vname)
{
# 090607
#  print(vname)
#  print(.transformType(vname))
#  print(.TRANSFORMS)
#  cat(sprintf("Option 1: %s\n", .transformType(vname) %in% .TRANSFORMS))
#  cat(sprintf("Option 2: %s\n", vname %in% names(crs$transforms)))
#  cat("Going with option 1\n")
  return(.transformType(vname) %in% .TRANSFORMS)
# 090612 Revert to the above since crs$transforms is not available in
# the pmml package! This is used in pmmltocibi and ends up wrapping
# transforms with FloatVal in the wrong place.
#
#  print(vname %in% names(crs$transforms))
#  return(vname %in% names(crs$transforms))

  # 090613 It seems that returning FALSE always is okay! Thus need to
  # remove this throughout the code.

#  return(FALSE)
}

.origVarType <- function(ty)
{
  if (ty %in% c(.TRANSFORMS.NORM.CONTINUOUS,
                .TRANSFORMS.IMPUTE, # 090617 But IZR and ICN are ambiguous.
                .TRANSFORMS.APPLY,
                .TRANSFORMS.BIN,
                .TRANSFORMS.OTHER.NUM))
    return("numeric")
  else
    return("factor")
}


.transformType <- function(tr)
{
  # 090607 Return the code of the type of transform representated by
  # the variable TR. The supplied parameter TR must be a transform
  # data structure (list) and this function simply returns the "type"
  # element of the list. 090609 However, pmmltocibi currently calls
  # this function often with a strin g variable, so let's keep that
  # working for now until migrated.
  
  # 090517 Should we be returning NULL if there is no type recorded?
  # Or perhaps an empty string so that %in% will still work. 09607
  # Note that there should always be a type.

  if (class(tr) == "character") # 090709 Remove this sometime - needed
                                # in pmmltocibi for now.
  {
    code <- sub("^([^_]*)_.*$", "\\1", tr)
    if (substr(code, 1, 2) %in% .TRANSFORMS.BIN)
      code <- substr(code, 1, 2)
    return(code)
  }
  else
    return(tr[[1]]$type)
}

.transformToDerived <- function(var)
{
  # 090607 REMOVE This function should not be needed any more with the
  # new transforms data structure (list). Remove it eventually, but
  # for now simply return the supplied var.

  return(var)
  
  if (.transformType(var) %in% .TRANSFORMS.APPLY)
    # No args
    return(var) 
  else if (.transformType(var) %in% .TRANSFORMS.IMPUTE)
    # One arg
    return(sub("^([^_]*_)(.*)(_[^_]*)$", "\\1\\2", var))
  else if (.transformType(var) %in% .TRANSFORMS.NORM.CONTINUOUS)
    # two args
    return(sub("^([^_]*_)(.*)(_[^_]*){2}$", "\\1\\2", var))
  else if (.transformType(var) %in% .TRANSFORMS.BIN)
  {
    nparms <- as.integer(sub("^..([^_]*)_.*$", "\\1", var)) + 1
    return(sub(sprintf("^([^_]*_)(.*)(_[^_]*){%s}$", nparms), "\\1\\2", var))
  }
  else if (.transformType(var) %in% .TRANSFORMS.INDICATOR)
    # 090603 Split the values out and return list of transform var
    # names
    return(paste(.transformToBasename(var), strsplit(var, "_")[[1]][-c(1,2)], sep="_"))
  else
    return(NULL)
}

.transformToBasename <- function(var)
{
  # 090607 REMOVE No longer used. Simply return var. Remove this eventually.

  return(var)
  
  # 090603 This is awaiting turning crs$transforms into a structure to
  # handle this properly, then we won't need this kludge. For now we
  # assume no variable has _ in its name, but limit this assumption
  # to indicator variables for now.

  # 090601 For TIN_Marital_Married_..., return TIN_Marital.

  if (substr(var, 1, 4) == "TIN_")
    # 090603 TIN now lists all levels in one go.
    return(sub("^([^_]+_[^_]+)_.*$", "\\1", var))
    # 090601 Bit of an assumption that the levels do not have
    # underscores in them, but assume so for now. Strip off the level
    # from the transformed variable's name and return that.
    ## return(sub("^(.*)_[^_]+$", "\\1", var))
  else
    return(var)
}

.transformToOriginal <- function(var)
{
  return(sub("^[^_]*_", "", .transformToDerived(var)))
}

.transformToParms <- function(var)
{
  if (.transformType(var) %in% .TRANSFORMS.BIN)
  {
    nparms <- as.integer(sub("^..([^_]*)_.*$", "\\1", var)) + 1
    return(unlist(strsplit(
           sub("^ ", "",
           gsub("_", " ",
           sub(sprintf("^([^_]*_)(.*)((_[^_]*){%s})$",
                       nparms), "\\3", var))), " ")))
  }
}

