# Copyright 2009-2014 Meik Michalke <meik.michalke@hhu.de>
#
# This file is part of the R package klausuR.
#
# klausuR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# klausuR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with klausuR.  If not, see <http://www.gnu.org/licenses/>.


#' Comparison of data sets
#'
#' The function \code{compare} will take two data.frames (or objects of class \code{\link[klausuR]{klausuR.answ-class}})
#' and compare them for equality. This is useful to check for typos before you calculate the results with
#' \code{\link[klausuR:klausur]{klausur}}. If you need to type in the given answers by hand, errors 
#' easily occur, so it is advisable to input all data at least twice (perhaps by different persons) and check for differences
#' with this function, which can then be corrected by looking up the original answer in the test.
#'
#' If you don't want to compare all variables but only a subset, you can use the \code{select} option (see examples below).
#' But be careful with this, at least all the values given in \code{id} are needed to produce the output table.
#'
#' If \code{new.set=TRUE}, a new data.frame will be returned, that is identical in both sets compared, but all dubious values
#' will be replaced by \code{NA}.
#' 
#' @param set1,set2 The data sets to be compared. Can be two data.frames or objects of class \code{\link[klausuR]{klausuR.answ-class}}.
#'    If the latter, their slots \code{id} and \code{items} will be compared.
#' @param select A vector with variables that should be compared, all others are omitted. At least all the values given in \code{id} are needed for the output!
#'    If \code{NULL}, all variables are examined.
#' @param ignore A vector with variables that should be dropped from both sets. See also \code{select}.
#' @param new.set Logical. If \code{TRUE}, a data.frame of the compared sets is returned, with all unequal cells set to NA.
#' @param trim Logical. Indicates wheter whitespace in character variables should be trimmed.
#' @param rename A named vector defining if variables in \code{set1} and \code{set2} need to be renamed into the klausuR name scheme. Accepts elements
#'    named \code{No}, \code{Name}, \code{FirstName}, \code{MatrNo}, \code{Pseudonym} and \code{Form}. The values of these elements
#'    represent the variable names of the input data.
#' @param id A named list of character vectors to help identify differing cases in the input data. The element names of this list will become column names in the
#'    generated output table, their values define the respective column names of the input data. If a value has more than one element, they will be collapsed into
#'    one string for the output.
#' @return If \code{new.set=FALSE}, a data.frame of the differences, if found (if not, just a message is returned). Otherwise returns a combined data.frame (see details).
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @keywords utilities
#' @seealso \code{\link[klausuR:klausur]{klausur}}
#' @export
#' @examples
#' \dontrun{
#' data(antworten)
#'
#' # create some differences
#' antworten2 <- antworten[-3, -7]
#' antworten2[4,6] <- NA
#' antworten2[8,8:10] <- antworten2[8,8:10] + 1
#'
#' # default comparison
#' compare(antworten, antworten2)
#'
#' # compare only variables 1 to 12
#' compare(antworten, antworten2, select=c(1:12))
#'
#' # omit variables 3 to 8 and create a new set called "antworten.comp"
#' # from the results
#' antworten.comp <- compare(antworten, antworten2, select=-c(3:8), new.set=TRUE)
#' }

compare <- function(set1, set2, select=NULL, ignore=NULL, new.set=FALSE, rename=c(), trim=FALSE, id=list(No="No", Name=c("FirstName", "Name"))){
  # to avoid NOTEs from R CMD check:
  MatrNo <- NULL

  # function to rename columns
  rename.vars <- function(rename=rename, set){
      vars.to.rename <- names(rename)
      id.names <- c("No", "Name", "FirstName", "MatrNo")
      id.possible.names <- c(id.names, "Pseudonym", "Form")
      if(any(vars.to.rename %in% names(set))){
        double.vars <- vars.to.rename[vars.to.rename %in% names(set)]
        warning(paste("Probably duplicate variable names found, please double check the outcome:\n ", paste(double.vars, collapse=", ")), call.=FALSE)
      } else {}
      id.invalid.names <- vars.to.rename[!vars.to.rename %in% id.possible.names]
      if(length(id.invalid.names) > 0){
        stop(simpleError(paste("Invalid variable names in 'rename':\n ",
          paste(id.invalid.names, collapse=", "))))
      } else {}
      # rename columns, if any
      for (ren.var in vars.to.rename){
        ren.from <- rename[ren.var]
        ren.to  <- ren.var
        colnames(set)[colnames(set) == rename[ren.var]] <- ren.var
      }
    return(set)
  }

  # get the names of the given sets, to better understand the outcome later
  set1.name <- deparse(substitute(set1))
  set2.name <- deparse(substitute(set2))

  # see if we have klausuR.answ class objects
  if(inherits(set1, "klausuR.answ")){
    set1 <- cbind(set1@id, subset(set1@items, select=-MatrNo))
  } else {}
  if(inherits(set2, "klausuR.answ")){
    set2 <- cbind(set2@id, subset(set2@items, select=-MatrNo))
  } else {}

  # see if we have factors; as of now, they are converted to character, with a warning
  data.factors1 <- names(set1)[sapply(names(set1), function(var){is.factor(set1[,var])})]
  data.factors2 <- names(set2)[sapply(names(set2), function(var){is.factor(set2[,var])})]
    if(length(data.factors1) > 0 | length(data.factors2) > 0){
      if(length(data.factors1) > 0){
        for (this.factor in data.factors1){
          set1[,this.factor] <- as.character(set1[,this.factor])
        }
        factor.message1 <- paste("Converted factors to character in ", set2.name, ":\n    ", paste(data.factors1, collapse=", "), "\n  ", sep="")
      } else {
        factor.message1 <- ""
      }
      if(length(data.factors2) > 0){
        for (this.factor in data.factors2){
          set2[,this.factor] <- as.character(set2[,this.factor])
        }
        factor.message2 <- paste("Converted factors to character in ", set1.name, ":\n    ", paste(data.factors2, collapse=", "), "\n  ", sep="")
      } else {
        factor.message2 <- ""
      }
      warning(paste("\nThe objects contained factors which were converted into character:\n  ", factor.message1, factor.message2, sep=""), call.=FALSE)
    }

  # probably trim whitespace
  if(isTRUE(trim)){
    data.char1 <- names(set1)[sapply(names(set1), function(var){is.character(set1[,var])})]
    data.char2 <- names(set2)[sapply(names(set2), function(var){is.character(set2[,var])})]
    if(length(data.char1) > 0){
      for (this.char in data.char1){
        set1[,this.char] <- gsub("^[[:space:]]+", "", gsub("[[:space:]]+$", "", set1[,this.char]))
      }
    } else {}
    if(length(data.char2) > 0){
      for (this.char in data.char2){
        set2[,this.char] <- gsub("^[[:space:]]+", "", gsub("[[:space:]]+$", "", set2[,this.char]))
      }
    } else {}
  } else {}

  # first thing, if both sets are indeed the same, we can quit immediately
  if(identical(set1, set2)) {
    message(paste("\nThe compared objects (",set1.name," & ",set2.name,") are identical!\n", sep=""))
    if(isTRUE(new.set)){
      # do we need renaming?
      if(length(rename) > 0){
        set1 <- rename.vars(rename=rename, set=set1)
      } else {}
      return(set1)
    } else {
      return(invisible(NULL))
    }
  } else {
    # before any values are even compared, check for equality of elements
    if(!identical(names(set1), names(set2))){
      missingVarsIn1 <- which(!names(set1) %in% names(set2))
      missingVarsIn2 <- which(!names(set2) %in% names(set1))
      missingColsIn1 <- names(set1)[missingVarsIn1]
      missingColsIn2 <- names(set2)[missingVarsIn2]
      if(length(missingVarsIn1)>0){
        set1 <- set1[,-missingVarsIn1]
        missing.message1 <- paste("Missing columns in ", set2.name, ":\n    ", paste(missingColsIn1, collapse=", "), "\n  ", sep="")
      } else {
        missing.message1 <- ""
      }
      if(length(missingVarsIn2)>0){
        set2 <- set2[,-missingVarsIn2]
        missing.message2 <- paste("Missing columns in ", set1.name, ":\n    ", paste(missingColsIn2, collapse=", "), "\n  ", sep="")
      } else {
        missing.message2 <- ""
      }
      warning(paste("\nThe objects do not include columns of the same name. Missing ones have been ignored in both:\n  ", missing.message1, missing.message2, sep=""), call.=FALSE)
    } else {}

    # do we need renaming?
    if(length(rename) > 0){
      set1 <- rename.vars(rename=rename, set=set1)
      set2 <- rename.vars(rename=rename, set=set2)
    } else {}

    # validate the demanded id values
    id.names <- names(id)
    id.values <- unlist(id)
    invalid.id.values <- id.values[!id.values %in% names(set1)]
    if(length(invalid.id.values)>0){
      stop(simpleError(paste("Invalid values for 'id', cannot be found in the given sets:\n  ", paste(invalid.id.values, collapse=", "), sep="")))
    } else {}

    # iterate though 'id' variables, see if both sets share the same cases
    id.differences1 <- unique(unlist(sapply(id.values, function(this.val){
        return(which(!set1[,this.val] %in% set2[,this.val]))
      })))
    id.differences2 <- unique(unlist(sapply(id.values, function(this.val){
        return(which(!set2[,this.val] %in% set1[,this.val]))
      })))
    if(length(id.differences1) > 0 | length(id.differences2) > 0){
      if(length(id.differences1) > 0){
        set1 <- set1[-id.differences1,]
        id.missings1 <- set1[id.differences1,id.values]
        id.miss.txt1 <- sapply(1:nrow(id.missings1), function(this.row){
            paste(id.missings1[this.row,], collapse=" ")
          })
        missing.message1 <- paste("Missing cases in ", set2.name, ":\n    ", paste(id.miss.txt1, collapse="\n    "), sep="")
      } else {
        missing.message1 <- ""
      }
      if(length(id.differences2) > 0){
        set2 <- set2[-id.differences2,]
        id.missings2 <- set2[id.differences2,id.values]
        id.miss.txt2 <- sapply(1:nrow(id.missings2), function(this.row){
            paste(id.missings2[this.row,], collapse=" ")
          })
        missing.message2 <- paste("Missing cases in ", set1.name, ":\n    ", paste(id.miss.txt2, collapse="\n    "), sep="")
      } else {
        missing.message2 <- ""
      }
      warning(paste("\nThe objects differ in the number of oservations.  Missing ones have been ignored in both:\n  ", missing.message1, missing.message2, sep=""), call.=FALSE)
    } else {}

    ## prepare sets for comparison
    # first create subsets, if necessary
    if(!is.null(select) | !is.null(ignore)){
      if(!is.null(select)){
        set1 <- subset(set1, select=select)
        set2 <- subset(set2, select=select)
        subset.varname <- paste(names(set1), collapse=", ")
        warning(paste("Only a subset was compared. Remaining columns:\n  ", paste(subset.varname, collapse=", "), sep=""), call.=FALSE)
      } else{
        subset.varname <- ""
      }
      if(!is.null(ignore)){
        set1 <- subset(set1, select=names(set1)[!names(set1) %in% ignore])
        set2 <- subset(set2, select=names(set2)[!names(set2) %in% ignore])
        warning(paste("Only a subset was compared. Dropped columns:\n  ", paste(ignore, collapse=", "), sep=""), call.=FALSE)
      } else{
        ignore.varname <- ""
      }
      # quit if the subset is identical
      if(identical(set1, set2)){
        message(paste("\nThe compared objects (",set1.name," & ",set2.name,") are identical!\n\n", sep=""))
        if(isTRUE(new.set)){
          # do we need renaming?
          if(length(rename) > 0){
            set1 <- rename.vars(rename=rename, set=set1)
          } else {}
          return(set1)
        } else {
          return(invisible(NULL))
        }
      } else {}
    } else {}
    # then sort the data according to "No"
    set1 <- set1[order(set1$No),]
    set2 <- set2[order(set2$No),]
    ## done with preparations here

    ## the actual comparison
    # create an array with indices of differences
    set.diff.array <- rbind(which(set1 != set2, arr.ind=TRUE), which(is.na(set1) != is.na(set2), arr.ind=TRUE))
    # total number of differences, needed for the sapply() call below
    uneq <- dim(set.diff.array)[1]
    # results are sorted by column, but we want them sorted by row
    set.diff.array <- set.diff.array[order(set.diff.array[,1]),]

    # the main act!
    differences <- t(
      sapply(1:uneq, function(x){
          if(uneq > 1) {
            idx <- set.diff.array[x,]
          } else {
            idx <- set.diff.array
          }
          varname <- names(set1)[idx[2]]
          # what values are present in each set?
          value1 <- set1[idx[1],idx[2]]
          value2 <- set2[idx[1],idx[2]]
          # find 'id' to that case
          result <- data.frame(paste(set1[idx[1],id[[id.names[1]]]], collapse=" "), stringsAsFactors=FALSE)
          colnames(result)[1] <- id.names[1]
          for(this.id in id.names[-1]){
            result[[this.id]] <- paste(set1[idx[1],id[[this.id]]], collapse=" ")
          }
          result$Item <- varname
          result$Set1 <- value1
          result$Set2 <- value2
          return(result)
        } # end function(x)
      ) # end sapply()
    ) # end t()

    # set the names of the given sets, to better understand the outcome
    colnames(differences)[4] <- set1.name
    colnames(differences)[5] <- set2.name
    rownames(differences) <- c(1:dim(differences)[1])

    # check if a new set should be created, with dubious values replaced by NA
    if(isTRUE(new.set)){
      new.df <- set1
      # there's only rows if there's also more than one found difference
      if(uneq > 1){
        for (this.val in 1:nrow(set.diff.array)){
          new.df[set.diff.array[this.val, 1], set.diff.array[this.val, 2]] <- NA
        }
      } else if(uneq == 1){
        new.df[set.diff.array[1, 1], set.diff.array[1, 2]] <- NA
      } else {
        # no differences, return the set
        return(new.df)
      }
      cat("\nThe objects differ.\n\n")
      print(differences)
      return(new.df)
    } else {
      cat("\nThe objects differ.\n\n")
      results <- data.frame(differences, stringsAsFactors=FALSE)
      # unlist all cols
      for (this.col in 1:ncol(results)){
        results[,this.col] <- unlist(results[,this.col])
      }
      return(results)
    }
  }
}
