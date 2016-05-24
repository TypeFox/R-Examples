# Copyright 2009-2015 Meik Michalke <meik.michalke@hhu.de>
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


# these are internal functions that are being called by some of the methods of klausuR
# they are not exported, hence not to be called by users themselves
# and are therefore only documented by the comments in this file.

## debug.message()
# small tool to communicate values of objects in test runs
debug.message <- function(x){
    obj.name <- deparse(substitute(x))
    message(paste0(obj.name, ": ", paste0(x, collapse=", ")))
    return(invisible(NULL))
}

## check.prefixes()
# used to define variable prefixes for items in data and correct answers vector
check.prefixes <- function(prefixes=c(), package="klausuR"){
  if("item" %in% names(prefixes)){
    item <- as.character(prefixes["item"])
  } else {
    if(identical(package, "SPSS")){
      item <- "item"
    } else {
      item <- "Item"
    }
  }
  if("corr" %in% names(prefixes)){
    corr <- as.character(prefixes["corr"])
  } else {
    if(identical(package, "SPSS")){
      corr <- "corr"
    } else {
      corr  <- "Item"
    }
  }
  return(c(item=item, corr=corr))
}
## end check.prefixes()


## data.check.klausur()
# this function is called by klausur() and klausur.mufo()
# for some sanity checks of the given data
data.check.klausur <- function(answ, corr, items, na.rm, prefixes=c(), keep.cases=NULL, recode.na=0){

    # check for var names to use
    prefix <- check.prefixes(prefixes=prefixes, package="klausuR")

    # in case no items were specified, take variables of names "Item##" as items
    if(is.null(items)){
      items <- grep(paste("^(item)([[:digit:]]{1,3})$", sep=""), names(answ), ignore.case=TRUE)
    } else{}

    # are all needed variables present in the answers data?
    needed.names <- c("No", "Name", "FirstName", "MatrNo")
    if(any(!needed.names %in% names(answ))){
        missing.vars <- needed.names[!needed.names %in% names(answ)]
        stop(simpleError(paste("Missing variables in observation data:\n ", paste(missing.vars, collapse=", "))))
    } else {}

    # we'll check for NAs only in variables we will use, so define them here
    relevant.items <- c(grep("^Name$|^FirstName$|^MatrNo$", names(answ)), items)

    # before we check for missing values, let's see if there are cases which should prevail
    if(!is.null(keep.cases)){
      missings.to.keep <- is.na(answ[answ[["MatrNo"]]==keep.cases,items])
      answ[answ[["MatrNo"]]==keep.cases,items[missings.to.keep]] <- recode.na
    } else {}
    # are there missing values in answ, at least in the relevant parts?
    if(sum(is.na(answ[, relevant.items]) > 0)){
      invalid.cases.row <- sort(unique(unlist(sapply(relevant.items, function(na.var){
          return(which(is.na(answ[,na.var])))
        }))))
      invalid.cases <- unlist(sapply(invalid.cases.row, function(this.case){
          case.summary <- paste("  MatrNo ", answ[this.case, "MatrNo"], " (", paste(names(answ)[is.na(answ[this.case,])], collapse=", "), ")", sep="")
        }))
      if(isTRUE(na.rm)){
        warning(paste("NAs were present in '",deparse(substitute(answ)),"' and cases have been removed:\n",
          paste(invalid.cases, collapse="\n"), sep=""), call.=FALSE)
        answ <- answ[-invalid.cases.row,]
#         for (na.var in relevant.items){
#           answ <- answ[!is.na(answ[, na.var]),]
#         }
      } else {
        warning(paste("NAs were present in '",deparse(substitute(answ)),"':\n", paste(invalid.cases, collapse=", "), sep=""), call.=FALSE)
      }
    } else{}
    if(sum(is.na(corr) > 0)){
      stop(simpleError(paste("NAs present in '",deparse(substitute(corr)),"'!\n", sep="")))
    } else{}

    # now let's check wheter all defined correct answers match the variables in answers data
    missing.items <- !names(corr) %in% names(answ[, items])
    missing.corr <- !names(answ[, items]) %in% names(corr)
    if(any(c(missing.items, missing.corr))){
      fehl.items.corr <- names(corr)[missing.corr]
      fehl.items.answ <- names(answ)[missing.items]
      stop(simpleError(paste("Please check:\n  ", paste(fehl.items.corr, fehl.items.answ, collapse=", "),
      "\n  The number of items differs between observed and correct answers!", sep="")))
    }

    # SPSS has a habbit of adding trailing spaces to character data. read.spss() cannot
    # always cope with this, so check items for factors and return stripped strings
    for (this.item in relevant.items){
      if(is.factor(answ[, this.item]) || is.character(answ[, this.item])){
        answ[, this.item] <- gsub("[[:space:]]+$", "", as.character(answ[, this.item]))
      } else {}
    }

  # return objects that have probably changed
  checked.data <- list(answ=answ, items=items)
  return(checked.data)
} ## end data.check.klausur()

## scoring.check.klausur()
# this function is called by klausur() and klausur.mufo()
# for some sanity checks of the desired scoring
scoring.check.klausur <- function(corr, marks, wght, score, maxp=NULL){

    # are there missing values?
    if(!is.null(marks)){
      if(sum(is.na(marks) > 0)){
        stop(simpleError(paste("NAs present in ",deparse(substitute(marks)),"!\n", sep="")))
      } else{}
    } else{}

    if(is.null(wght)){
      message("No weight vector (wght) given. Maximum score is determined from items and scoring method.")
    } else{
      if(length(wght) != length(corr)){
      stop(simpleError("The number of weights differs from the number if items!"))
      } else{}
    }
    if(!identical(score, "solved")){
      if(!score %in% c("partial", "liberal", "pick-n", "NR", "ET", "NRET", "NRET+")){
        stop(simpleError("Invalid value for score, must be either \"solved\", \"partial\", \"liberal\", \"pick-n\", \"NR\", \"ET\", \"NRET\", or \"NRET+\"!"))
      } else if(identical(score, "partial")){
        warning("Partially answered items were allowed (but only if no wrong alternative was checked).", call.=FALSE)
      } else if(score %in% c("liberal", "pick-n")){
        warning("Partially answered items were allowed (wrong alternatives were ignored but didn't invalidate a whole answer).", call.=FALSE)
      } else if(identical(score, "ET")){
        warning("Partially answered items were allowed (scored according to Elimination Testing).", call.=FALSE)
      } else if(score %in% c("NRET", "NRET+")){
        warning("Partially answered items were allowed (scored according to Number Right Elimination Testing).", call.=FALSE)
      } else {}
    } else {}

    if(!is.null(maxp) && !is.numeric(maxp)){
      stop(simpleError("'maxp' must either be numeric or NULL!"))
    } else {}

  return(TRUE)
} ## end scoring.check.klausur()

## gen.item.names()
# an internal function to generate item names
# takes the number of items to be generated
gen.item.names <- function(num, prefix="Item"){
  if(length(num) == 1){
    # NULL if no valid item number given
    if(num < 1){
      return(NULL)
    } else{}
    # currently, 999 items are the theoretical limit
    if(num >= 1000){
      return(NULL)
    } else{}
    item.num <- c(1:num)
  } else {
    # if num is a vector, create item names from that vector
    item.num <- num
    num <- length(item.num)
  }

    items <- if(num < 10){
      paste(prefix, item.num, sep="")
    } else {
      if(num < 100){
        paste(prefix, sprintf("%02d", item.num), sep="")
      } else {
        paste(prefix, sprintf("%03d", item.num), sep="")
      }
    }
    return(items)
  } ## end gen.item.names()

## calc.cronbach.alpha()
# calculates cronbachs alpha, needs a matrix with dichotomous data
# this function uses the package "psychometric"!
calc.cronbach.alpha <- function(dichot.matrix){
  try.cron.alpha <- function(){
    cr.alpha <- alpha(dichot.matrix)
    cr.alpha.ci <- alpha.CI(cr.alpha, k=dim(dichot.matrix)[2], N=dim(dichot.matrix)[1], level=0.95)
    cr.deleted <- as.matrix(sapply(1:dim(dichot.matrix)[2], function(x){alpha(dichot.matrix[-x])}))
  rownames(cr.deleted) <- names(dichot.matrix)
  colnames(cr.deleted) <- "alphaIfDeleted"
    return(list(alpha=cr.alpha, ci=cr.alpha.ci, deleted=cr.deleted))
  }
  cron.alpha.list <- tryCatch(try.cron.alpha(), error=function(e){
  warning("Cronbach's alpha calculation failed!\n", e, call.=FALSE)
  return(list(alpha=NA, ci=NA, deleted=NA))})
  } ## end function calc.cronbach.alpha()

## calc.item.analysis()
# performes basic item analysis, like discrimatory power etc.
# this function uses the package "psychometric"!
calc.item.analysis <- function(dichot.matrix, cron.alpha.list){
  try.item.analysis <- function(){
    item.anal <- cbind(item.exam(dichot.matrix, discrim=TRUE), alphaIfDeleted=cron.alpha.list$deleted)
    # add selection index (selektionskennwert) as suggested by lienert
    item.anal <- cbind(item.anal, selIdx=item.anal[["Item.Tot.woi"]] / (2 * sqrt(item.anal[["Difficulty"]] * (1 - item.anal[["Difficulty"]]))))
    return(item.anal)
  }
  item.analyse <- tryCatch(try.item.analysis(), error=function(e){
  warning("Item analysis failed!\n", e, call.=FALSE)
  return(data.frame(NA))})
  } ## end function calc.item.analysis()

## global.results()
# glues together parts of results into an object
global.results <- function(answ, points, maxp, mark, minp=0){
    if(!is.null(answ$No)){
      results <- data.frame(  No=answ$No,
          Name=answ$Name)
    } else {
      results <- data.frame(Name=answ$Name)
    }
    # write all desired information into the data object
    results$FirstName <- answ$FirstName
    results$MatrNo    <- answ$MatrNo
    results$Points    <- points
    results$Percent   <- round(100*((points-minp)/(maxp-minp)), digits=1)
    results$Mark      <- mark
    # if pseudonyms were given, include them for anonymous feedback
    if(!is.null(answ$Pseudonym)){
      results$Pseudonym <- answ$Pseudonym
    } else{}
  return(results)
} ## end global.results()

## anon.results()
anon.results <- function(glob.res){
    if(!is.null(glob.res$Pseudonym)){
      results <- data.frame(  Pseudonym=glob.res$Pseudonym,
         Points=glob.res$Points,
        Percent=glob.res$Percent,
           Mark=glob.res$Mark)
    } else {
      results <- data.frame(Anonymous="(no pseudonyms for anonymous feedback available)")
    }
  return(results)
} ## end anon.results()

## answ.alternatives()
# splits items with multiple answer alternatives into its parts
# the result can be a list with every alternative,
# or, if latex=TRUE, one object with all alternatives separated by commas,
# which is used for nicer looking output in klausur.report()
answ.alternatives <- function(answ, latex=FALSE){
  pasteAlt <- function(x){
    if(identical(x, character()) | isTRUE(is.na(x))){
      return("{}")
    } else {}
    return(paste(x, collapse=", "))
  }
  # divide all answers into its parts
  answ.parts <- strsplit(as.character(answ), "")
  names(answ.parts) <- names(answ)
  if(isTRUE(latex)){
    if(length(answ.parts) > 1){
      answ.parts <- lapply(answ.parts, pasteAlt)
    } else {
      answ.parts <- pasteAlt(answ.parts[[1]])
    }
  } else{}

  return(answ.parts)
} ## end answ.alternatives()


## find.partial()
# is being called by partial(), see below
find.partial <- function(item.answ, corr, item, wght=NULL, mode="absolute", strict=TRUE, pickN=FALSE, digits=NULL){
  # divide all correct answers into their parts
  corr.parts <- answ.alternatives(corr)
  answ.parts <- lapply(item.answ, answ.alternatives)
#debug.message(corr.parts)
#debug.message(answ.parts)

  # how many correct answers are there?
  corr.length <- length(corr.parts[item][[1]])
  if(corr.length > 1 | (isTRUE(pickN) & length(unlist(answ.parts)) > 1)){
    result.list <- lapply(answ.parts[[1]],
        function(x){
          # count only if no more answers were checked than correct answers available
          if(!isTRUE(pickN) && length(x) > corr.length){
            return(0)
          } else {
            abs.correct <- !is.na(pmatch(x, corr.parts[item][[1]]))
            abs.false   <- sum(is.na(pmatch(x, corr.parts[item][[1]])))
            # if in strict mode, discard if more answers were checked then correct one,
            # that is, if at least one wrong answer was given, return no points at all
            if(isTRUE(strict) && abs.false > 0){
              return(0)
            } else {
              return(abs.correct)
            }
          }
        })
    # this corresponds to the "absolute" value
    result <- unlist(lapply(result.list, sum))
    # if we want the percentage instead:
    if(identical(mode, "percent")){
      result <- (result * wght)/corr.length
      if(!is.null(digits)){
        result <- round(result, digits=digits)
      }
    } else{}
  } else {
    # this corresponds to the "absolute" value
    # one exception: in pickN items with missing correct answers (i.e., all
    # alternatives are in fact wrong) must *not* be counted here!
    if(isTRUE(pickN) & (identical(as.character(corr[item]), "") | is.na(corr[item]) | is.null(corr[item]))){
      result <- rep(0, length(unlist(item.answ)))
    } else {
      result <- as.numeric(item.answ == corr[item])
#debug.message(item.answ)
#debug.message(corr[item])
    }
    # percentage is irrelevant for dichotomous items,
    # but there might be a weight vector
    if(identical(mode, "percent")){
      result <- result * wght
    } else{}
  }
  return(result)
} ## end find.partial()


## partial()
# this function computes results of one item including partially correct answers
# if multiple correct answer alternatives are possible.
# item.answ must be a vector of given answers to one item with the item name given
#
# mode can be
# - "absolute" (number of correct alternatives, ignores weights -- default)
# - "percent"  (percent of correct alternatives, can be combined with wght to weight the result)
#
# if strict=TRUE, only answers are counted if *no* wrong alternative was checked at all
partial <- function(item.answ, corr, wght=NULL, mode="absolute", strict=TRUE, pickN=FALSE, digits=NULL){
  # check for partially correct answers
  # firstly, extract the item name from the answ vector
  item <- names(item.answ)
  if(length(item) == 0){
    stop(simpleError("Partial results wanted, but incorrect item answers given: no item name defined!\n"))
  } else{}

  if(is.null(wght)){
    wght <- 1
  } else{}

  part.results <- find.partial(item.answ=item.answ, corr=corr, item=item, wght=wght, mode=mode,
    strict=strict, pickN=pickN, digits=digits)

  return(part.results)
} ## end partial()


## pickN()
# participants get 1/k points for each correctly checked answer as well as for each
# correctly *un*checked distractor. if a distractor is falsely checked, no points are given.
# calls partial() for both correct and wrong answers in absolute mode
pickN <- function(item.answ, corr, wrong, wght=NULL, mode="percent", digits=NULL){
  right.alternatives <- sapply(answ.alternatives(corr), length)
  wrong.alternatives <- sapply(answ.alternatives(wrong), length)
  all.alternatives <- right.alternatives + wrong.alternatives

  partial.right <- partial(item.answ=item.answ, corr=corr, wght=1, mode="absolute", strict=FALSE, pickN=TRUE)
#debug.message(partial.right)
  # count how many wrong alternatives were falsely chosen
  partial.wrong <- partial(item.answ=item.answ, corr=wrong, wght=1, mode="absolute", strict=FALSE, pickN=TRUE)
  item <- names(item.answ)
#debug.message(partial.wrong)

  points.absolute <- partial.right + wrong.alternatives[item] - partial.wrong
  
  if(is.null(wght)){
    wght <- 1
  } else{}

  if(identical(mode, "absolute")){
    results <- points.absolute * wght
  } else {
    results <- points.absolute * wght / all.alternatives[item]
    if(!is.null(digits)){
      results <- round(results, digits=digits)
    }
  }
  
  return(results)
} ## end pickN()


## function nret.score()
# as alternative scoring functions, nret.score() implements three modes:
#  - NR: traditional number right
#      c(true.pos=1, false.pos=0, true.neg=0, false.neg=0, miss=0)
#  - ET: elimination testing (strike wrong alternatives)
#      c(true.pos=0, false.pos=0, true.neg=1, false.neg=1-num.alt, miss=0)
#  - NRET: number right elimination testing (strike wrong alternatives, mark one as true)
#      c(true.pos=1, false.pos=0, true.neg=1, false.neg=1-num.alt, miss=0)
#  - NRET+: number right elimination testing (strike wrong alternatives, mark one as true; 0 points if too many "+")
#      c(true.pos=1, false.pos=0, true.neg=1, false.neg=1-num.alt, miss=0)
#
# errors (like both alternatives marked) will be evaluated like missings, pointwise
#
# answ: given answer to one item (e.g. "-+--")
# corr: the correct pattern (like answ)
# num.alt: number of answer alternatives, counted automatically if NULL
nret.score <- function(answ, corr, score="NRET", is.true="+", is.false="-", missing="0", err="*",
  num.alt=NULL, true.false=FALSE) {

  # count answer alternatives
  if(is.null(num.alt)){
    num.alt <- nchar(corr)
  } else if(!is.numeric(num.alt)){
    stop(simpleError("Value of \"num.alt\" must be NULL or a number!"))
  }

  # in which mode will be scored?
  if(identical(score, "NR")){
    mtx <- c(true.pos=1, false.pos=0, true.neg=0, false.neg=0, miss=0, err.true=1, err.false=0)
  } else if(identical(score, "ET")){
    mtx <- c(true.pos=0, false.pos=0, true.neg=1, false.neg=as.numeric(1-num.alt), miss=0, err.true=as.numeric(1-num.alt), err.false=1)
  } else if(identical(score, "NRET")){
    mtx <- c(true.pos=1, false.pos=0, true.neg=1, false.neg=as.numeric(1-num.alt), miss=0, err.true=as.numeric(2-num.alt), err.false=1)
  } else if(identical(score, "NRET+")){
    mtx <- c(true.pos=1, false.pos=0, true.neg=1, false.neg=as.numeric(1-num.alt), miss=0)
  } else {
    stop(simpleError(paste("Unknown scoring mode:", score)))
  }

  all.results <- lapply(answ, function(curr.answ){
    # first split character vectors into atomic vectors
    answ.split <- unlist(strsplit(curr.answ, split=""))
    corr.split <- unlist(strsplit(corr, split=""))
    # check for equal length
    if(length(answ.split) != length(corr.split)){
      stop(simpleError("Given and correct answers are of unequal length!"))
    } else {}
    # check for correct input
    if(sum(!answ.split %in% c(is.true, is.false, missing, err)) > 0){
      stop(simpleError("Given answer vector includes invalid characters!"))
    } else {}
    if(sum(!corr.split %in% c(is.true, is.false, missing, err)) > 0){
      stop(simpleError("Correct answer vector includes invalid characters!"))
    } else {}

    ## plausibility checks
    # too many "correct" answers?
    num.yeses <- sum(answ.split %in% is.true)
    num.trues <- sum(corr.split %in% is.true)
      # if more than one correct answer is defined, the penalty for
      # marking one "wrong" must be aligned
      if(num.trues > 1 && score %in% c("ET","NRET","NRET+")){
        mtx["false.neg"] <- (1-num.alt)/num.trues
      } else {}
    if(identical(score, "NRET+") && num.yeses > num.trues && !isTRUE(true.false)){
      result <- 0
    } else {
      # then compare answer by answer
      points <- sapply(1:length(corr.split), function(idx){
          answ.given <- as.character(answ.split[idx])
          answ.crrct <- as.character(corr.split[idx])
          # score NR will only return TRUE or FALSE, let's take care of that first
          if(identical(score, "NR") & isTRUE(true.false) & identical(answ.crrct, is.true)){
            if(identical(answ.given, is.true)){
              return(TRUE)
            } else {
              return(FALSE)
            }
          } else if(identical(score, "NR") & isTRUE(true.false) & !identical(answ.crrct, is.true)){
            return()
          } else {}
          
          if(identical(answ.given, answ.crrct)){
            if(identical(answ.given, is.true)){
              # this is a true positive
              if(isTRUE(true.false)){
                  return("P")
              } else {
                return(mtx["true.pos"])
              }
            } else if(identical(answ.given, is.false)){
              # this is a true negative
              if(isTRUE(true.false)){
                  return("N")
              } else {
                return(mtx["true.neg"])
              }
            } else {
              # this is impossible...
              stop(simpleError("Are you sure your answer vector is correct?!"))
            }
          } else {
            if(identical(answ.given, is.true)){
              # this is a false positive
              if(isTRUE(true.false)){
                  return("p")
              } else {
                return(mtx["false.pos"])
              }
            } else if(identical(answ.given, is.false)){
              # this is a false negative
              if(isTRUE(true.false)){
                  return("n")
              } else {
                return(mtx["false.neg"])
              }
            } else if(answ.given %in% c(missing, err)){
              if(isTRUE(true.false)){
                  return(as.character(answ.given))
              } else {
                if(score %in% c("NR","ET","NRET")){
                  # the authors didn't discuss failed answers, so this is by the book
                  if(identical(answ.given, missing)){
                    return(mtx["miss"])
                  } else {
                    # ok, it's an error. did it happen to a wrong or right alternative?
                    if(identical(answ.crrct, is.true)){
                      return(mtx["err.true"])
                    } else {
                      return(mtx["err.false"])
                    }
                  }
                } else {
                  return(mtx["miss"])
                }
              }
            }
          }
        }
      )

      if(isTRUE(true.false)){
        # return true/false indicators, e.g. "PNNn0"
        if(identical(score, "NR")){
          result <- points
        } else {
          result <- paste(points, collapse="")
        }
      } else {
        result <- sum(points)
      }

    }

    return(result)
    })
    
  ## currently, no negative points are valid for mark assignments
  # so to be sure, we'll globally add num.alt-1 points, so 0 is the minimum
  ## if this gets changed, take care of the calculation of max. points as well!!!
  if(isTRUE(true.false) || identical(score, "NR")){
    all.results <- unlist(all.results)
  } else {
    all.results <- unlist(all.results) + (num.alt-1)
  }

  return(all.results)
} ## end function nret.score()


## function nret.minmax()
# compute minimum/maximum points, number of alternatives and baseline
nret.minmax <- function(corr, score="NRET", is.true="+", is.false="-", quiet=FALSE){
  # initial value for minimum score
  min.score <- 0
  # split whole answer vector
  corr.split <- unlist(strsplit(as.character(corr), split=""))
  num.trues <- sum(corr.split %in% is.true)
  num.false <- sum(corr.split %in% is.false)

  if(score %in% c("ET", "NRET", "NRET+")){
    # these need some special treatment, because there can be more points than items
    ## currently, no negative points are valid for mark assignments
    # so to be sure, we'll globally add num.alt-1 points, so 0 is the minimum
    # see also the klausur.gen.marks function below, since this had to be cosidered there, too!
    #
    # get all alternatives
    num.alt.all <- nchar(corr)
    # check if they're all of equal length
    if(all(num.alt.all == num.alt.all[1])){
      num.alt <- as.numeric(num.alt.all[1])
    } else {
      num.alt <- max(num.alt.all)
      min.score <- sum(num.alt - num.alt.all)
      if(!isTRUE(quiet)){
        warning(paste("Items differ in number of answer alternatives: ", min(num.alt.all), "-", num.alt,
          "\n  Took the maximum (", num.alt, ") to determine additive constant to avoid negative points.",
          "\n  In effect, the lowest achievable score is ", min.score ," points.", sep=""), call.=FALSE)
      } else {}
    }
    # compute baseline, that is, what do you get with all missings?
    baseline <- length(num.alt.all) * (num.alt-1)
    if(!isTRUE(quiet)){
      warning(paste("The baseline (all missings) used for solved percentage is ", baseline ," points.", sep=""), call.=FALSE)
    } else {}

    if(score %in% c("NRET", "NRET+")){
      maxp <- num.trues + num.false + (length(corr) * (num.alt-1))
    }  else {
      maxp <- num.false + (length(corr) * (num.alt-1))
    }
  } else if(identical(score, "NR")){
    # in case no weights were given, count each item as one point
    maxp <- num.trues
    num.alt <- NULL
    baseline <- 0
    min.score <- 0
  } else {
    # in case this is no (NR)ET data at all
    maxp <- length(corr)
    num.alt <- NULL
    baseline <- 0
    min.score <- 0
  }

  results <- c(maxp=maxp, minp=min.score, baseline=baseline, num.alt=num.alt)
  return(results)
} ## end function nret.minmax()

## marks.summary()
# this function takes a vector with marks and returns a summarising matrix
# with the effective ranges of points and percentage for each mark defined
marks.summary <- function(marks, minp=0, add.const=0){
  # since 0 doesn't get counted, adjust add.const values below 0
  marks.levels <- levels(as.factor(marks))
  maxp <- length(marks) + add.const
  marks.matrix <- sapply(marks.levels, function(x){
    mark.min <- max(c(minp, min(which(marks == x)))) + add.const
    mark.max <- max(which(marks == x)) + add.const
    # avoid strange starting points above zero
    if((mark.min-minp) <= 1 & ((add.const >= 0 & minp >= 0) | add.const == minp)){
      mark.min.pct <- ""
    } else {
      mark.min.pct <- ceiling((mark.min-minp) / (maxp-minp) * 100)
    }
    mark.max.pct <- ceiling((mark.max-minp) / (maxp-minp) * 100)
    # if it's only one point value, don't display a range
    if(mark.min == mark.max){
      m.f.points <- sprintf("%3s",mark.min)
      m.f.pct <- sprintf("%3s",mark.min.pct)
    } else {
      m.f.points <- paste(sprintf("%3s", mark.min), " -- ", sprintf("%3s",mark.max), sep="")
      m.f.pct <- paste(sprintf("%3s", mark.min.pct), " < ", sprintf("%3s",mark.max.pct), sep="")
    }

    mark.frame <- c(
      Points=m.f.points,
      Percent=m.f.pct
      )
    return(mark.frame)
    }
  )
  return(t(marks.matrix))
} ## end marks.summary()


## function distrct.analysis()
# answ: a data.frame containing all items in columns and all answers in its rows
# corr: named vector with the correct answers (names indicate items); ignored if NULL
# points: a matrix with columns "MatrNo" and one for each item, listing the results in points, already weighted!
# results: a data.frame with two columns, "MatrNo" and "Points"; ignored if NULL
# partWHole: logical, wheter part-whole correction should be applied; since this is tricky to implement and
#   interpret for scoring functions other than "NR" and "solved", it can be turned off
distrct.analysis <- function(answ, corr=NULL, points=NULL, results=NULL, partWhole=FALSE){
  # set MatrNo NULL to fulfill CRAN check's needs
  MatrNo <- NULL
  # check if MatrNo needs to be stripped off
  if("MatrNo" %in% names(answ)){
    if(!is.null(results)){
      if(any(!identical(sort(answ[["MatrNo"]]), sort(results[["MatrNo"]])), !identical(sort(answ[["MatrNo"]]), sort(points[["MatrNo"]])))){
        stop(simpleError("Distractor analysis: MatrNo in given items doesn't match those in the given points/results!"))
      } else {}
      if(!"Points" %in% names(results)){
        stop(simpleError("Distractor analysis: Missing column \"Points\" in given object!"))
      } else {}
      # "MatrNo" in all data sets are identical, so we'll sort both by it
      answ <- answ[order(answ[["MatrNo"]]),]
      results.partWhole <- results <- results[order(results[["MatrNo"]]),]
      points <- points[order(points[["MatrNo"]]),]
    } else {}
    # create an object without "MatrNo"
    answNoMatr <- subset(answ, select=-MatrNo)
  } else {
    results <- points <- NULL
  }

  if(!is.null(corr)){
    if(!identical(names(answNoMatr), names(corr))){
      stop(simpleError("Distractor analysis: Item names (answ) don't match those in the given vector with correct answers!"))
    } else {}
    if(!identical(names(subset(points, select=-MatrNo)), names(corr))){
      stop(simpleError("Distractor analysis: Item names (points) don't match those in the given vector with correct answers!"))
    } else {}
  } else {}

  results <- lapply(names(answNoMatr), function(thisItem){
      selected.absolute <- summary(as.factor(answNoMatr[[thisItem]]))
      selected.percent <- selected.absolute / sum(selected.absolute) * 100
      selected.all <- data.frame(
        answer=names(selected.absolute),
        absolute=selected.absolute,
        percent=selected.percent,
        correct="",
        discrim=NA,
        points=NA,
        stringsAsFactors=FALSE)
      # add a cloumn as indicator for correct answer
      if(!is.null(corr)){
        correct <- as.character(corr[[thisItem]])
        # was the correct alternative checked at all?
        if(correct %in% selected.all[["answer"]]){
          selected.all[selected.all[["answer"]] == correct, "correct"] <- "*"
        } else {}
      } else {}
      # compute correlations with the end result, if available,
      # as an indicator for discriminate power
      if(!is.null(results) && !is.null(points) && nrow(selected.all) > 1){
        selected.all[["discrim"]] <- sapply(selected.all[["answer"]], function(thisAnswer){
            # return(suppressWarnings(polyserial(results[["Points"]], answ[[thisItem]] == thisAnswer)))
            # point-biserial correlations seems to be the proper one here, which is equivalent to pearson
            return(suppressWarnings(cor(results[["Points"]], answ[[thisItem]] == thisAnswer, method="pearson")))
          })
        selected.all[["points"]] <- sapply(selected.all[["answer"]], function(thisAnswer){
            return(suppressWarnings(mean(results[answ[[thisItem]] == thisAnswer, "Points"])))
          })
        if(isTRUE(partWhole)){
          indiv.points <- points[[thisItem]]
          results.partWhole[["Points"]] <- results.partWhole[["Points"]] - indiv.points
          selected.all[["discrim.partWhole"]] <- sapply(selected.all[["answer"]], function(thisAnswer){
              return(suppressWarnings(cor(results.partWhole[["Points"]], answ[[thisItem]] == thisAnswer, method="pearson")))
            })
          selected.all[["points.partWhole"]] <- sapply(selected.all[["answer"]], function(thisAnswer){
              return(suppressWarnings(mean(results.partWhole[answ[[thisItem]] == thisAnswer, "Points"])))
            })
        } else {}
      } else {}
      rownames(selected.all) <- NULL
      return(selected.all)
    })

  names(results) <- names(answNoMatr)

  return(results)
} ## end function distrct.analysis()


## klausur.reorderItems()
# put items in correct order (multiple test forms)
klausur.reorderItems <- function(slot, order){
  # to avoid NOTEs from R CMD check:
  MatrNo <- NULL

  part.slot <- subset(slot, select=-MatrNo)
  slot.names <- names(part.slot)
  matn.slot <- subset(slot, select=MatrNo)
  # now let's reorder the stuff
  part.slot.reordered <- part.slot[,order]
  names(part.slot.reordered) <- slot.names
  # finally glue MatrNo back
  reordered.items <- cbind(matn.slot, part.slot.reordered)
  return(reordered.items)
} ## end klausur.reorderItems()

## function plot.merger()
# a hack to produce combined histograms from several objects of class klausuR
# takes a list of klausuR objects and returns a combined results slot
plot.merger <- function(klsr=list()){
  # first combine all result slots into one data.frame
  k.merged.results <- data.frame()
  for(this.klsr in klsr){
    stopifnot(inherits(this.klsr, "klausuR"))
    k.merged.results <- rbind(k.merged.results, this.klsr@results)
  }
  return(k.merged.results)
} ## end function plot.merger()


## function latex.umlaute()
# this function will replace German umlauts with LaTeX equivalents
# some sanitizing is also done
# it's used in tabellenbau() of klausur.report()
latex.umlaute <- function(input){
  output <- gsub("(\\\\?)&", '\\\\&', as.character(input), perl=TRUE)
  output <- gsub("(\\\\?)_", '\\\\_{}', as.character(output), perl=TRUE)
  output <- gsub("(\\\\?)#", '\\\\#', as.character(output), perl=TRUE)
  output <- gsub("(\\\\?)\\^", '\\\\^{}', as.character(output), perl=TRUE)
  output <- encoded_text_to_latex(enc2utf8(output), "utf8")
  return(output)
} ## end function latex.umlaute()


## function file.umlaute()
# this function will replace German umlauts and other special chars
# for filenames it's used in tabellenbau() of klausur.report()
file.umlaute <- function(input){
  output <- gsub("\u00C0|\u00C1|\u00C2|\u00C3|\u00C5","A",as.character(input)) # À Á Â Ã Å
  output <- gsub("\u00C4|\u00C6","Ae",as.character(output)) # Ä Æ
  output <- gsub("\u00C7","C",as.character(output)) # Ç  c3 87  LATIN CAPITAL LETTER C WITH CEDILLA
  output <- gsub("\u00C8|\u00C9|\u00CA|\u00CB|\u00D0","E",as.character(output)) # È É Ê Ë Ð (Eth)
  output <- gsub("\u00CC|\u00CD|\u00CE|\u00CF","I",as.character(output)) # Ì Í Î Ï
  output <- gsub("\u00D1","N",as.character(output)) # Ñ
  output <- gsub("\u00D2|\u00D3|\u00D4|\u00D5","",as.character(output)) # Ò Ó Ô Õ
  output <- gsub("\u00D6|\u00D8","Oe",as.character(output)) # Ö Ø
  output <- gsub("\u00D9|\u00DA|\u00DB","U",as.character(output)) # Ù Ú Û
  output <- gsub("\u00DC","Ue",as.character(output)) # Ü
  output <- gsub("\u00DD","Y",as.character(output)) # Ý
  output <- gsub("\u00DE","Th",as.character(output)) # Þ
  output <- gsub("\u00DF","ss",as.character(output)) # ß
  output <- gsub("\u00E0|\u00E1|\u00E2|\u00E3|\u00E5","a",as.character(output)) # à á â ã å
  output <- gsub("\u00E4|\u00E6","ae",as.character(output)) # ä æ
  output <- gsub("\u00E7","c",as.character(output)) # ç
  output <- gsub("\u00E8|\u00E9|\u00EA|\u00EB|\u00F0","e",as.character(output)) # è é ê ë ð (eth)
  output <- gsub("\u00EC|\u00ED|\u00EE|\u00EF","i",as.character(output)) # ì í î ï
  output <- gsub("\u00F1","n",as.character(output)) # ñ  c3 b1  LATIN SMALL LETTER N WITH TILDE
  output <- gsub("\u00F2|\u00F3|\u00F4|\u00F5","o",as.character(output)) # ò ó ô õ
  output <- gsub("\u00F6|\u00F8","oe",as.character(output)) # ö ø
  output <- gsub("\u00F9|\u00FA|\u00FB","u",as.character(output)) # ù ú û
  output <- gsub("\u00FC","ue",as.character(output)) # ü
  output <- gsub("\u00FD|\u00FF","y",as.character(output)) # ý ÿ
  output <- gsub("\u00FE","th",as.character(output)) # þ
  output <- gsub("\\(|\\)|\\[|\\]|\\.|\\*|\\#|\\'", "", as.character(output))
  return(output)
} ## end function file.umlaute()


## function axis.breaks()
# makes sure that the range of breaks for histograms' labels captures all values
axis.breaks <- function(min, max){
  min <- floor(min) - 1
  max <- ceiling(max)
  breaks <- c(min:max)
  if(min(breaks) > min){
    min <- min - 1
    breaks <- axis.breaks(min, max)
  } else {}
  
  if(max(breaks) < max){
    max <- max + 1
    breaks <- axis.breaks(min, max)
  } else {}
  return(list(breaks=breaks, min=min, max=max))
}
## end function axis.breaks()
