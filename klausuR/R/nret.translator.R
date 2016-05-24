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


#' Convert NRET/ET data between klausuR and other software
#'
#' This function should help to interchange answer data between R and other statistical software packages -- especially
#' SPSS, but it's probably useful for other products as well.
#'
#' \code{\link[klausuR:klausur]{klausur}} expects data in a special format if it should be evaluated according to (Number Right) Elimination
#' Testing (NRET/ET), only one variable per item. Other software products might not be able to process this rather
#' condensed format. In that case, you will most likely need several variables for each item, i.e. one per answer alternative.
#' Adding to that, the coding of answers is by default done with "+", "-", "0" and "*" in \code{klausuR}, again a
#' solution that might confuse other products.
#'
#' This function translates data in both directions, and does also convert vectors giving the correct answer. The latter
#' will turn a \code{klausuR} type answer string into a number indicating the correct alternative (and the other way round).
#' This means that it will only work if there's exactly one valid answer to each item. If you convert towards SPSS,
#' the resulting list will also include SPSS syntax to define variables respectively.
#'
#' @note The conversion is done on an object basis, that is, \code{nret.translator} will not open or write files,
#' but take and return R objects. The function should ignore any other columns/variables in the object.
#'
#' @param dat A data.frame, the object to convert.
#' @param items Optional vector defining the columns to convert. If \code{NULL}, the function will try to autodetect
#'    Items: \code{klausuR} type items are expected to be named \code{"ItemXXX"}, with XXX indicating the item number,
#'    SPSS type items \code{"itemXXXaYY"},  with XXX indicating the item number and YY the number of the answer alternative.
#' @param spss Either "in" or "out", depending on the direction of conversion.
#' @param corr Logical. Set to \code{TRUE} if \code{dat} is a vector with the correct answers. If \code{corr=TRUE} and
#'    \code{spss="in"}, you must also set \code{num.alt} accordingly!
#' @param num.alt A numeric value definig the number of answer alternatives for each item. Can be a vector, if items
#'    have different numbers of options. If it is shorter than the number of items, it will be repeated for all items.
#' @param klausuR.alt A named vector defining the codes for \code{klausuR} type of answers.
#' @param spss.alt A named vector defining the codes for SPSS type of answers.
#' @param rm.old.vars Logical. If \code{TRUE}, the converted columns will not be returned. Only relevant if \code{corr=FALSE}.
#' @param items.only Logical. If \code{TRUE}, only the converted columns will be returned. Only relevant if \code{corr=FALSE}.
#' @param klausuR.prefix A named character vector with two optional elements, \code{item} and \code{corr}, defining the name prefix
#'    used for the items in the test data and the vector with correct answers, respectively. Defaults to \code{item="Item"} and \code{corr="Item"}.
#' @param spss.prefix Like \code{klausuR.prefix}, but for the SPSS data. Defaults to \code{item="item"} and \code{corr="corr"}.
#' @return If \code{corr=FALSE}, a data.frame with more or less columns (depending on \code{rm.old.vars} and \code{items.only}).
#'    If \code{corr=TRUE}, returns a named vector if \code{spss="in"} and a list if \code{spss="out"} (containing SPSS syntax
#'    in the element \code{syntax} and also a named vector, called \code{answ}).
#' @author m.eik michalke \email{meik.michalke@@uni-duesseldorf.de}
#' @seealso \code{\link[klausuR:klausur]{klausur}}
#' @keywords misc
#' @export
#' @examples
#' \dontrun{
#' # from SPSS to R
#' data(spss.data)
#' klausuR.data <- nret.translator(spss.data, spss="in")
#' spss.corr <- c(
#'    item01=2, item02=3, item03=3, item04=3, item05=2,
#'    item06=2, item07=3, item08=1, item09=1, item10=2)
#' klausuR.corr <- nret.translator(spss.corr, spss="in", corr=TRUE, num.alt=3)
#'
#' # from R to SPSS
#' spss.data <- nret.translator(klausuR.data)
#' spss.corr <- nret.translator(klausuR.corr, corr=TRUE, num.alt=3)
#' # if you find the syntax useful
#' cat(spss.corr$syntax, file="~/somewhere/NRET.sps")
#' }

nret.translator <- function(dat, items=NULL, spss="out", corr=FALSE, num.alt=NULL,
  klausuR.alt=c(is.true="+", is.false="-", missing="0", err="*"),
  spss.alt=c(is.true="2", is.false="1", missing="0", err="3"),
  rm.old.vars=TRUE, items.only=FALSE, klausuR.prefix=c(), spss.prefix=c()){

  klausuR.prefix <- check.prefixes(prefixes=klausuR.prefix, package="klausuR")
  spss.prefix    <- check.prefixes(prefixes=spss.prefix, package="SPSS")

  # local copy of data:
  dat.orig <- dat

  # function to translate alternatives
  trans.alt <- function(data, items, alt.in, alt.out){
    data[,items] <- sapply(items, function(item){
        for (check in c("is.true", "is.false", "missing", "err")){
          data[data[,item] == alt.in[check],item] <- alt.out[check]
        }
        return(data[,item])
      })
    return(data)
  }

  if(identical(spss, "out")){
    ## spss="out"
    if(isTRUE(corr)){
      if(!is.vector(dat)){
        stop(simpleError("If 'corr=TRUE', 'dat' must be a vector!"))
      } else {}
      num.items <- length(dat)
      # transform correct answers
      # our SPSS syntax knows only one information: number of the true alternative
      # this only works if there's only one true answer!
      new.answers <- sapply(dat, function(item){
          answ.parts <- unlist(strsplit(item, split=""))
          # find correct one
          true.one <- which(answ.parts == klausuR.alt["is.true"])
          # more than one true answer? we can't handle that
          if(length(true.one) > 1){
            stop(simpleError("Sorry, can't handle items with more than one correct answer :-("))
          } else {
            return(true.one)
          }
        })
      old.digits <- as.numeric(gsub(paste("^(", klausuR.prefix[["corr"]], ")([[:digit:]]{1,3})$", sep=""), "\\2", names(dat), perl=TRUE))
      new.names <- paste(spss.prefix[["corr"]], old.digits, sep="")
      names(new.answers) <- new.names
      # for convenience, create some SPSS syntax
      SPSS.compute <- sapply(1:num.items, function(item.idx){
          cmp.line <- paste("COMPUTE ", new.names[item.idx], "=", new.answers[item.idx], ".\n", sep="")
        })
      SPSS.syntax <- paste("NUMERIC ", new.names[1]," TO ", new.names[num.items]," (F2.0).\n",
        paste(SPSS.compute, collapse=""), "EXECUTE.\n", sep="")
      results <- list(syntax=SPSS.syntax, answ=new.answers)
    } else {
      # transform data
      # extract items (ItemXXX)
      # get columns with items
      if(is.null(items)){
        items.idx <- grep(paste("^(", klausuR.prefix[["item"]], ")([[:digit:]]{1,3})$", sep=""), names(dat))
      } else {
        items.idx <- items
      }
      # just in case item names are in the wrong order
      items.idx <- items.idx[order(names(dat[, items.idx]))]

      # as a precaution, replace NAs with valid "missing" values
      dat[is.na(dat)] <- klausuR.alt["missing"]

      # split answers into alternatives
      old.item.names <- names(dat[, items.idx])
      num.items <- length(old.item.names)
      old.digits <- as.numeric(gsub(paste("^(", klausuR.prefix[["item"]], ")([[:digit:]]{1,3})$", sep=""), "\\2", old.item.names, perl=TRUE))
      # if the data is ok, the numbers of alternatives in the first row
      # should suffice to get them for all subjects
      num.alternatives <- nchar(dat[1,items.idx])
      new.item.names <- as.vector(sapply(1:num.items, function(item.idx){
          item.pre <- if(num.items < 10){
            paste(spss.prefix[["item"]], old.digits[item.idx], sep="")
            } else if(num.items < 100){
              paste(spss.prefix[["item"]], sprintf("%02d", old.digits[item.idx]), sep="")
            } else {
              paste(spss.prefix[["item"]], sprintf("%03d", old.digits[item.idx]), sep="")
            }

          n.items <- if(max(num.alternatives) < 10){
            paste(item.pre, "a", c(1:num.alternatives[item.idx]), sep="")
            } else {
              paste(item.pre, "a", sprintf("%02d", c(1:num.alternatives[item.idx])), sep="")
            }
          return(n.items)
        }))
      # split data, item-wise
      new.dat <- matrix(nrow=dim(dat)[[1]])[,-1]
      for(item in old.item.names){
          new.dat <- cbind(new.dat, t(as.data.frame(strsplit(dat[,item], split=""), stringsAsFactors=FALSE)))
        }

      # rename the results (itemXaY (X=Item, Y=Antwortalternative))
      dimnames(new.dat) <- list(NULL, new.item.names)

      # translate alternatives
      new.dat <- trans.alt(data=new.dat, items=new.item.names, alt.in=klausuR.alt, alt.out=spss.alt)
      # re-order by new names
      new.dat <- new.dat[,order(names(new.dat))]
      
      # output R object
      if(isTRUE(items.only)){
        results <- as.data.frame(new.dat, stringsAsFactors=FALSE)
      } else if(isTRUE(rm.old.vars)){
        results <- cbind(dat[,-items.idx], new.dat)
      } else {
        results <- cbind(dat.orig, new.dat)
      }
    } # end !corr
  } else if(identical(spss, "in")){
  
    ## spss="in"
    if(isTRUE(corr)){
      if(!is.vector(dat)){
        stop(simpleError("If 'corr=TRUE', 'dat' must be a vector!"))
      } else {}
      num.items <- length(dat)
      # transform correct answers
      if(is.null(num.alt)){
        stop(simpleError("You must specify the number of answer alternatives!"))
      } else {
        if(!identical(num.items, length(num.alt))){
          # check if num.alt can be replicated sanely
          stopifnot(num.items %% length(num.alt) == 0)
          num.fctr <- num.items/length(num.alt)
          num.alt <- rep(num.alt, num.fctr)
        } else {}
      }
      # create the new answer vector
      new.answers <- sapply(1:num.items, function(item.idx){
          answ <- rep(klausuR.alt["is.false"], num.alt[item.idx])
          answ[dat[item.idx]] <- klausuR.alt["is.true"]
          answ <- paste(answ, collapse="")
          return(answ)
        })
      old.digits <- as.numeric(gsub(paste("^(", spss.prefix[["corr"]], ")([[:digit:]]{1,3})$", sep=""), "\\2", names(dat), perl=TRUE))
      new.item.names <- gen.item.names(old.digits, prefix=klausuR.prefix[["corr"]])
      names(new.answers) <- new.item.names
      results <- new.answers
    } else {
      # extract items (itemXaY (X=Item, Y=Antwortalternative))
      # get columns with items
      if(is.null(items)){
        items.idx <- grep(paste("^(", spss.prefix[["item"]], ")([[:digit:]]{1,3}a([[:digit:]]{1,2}))$", sep=""), names(dat))
      } else {
        items.idx <- items
      }
      # just in case item names are in the wrong order
      items.idx <- items.idx[order(names(dat[, items.idx]))]

      # this holds just the "itemXX" prefix
      items.pre <- gsub(paste("^((", spss.prefix[["item"]], ")[[:digit:]]{1,3})(a[[:digit:]]{1,2})$", sep=""), "\\1", names(dat[, items.idx]), perl=TRUE)
      # unique item names
      item.names <- unique(items.pre)

      # as a precaution, replace NAs with valid "missing" values
      dat[is.na(dat)] <- spss.alt["missing"]

      # translate alternatives
      dat <- trans.alt(data=dat, items=items.idx, alt.in=spss.alt, alt.out=klausuR.alt)

      # combine alternatives into single items
      new.items <- sapply(item.names, function(item){
          alts.of.item <- grep(paste(item,"a", sep=""), names(dat[,items.idx]), perl=TRUE)
          # combine values row by row
          n.item.rows <- sapply(1:dim(dat)[[1]], function(row.num){
              item.combi <- paste(dat[,items.idx][row.num,alts.of.item], collapse="")
              return(item.combi)
            })
          return(n.item.rows)
        })

      # rename the results (ItemXXX)
      old.digits <- as.numeric(gsub(paste("^(", spss.prefix[["item"]], ")([[:digit:]]{1,3})$", sep=""), "\\2", item.names, perl=TRUE))
      new.item.names <- gen.item.names(old.digits, prefix=klausuR.prefix[["item"]])
      dimnames(new.items) <- list(NULL,new.item.names)
      new.items <- as.data.frame(new.items, stringsAsFactors=FALSE)
      # re-order by new names
      new.items <- new.items[,order(names(new.items))]

      # output R object
      if(isTRUE(items.only)){
        results <- new.items
      } else if(isTRUE(rm.old.vars)){
        results <- cbind(dat[,-items.idx], new.items)
      } else {
        results <- cbind(dat.orig, new.items)
      }
    } #end !corr
  } else {
    stop(simpleError(paste("Wrong value for 'spss':", spss)))
  }

  return(results)
}
