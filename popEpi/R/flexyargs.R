


#' @title Flexible Variable Usage in \pkg{popEpi} Functions
#' @author Joonas Miettinen
#' @name flexible_argument
#' @aliases flexible_popEpi_arguments, flexible_arg,
#' flexible_args, flexible_arguments
#' @description Certain arguments in \pkg{popEpi} can be passed in multiple 
#' ways. This document shows the usage and a pitfall in the
#' usage of such flexible arguments.
#' 
#' @details 
#' 
#' Flexible arguments in \pkg{popEpi} are used to pass variables existing
#' in your data or in the environment where the function is used 
#' (for everyday users this is the global environment - in simple terms,
#' where your data is / your work space). The flexible arguments
#' are modeled after the \code{by} argument in \code{data.tables} - 
#' see \code{?data.table}. There are many ways to supply the same information
#' to certain functions in \pkg{popEpi}, but the possible ways listed below
#' may be limited in some of them to only allow for using only a part of them.
#' 
#' @section Everyday usage:
#' 
#' Most commonly you may pass
#' variable names as character strings, e.g.
#' 
#' \code{FUN(arg = c("V1", "V2"), data = x)}
#' 
#' which may be stored in advance:
#' 
#' \code{vars <- c("V1", "V2")}
#' 
#' \code{FUN(arg = vars, data = x)}
#' 
#' where \code{x} contains those variables. You may also supply variable
#' names as symbols:
#' 
#' \code{FUN(arg = V1, data = x)}
#' 
#' Or as a list of symbols (similarly to as in \code{\link{aggregate}}):
#' 
#' \code{FUN(arg = list(V1, V2), data = x)}
#' 
#' Or as a list of expressions:
#' 
#' \code{FUN(arg = list(V1 + 1, factor(V2)), data = x)}
#' 
#' A formula without a left-hand-side specified is sometimes allowed as well:
#' 
#' \code{FUN(arg = ~ I(V1 + 1) + factor(V2), data = x)}
#' 
#' Using a symbol or a list of symbols/expressions typically
#' causes the function to look for the variable(s)
#' first in the supplied data (if any) and then where the function was called.
#' For everyday users this means you might define e.g.
#' 
#' \code{V3 <- factor(letters)}
#' 
#' and do e.g.
#' 
#' \code{FUN(arg = list(V1 + 1, factor(V2), V3), data = x)}
#' 
#' provided \code{V1} and \code{V2} exist in \code{x} or in the function calling
#' environment.
#' 
#' @section A pitfall:
#' 
#' There is one way to use flexible arguments incorrectly: By supplying
#' the name of a variable which exists both in the supplied data
#' and the calling environment, and intending the latter to be used. E.g.
#' 
#' \code{vars <- c("V2")}
#' 
#' \code{FUN(arg = V3, data = x)}
#' 
#' where \code{x} has a column named \code{vars}. This causes the function to
#' use \code{x$vars} and NOT \code{x$V2}.
#' 
#' @section Advanced:
#' 
#' Function programmers are adviced to pass character strings
#' whenever possible. To fool-proof against conflicts as described in the
#' section above, refer to the calling environment explicitly when
#' passing the variable containing the character strings:
#' 
#' \code{TF <- environment() ## current env to refer to}
#' 
#' \code{vars <- c("V1", "V2")}
#' 
#' \code{FUN(arg = TF$vars, data = x)}
#' 
#' Even if \code{x} has columns named \code{vars} and \code{TF}, 
#' using \code{TF$vars} does not use those columns but only evaluates
#' \code{TF$vars}
#' in the calling environment. This is made possible by the fact
#' that data is always passed as a \code{data.frame}, within which evaluation
#' of expressions using the dollar operator is not possible. Therefore
#' it is safe to assume the data should not be used. However, lists of 
#' expressions will not be checked for dollar use and will fail in conflict
#' situations:
#' 
#' \code{TF <- environment() ## current env to refer to}
#' 
#' \code{vars <- letters[1:5]}
#' 
#' \code{x <- data.frame(vars = 1:5, TF = 5:1, V1 = 10:6)}
#' 
#' \code{FUN(arg = list(TF$vars, V1), data = x)}
#' 
#' On the other hand you may typically also pass quoted (\code{\link{quote}})
#' or substituted \code{\link{substitute}} expressions etc., where
#' the \code{env$object} trick will work as well:
#' 
#' \code{q <- quote(list(vars, V1))}
#' 
#' \code{FUN(arg = TF$q, data = x)}
#' 
#' This works even with
#' 
#' \code{a <- 1:5}
#' 
#' \code{V1 <- quote(TF$a)}
#' 
#' \code{FUN(arg = TF$V1, data = x)}
#' 
#' So no conflicts should occur.
#' @examples 
#' 
#' data(sire)
#' ## prepare data for e.g. 5-year "period analysis" for 2008-2012
#' ## note: sire is a simulated cohort integrated into popEpi.
#' BL <- list(fot=seq(0, 5, by = 1/12))
#' x <- lexpand(sire, birth = bi_date, entry = dg_date, exit = ex_date,
#'              status = status %in% 1:2,
#'              breaks = BL)
#'               
#' x <- aggre(x, by = fot)
#'
#' ## silly example of referring to pyrs data by fixed character string;
#' ## its possible that the real name wont be fixed in a real-life application.
#' pyrs <- "actual_pyrs"  
#' TF <- environment()
#' x$actual_pyrs <- as.numeric(x$pyrs)
#' x$pyrs <- 1
#'              
#' ## this works (uses actual_pyrs eventually)
#' st <- survtab_ag(fot ~ 1, data = x, surv.type = "surv.obs",
#'                  pyrs = TF$pyrs, d = from0to1, 
#'                  surv.method = "hazard")
#' ## this would be wrong (sees expression 'pyrs' and uses that column,
#' ## which is not what is intended here)
#' st <- survtab_ag(fot ~ 1, data = x, surv.type = "surv.obs",
#'                  pyrs = pyrs, d = from0to1,
#'                  surv.method = "hazard")

NULL









