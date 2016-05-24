{############################################################################### 
# base.patterns.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 7/30/2012
# 
# DESCRIPTION
# ===========
# Base Patterns for building off.
# 
# LICENSE
# ========
# lint is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# lint is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' @name base-patterns
#' @rdname base-patterns
#' @title Base Patterns
#' @description
#'  Use these Perl regular expressions to help build pattern based styles.
#' @exportPattern .*\\.characters
#' @exportPattern .*\\.constant
NULL

if(getRversion() < "2.15.0") {
    paste0 <- function(..., collapse = NULL) {
        paste(..., collapse=collapse, sep = '')
    }
}


#' @rdname base-patterns
start.characters <- "[\\p{L}\\.]"

#' @rdname base-patterns
following.characters <- "[\\p{L}\\p{N}\\._\\d]"

#' @rdname base-patterns
name.pattern <- sprintf("%s%s*", start.characters, following.characters)

#' @rdname base-patterns
real.constant <- c("\\d+\\.?\\d*", "\\d*\\.?\\d+")

#' @rdname base-patterns
exp.constant <- "\\d+\\.?\\d*[eE][+-]?\\d+"

#' @rdname base-patterns
int.constant <- "\\dL"

#' @rdname base-patterns
complex.constant <- sprintf("%si", real.constant)

#' @rdname base-patterns
numeric.all.constant <- c(real.constant, exp.constant, int.constant
                          , complex.constant)

#' @rdname base-patterns
numeric.constant <- paste('(', paste(numeric.all.constant, collapse='|'), ')')


.no.exclude <- character(0)
rx.lb.no.opp <- '(?<![\\-+*<>=!%])'
rx.la.no.opp <- '(?![\\-+*<>=!%]'
arith.opp <- c( '+'  = '\\+'
              , '*'  = '(?<![*])\\*(?![*])'
              , '/'  = '\\/'
              , '^'  = '\\^'
              , '-'  = '(?<![<])(-)(?![>])'
              , '**' = '\\*\\*'
              )
logical.opp <- c(  '|' = '(?<![|])\\|(?![|])'
                ,  '<' = '(?<![<-])(>)(?![>=])'
                ,  '>' = '(?<![<])(<)(?![>=-])'
                ,  '&' = '(?<![&])&(?![&])'
                , '||' = '\\|\\|'
                , '&&' = '&&'
                , '<=' = '<='
                , '==' = '=='
                , '!=' = '!='
                , '>=' = '>='
                )
assign.opp  <- c( '='   = '(?<![<>=!])(=)(?!=)'
                , '<-'  = '(?<![<])(<-)'
                , '->'  = '(->)(?![>])'
                , '<<-' = '<<-'
                , '->>' = '->>'
                )
special.opp <- c('%[^%]*%')
all.opp    <- c(arith.opp, logical.opp, assign.opp, special.opp)
infix.noeq <- setdiff(all.opp, all.opp['='])
no.lead.rx <- "(&<[^\\s\\^\\-!%+*/<>=\\|&)"
no.preceeding.space.rx <- "(?<=[^\\s])"
no.trailing.space.rx   <- "(?=[^\\s])"
no.preceeding.percent  <- "(?<!%)"
no.trailing.percent    <- "(?!%)"

