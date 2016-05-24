## Copyright (C) 2012 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 2 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


#### The real plotmath rendering engine is
#### ~/R/D/r-devel/R/src/main/plotmath.c (3237 lines as of 2013-01 !!)

## From ~/R/D/r-devel/R/src/main/plotmath.c  [and edited]
Ugreek <- lapply(c.Ugreek <- c(
     "Alpha",
     "Beta",
     "Chi",
     "Delta",
     "Epsilon",
     "Phi",
     "Gamma",
     "Eta",
     "Iota",
     "theta1",
     "vartheta",
     "Kappa",
     "Lambda",
     "Mu",
     "Nu",
     "Omicron",
     "Pi",
     "Theta",
     "Rho",
     "Sigma",
     "Tau",
     "Upsilon",
     "sigma1",
     "varsigma",
     "stigma",
     "Omega",
     "Xi",
     "Psi",
     "Zeta"
    ), as.symbol)

Lgreek <- lapply(c.Lgreek <- c(
     "alpha",
     "beta",
     "chi",
     "delta",
     "epsilon",
     "phi",
     "gamma",
     "eta",
     "iota",
     "phi1",
     "varphi",
     "kappa",
     "lambda",
     "mu",
     "nu",
     "omicron",
     "pi",
     "theta",
     "rho",
     "sigma",
     "tau",
     "upsilon",
     "omega1",
     "omega",
     "xi",
     "psi",
     "zeta"
     ), as.symbol)

## BinTable <- lapply(c.BinTable <- c(
BinTable <- c(
    "*" ~ " ", # ! math multiplication : just "a" space
    "^" ~ "^", ## <- added {not here in plotmath}
    "+" ~ "+",
    "-" ~ "-",
    "/" ~ "/",
    ":" ~ ":",
    "%+-%" ~ "\\pm",
    "%*%" ~ "\\times",
    "%/%" ~ "\\div",
    "%intersection%" ~ "\\cap",
    "%union%" ~ "\\cup",
    "%.%" ~ "\\cdot" ##  cdot or dotmath
    )

RelTable <- c(
    "<" ~ "<",
    "==" ~ "=",
    "=" ~ "=",
    ">" ~ ">",
    "%=~%" ~ "\\cong",    	## congruent
    "!=" ~ "\\ne",
    "<=" ~ "\\le",
    ">=" ~ "\\ge",
    "%==%" ~ "\\equiv",    	## equivalence
    "%~~%" ~ "\\approx",    	## approxequal
    "%prop%" ~ "\\propto", ## proportional to
    "%~%" ~ "\\sim",       ## distributed as {only from R 3.0.0 !!}
    ## Arrows :
    "%<->%" ~ "\\leftrightarrow",
    "%<-%" ~ "\\leftarrow",
    "%up%" ~ "\\uparrow",
    "%->%" ~ "\\rightarrow",
    "%down%" ~ "\\downarrow",
    "%<=>%" ~ "\\Leftrightarrow",
    "%<=%"  ~ "\\Leftarrow",
    "%dblup%" ~ "\\Uparrow",
    "%=>%" ~ "\\Rightarrow",
    "%dbldown%" ~ "\\Downarrow",
    ## TeX set symbols
    "%supset%" ~ "\\supset",
    "%supseteq%" ~ "\\supseteq",
    "%notsubset%" ~ "\\notsubset",
    "%subset%" ~ "\\subset",
    "%subseteq%" ~ "\\subseteq",
    "%in%" ~ "\\in",
    "%notin%" ~ "\\not\\in" #not standard LaTeX: "\\nin"
    )

c.BinTable <- vapply(BinTable, `[[`, "", 2)
c.RelTable <- vapply(RelTable, `[[`, "", 2)

AccentTable <- lapply(c.AccentTable <- c(
    "hat",
    "ring",
    "tilde",
    "dot"
    ), as.symbol)

c.delimTab <- c("(","[","{", "[[") ## the right-version does not appear in plotmath {I think}



isTab <- function(tab2)
    is.list(tab2) && all(vapply(tab2, class, "") == "formula")

if(FALSE) ## unused here
tab2sym <- function(tab2) {
    stopifnot(isTab(tab2))
    lapply(lapply(tab2, `[[`, 2L), as.symbol)
}

getTab <- function(d.a, tab2) {
    stopifnot(isTab(tab2))
    pos <- which(d.a == vapply(tab2, `[[`, "", 2L))
    if(identical(r <- tab2[[pos]][[3L]], quote(NY))) "NOT YET" else r
}


##' deparse() with
mDeparse <- function(ex, ..., collapse="")
    paste(sub("^ +", '', sub(" +$", "", deparse(ex, width.cutoff=500, ...))),
          collapse=collapse)


##' Is 'ex' *not* a "regular R name":
##' @param ex language, symbol, expression, or character
isOp <- function(ex)
    is.na(match(substr(ex,1,1), c(".", letters,LETTERS, 0:9)))

renderAtom <- function(a, Len, d.a = mDeparse(a))
{
### TODO MM Len is not used anymore
    stopifnot(is.numeric(Len), Len >= 1)
    switch(d.a,
           "~" = "\\ ", # space
           ## otherwise:
           if(d.a %in% c(c.Ugreek,c.Lgreek)) paste0("\\", d.a)
           ## these may not be needed here:
           ## else if(Len > 1 && d.a %in% c.BinTable) getTab(d.a, BinTable)
           ## else if(Len > 1 && d.a %in% c.RelTable) getTab(d.a, RelTable)
           ## these may not be needed here:
           else if(d.a %in% c.BinTable) getTab(d.a, BinTable)
           else if(d.a %in% c.RelTable) getTab(d.a, RelTable)
           else d.a)
}

if(FALSE) ## unused here
isBinary <- function(a, d.a = mDeparse(a))
{
    d.a %in% c.BinTable ||
    d.a %in% c.RelTable
}

expr2latex <- function(expr) {
    L <- length(expr)
    if(!L) "" else {
        Symb <- is.symbol(expr)
        F <- if(Symb) expr else expr[[1]]
        cF <- mDeparse(F)
        FF <- renderAtom(F, Len=L, d.a = cF)
        if(Symb && L != 1)
            stop("is.symbol(.), but length(.) = ", L, " != 1")
        else if(!Symb && typeof(expr) != "language" && L != 1)
            stop("is not language nor symbol), but length(.) = ", L, " != 1")
        switch(L,
               ## length 1:
               FF,

           { ## length 2: e.g.  "- 1", "+ x", "!TRUE",  "~ ff",
               rhs <- expr2latex(expr[[2]])
               if       (cF == "bold") paste0("\\mathbf{", rhs, "}")
               else if(cF == "italic") paste0("\\mathit{", rhs, "}")
               else if(!isOp(cF)) # not a binary operator ==> "function call":
                   paste0(FF,"(",rhs,")") ## e.g. "O(n)"
               else if(cF == "{") paste0("{", rhs, "}")
               else if(cF == "(") paste0("(", rhs, ")")
               else paste(FF, rhs)
           },

           { ## length 3:
               lhs <- expr2latex(expr[[2]])
               rhs <- expr2latex(expr[[3]])
               if(cF == "[") ## subscript
                   paste0(lhs, "_{", rhs, "}")
               else if(cF == "~") ## space
                   paste(lhs, "\\", rhs)
               ## not treated, as plotmath() does neither :
               ## else if(cF == "[[")
               ##     paste0(lhs, "[[", rhs, "]]")
               else if(cF %in% c.BinTable)
                   paste(lhs, getTab(cF, BinTable), rhs)
               else if(cF %in% c.RelTable)
                   paste(lhs, getTab(cF, RelTable), rhs)
	       else if(isOp(cF)) ## e.g.   U + x
                   paste(lhs, FF, rhs)
               else ## log(x, 2)
                   paste0(FF, "(", lhs, ",", rhs, ")")
           },

               ## length >=4 : F(a, b, c, ...)
               stop("length(expr) = ",L," (>= 4);  not yet implemented") # TODO MM
               )## end{switch}
    }
}## { end } expr2latex


## original escape_latex from fortunes:::toLatex.fortune
escapeLatex <- function(x) {
    x <- gsub("\\\\ ", "\\textbackslash\\ ", x, fixed = TRUE)
    x <- gsub("\\\\", "\\textbackslash ", x, fixed = TRUE)
    x <- gsub("\\n", "\\textbackslash n", x, fixed = TRUE)
    x <- gsub("| ", "\\textbar\\ ", x, fixed = TRUE)
    x <- gsub("|", "\\textbar", x, fixed = TRUE)
    x <- gsub("#", "\\#", x, fixed = TRUE)
    x <- gsub("$", "\\$", x, fixed = TRUE)
    x <- gsub("&", "\\&", x, fixed = TRUE)
    x <- gsub("~ ", "\\textasciitilde\\ ", x, fixed = TRUE)
    x <- gsub("~", "\\textasciitilde ", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("^", "\\verb|^|", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x <- gsub("{", "\\{", x, fixed = TRUE)
    x <- gsub("}", "\\}", x, fixed = TRUE)
    x <- gsub(" '", " `", x, fixed = TRUE)
    x <- gsub(" \"", " ``", x, fixed = TRUE)
    x <- gsub("... ", "\\dots\\ ", x, fixed = TRUE)
    x <- gsub("...",  "\\dots", x, fixed = TRUE)
    x <- gsub(" - ", " -- ", x, fixed = TRUE)
    x
} ## {escape_latex}
