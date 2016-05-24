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


### 'Poor-man's approach' how to create lines for a LaTeX table

## Now relies on R (>= 3.0.0)'s  format.ftable():

##' Just a version of format.ftable(), returning our "matrix + attr" object:
fftable <- function(x, lsep = " | ", quote = FALSE, method = "compact", ...)
{
    stopifnot(inherits(x, "ftable"))
    structure(format(x, method=method, quote=quote, lsep=lsep, ...),
	      ncv = length(attr(x, "col.vars")),
              ## (method=="non.compact" || method=="col.compact"),
              ## => this would ensure the midrule to come after the names of
              ##    the columns containing row names, but complicates other
              ##    issues such as introducing cmidrules of length 1 as the
              ##    head is then one line longer...
	      nrv = length(attr(x, "row.vars")))
}

##' @title Create and cat rows of a LaTeX table from a given matrix
##' @param x character or numeric matrix
##' @param rsep character string inserted at the end of each row
##' @param csep character string for separating different cells in a row
##' @param include.rownames logical indicating whether row names are included
##'        in the first column
##' @return nothing; LaTeX table rows are cat'ed
##' @author Marius Hofert
cattablines <- function(x, rsep = "\\\\", csep = " & ", include.rownames = TRUE) {
    stopifnot(is.matrix(x))
    if(include.rownames) {
        x <- cbind(rownames(x), x)
        rownames(x) <- NULL
    }
    cat(paste0(apply(x, 1L, function(row) paste0(row, collapse=csep)), rsep),
	sep="\n") # table content
}

##' @title Ingredients for Converting an ftable to a LaTeX Table
##' @param x character matrix with attributes 'nrv' (number or row variables) and
##'	   'ncv' (number of column variables)
##' @param align either a character vector...
##'	  - ... of length > 1: e.g., c("c", "c", "c", "S[table-format=1.2]")
##'	  - ... of length 1: e.g., c("*{3}{c} S[table-format=1.2]")
##'	  or NULL in which case a useful default is constructed ('r' for all table
##'       entries, 'l' for all columns of row names)
##' @param booktabs logical indicating whether a booktabs-compliant LaTeX table
##'	   is created (requires \usepackage{booktabs} in the preamble)
##' @param head either
##'	   - a character vector containing the lines of the header
##'	   - NA (= no header)
##'	   - NULL (in which case a default is constructed)
##' @param rsep character string inserted at the end of each row
##' @param sp numeric scaling factor for separating blocks of rows
##' @param rsep.sp numeric vector of length equal to the number of different
##'	   groups or row variables minus 1, giving the spaces (interpreted as pt)
##'	   between the groups
##' @param csep character string for separating different cells in a row
##' @param quote see ?format.ftable
##' @return list with components
##'	    body     : character vector of lines of the table body
##'	    body.raw : character matrix of cells of the table body
##'	    head     : character vector of lines of the table header
##'	    head.raw : character matrix of cells of the table header
##'	    align    : alignment string
##'	    rsepcol  : character vector of row separators (last entries of each row)
##' @author Marius Hofert
##' { tablines }
tablines <- function(x, align = NULL, booktabs = TRUE,
		     head = NULL,
		     rsep = "\\\\", sp = if(booktabs) 3 else 1.25, rsep.sp = NULL,
		     csep = " & ", quote = FALSE)
{
    ## checks
    nrv <- attr(x, "nrv")
    ncv <- attr(x, "ncv")
    stopifnot(is.matrix(x), is.character(x),
	      is.numeric(nrv), length(nrv) == 1,
	      is.numeric(ncv), length(ncv) == 1,
	      is.null(align) || is.character(align))
    nr <- max(1, nrv)
    nc <- max(1, ncv)

    if(is.null(rsep.sp)) ## default for rsep.sp
	rsep.sp <- sp*rev(seq_len(nr-1))
    else stopifnot(length(rsep.sp) == nr-1)

    ## trim leading and trailing whitespace:
    x2 <- gsub("^\\s+|\\s+$", "", x)

    ## remove header
    head.raw <- x2[ seq_len(nc), ,drop=FALSE] # original head
    x3	     <- x2[-seq_len(nc), ,drop=FALSE] # body.raw; table without head
    if(all(!nzchar(x3[1,]))) x3 <- x3[-1, ,drop=FALSE] # remove possible empty first row (for the case where nrv=0 and booktabs=FALSE)

    ## determine row separators (easy with stripped whitespace and removed col names)
    x3. <- x3 # x3. will just be a convenience dummy here
    if(nrow(x3.) > 0 && nrv > 0) x3.[1, 1:nrv] <- "" # remove names of columns containing row names (if there are any such columns); this is for the algorithm determining the row separators to work correctly
    rsepcol <- rep(rsep, nrow(x3.)) # last column will contain the row
				      # separators (partly overwritten below)
    if(length(rsep.sp)) {
	rsepvec <- paste0(rsep,
			  c(paste0(if(booktabs) " \\addlinespace[" else "\\\\[",
				   rsep.sp, "pt]"), ""))
	for(j in rev(seq_len(nr))) { # walk from right to left through all (except the
			  # last) column containing row labels
	    labs <- x3.[,j] != "" # determine in which row the formatted
				    # ftable has labels
	    stopifnot((il <- which(labs) - 1L) > 0)
	    rsepcol[il] <- rsepvec[j] # set the corresponding
					# separator for all rows except those
					# containing column header;
					# overwritten by subsequent j's
	}
    } # don't append the column of row separators yet

    ## determine alignment of table entries if not given
    ncl <- ncol(x3) # align has to be a character vector of this length
    align <- if(is.null(align)) {
        if(nrv == 0) paste0("l*{", ncl-1, "}{r}") # takes care of the 'wide' format with no row names
                                        # (1st col [containing col names] should be left aligned then)
        else paste0("*{", nrv, "}{l}", # cols with row names will be left aligned
                    "*{", ncl-nrv, "}{r}") # all other cols will be right aligned
    } else {
        ## OLD version:
        ## if(length(align) == 1)
        ##     paste0("*{", ncl, "}{", align, "}")
        ## else {
        ##     stopifnot(length(align) == ncl)
        ##     paste(align, collapse="")
        ## }
        stopifnot(length(align) == 1)
        paste(align, collapse="") # => can be a vector or already a string (typically the latter)
    }
    nrh <- nrow(head.raw) # number of rows of the head
    if(!booktabs && nrh > 1) align <- paste0("@{\\extracolsep{0.6em}}", align) # for \cline's to be separated

    ## determine head if not given (in booktabs style: everything strictly between '\toprule' and '\midrule')
    if(is.null(head)) { # construct default header (lines)
	nch <- ncol(head.raw) # number of columns of the head
	## function to determine title rules
	titlerule <- function(x, booktabs) {
	    stopifnot(nr >= 1)
	    ok <- (nchar(x) > 0)[-seq_len(nr)] # which columns without labels contain non-empty strings
	    beg <- nr + which(ok) # start positions
	    end <- c(beg[-1]-1, nch) # end positions
	    rr <- paste0(if(booktabs) "\\cmidrule(lr){" else "\\cline{",
			 beg, "-", end, "}", collapse=" ")
	    if(booktabs) rr else paste(rr, "\\noalign{\\smallskip}") # put in more vertical space
	}
	## function to determine multicolumns
	## (precisely for all headers of columns containing no row variables)
	multicolumn <- function(x, m) { ##' m: #{initial columns with (ftable) labels}
	    ok <- nchar(x) > 0 # columns containing non-empty strings
	    ii <- seq_len(m)
	    if(m) ok <- ok[-ii] # drop columns without labels
	    T <- if(!any(ok)) rsep else {
		beg <- m + which(ok) # start positions
		d <- diff(c(beg, nch+1))
		r <- paste0("\\multicolumn{", d, "}{c}{", x[beg], "}") # vector
		paste(paste0(r, collapse=csep), rsep) # attach line ending to the body of the row
	    }
	    if(m) T <- c(paste(x[ii], collapse=csep), T)
	    paste(T, collapse=csep)
	}
	## determine head
	head <- character(2*nrh-1)
	for(i in seq_len(nrh-1)) {
	    head[2*i-1] <- multicolumn(head.raw[i,], m=nrv) # title
	    head[2*i]	<- titlerule  (head.raw[i,], booktabs) # rule
	}
	head[2*nrh-1] <- multicolumn(head.raw[nrh,], m=nrv)
    } else if(is.na(head)) head <- NULL # in case head==NA => no header

    ## put in column separators (=> character vector)
    x4 <- apply(x3, 1, paste, collapse=csep) # put in cell separators
    x5 <- gsub("^\\s+", "", x4) # beautification: trim leading whitespace (for indent)
    x6 <- paste(x5, rsepcol) # put in row separators

    ## return list
    list(body = x6, body.raw = x3, head = head, head.raw = head.raw,
	 align = align, rsepcol = rsepcol)
}
##' { end } tablines

##' @title Wrapper for a floating LaTeX Table
##' @param x a character vector containing the lines of the body of the table
##'	   (everything strictly between '\midrule' and '\bottomrule' [in case of
##'	   booktabs=TRUE]); a table header can be passed via attributes
##' @param align table columns alignment string (e.g., "cccS[table-format=1.2]")
##' @param do.table logical indicating whether \begin{table} ... \end{table}
##'        is used
##' @param placement table placement string
##' @param center logical indicating whether to center the table via \centering
##' @param fontsize possible fontsize adjustment string: "tiny",
##'	   "scriptsize", "footnotesize", "small", "normalsize", "large", "Large",
##'	   "LARGE", "huge", or "Huge"
##' @param booktabs logical indicating whether a booktabs-compliant LaTeX table
##'	   is created (requires \usepackage{booktabs} in the preamble)
##' @param caption table caption or NULL (for no caption)
##' @param label table label or NULL (for no label)
##' @return return value of writeLines()
##' @author Marius Hofert
##' @note - LaTeX requirement: 'tabularx', 'booktabs' (if booktabs=TRUE), and
##'	    siunitx (depends on the table header)
##'	  - For appending to a file, use:
##'	       myfile <- file("foo.txt", "w")
##'	       writeLines(wrapLaTable(...), con=myfile)
##'	       close(myfile)
##'	  - wrapLaTable can't determine 'align' in a sensible way (from strings)
##'	  - both 'xtable' and 'tables' use \hline in the case booktabs=FALSE
##' { wrapLaTable }
wrapLaTable <- function(x, align, do.table = TRUE, placement = "htbp",
			center = TRUE, fontsize = "normalsize", booktabs = TRUE,
			caption = NULL, label = NULL)
{
    stopifnot(is.character(x))
    head <- attr(x, "head") # NA or character
    structure(c(if(do.table)
		paste0("\\begin{table}[", placement, "]"), # float with placement
		paste0(if(center) "  \\centering", # centering and fontsize
		       if(fontsize != "normalsize") paste0("\\", fontsize)),
		paste0("  \\begin{tabular}{", align, "}"), # tabular with alignment
		if(booktabs) "    \\toprule" else "    \\hline\\noalign{\\smallskip}",
		if(is.character(head)) # put in header lines (essentially if head is not NULL)
		paste0("    ",
		       c(head, if(booktabs) "\\midrule"
			 else "\\hline\\noalign{\\smallskip}")),
		## actual entries
		paste0("    ", x),
		##             ==
		## footer
		if(booktabs) "    \\bottomrule" else "    \\hline",
		"  \\end{tabular}",
		if(do.table)
		c(if(!is.null(caption)) paste0("  \\caption{", caption, "}"),
		  if(!is.null(label))   paste0("  \\label{", label, "}"),
		  "\\end{table}")),
	      usepackage = if(booktabs) "booktabs", # <- another atttribute for 'Latex" class
	      class = "Latex")
}
##' { end } wrapLaTable

##' check whether a character string is 'math' (for putting it in $ $)
##' one.is.math = TRUE => single letters count as 'math' (variables)
is.probably.latex.math <- function(ch, one.is.math = FALSE) {
    (if(one.is.math) grepl("^[A-Za-z]$", ch) else FALSE) |
    grepl("\\\\", ch) |
    grepl("\\{.*\\}", ch) |
    grepl("[+^]", ch) ## more "math operators" (which ?)
}

##' do *not* escape the math latex expression, but "mathify"
escapeORmath <- function(x, exprFUN, escapeFUN, one.is.math = TRUE) {
    r <- vapply(x, exprFUN, "")
    isL <- is.probably.latex.math(r, one.is.math=one.is.math)
    r[ isL] <- paste("\\(", r[isL], "\\)") # ~= $ ... $
    r[!isL] <- escapeFUN(r[!isL])
    r
}

##' @title Create an ftable with Expressions Converted and Symbols Escaped
##' @param x ftable in 'character' format
##' @param vList variable specification list
##' @param x.escape logical indicating whether the actual table entries are escaped
##' @param exprFUN function to convert R expressions to LaTeX
##' @param escapeFUN function to escape symbols for LaTeX
##' @return ftable with "expressions", col.vars, and row.vars converted to LaTeX
##'         and possibly escaped.
##' @author Martin Maechler
##' { ftable2latex }
ftable2latex <- function(x, vList = NULL, x.escape,
			 exprFUN = expr2latex, escapeFUN = escapeLatex)
{
    ## checks
    stopifnot(is.function(exprFUN), is.function(escapeFUN))
    cl <- class(x)
    dn <- c(r.v <- attr(x, "row.vars"),
	    c.v <- attr(x, "col.vars"))
    if(is.null(vList)) {
	nvl <- names(vList <- dimnames2varlist(dn))
    } else {
	stopifnot(names(dn) %in% (nvl <- names(vList)))
    }
    vl <- .vl.as.list(vList)
    ## apply escapeORmath() to expressions of column and row variables
    names(c.v) <- lapply(lapply(vl[match(names(c.v), nvl)], `[[`, "expr"),
			 escapeORmath, exprFUN=exprFUN, escapeFUN=escapeFUN)
    names(r.v) <- lapply(lapply(vl[match(names(r.v), nvl)], `[[`, "expr"),
			 escapeORmath, exprFUN=exprFUN, escapeFUN=escapeFUN)
    ## for the entries of 'x' itself, we cannot apply exprFUN(.) everywhere,
    ## only ``where expr''
    exprORchar <- function(u) {
	lang <- vapply(u, is.language, NA) # TRUE if 'name', 'call' or 'expression'
	u[ lang] <- exprFUN	(u[ lang]) # apply (per default) expr2latex()
	u[!lang] <- as.character(u[!lang]) # or format()?
	u
    }
    x <- exprORchar(x) # converts expressions (and only those) to LaTeX
    if(x.escape) x <- escapeFUN(x) # escapes LaTeX expressions
    ## now the transformed row and col names
    attr(x, "row.vars") <- lapply(r.v, escapeFUN)
    attr(x, "col.vars") <- lapply(c.v, escapeFUN)
    class(x) <- cl
    x
}
##' { end } ftable2latex

##' @title Converting an ftable to a LaTeX Table
##' @param object an ftable
##' @param vList variable specification list
##' @param x.escape logical indicating whether the actual table entries are escaped
##'        (FALSE by default to allow colors for entries and other user-crazy stuff)
##' @param exprFUN function to convert R expressions to LaTeX
##' @param escapeFUN function to escape symbols for LaTeX
##' Arguments from tablines():
##' @param align either a character vector...
##'	  - ... of length > 1: e.g., c("c", "c", "c", "S[table-format=1.2]")
##'	  - ... of length 1: e.g., c("*{3}{c} S[table-format=1.2]")
##'	  or NULL in which case a useful default is constructed ('r' for all table
##'       entries, 'l' for all columns of row names)
##' @param booktabs logical indicating whether a booktabs-compliant LaTeX table
##'	   is created (requires \usepackage{booktabs} in the preamble)
##' @param head either
##'	   - a character vector containing the lines of the header
##'	   - NA (= no header)
##'	   - NULL (in which case a default is constructed)
##' @param rsep character string inserted at the end of each row
##' @param sp numeric scaling factor for separating blocks of rows
##' @param rsep.sp numeric vector of length equal to the number of different
##'	   groups or row variables minus 1, giving the spaces (interpreted as pt)
##'	   between the groups
##' @param csep character string for separating different cells in a row
##' Arguments from format.ftable():
##' @param quote see ?format.ftable
##' @param lsep see ?format.ftable
##' Arguments from wrapLaTable():
##' @param do.table logical indicating whether \begin{table} ... \end{table}
##'        is used
##' @param placement table placement string
##' @param center logical indicating whether to center the table via \centering
##' @param fontsize possible fontsize adjustment string: "tiny",
##'	   "scriptsize", "footnotesize", "small", "normalsize", "large", "Large",
##'	   "LARGE", "huge", or "Huge"
##' @param caption table caption or NULL (for no caption)
##' @param label table label or NULL (for no label)
##' @param ... additional arguments passed to fftable()
##' @return return object of wrapLaTable()
##' @author Marius Hofert and Martin Maechler
##' { toLatex.ftable }
toLatex.ftable <- function(object, vList = NULL, x.escape = FALSE,
			   exprFUN = expr2latex, escapeFUN = escapeLatex,
			   align = NULL, booktabs = TRUE, head = NULL,
			   rsep = "\\\\", sp = if(booktabs) 3 else 1.25,
                           rsep.sp = NULL, csep = " & ", quote = FALSE,
                           lsep=" \\textbar\\ ", do.table = TRUE,
                           placement = "htbp", center = TRUE,
                           fontsize = "normalsize", caption = NULL, label = NULL,
			   ...)
{
    ## convert expressions, leave rest:
    ft <- ftable2latex(object, vList, x.escape=x.escape,
		       exprFUN=exprFUN, escapeFUN=escapeFUN)
    ## ftable -> character matrix (formatted ftable) with attributes 'ncv' and 'nrv'
    ft <- fftable(ft, quote=quote, lsep=lsep, ...)
    ## character matrix -> latex {head + body}:
    tlist <- tablines(ft, align=align, booktabs=booktabs,
		      head=head, rsep=rsep, sp=sp, rsep.sp=rsep.sp, csep=csep)
    ## wrap table and return 'Latex' object:
    wrapLaTable(structure(tlist$body, head = tlist$head),
		do.table = do.table, align = tlist$align,
		placement = placement, center = center, booktabs = booktabs,
		fontsize = fontsize, caption = caption, label = label)
}
##' { end } toLatex.ftable

##' arguments similar to toLatex.ftable()
##' @author Martin Maechler
##' { toLatex.varlist }
toLatex.varlist <-
    function(object,
	     col.vars = c("Variable", "expression", "type", "value"),
	     exprFUN = expr2latex, escapeFUN = escapeLatex,
	     align = NULL, booktabs = TRUE, head = NULL,
	     rsep = "\\\\", sp = if(booktabs) 3 else 1.25, rsep.sp = NULL, csep = " & ",
	     do.table = TRUE, placement = "htbp", center = TRUE,
	     fontsize = "normalsize",
	     caption = NULL, label = NULL, ...)
{
    stopifnot((lc <- length(col.vars)) >= 3) # col.vars must at least contain 3 columns (variable, type, value)
    ## fmat := character matrix w/ "nrv" and "ncv" attr()ibutes
    fmat <- cbind(
	Variable = {
	    paste0("\\texttt{", escapeFUN(names(object)), "}")
	},
	expression = {
	    escapeORmath(lapply(object, `[[`, "expr"),
			 exprFUN=exprFUN, escapeFUN=escapeFUN)
	},
	type = {
	    escapeFUN(lapply(object, `[[`, "type"))
	},
	value = {
	    vlis <- mkNms(object) # TODO MM? improve mkNms() allowing drop0.. ?
	    vapply(vlis, function(ch) paste(ch, collapse=", "), "")
	})[, col.vars ]
    fmat <- rbind(col.vars, fmat, deparse.level=0L)
    attr(fmat, "nrv") <- 0
    attr(fmat, "ncv") <- 1
    ## char matrix -> latex {head + body}:
    if(is.null(align)) align <- paste0("l*{",lc-2,"}{c}r")
    tlist <- tablines(fmat, align = align, booktabs = booktabs,
		      head = head, rsep = rsep, sp = sp, rsep.sp = rsep.sp,
                      csep = csep)
    ## wrap table and return "Latex" object:
    wrapLaTable(structure(tlist$body, head = tlist$head),
		do.table = do.table, align = tlist$align,
		placement = placement, center = center, booktabs = booktabs,
		fontsize = fontsize, caption = caption, label = label)
}
##' { end } toLatex.varlist
