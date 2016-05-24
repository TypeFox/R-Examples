####  R --> LaTeX
####  ----------- Exemplifying	 toLatex(), expr2latex() and similar utilities

require("simsalapar")

###---------- Part I --------- Plotmath Expr --> LaTeX -----------------------

## Some "plotmath" formulas we'd want to see:
Lex <- list(quote(1), quote(-1), quote(x), quote(foo.3)
	    , quote(3.1415e100), quote(+3.45), quote(alpha)
	    , quote( N[sim] )
	    , quote( N[sim] ~ O(n) )
	    , quote(x %notin% N)
	    , quote(x %+-% epsilon)
	    , quote(N[s*m^2])
	    , quote( 2^{N[sim] - 3} ~~~ O(n^{n^2}) )
	    )

tab2sym <- function(tab2) {
    stopifnot(isTab(tab2))
    lapply(lapply(tab2, `[[`, 2L), as.symbol)
}

mDeparse <- simsalapar:::mDeparse

showLexpr <- function(lex) {
    data.frame(expr = vapply(lex, mDeparse, ""),
	  length = vapply(lex, length, 1),
	  class = vapply(lex, class, ""),
	  typeof= vapply(lex, typeof, ""),
	  is.atomic= vapply(lex, is.atomic, NA))
}

showLexpr(Lex)

viewPDF <- function(file, pdfviewer = getOption("pdfviewer"))
	system2(pdfviewer, file, wait=FALSE)

##' @title LaTeX "text" -> PDF
##' @param ch character vector of valid LaTeX "text"
##' @param file LaTeX file name to be used
##' @param clean logical indicating if *.log, *.aux, etc are to be removed at the end
##' @param view logical indicating if the PDF viewer (\code{pdfviewer})
##'   should be called to display the resulting pdf file.
##' @param pdfviewer character string specifying a PDF viewer program (to be on the PATH)
##' @return invisibly, a list with components \code{texfile} and \code{pdf}
##'  each with the filename of the corresponding file.
##' @author Martin Maechler
latex2pdf <- function(ch, file = tempfile("Rlatex2pdf", fileext = ".tex"),
                      docclass = c("standalone", "article"),
                      math = FALSE, clean = TRUE,
                      view=TRUE, pdfviewer = getOption("pdfviewer"))
{
    if(!is.character(ch)) ch <- as.character(ch) # {checking "character"}
    stopifnot(require(tools),# -> texi2pdf()
	      length(ch) > 0)
    docclass <- match.arg(docclass)
    usePkgs <- attr(ch, "usepackage")# from 'Latex' object
    is.art <- docclass %in% c("article")# .. more to come ?
    ## Without the '[math]', the cropping is cutting off too much (for some 'expr')
    writeLines(c(paste0("\\documentclass", if(math) "[math]",
			"{", docclass, "}"),
		 if(length(usePkgs)) paste0("\\usepackage{",usePkgs,"}"),
		 "\\begin{document}", ch, "\\end{document}"),
	       con=file)
    message("file: ", file)
    owd <- getwd(); on.exit(setwd(owd))
    setwd(dirname(file))# -> so the log, pdf, ... files are there, too!
    texi2pdf(file, index=FALSE, clean=clean)
    pfile <- sub("\\.tex", ".pdf", file)
    if(!file.exists(pfile))
	warning("something went wrong: PDF (", pfile,") is missing")
    else if(view) viewPDF(pfile, pdfviewer)
    invisible(list(texfile = file, pdf = pfile))
} ## { end } latex2pdf

wrapExpr <- function(ex, LHS="Alpha", RHS="Omega")
    parse(text=paste(LHS, ex, RHS))[[1L]]

##' @title plotmath expression -> LaTeX -> PDF
##' @param ex plotmath "expression"
##' @param make2op logical indicating if wrapExpr() should be applied when isOp()
##' @param file LaTeX file name to be used
##' @param clean logical indicating if *.log, *.aux, etc are to be removed at the end
##' @param view logical indicating if the PDF viewer (\code{pdfviewer})
##'   should be called to display the resulting pdf file.
##' @param pdfviewer character string specifying a PDF viewer program (to be on the PATH)
##' @return invisibly, a list with components \code{texfile} and \code{pdf}
##'  each with the filename of the corresponding file.
##' @author Martin Maechler
pdfExpr <- function(ex, make2op=TRUE, LHS="Alpha", RHS="Omega",
		    file=tempfile("expr2latex", fileext = ".tex"), clean=TRUE,
		    view=TRUE, pdfviewer = getOption("pdfviewer"))
{
    if(!is.language(ex) && length(ex) != 1)
	stop("'ex' must be 'language' or of length 1")
    if(make2op && length(ex) == 1 && isOp(ex))
	ex <- wrapExpr(ex, LHS=LHS, RHS=RHS)
    latex2pdf(c("\\(", expr2latex(ex), "\\)"), math=TRUE,
	      file=file, clean=clean, view=view, pdfviewer=pdfviewer)
}

##' @title plotmath expression -> Simple plot showing it
##' @param ex plotmath "expression"
##' @param x  x- and
##' @param y  y- coordinate of where the graphic should be put.
##' @param add logical specifying if the plot should be added to an existing one.
##' @param make2op logical indicating if binary operators should be
##'   \dQuote{wrapped with} two arguments, \code{LHS} and \code{RHS},
##'   respectively.
##' @return (invisibly, the result of \code{\link{text}})
##' @author Martin Maechler
plExpr <- function(ex, x=.5, y=0, add=FALSE, cex = 4, d.lab= 1.25,
                   make2op=TRUE, LHS="Alpha", RHS="Omega", verbose=FALSE) {
    if(!is.language(ex) && length(ex) != 1)
	stop("'ex' must be 'language' or of length 1")
    ## else
    ## (is.symbol(ex) || is.atomic(ex))))

    Tex <- paste(mDeparse(ex),":")
    if(make2op && length(ex) == 1 && isOp(ex))
	ex <- wrapExpr(ex, LHS=LHS, RHS=RHS)

    ## The plot :  nchar(Tex)  character wide,  (1 + d.lab + cex) character heigh:

    ## Better design (to do in future):
    ## require("grid") and return a GROB --> can have size depend on device!

    if(!add) {
        ## Workaround "bug" in par() functioning:
        ## dev.new(width=2, height=1); op <- par(no.readonly=TRUE) # negative op$pin !
        ## par(mar=rep(0,4))
        o1 <- par(mar = rep(0,4)); op <- par(no.readonly=TRUE); op$mar <- o1$mar
        on.exit(par(op))
        plot.new()
        plot.window(0:1, 0:1, xaxs="i", yaxs="i") ## --> par("usr") =^=  [0,1]^2
    }

    wid <- max(strwidth(ex)*cex, strwidth (Tex))# grid has 'stringWidth(.)' ..
    hei <- (h.T <- strheight(Tex))*(1 + d.lab) + (h.ex <- strheight(ex))*cex

    if(verbose) message(sprintf("(wid = %.7g, hei = %.7g)", wid,hei))
    f.cex <- 1/max(wid,hei)
    opc <- par(cex= f.cex)
    if(add) on.exit(par(opc)) # otherwise, reset much more
    if(verbose) message("par('cex'): ", format(par("cex")))

    d.y <- f.cex * (d.lab*h.T +h.ex*cex)
    text(x, y,	    adj = c(NA,0), labels =  ex, cex=cex, col=2)
    text(x, y+ d.y, adj = c(NA,0), labels = Tex)
}

plExpr(quote(x %notin% N))
plExpr(quote(x %+-% epsilon))

plExpr(quote(N[s*m^2]))

plExpr (quote(x %notin% N))
pdfExpr(quote(x %notin% N))

## Now, try seeing *both*, side by side :

pl.both <- function(ex,
                    pfil = tempfile("expr2L_", fileext = ".pdf"),
                    width = 4, height = width/2,
                    file = tempfile("pl_both2pdf_", fileext = ".tex"),
                    frame2 = TRUE, clean = TRUE,
                    view = TRUE, pdfviewer = getOption("pdfviewer"))
{
    r.pdf <- pdfExpr(ex, view=FALSE)

    pdf(pfil, width=width, height=height)
    plExpr(ex)
    dev.off()

    writeLines(c("\\documentclass{standalone}",
                 "\\usepackage{graphicx}",
                 "\\pagestyle{empty}",
		 "\\begin{document}",

                 paste0("\\includegraphics[width=5in]{",pfil,"}"),

                 if(frame2) "\\framebox{",
                 paste0("\\includegraphics[width=5in]{",r.pdf$pdf,"}"),
                 if(frame2) "}",

                 "\\end{document}"),
	       con= file)
    message("file: ", file)
    owd <- getwd(); on.exit(setwd(owd))
    setwd(dirname(file))# -> so the log, pdf, ... files are there, too!
    texi2pdf(file, index=FALSE, clean=clean)
    pfile <- sub("\\.tex", ".pdf", file)
    if(!file.exists(pfile))
	warning("something went wrong: PDF (", pfile,") is missing")
    else if(view) viewPDF(pfile, pdfviewer)
    invisible(list(texfile = file, pdf = pfile))
}

pl.both(quote(x %notin% N))
pl.both(quote(N[s*m^2]))


###----------- Part II --------- VarList -> Ftable --> LaTeX Tables --------------

### a)   ftable --> LaTeX : -------------------------------------

ft1 <- ftable(Titanic, col.vars = 1:4)
ft2 <- ftable(Titanic, row.vars = 1)
(ft3 <- ftable(Titanic, row.vars = 1:2))
(ft4 <- ftable(Titanic, row.vars = 1:3))
ft5 <- ftable(Titanic, row.vars = 1:4)

## doc.class "standalone"  does *not* work with begin{table} .. end{table}
latex2pdf(toLatex(ft4), docclass="article")

## However it *does* very nicely (!) tabular directly (and with*out* booktabs):
latex2pdf(toLatex(ft4, do.table=FALSE, booktabs=FALSE))# "standalone"
## quite nice ... but "ugly":
## "Age -- Survived" header, then Child / Adult is below "Survived" instead of "Age"

## and 'standalone' now even works with  'booktabs':
latex2pdf(toLatex(ft4, do.table=FALSE))


latex2pdf(toLatex(ft3, do.table=FALSE, method="compact"))
## -->  MM:"really ugly -- format.ftable() is much nicer"
latex2pdf(toLatex(ft3, do.table=FALSE, method="col.compact"))
latex2pdf(toLatex(ft3, do.table=FALSE, method="row.compact"))# MH: "ugly"
latex2pdf(toLatex(ft3, do.table=FALSE, method="non.compact"))


## What tablines() returns
tablines(fftable(ft2))

## LaTeX (booktabs/non-booktabs) versions
for(ft in list(ft1,ft2,ft3,ft4,ft5))
    for(bookt in c(FALSE, TRUE))
        for(meth in c("non.compact", "row.compact", "col.compact", "compact")) {
            lt <- toLatex(ft, booktabs=bookt, method=meth, do.table=FALSE)
            latex2pdf(lt, file = paste0("R_ft", nrow(ft),"-",if(bookt)"T","-",
                          abbreviate(meth),".tex"))
        }


### b)   varlist --> LaTeX : -------------------------------------

vList <- varlist(
    n.sim = list(value = 1000, expr = quote(N[sim])),
    n     = list(type="grid", value = c(20, 100, 500)), # sample sizes
    meth  = list(type="grid", expr = quote(italic(method)),
                 value = c("classical", "robust")),
    alpha = list(value = 0.95))


latex2pdf(toLatex(vList, do.table=FALSE))

arrTitan <- as.table(Titanic)
vlTit <- dimnames2varlist(dimnames(arrTitan))

latex2pdf(toLatex(vlTit, do.table=FALSE))
