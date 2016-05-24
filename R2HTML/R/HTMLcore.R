#     $Id: HTMLcore.R 47 2008-05-23 17:29:31Z mentus $
#     R2HTML - Library of exportation to HTML for R
#     Copyright (C) 2002-2004 - Eric Lecoutre

#     R2HTML Package

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#----------------------------------------------------------------------------------------------------#
#
#     Contact:
#
#     Eric Lecoutre
#     <lecoutre@stat.ucl.ac.be>
#
#     Institut de statistique
#     Voie du Roman Pays, 20
#     1348 Louvain-la-Neuve
#     BELGIQUE
#
#----------------------------------------------------------------------------------------------------#
".HTMLEnv" <- new.env(parent=emptyenv())

#----------------------------------------------------------------------------------------------------#
"HTMLSetFile" <- function(file) {
    assign(".HTML.file", file, .HTMLEnv)
    file
}

#----------------------------------------------------------------------------------------------------#
"HTMLGetFile" <- function() {
    if(exists(".HTML.file", .HTMLEnv))
         get(".HTML.file", .HTMLEnv)
    else
         stop("not default HTML output file defined; please call HTMLSetFile() to set it")
}

#----------------------------------------------------------------------------------------------------#

"HTML"<- function(x,...) {
	UseMethod("HTML")
  	}

#----------------------------------------------------------------------------------------------------#

"HTML.default"<-
function(x, file=HTMLGetFile(), append=TRUE,...)
{
	HTML(paste(capture.output(x),collapse="\n<br>\n"),file=file,append=append,...)
	invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.atomic"<- function(x, file=HTMLGetFile(), append=TRUE, ...){
	cat(paste("\n<p class='atomic'>",paste(x,collapse="&nbsp; "),"</p>\n",sep="",collapse=""), file= file, append = append, sep = " ")
}

#----------------------------------------------------------------------------------------------------#

"HTML.complex"<- function(x, file=HTMLGetFile(), append=TRUE,...){
	cat(paste("\n<p><font class='complexRe'>",Re(x),"</font>",ifelse(sign(Im(x))<0,"-","+"),"<font class='complexIm'>",Im(x),"</font><font class='complexI'>i</font>","</p>\n",sep="",collapse=""), file= file, append = append, sep = " ")
	}

#----------------------------------------------------------------------------------------------------#

"HTML.numeric"<- function(x, file=HTMLGetFile(), append=TRUE, ...){
	if(!is.null(names(x))) {
		HTML(as.table(x),file=file,append=append,...)
		}
	else {
		cat(paste("\n<p class='numeric'>",paste(x,collapse="&nbsp; "),"</p>\n",sep="",collapse=""), file= file, append = append, sep = " ")
		}
	}
#----------------------------------------------------------------------------------------------------#

"HTML.integer"<- function(x, file=HTMLGetFile(), append=TRUE, ...){
	cat(paste("\n<p class='integer'>",paste(x,collapse="&nbsp; "),"</p>\n",sep="",collapse=""), file= file, append = append, sep = " ")
	}

#----------------------------------------------------------------------------------------------------#

"HTML.logical"<- function(x, file=HTMLGetFile(), append=TRUE,...){
	cat(paste("\n<p class='logical'>",paste(x,collapse="&nbsp; "),"</p>\n",sep="",collapse=""), file= file, append = append, sep = " ")
	}

#----------------------------------------------------------------------------------------------------#

"HTML.character"<- function(x, file=HTMLGetFile(), append=TRUE, ...){
	cat(paste("\n<p class='character'>",paste(x,collapse="&nbsp; "),"</p>\n",sep="",collapse=""), file= file, append = append, sep = " ")
	}

#----------------------------------------------------------------------------------------------------#

"HTML.call"<- function(x, file=HTMLGetFile(), append=TRUE, ...){
	cat(paste("<font class='call'>",deparse(x),"</font>",sep="",collapse=""), file= file, append = append, sep = " ")
	}

#----------------------------------------------------------------------------------------------------#

"HTML.function"<-function(x,file=HTMLGetFile(), append=TRUE,...){
	 cat(paste("\n<br>\n<xmp class=function>",
	 paste(capture.output(x),collapse="\n"),"\n</xmp><br>\n",sep=""),
	file=file,append=append,sep="\n<br>\n")
	invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.environment"<-function(x,file=HTMLGetFile(), append=TRUE,...){
	cat(paste("\n<br>environment: <font class='environment'>",attributes(x)$name,"</font><br>\n",sep=""),
	file=file,append=append)
	invisible(x)
}
#----------------------------------------------------------------------------------------------------#

"HTML.formula"<-function(x,file=HTMLGetFile(), append=TRUE,...) {
	HTML(paste("<font class='formula'>",deparse(unclass(x)),"</font>",collapse=""),file=file,append=append,...)
	}

#----------------------------------------------------------------------------------------------------#

"HTML.array"<- function(x, file=HTMLGetFile(), append=TRUE, ...)
{
	odometer <- function(current, radix)
	{
		if(any(c(current, radix) < 0))
			stop("arguments must be non-negative")
		lc <- length(current)
		if(length(radix) != lc)
			radix <- rep(radix, length = lc)
		radix <- radix - 1
		for(i in 1:lc) {
			if((ii <- current[i]) < radix[i]) {
				current[i] <- ii + 1
				return(current)
			}
			else current[i] <- 0
		}
		current
	}


	d <- dim(x)
	ndim <- length(d)
	dn <- dimnames(x)
	if(ndim == 1)
		HTML.matrix(matrix(x, 1, dimnames = list("", if(is.null(
			dn)) paste("[", 1:d[1], "]", sep = "") else dn[[1]])),
			file = file, append=append,...)
	else if(ndim == 2)
		HTML.matrix(x, Border = 0, file = file, append=append,...)
	else {
		if(length(dn) < ndim)
			dn <- vector("list", ndim)
		for(i in 3:ndim)
			if(length(dn[[i]]) < d[i]) dn[[i]] <- paste(1:d[i])
		xm <- array(x[1], d[1:2])
		dimnames(xm) <- dn[1:2]
		d <- d[ - (1:2)]
		nm <- length(xm)
		which <- 1:nm
		dn <- dn[ - (1:2)]
		ndim <- ndim - 2
		counter <- rep(0, length(d))
		for(i in 1:(length(x)/nm)) {
			cat("<br>, , ", file = file, append = TRUE)
			for(j in 1:ndim)
				cat(dn[[j]][counter[j] + 1], if(j < ndim) ", "
				   else "<br>", sep = "", file = file, append
				   = TRUE)
			xm[1:nm] <- x[which]
			HTML.matrix(xm, Border = 0, file = file, append=TRUE,...)
			counter <- odometer(counter, d)
			which <- which + nm
		}
	}
	invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.by"<- function (x, file=HTMLGetFile(),vsep="\n<hr size=1 width=100%>\n",append=TRUE,...)
{

    HTML("\n",file=file,append=append,...)
    d <- dim(x)
    dn <- dimnames(x)
    dnn <- names(dn)
    if (missing(vsep))
        vsep <- "\n<hr size=1 width=100%>\n"
    lapply(seq(along = x), function(i, x, vsep, ...) {
        if (i != 1 && !is.null(vsep))
            HTML(vsep, file=file,append=TRUE)
        ii <- i - 1
        for (j in seq(along = dn)) {
            iii <- ii%%d[j] + 1
            ii <- ii%/%d[j]
            HTML(paste(dnn[j], ": ", dn[[j]][iii], "\n<br>", sep = ""),file=file,append=TRUE,...)
        }
        HTML(x[[i]], file=file,append=TRUE)
    }, x, vsep, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.family" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    HTML(paste("\n<br><b>Family</b>:<font class='family'>", x$family, "\n</font><br>",sep=""),file=get(".HTML.file",pos=1),append=append,...)
    HTML(paste("\n<b>Link function</b>:<font class='link'>", x$link, "\n</font><br>\n<br>",sep=""),file=get(".HTML.file",pos=1),append=TRUE,...)
}

#----------------------------------------------------------------------------------------------------#

"HTML.terms" <- function (x, file=HTMLGetFile(), append=TRUE,...)	HTML.default(paste("<font class='terms'>",unclass(x),"</font>",sep=""),file=file,append=append,...)

#----------------------------------------------------------------------------------------------------#

"HTML.factor" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    HTML("\n\n<font class='factor'>",file=file,append=append,...)
    if (length(x) <= 0)
        HTML("factor(0)\n<br>\n",file=file,append=TRUE,...)
    else HTML(as.character(x), file=file,append=TRUE, ...)
    HTML("</font>\n",file=file,append=TRUE,...)
    HTMLbr(file=file,append=TRUE,...)
    HTML(paste("Levels:<font class='factorlevels'> ", paste(levels(x), collapse = " "), "</font>\n<br>",sep=""),file=file,append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
"HTML.density" <- function (x,file=HTMLGetFile(),  digits=4,append=TRUE,...)
{

    HTML(paste("\n<br><b>Call</b>:<font class='call'>\n      ", deparse(x$call), "</font><br><br>\n\n<b>Data</b><font class='dataname'>: ", x$data.name,
        "</font> (", x$n, " obs.);", " <b>Bandwidth</b> 'bw' = ", round(x$bw, digits), "\n<br>\n<br>", sep = ""),append=append,file=file)
    HTML(summary(as.data.frame(x[c("x", "y")])),append=TRUE, ...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#
"HTML.infl" <- function (x,  file=HTMLGetFile(),digits = max(3, getOption("digits") - 4),append=TRUE,...)
{
    HTML(paste("\n<br>Influence measures of\n<br>      <font class='call'>  ", deparse(x$call), ":</font>\n<br>\n<br>",sep=""),file=file,append=append,...)
    is.star <- apply(x$is.inf, 1, any, na.rm = TRUE)
    HTML(data.frame(round(x$infmat,digits), inf = ifelse(is.star, "*", " ")),file=file, append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.lm"<-function(x,file=HTMLGetFile(),digits= max(3, getOption("digits") - 3),append=TRUE,...)
{
	HTMLli(paste("Call: <font class='call'>",deparse(x$call),"</font>",sep=""),file=file,append=append,...)
	HTMLli("Coefficients<br>",file=file,append=TRUE,...)
	HTML(round(x$coeff,3),file=file,append=TRUE,...)

}

#----------------------------------------------------------------------------------------------------#
"HTML.lm.null" <- function (x, file=HTMLGetFile(),digits = max(3, getOption("digits") - 3),append=TRUE,...)
{
    HTMLli(paste("Call: <font class='call'>", deparse(x$call),"</font>", "\n<br>", sep = ""),file=file,append=append,...)
    HTMLli("No coefficients<br>\n",append=TRUE,...)
    invisible(x)
}
#----------------------------------------------------------------------------------------------------#


"HTML.ftable" <- function (x,  file=HTMLGetFile(),digits = getOption("digits"),append=TRUE,...)
{
 if (!inherits(x, "ftable"))
        stop("x must be an `ftable'")
    ox <- x
    makeLabels <- function(lst) {
        lens <- sapply(lst, length)
        cplensU <- c(1, cumprod(lens))
        cplensD <- rev(c(1, cumprod(rev(lens))))
        y <- NULL
        for (i in rev(seq(along = lst))) {
            ind <- 1 + seq(from = 0, to = lens[i] - 1) * cplensD[i +
                1]
            tmp <- character(length = cplensD[i])
            tmp[ind] <- lst[[i]]
            y <- cbind(rep(tmp, times = cplensU[i]), y)
        }
        y
    }
    makeNames <- function(x) {
        nmx <- names(x)
        if (is.null(nmx)) nmx <- rep("", length = length(x))
        nmx
    }
    xrv <- attr(x, "row.vars")
    xcv <- attr(x, "col.vars")
    LABS <- cbind(rbind(matrix("", nrow = length(xcv), ncol = length(xrv)), makeNames(xrv), makeLabels(xrv)), c(makeNames(xcv),rep("", times = nrow(x) + 1)))
    DATA <- rbind(t(makeLabels(xcv)), rep("", times = ncol(x)), format(unclass(x), digits = digits))
    x <- cbind(apply(LABS, 2, format, justify = "left"), apply(DATA, 2, format, justify = "right"))
    HTML(x,file=file,append=append,...)
    invisible(ox)
}

#----------------------------------------------------------------------------------------------------#

"HTML.POSIXlt" <- function (x, file=HTMLGetFile(), append=TRUE,...) HTML(paste("<P class='POSIXlt'>",format(x, usetz = TRUE),"</p>",sep=""), file=file,append=append,...)

#----------------------------------------------------------------------------------------------------#

"HTML.POSIXct" <- function (x, file=HTMLGetFile(), append=TRUE,...) HTML(paste("<P class='POSIXct'>",format(x, usetz = TRUE),"</p>",sep=""), file=file,append=append,...)


#----------------------------------------------------------------------------------------------------#

"HTML.octmode" <- function (x, file=HTMLGetFile(), append=TRUE,...)  HTML(paste("<P class='octmode'>",format(x),"</p>",sep=""), file=file,append=append,...)

#----------------------------------------------------------------------------------------------------#

"HTML.rle" <- function (x, digits = getOption("digits"), file=HTMLGetFile(), append=TRUE,...)
{
    HTML("<b><center>Run Length Encoding</center></b>\n<br>\n",file=file,append=append,...)
	tab<-rbind(x$length,x$values)
	tab<-cbind(c("Length","Values"),tab)
    HTML(tab,file=file,append=TRUE,...)
}

#----------------------------------------------------------------------------------------------------#

"HTML.logLik" <- function (x, file=HTMLGetFile(),digits = getOption("digits"),append=TRUE,...)    HTML(paste("<p>`log Lik.' ", format(c(x), digits = digits), " (df=",  format(attr(x, "df")), ")\n</p>", sep = ""),file=file,append=append,...)

#----------------------------------------------------------------------------------------------------#

 "HTML.xtabs" <- function (x,file=HTMLGetFile(), append=TRUE,...)
{
    ox <- x
    attr(x, "call") <- NULL
    HTML.table(x,file=file, append=append,...)
    invisible(ox)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.lm"<-function (x, file=HTMLGetFile(),digits = max(3, getOption("digits") - 3), symbolic.cor = p >   4, signif.stars = getOption("show.signif.stars"),append=TRUE,...)
{

	HTML("\n",file=file,append=append)
	HTMLli(paste("Call:<font class='call'> ",deparse(x$call),"</font>","\n", sep = "", collapse = ""),file=file,append=TRUE)

	resid <- x$residuals
	df <- x$df
	rdf <- df[2]

	HTMLli(paste(if (!is.null(x$w) && diff(range(x$w))) "Weighted "," Residuals<br>\n"),file=file,append=TRUE)
	if (rdf > 5) {
	    nam <- c("Min", "1Q", "Median", "3Q", "Max")
	    rq <- if (length(dim(resid)) == 2)
		structure(apply(t(resid), 1, quantile), dimnames = list(nam,   dimnames(resid)[[2]]))
	    else structure(quantile(resid), names = nam)
	    HTML(rq,  file=file,append=TRUE,...)
	}
	else if (rdf > 0) {
	    HTML(resid,file=file,append=TRUE,...)
	}
	else {
	    HTML(paste("ALL", df[1], "residuals are 0: no residual degrees of freedom!<br>\n",sep=""),file=file,append=TRUE,...)
	}
	if (nsingular <- df[3] - df[1])

		HTMLli(paste("Coefficients (",nsingular, "not defined because of singularities)<br>\n",sep=""),file=file,append=TRUE)
	else HTMLli("Coefficients\n",file=file,append=TRUE)


	HTML.coefmat(x$coef, digits = digits, signif.stars = signif.stars, file=file,append=TRUE,...)

	HTMLli(paste("Residuals standard error: ",round(x$sigma,digits)," on ",rdf," degrees of freedom\n",sep=""),file=file,append=TRUE)



	if (!is.null(x$fstatistic)) {
		HTMLli(paste("Multiple R-Squared:<b>",round(x$r.squared,digits),"</b>",sep=""),file=file,append=TRUE)
		HTMLli(paste("Adjusted R-Squared:<b>",round(x$adj.r.squared,digits),"</b>",sep=""),file=file,append=TRUE)
	    	HTMLli(paste("F-statistics: <b>", round(x$fstatistic[1],digits), "</b> on ",x$fstatistic[2], " and ", x$fstatistic[3], " DF. P-value:<b>",round(1-pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3]),digits),"</b>." ,sep=""),file=file,append=TRUE)
	 	}
	correl <- x$correlation
	if (!is.null(correl)) {
	    p <- NCOL(correl)
	    if (p > 1) {
		HTMLli("Correlation of Coefficients:\n",file=file,append=TRUE,...)
		if (symbolic.cor)
		    HTML(symnum(correl)[-1, -p],file=file,append=TRUE,...)
		else {
		    correl[!lower.tri(correl)] <- NA
		    HTML(correl[-1, -p, drop = FALSE],file=file,append=TRUE,...)
		}
	    }
	}
	invisible(x)
}


#----------------------------------------------------------------------------------------------------#
"HTML.coefmat"<- function (x, digits = max(3, getOption("digits") - 2), signif.stars = getOption("show.signif.stars"),
    dig.tst = max(1, min(5, digits - 1)), cs.ind = 1:k, tst.ind = k +
        1, zap.ind = integer(0), P.values = NULL, has.Pvalue = nc >=
        4 && substr(colnames(x)[nc], 1, 3) == "Pr(", na.print = "",file=HTMLGetFile(), append=TRUE,...)
{
   cat("\n",file=file,append=append,...)
    if (is.null(d <- dim(x)) || length(d) != 2)
        stop("1st arg. 'x' must be coefficient matrix/d.f./...")
    nc <- d[2]
    if (is.null(P.values)) {
        scp <- getOption("show.coef.Pvalues")
        if (!is.logical(scp) || is.na(scp)) {
            warning("option `show.coef.Pvalues' is invalid: assuming TRUE")
            scp <- TRUE
        }
        P.values <- has.Pvalue && scp
    }
    else if (P.values && !has.Pvalue)
        stop("'P.values is TRUE, but has.Pvalue not!")
    if (has.Pvalue && !P.values) {
        d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
        nc <- nc - 1
        has.Pvalue <- FALSE
    }
    else xm <- data.matrix(x)
    k <- nc - has.Pvalue - (if (missing(tst.ind))
        1
    else length(tst.ind))
    if (!missing(cs.ind) && length(cs.ind) > k)
        stop("wrong k / cs.ind")
    Cf <- array("", dim = d, dimnames = dimnames(xm))
    ok <- !(ina <- is.na(xm))
    if (length(cs.ind) > 0) {
        acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
        digmin <- 1 + floor(log10(range(acs[acs != 0], na.rm = TRUE)))
        Cf[, cs.ind] <- format(round(coef.se, max(1, digits -
            digmin)), digits = digits)
    }
    if (length(tst.ind) > 0)
        Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
            digits = digits)
    if (length(zap.ind) > 0)
        Cf[, zap.ind] <- format(zapsmall(xm[, zap.ind], digits = digits),
            digits = digits)
    if (any(r.ind <- !((1:nc) %in% c(cs.ind, tst.ind, zap.ind,
        if (has.Pvalue) nc))))
        Cf[, r.ind] <- format(xm[, r.ind], digits = digits)
    okP <- if (has.Pvalue)
        ok[, -nc]
    else ok
    x0 <- (xm[okP] == 0) != (as.numeric(Cf[okP]) == 0)
    if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
        Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1,
            digits - 1))
    }
    if (any(ina))
        Cf[ina] <- na.print
    if (P.values) {
        if (!is.logical(signif.stars) || is.na(signif.stars)) {
            warning("option `show.signif.stars' is invalid: assuming TRUE")
            signif.stars <- TRUE
        }
        pv <- xm[, nc]
        if (any(okP <- ok[, nc])) {
            Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst)
            signif.stars <- signif.stars && any(pv[okP] < 0.1)
            if (signif.stars) {
                Signif <- symnum(pv, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
                Cf <- cbind(Cf, formatC(unclass(Signif)))
            }
        }
        else signif.stars <- FALSE
    }
    else signif.stars <- FALSE

    HTML.matrix(Cf, file=file,  ...)
    if (signif.stars)     HTML(paste("\n<p>--- Signif. codes: ", attr(Signif, "legend"), "</p>\n",sep=""),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.table"<- function(x, file=HTMLGetFile(), append=TRUE,digits=4,...)
{
	cat("\n",file=file,append=append)
	if (!is.null(digits) && is.numeric(x)) x <- round(x,digits) # PhG, because summary(iris) returns a table, but it is not numeric!
	if (is.null(dim(x))) HTML(t(as.matrix(x)),file=file,append=TRUE,digits=NULL,...)
	else HTML(unclass(x),file=file,append=TRUE,...)
}


#----------------------------------------------------------------------------------------------------#

"HTML.listof" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
   cat("\n",file=file,append=append,...)
    nn <- names(x)
    ll <- length(x)
    if (length(nn) != ll)
        nn <- paste("Component ", seq(ll))
    for (i in seq(length = ll)) {
        HTMLli(paste(nn[i],":\n<br>",sep=""),file=file)
        HTML(x[[i]], file=file)
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.ts" <- function (x, calendar=NULL, file=HTMLGetFile(), append=TRUE,...)
{
   cat("\n", file=file,append=append,...)
    x.orig <- x
    x <- as.ts(x)
    fr.x <- frequency(x)
    if (missing(calendar))
        calendar <- any(fr.x == c(4, 12))
    if (!calendar)
        header <- function(x) {
            if ((fr.x <- frequency(x)) != 1)
		HTML(paste("\n<br><b>Time series</b>:\n<br><li>Start=",deparse(start(x)),"\n<br><li>End=",deparse(end(x)),"\n<br><li>Frequency=",deparse(fr.x),"\n<br>",sep=""),file=file)
            else
            HTML(paste("\n<br><b>Time series</b>:\n<br><li>Start=",format(tsp(x)[1]),"\n<br><li>End=",format(tsp(x)[2]),"\n<br><li>Frequency=",deparse(fr.x),"\n<br>",sep=""),file=file)
	        }
    if (NCOL(x) == 1) {
        if (calendar) {
            if (fr.x > 1) {
                dn2 <- if (fr.x == 12)
                  month.abb
                else if (fr.x == 4) {
                  c("Qtr1", "Qtr2", "Qtr3", "Qtr4")
                }
                else paste("p", 1:fr.x, sep = "")
                if (NROW(x) <= fr.x && start(x)[1] == end(x)[1]) {
                  dn1 <- start(x)[1]
                  dn2 <- dn2[1 + (start(x)[2] - 2 + seq(along = x))%%fr.x]
                  x <- matrix(format(x, ...), nrow = 1, byrow = TRUE,
                    dimnames = list(dn1, dn2))
                }
                else {
                  start.pad <- start(x)[2] - 1
                  end.pad <- fr.x - end(x)[2]
                  dn1 <- start(x)[1]:end(x)[1]
                  x <- matrix(c(rep("", start.pad), format(x,
                    ...), rep("", end.pad)), ncol = fr.x, byrow = TRUE,
                    dimnames = list(dn1, dn2))
                }
            }
            else {
                tx <- time(x)
                attributes(x) <- NULL
                names(x) <- tx
            }
        }
        else {
            header(x)
            attr(x, "class") <- attr(x, "tsp") <- attr(x, "na.action") <- NULL
        }
    }
    else {
        if (calendar && fr.x > 1) {
            tm <- time(x)
            t2 <- 1 + round(fr.x * ((tm + 0.001)%%1))
            p1 <- format(floor(tm))
            rownames(x) <- if (fr.x == 12)
                paste(month.abb[t2], p1, sep = " ")
            else paste(p1, if (fr.x == 4)
                c("Q1", "Q2", "Q3", "Q4")[t2]
            else format(t2), sep = " ")
        }
        else {
            if (!calendar)
                header(x)
            rownames(x) <- format(time(x))
        }
        attr(x, "class") <- attr(x, "tsp") <- attr(x, "na.action") <- NULL
    }
    NextMethod("HTML", x, file=file, ...)
    invisible(x.orig)
}

#----------------------------------------------------------------------------------------------------#


"HTML.list" <- function(x,file=HTMLGetFile(),first=TRUE,append=TRUE,...)
{
	cat("\n", file=file,append=append,...)
	if (first) {HTML("<hr class='hr'>",file=file,append=TRUE,sep="\n")}
	for (i in 1:length(x))  {
		cat("<ul>",file=file,append=TRUE,sep="\n")
		cat("</center><li>",file=file,append=TRUE,sep="\n")
		HTML(x[[i]],file=file,first=FALSE,...)
		cat("</ul>",file=file,append=TRUE,sep="\n")

	}
	cat("\n<br><hr class='hr'>",file=file,append=TRUE,sep="\n")
}
#----------------------------------------------------------------------------------------------------#

"HTML.pairlist" <- function(x,file=HTMLGetFile(),first=TRUE,append=TRUE,...)
{
	cat("\n", file=file,append=append,...)
	if (first) {HTML("<hr class='hr'>",file=file,append=TRUE,sep="\n")}
	for (i in 1:length(x))  {
		cat("<ul>",file=file,append=TRUE,sep="\n")
		cat("</center><li>",file=file,append=TRUE,sep="\n")
		HTML(x[[i]],file=file,first=FALSE,...)
		cat("</ul>",file=file,append=TRUE,sep="\n")

	}
	cat("\n<br><hr class='hr'>",file=file,append=TRUE,sep="\n")
}



#----------------------------------------------------------------------------------------------------#

# row.names option contributed by
# Tobias Verbeke on 2006-05-27
#
# Fixed bug of invalid HTML output when using
# row.names = FALSE, as patch contributed
# by Michael Irskens on 2006-11-04
#

"HTML.data.frame" <- function(
            x, file=HTMLGetFile(),
            Border = 1, innerBorder = 0,
            classfirstline = "firstline",
            classfirstcolumn = "firstcolumn",
            classcellinside = "cellinside",
            append = TRUE,
            align = "center",
            caption = "",
            captionalign = "bottom",
            classcaption = "captiondataframe",
            classtable = "dataframe",
            digits = getOption("R2HTML.format.digits"),
            nsmall = getOption("R2HTML.format.nsmall"),
            big.mark = getOption("R2HTML.format.big.mark"),
            big.interval = getOption("R2HTML.format.big.interval"),
            decimal.mark = getOption("R2HTML.format.decimal.mark"),
            sortableDF = getOption("R2HTML.sortableDF"),
            row.names = TRUE,
            ...)
{
   cat("\n", file = file, append = append)

    # Handle sortableDF argument
    if (is.null(sortableDF)) sortableDF = FALSE
    if (sortableDF)
      cat(paste(c("\n<style>", ".tablesort  {",
                  "cursor: pointer ;",
                  " behavior:url(tablesort.htc);",
                  " -moz-binding: url(moz-behaviors.xml#tablesort.htc);",
                  "}",
                  "</style>\n"),
                  collapse="\n"),
          file = file, append = TRUE)


   # if (!is.null(digits)) x[] = lapply(x, FUN = function(vec) if (is.numeric(vec)) round(vec, digits) else vec)

   txt <- paste("\n<p align=",align,">")
   txtcaption <- ifelse(is.null(caption),
                        "",
                        paste("\n<caption align=", captionalign,
                              " class=", classcaption, ">",
                              caption,
                              "</caption>\n", sep=""))

   if (!is.null(Border))
     txt <- paste(txt, "\n<table cellspacing=0 border=", Border, ">",
                  txtcaption,"<tr><td>",
                  "\n\t<table border=", innerBorder, " class=",classtable,">",
                  sep = "")
   else txt <- paste(txt, "\n<table border=", innerBorder,
                     " class=",classtable," cellspacing=0>",
                     txtcaption, sep = "")
   txt <- paste(txt,"\t<tbody>",sep="\n")

   VecDebut <- c(
        if(row.names)
          paste("\n\t\t<th>",
                if(sortableDF) '<b class="tablesort">',
                sep = "", collapse = ""),
        rep(paste("\n\t\t<th>",
                  if(sortableDF) '<b class="tablesort">',
                  sep = "", collapse = ""), ncol(x) - 1)
                )
   VecMilieu <- c(
                 if(row.names) "&nbsp;",
                 as.character(dimnames(x)[[2]])
                 )
   VecFin <- c(
              if(row.names)
                paste(if(sortableDF) '</b>', "", "</th>", collapse = ""),
              rep(
                  paste(if(sortableDF) '</b>',"", "</th>", collapse = ""), ncol(x) - 1
                 ),
              "</th>"
              )
   txt <- paste(txt, "\n\t<tr class=", classfirstline, ">",
                paste(VecDebut, VecMilieu, VecFin, sep = "", collapse = ""),
                "\n\t</tr>"
                )

   x.formatted <- format(x, digits = digits, nsmall = nsmall,
                         big.mark = big.mark, big.interval = big.interval,
                         decimal.mark = decimal.mark)
   x.formatted <- as.matrix(x.formatted)
   x.formatted[is.na(x.formatted)] <- " "
   x.formatted[is.nan(x.formatted)] <- " "

   for(i in 1:dim(x)[1]) {
      if(i == 1) {
         VecDebut <- c(if(row.names)
                         paste("\n<td class=", classfirstcolumn, ">",
                               sep = ""),
                       paste("\n<td class=", classcellinside, ">", sep = ""),
                       rep(paste("\n<td class=", classcellinside, ">",
                                 sep = ""),
                           dim(x)[2] - 1)
                      )
         VecMilieu <- c(if(row.names)
                          dimnames(x)[[1]][i],
                        HTMLReplaceNA(x.formatted[i,])
                       )
         VecFin <- c(if(row.names) "\n</td>",
                     rep("\n</td>", dim(x)[2] - 1),
                     "\n</td></tr>\n"
                    )
      }
      else {
         VecDebut <- c(if(row.names)
                         paste("\n<td class=", classfirstcolumn, ">",
                               sep = ""),
                       paste(rep(paste("\n<td class=", classcellinside, ">",
                                       sep = ""),
                                 dim(x)[2])
                            )
                      )
         VecMilieu <- c(if(row.names)
                          dimnames(x)[[1]][i],
                        HTMLReplaceNA(x.formatted[i,]))
         VecFin <- c(if(row.names) "\n</td>",
                     rep("\n</td>", dim(x)[2] - 1),
                     "\n</td></tr>\n")
      }
      txt <- paste(txt,  "\n<tr>",
                   paste(VecDebut, VecMilieu, VecFin, sep = "", collapse = ""))
   }
   txt <- paste(txt, "\n\t</tbody>\n</table>\n",
                if (!is.null(Border)) "</td></table>\n","<br>")
   cat(txt, "\n", file = file, sep = "", append = TRUE)

}

#----------------------------------------------------------------------------------------------------#

"HTML.matrix" <- function(x, file=HTMLGetFile(), Border = 1, innerBorder = 0, classfirstline = "firstline", classfirstcolumn = "firstcolumn", classcellinside = "cellinside",  append=TRUE,align="center",caption="",captionalign="bottom",classcaption="captiondataframe",classtable="dataframe",digits=getOption("R2HTML.format.digits"),nsmall = getOption("R2HTML.format.nsmall"), big.mark = getOption("R2HTML.format.big.mark"), big.interval = getOption("R2HTML.format.big.interval"), decimal.mark = getOption("R2HTML.format.decimal.mark"),...)
{
   cat("\n", file=file,append=append)

   # if (is.numeric(x) & !is.null(digits)) x<-round(x,digits=digits)

   txt <- paste("\n<p align=",align,">")
   txtcaption <- ifelse(is.null(caption),"",paste("<caption align=",captionalign," class=",classcaption,">",caption,"</caption>\n",sep=""))

   if (!is.null(Border)) txt <- paste(txt, "\n<table cellspacing=0 border=",Border,">",txtcaption,"<tr><td>","\n\t<table border=", innerBorder,  " class=",classtable,">", sep = "")
   else txt <- paste(txt, "\n\t<table border=", innerBorder, " class=", classtable," cellspacing=0>", txtcaption, sep = "")


   txt <- paste(txt,"\t<tbody>",sep="\n")


   if(is.null(dimnames(x)[[2]]) == FALSE) {
      VecDebut <- c(if(is.null(dimnames(x)[[1]]) == FALSE) paste(
            "<th>", sep = ""),
         rep(paste("<th>", sep = ""), dim(
         x)[2] - 1))
      VecMilieu <- c(if(is.null(dimnames(x)[[1]]) == FALSE) "",
         as.character(dimnames(x)[[2]]))
      VecFin <- c(if(is.null(dimnames(x)[[1]]) == FALSE) "</th>", rep(
         "</th>", dim(x)[2] - 1), "</th>")
      txt <- paste(txt,"<tr class=",classfirstline,">", paste(VecDebut, VecMilieu, VecFin, sep = "",collapse = ""),"</tr>\n")
   }

     x.formatted <- format(x, digits=digits, nsmall=nsmall, big.mark=big.mark, big.interval=big.interval, decimal.mark=decimal.mark)
   x.formatted <- as.matrix(x.formatted)
   x.formatted[is.na(x.formatted)] <- " "
   x.formatted[is.nan(x.formatted)] <- " "

   for(i in 1:dim(x)[1]) {
      if(i == 1) {
         VecDebut <- c(if(is.null(dimnames(x)[[1]]) == FALSE) paste(
              "\n<tr><td class=", classfirstcolumn, ">", sep = ""),
            paste("\n<td class=", classcellinside, ">", sep = ""),
            rep(paste("\n<td class=", classcellinside, ">", sep =
            ""), dim(x)[2] - 1))
         VecMilieu <- c(if(is.null(dimnames(x)[[1]]) == FALSE)
              dimnames(x)[[1]][i],
              HTMLReplaceNA(x.formatted[i,]))
         VecFin <- c(if(is.null(dimnames(x)[[1]]) == FALSE) "</td>",
            rep("</td>", dim(x)[2] - 1), "</td></tr>\n")
      }
      else {
         VecDebut <- c(if(is.null(dimnames(x)[[1]]) == FALSE) paste(
              "\n<tr><td class=", classfirstcolumn, ">", sep = ""),
            paste(rep(paste("\n<td class=", classcellinside, ">", sep
             = ""), dim(x)[2])))
         VecMilieu <- c(if(is.null(dimnames(x)[[1]]) == FALSE)
              dimnames(x)[[1]][i],
              HTMLReplaceNA(x.formatted[i,]))
         VecFin <- c(if(is.null(dimnames(x)[[1]]) == FALSE) "</td>",
            rep("</td>", dim(x)[2] - 1), "</td></tr>\n")
      }
      txt <- paste(txt, paste(VecDebut, VecMilieu, VecFin, sep = "",collapse = ""))
   }
   txt <- paste(txt, "\n\t</tbody>\n</table>\n",if (!is.null(Border)) "</td></table>\n","<br>")
   cat(txt, "\n", file = file, sep = "", append=TRUE)
   }

#----------------------------------------------------------------------------------------------------#

"HTML.structure"<-
function(x, a = attributes(x), prefix = "", file=HTMLGetFile(), append=TRUE, ...)
{
	cat("\n",file=file,append=append,...)
	n <- length(dim(x))
	nn <- names(a)
	ate <- character(0)
	if(n > 0) {
		if(n == 2)
			HTML.matrix(x, file = file,append=TRUE, ...)
		else HTML.array(x, file = file,append=TRUE, ...)
		ate <- c("dim", "dimnames")
		if(n == 1)
			ate <- c(ate, "names")
	}
	else if(!is.atomic(x)) {
		HTML(as.vector(x), file = file,append=TRUE, ...)
		ate <- "names"
	}
	else if(length(tsp(x))) {
		HTML.ts(x, file = file,append=TRUE, ...)
		ate <- "tsp"
	}
	else if(length(names(x))) {
		HTML.matrix(matrix(x, 1, dimnames = list("", names(x))),
			file = file,append=TRUE, ...)
		ate <- "names"
	}
	else HTML(as.vector(x), file = file,append=TRUE, ...)
	ii <- !match(nn, ate, nomatch = FALSE)
	nn <- nn[ii]
	a <- a[ii]
	for(i in seq(nn)) {
		this <- paste("attr(", prefix, ", \"", nn[i], "\")", sep = "")
		HTML(this, file=file,append=TRUE)
		HTML(a[[i]], file = file, append=TRUE, ...)
	}
	invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.connection" <- function(x,file=HTMLGetFile(), append=TRUE,...) HTML(paste("<font class='connection'>",unlist(summary(x)),"</font>",sep=""),file=file,append=append,...)

#----------------------------------------------------------------------------------------------------#

"HTML.socket" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    if (length(port <- as.integer(x$socket)) != 1)
        stop("invalid `socket' argument")
    HTML(paste("Socket connection #", x$socket, "to", x$host, "on port",
        x$port, "\n<br>",sep=""),file=file,append=append,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
"HTML.htest" <- function (x, digits = 4, quote = TRUE, prefix = "",file=HTMLGetFile(), append=TRUE, ...)
{
            HTML("\n", file=file,append=append)
            HTML(as.title(paste("&nbsp;",x$method,sep="")),file=file,append=TRUE,...)
            HTMLli(paste("\n data:<font class=dataname>",x$data.name,"</font>\n",sep=""),file=file,append=TRUE,...)
           out <- character()
            if (!is.null(x$statistic))
                        out <- c(out, paste(names(x$statistic), "=<b>", format(round(x$statistic,4)),"</b>"))
            if (!is.null(x$parameter))
                        out <- c(out, paste(names(x$parameter), "=<b>", format(round(x$parameter,3)),"</b>"))
            if (!is.null(x$p.value))
                        out <- c(out, paste("p-value =<font class='pvalue'>", format.pval(x$p.value,digits = digits),"</font>"))
            HTMLli(paste(out,collapse=" , "),file=file,append=TRUE,...)
    if (!is.null(x$alternative)) {
        HTMLli("alternative hypothesis: ",file=file)
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1) {
               alt.char <- switch(x$alternative, two.sided = "not equal to",
                  less = "less than", greater = "greater than")
                HTML(paste("true", names(x$null.value), "is", alt.char,
                 x$null.value, "\n"),file=file,append=TRUE,...)
            }
            else {
               HTMLli(paste(x$alternative, "\nnull values:\n<br>"),file=file,append=TRUE,...)
               HTML(x$null.value, file=file,append=TRUE,...)
            }
        }
        else HTML(paste(x$alternative, "\n<br>"),file=file,append=TRUE,...)
    }
    if (!is.null(x$conf.int)) {
        HTMLli(paste("<b>",format(100 * attr(x$conf.int, "conf.level")), "</b> percent confidence interval:\n",
         "<b>[", paste(format(c(x$conf.int[1], x$conf.int[2])),sep="",collapse=" ;"),"]</b>",sep=""),file=file,append=TRUE,...)
    }
    if (!is.null(x$estimate)) {
        HTMLli("sample estimates:\n",file=file,...)
        HTML(t(as.matrix(x$estimate)),file=file,...)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

 "HTML.aov" <- function (x, intercept = FALSE, tol = .Machine$double.eps^0.5, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file,append=append,...)
    if (!is.null(cl <- x$call))  HTMLli(paste("Call:\n<br><font class='call'>", deparse(cl)),"</font>",file=file)
    asgn <- x$assign[x$qr$pivot[1:x$rank]]
    effects <- x$effects
    if (!is.null(effects))
        effects <- as.matrix(effects)[seq(along = asgn), , drop = FALSE]
    rdf <- x$df.resid
    uasgn <- unique(asgn)
    nmeffect <- c("(Intercept)", attr(x$terms, "term.labels"))[1 + uasgn]
    nterms <- length(uasgn)
    nresp <- NCOL(effects)
    df <- numeric(nterms)
    ss <- matrix(NA, nterms, nresp)
    if (nterms) {
        for (i in seq(nterms)) {
            ai <- asgn == uasgn[i]
           df[i] <- sum(ai)
            ef <- effects[ai, , drop = FALSE]
            ss[i, ] <- if (sum(ai) > 1)
                colSums(ef^2)
            else ef^2       }
        keep <- df > 0
        if (!intercept && uasgn[1] == 0)
            keep[1] <- FALSE
        nmeffect <- nmeffect[keep]
        df <- df[keep]
        ss <- ss[keep, , drop = FALSE]
        nterms <- length(df)    }
    HTMLli("Terms:\n<br>",file=file)
    if (nterms == 0) {
        if (rdf > 0) {
            ss <- colSums(as.matrix(x$residuals)^2)
            ssp <- sapply(ss, format)
            if (!is.matrix(ssp))
                ssp <- t(ssp)
            tmp <- as.matrix(c(ssp, format(rdf)))
            if (length(ss) > 1) {
                rn <- colnames(x$fitted)
                if (is.null(rn))
                  rn <- paste("resp", 1:length(ss))
            }
            else rn <- "Sum of Squares"
            dimnames(tmp) <- list(c(rn, "Deg. of Freedom"), "Residuals")
            HTML(as.data.frame(tmp), file=file,..)
            HTMLli(paste("Residual standard error:", paste(sapply(sqrt(ss/rdf),format),collapse=" "), "\n"),file=file)
        }
        else HTML.matrix(matrix(0, 2, 1, dimnames = list(c("Sum of Squares","Deg. of Freedom"), "<empty>")),file=file)
    }
    else {
        if (rdf > 0) {
            resid <- as.matrix(x$residuals)
            nterms <- nterms + 1
            df <- c(df, rdf)
            ss <- rbind(ss, colSums(resid^2))
            nmeffect <- c(nmeffect, "Residuals")        }
        ssp <- apply(zapsmall(ss), 2, format)
        tmp <- t(cbind(ssp, format(df)))
        if (ncol(effects) > 1) {
            rn <- colnames(x$coef)
            if (is.null(rn))
                rn <- paste("resp", seq(ncol(effects)))        }
        else rn <- "Sum of Squares"
        dimnames(tmp) <- list(c(rn, "Deg. of Freedom"), nmeffect)
        HTML(as.data.frame(tmp), file=file)
       rank <- x$rank
        int <- attr(x$terms, "intercept")
        nobs <- NROW(x$residuals) - !(is.null(int) || int ==      0)
        if (rdf > 0) {
            rs <- sqrt(colSums(as.matrix(x$residuals)^2)/rdf)
            HTMLli(paste("Residual standard error:", paste(sapply(rs,format),collapse=" "), "\n"),file=file)       }
        coef <- as.matrix(x$coef)[, 1]
        R <- x$qr$qr
       R <- R[1:min(dim(R)), , drop = FALSE]
        R[lower.tri(R)] <- 0
        if (rank < (nc <- length(coef))) {
            HTMLli(paste(nc - rank, "out of", nc, "effects not estimable\n"),file=file)
            R <- R[, 1:rank, drop = FALSE]        }
        d2 <- sum(abs(diag(R)))
        diag(R) <- 0
        if (sum(abs(R))/d2 > tol)
            HTMLli("Estimated effects may be unbalanced\n",file=file)
        else HTMLli("Estimated effects are balanced\n",file=file)
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.anova" <- function (x, digits = max(getOption("digits") - 2, 3), signif.stars = getOption("show.signif.stars"),file=HTMLGetFile(), append=TRUE,...)
{
   cat("\n", file=file,append=append,...)
    if (!is.null(heading <- attr(x, "heading")))
        HTML(paste("<p><b>",heading, "</b></p>"),file=file)
   nc <- (d <- dim(x))[2]
    if (is.null(cn <- colnames(x)))
        stop("anova object must have colnames(.)!")
   ncn <- nchar(cn)
    has.P <- substr(cn[nc], 1, 3) == "Pr("
    zap.i <- 1:(if (has.P) nc - 1 else nc)
    i <- which(substr(cn, 2, 7) == " value")
    i <- c(i, which(!is.na(match(cn, c("FALSE", "Cp", "Chisq")))))
    if (length(i))
        zap.i <- zap.i[!(zap.i %in% i)]
    tst.i <- i
    if (length(i <- which(substr(cn, ncn - 1, ncn) == "Df")))
        zap.i <- zap.i[!(zap.i %in% i)]
    HTML.coefmat(x, digits = digits, signif.stars = signif.stars,
        has.Pvalue = has.P, P.values = has.P, cs.ind = NULL,
        zap.ind = zap.i, tst.ind = tst.i, na.print = "", file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.glm" <- function (x, digits = max(3, getOption("digits") - 3), na.print = "", file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file,append=append,...)
    HTMLli(paste("Call: <font class='call'>", deparse(x$call),"</font>", "\n<br>\n<br>"),file=file)
    HTMLli("Coefficients",file=file)
    if (is.character(co <- x$contrasts))
        HTML(paste("  [contrasts: ", apply(cbind(names(co), co), 1,
            paste, collapse = "="), "]"),file=file)
    HTMLbr(file=file)
    HTML(format(x$coefficients, digits = digits),file=file)
    HTMLli(paste("\nDegrees of Freedom:<b>", x$df.null, "</b>Total (i.e. Null);<b> ",
        x$df.residual, "</b> Residual\n"),file=file)
    HTMLli(paste("Null Deviance:<b>    ", format(signif(x$null.deviance,
        digits)), "</b> &nbsp;&nbsp; Residual Deviance:<b>", format(signif(x$deviance,
        digits)), " </b>&nbsp;&nbsp;    AIC:<b>  ", format(signif(x$aic, digits)), "</b>\n<br>"),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

 "HTML.tables.aov" <-  function (x, digits = 4, file=HTMLGetFile(),...)
 {
HTML("<center>",file=file)
     tables.aov <- x$tables
     n.aov <- x$n
     se.aov <- if (se <- !is.na(match("se", names(x))))
         x$se
     type <- attr(x, "type")
     switch(type, effects = HTML("<p class=partitle>Tables of effects\n</p>",file=file), means = HTML("<P CLASS=partitle>Tables of means\n</p>",file=file),
         residuals = if (length(tables.aov) > 1)
             HTML("<p class=partitle>Table of residuals from each stratum\n</p>",file=file))
     if (!is.na(ii <- match("Grand mean", names(tables.aov)))) {
         HTML("<p>Grand mean\n</p>",file=file)
         gmtable <- tables.aov[[ii]]
         HTML.mtable(gmtable, digits = digits, file=file)
     }
     for (i in names(tables.aov)) {
         if (i == "Grand mean")
             next
         table <- tables.aov[[i]]
         HTML(paste("\n<p>", i, "\n</p>"),file=file)
         if (!is.list(n.aov))
             HTML.mtable(table, digits = digits,file=file,append=TRUE, ...)
         else {
             n <- n.aov[[i]]
             if (length(dim(table)) < 2) {
                 table <- rbind(table, n)
                 rownames(table) <- c("", "rep")
                 HTML(table, digits = digits, file=file)
             }
             else {
                 ctable <- array(c(table, n), dim = c(dim(table),
                   2))
                 dim.t <- dim(ctable)
                 d <- length(dim.t)
                 ctable <- aperm(ctable, c(1, d, 2:(d - 1)))
                 dim(ctable) <- c(dim.t[1] * dim.t[d], dim.t[-c(1,
                   d)])
                 dimnames(ctable) <- c(list(format(c(rownames(table),
                   rep("rep", dim.t[1])))), dimnames(table)[-1])
                 ctable <- eval(parse(text = paste("ctable[as.numeric(t(matrix(seq(nrow(ctable)),ncol=2)))",
                   paste(rep(", ", d - 2), collapse = " "), "]")))
                 names(dimnames(ctable)) <- names(dimnames(table))
                 class(ctable) <- "mtable"
                 HTML.mtable(ctable, digits = digits,file=file, append=TRUE,...)
             }
         }
     }
     if (se) {
         if (type == "residuals")
             rn <- "df"
         else rn <- "replic."
         switch(attr(se.aov, "type"), effects = HTML("\n<p class=partitle>Standard errors of effects\n</p>",file=file),
             means = HTML("\n<p class=partitle>Standard errors for differences of means\n</p>",file=file),
             residuals = HTML("\n<p class=partitle>Standard errors of residuals\n</p>",file=file))
         if (length(unlist(se.aov)) == length(se.aov)) {
             n.aov <- n.aov[!is.na(n.aov)]
             se.aov <- unlist(se.aov)
             cn <- names(se.aov)
             se.aov <- rbind(format(se.aov, digits = digits),
                 format(n.aov))
             dimnames(se.aov) <- list(c(" ", rn), cn)
             HTML.matrix(se.aov,file=file)
         }
         else for (i in names(se.aov)) {
             se <- se.aov[[i]]
             if (length(se) == 1) {
                 se <- rbind(se, n.aov[i])
                 dimnames(se) <- list(c(i, rn), "")
                 HTML(se, file=file)
             }
             else {
                 dimnames(se)[[1]] <- ""
                 HTML(paste("\n<p>", i, "\n</p>"),file=file)
                 HTML("When comparing means with same levels of:\n<br>",file=file)
                 HTML(se, file=file, ...)
                 HTML(paste("replic.", n.aov[i], "\n<br>"),file=file)
             }
         }
     }
	HTML("</center>",file=file)
     invisible(x)
 }


#----------------------------------------------------------------------------------------------------#

"HTML.mtable" <- function (x, digits = getOption("digits"),file=HTMLGetFile(), append=TRUE,...)
{
   cat("\n", file=file,append=append,...)
    xxx <- x
    xx <- attr(x, "Notes")
    nn <- names(dimnames(x))
    a.ind <- match(names(a <- attributes(x)), c("dim", "dimnames",
        "names"))
    a <- a[!is.na(a.ind)]
    class(x) <- attributes(x) <- NULL
    attributes(x) <- a
    if (length(x) == 1 && is.null(names(x)) && is.null(dimnames(x)))
        names(x) <- rep("", length(x))
    if (length(dim(x)) && is.numeric(x)) {
        xna <- is.na(x)
        x <- format(zapsmall(x, digits))
        x[xna] <- "  "
    }
    HTML(x, file=file, ...)
    if (length(xx)) {
        HTML("\n<br>Notes:\n<br>",file=file)
        HTML(xx,file=file)
    }
    invisible(xxx)
}

#----------------------------------------------------------------------------------------------------#

"HTML.integrate" <- function (x, digits = getOption("digits"), file=HTMLGetFile(), append=TRUE,...)
{
   cat("\"n", file=file,append=append,...)
    if (x$message == "OK")
        HTML(paste("<p>",format(x$value, digits = digits), " with absolute error < ",
            format(x$abs.error, digits = 2), "\n</p>", sep = ""),file=file)
    else HTML(paste("<p>failed with message `", x$message, "'\n</p>", sep = ""),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.lm.null" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{

    cat("\"n", file=file,append=append,...)
    HTMLli(paste("<br><p>Call:<font class=call> ", paste(deparse(x$call), sep = "\n<br>", collapse = "\n<br>"), "</font></p>" ),file=file)
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
    if (rdf > 5) {
        HTMLli("Residuals:\n<br>",file=file)
        if (length(dim(resid)) == 2) {
            rq <- apply(t(resid), 1, quantile)
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q",
                "Max"), dimnames(resid)[[2]])
        }
        else {
            rq <- quantile(resid)
            names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        }
        HTML(round(rq, digits) ,file=file)
    }
    else if (rdf > 0) {
        HTMLli("Residuals:\n<br>",file=file)
        HTML(round(resid, digits ), file=file)
    }
    else HTMLli("\n<br>No Coefficients:\n<br>",file=file)
    HTMLli(paste("\n<br>Residual standard error:<b> ", format(signif(x$sigma,
        digits)), "on <b> ", rdf, " </b>degrees of freedom\n<br><br>",sep=""),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.glm" <- function (x, digits = max(3, getOption("digits") - 3), na.print = "",
    symbolic.cor = p > 4, signif.stars = getOption("show.signif.stars"), file=HTMLGetFile(), append=TRUE,
    ...)
{
    cat("\n", file=file,append=append,...)
    HTMLli(paste("\n<p>Call: <font class=call>",paste(deparse(x$call),collapse=" "),"</font>"),file=file)

    HTML("<p>Deviance Residuals: \n</p>",file=file)
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q",
            "Max")
    }
    HTML(t(round(x$deviance.resid,digits)) , file=file)
    HTML("\n<p>Coefficients:\n</p>",file=file)
    HTML.coefmat(x$coef, signif.stars = signif.stars, file=file)

    HTML(paste("\n<p>(Dispersion parameter for ", x$family$family, " family taken to be ",
        format(x$dispersion), ")\n</p>\n"),file=file)

       HTML(paste("<li>Null deviance:<b>", round(x$null.deviance,digits), "</b> on <b>", x[c("df.null")],"</b> degrees of freedom."),file=file)

       HTML(paste("<li>Residual deviance:<b>", round(x$deviance,digits), "</b> on <b>", x[c("df.residual")],"</b> degrees of freedom."),file=file)


       HTML(paste("<p>AIC:<b> ", format(x$aic, digits = max(4, digits + 1)), "</b>\n</p>\n<p>Number of Fisher Scoring iterations: <b>",     x$iter, "</b>\n</p>", sep = ""),file=file)
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            HTML("\n<p>Correlation of Coefficients:\n</p>")
            if (symbolic.cor)
                HTML(symnum(correl)[-1, -p],file=file)
            else {
                correl[!lower.tri(correl)] <- NA
                HTML(correl[-1, -p, drop = FALSE], file=file)
            }
        }
    }
    HTMLbr(file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.hsearch" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
   cat("\"n", file=file,append=append,...)
    fields <- paste(x$fields, collapse = " or ")
    db <- x$matches
    if (NROW(db) > 0) {
        HTML(paste("<p>Help files with ", fields, " matching `",
            x$pattern, "',\n", "type `help(FOO, package = PKG)' to inspect ",
            "entry `FOO(PKG) TITLE':", "\n</p>", sep = ""), file=file)
        dbnam <- paste(db[, "name"], "(", db[, "Package"], ")",sep = "")
        dbtit <- paste(db[, "title"], sep = "")
        HTML(cbind(dbnam, dbtit), file=file)
    }
    else HTML(paste("<p>No help files found with ", fields, " matching `", x$pattern, "'\n</p>", sep = ""),file=file)
}

#----------------------------------------------------------------------------------------------------#

"HTML.aov" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
NextMethod("HTML")
}

"HTML.aovlist" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
   cat("\"n", file=file,append=append,...)
    cl <- attr(x, "call")
    if (!is.null(cl)) {
        cat("\nCall:\n<font class=call>",file=file,append=TRUE,...)
        dput(cl,file=file)
        cat("\n</font>",file=file,append=TRUE,...)
    }
    if (!is.null(attr(x, "weights")))
        cat("Note: The results below are on the weighted scale\n",file=file,append=TRUE,...)
    nx <- names(x)
    if (nx[1] == "(Intercept)") {
        mn <- x[[1]]$coef
        if (is.matrix(mn)) {
            cat("\nGrand Means:\n",file=file,append=TRUE,...)
            cat(format(mn[1, ]), file=file,append=TRUE,...)
        }
        else cat("\nGrand Mean:", format(mn[1]), "\n",file=file,append=TRUE,...)
        nx <- nx[-1]
    }
    for (ii in seq(along = nx)) {
        i <- nx[ii]
        cat("\nStratum ", ii, ": ", i, "\n", sep = "",file=file,append=TRUE,...)
        xi <- x[[i]]
        cat(xi,file=file,append=TRUE, ...)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.SavedPlots" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
	cat("\"n",file=file,append=append,...)
    if (x[[1]] != 31416) {
        HTML("<p>object is not of class `SavedPlots'</p>\n<br>",file=file)
        return()
    }
    HTML("<p>Saved Plots from R version 1.4.0 or later</p>\n<br>\n<br>",file=file,append=TRUE,...)
    HTML("  Contains", x[[2]], "out of a maximum", x[[3]], "plots\n",file=file,append=TRUE,...)
    lens <- sapply(x[[5]], length)[1:x[[2]]]
    cat("  #plot calls are", paste(lens, collapse = ", "), "\n",file=file,append=TRUE,...)
    cat("  Current position is plot", 1 + x[[4]], "\n",file=file,append=TRUE,...)
}

#----------------------------------------------------------------------------------------------------#

"HTML.ordered" <- function (x, quote = FALSE,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (length(x) <= 0)
        HTML("\n<p>ordered(0)\n</p>",file=file,append=TRUE,...)
    else HTML(as.character(x), file,file, append=TRUE,...)
    HTML(paste("\n<p>Levels<font class=factorlevels>: ", paste(levels(x), collapse = " < "), "</font>\n</p>"),file=file,append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.difftime" <- function (x, digits = getOption("digits"),file=HTMLGetFile(), append=TRUE, ...)
{
    cat("\n",file=file,append=append,...)
    if (length(x) > 1)
        HTML(paste("<p>Time differences of ", paste(format(unclass(x),
            digits = digits), collapse = ", "), " ", attr(x,
            "units"), "\n</p>", sep = ""),file=file,append=TRUE,...)
    else HTML(paste("<p>Time difference of ", format(unclass(x), digits = digits),
        " ", attr(x, "units"), "\n", sep = ""),file=file,append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.dummy.coef" <- function (x, file=HTMLGetFile(), append=TRUE,title="",...)
{
    cat("\n",file=file,append=append,...)
    terms <- names(x)
    n <- length(x)
    nm <- max(sapply(x, length))
    ans <- matrix("", 2 * n, nm)
    rn <- rep("", 2 * n)
    line <- 0
    for (j in seq(n)) {
        this <- x[[j]]
        n1 <- length(this)
        if (n1 > 1) {
            line <- line + 2
            ans[line - 1, 1:n1] <- names(this)
            ans[line, 1:n1] <- format(this, ...)
            rn[line - 1] <- paste(terms[j], ":   ", sep = "")
        }
        else {
            line <- line + 1
            ans[line, 1:n1] <- format(this, ...)
            rn[line] <- paste(terms[j], ":   ", sep = "")
        }
    }
    rownames(ans) <- rn
    colnames(ans) <- rep("", nm)
    HTML(paste("\n<p>",if (title=="")
        "Full coefficients are"
    else title, "\n</p>"),file=file,append=TRUE,...)
    HTML.matrix(ans[1:line, , drop = FALSE],file=file,append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.dummy.coef.list" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    for (strata in names(x)) HTML.dummy.coef(x[[strata]], file=file, title = paste("\n<p>     Error:", strata,"</p>"),append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

 "HTML.glm.null" <- function (x, digits = max(3, getOption("digits") - 3), na.print = "",
    file=HTMLGetFile(), append=TRUE,...)
{

      cat("\n",file=file,append=append,...)
    HTMLli(paste(" Call: <font class='call'>", deparse(x$call),"</font>", "\n<br>\n"),file=file)
    HTMLli("No coefficients\n<br>")
    HTMLli(paste("Degrees of Freedom:<b>", length(x$residuals), "</b> Total; <b>",
        x$df.residual, " </b>Residual\n<br>"),file=file)
    HTMLli(paste("Null Deviance:<b>", format(signif(x$null.deviance, digits)),
        "</b>\n<br>"),file=file)
    HTMLli(paste("Residual Deviance: <b>", format(signif(x$deviance, digits)),
        " </b><br>\n"),file=file)
    HTMLli(paste("AIC:<b>", format(signif(x$aic, digits)), "</b><br>\n"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.MethodsFunction"<- function (x,file=HTMLGetFile(), append=TRUE, ...)
{
    	cat("\n",file=file,append=append,...)
	info=attr(x,"info")
	if (dim(info)[1]==0) HTML("<p>No available generic function for the class",file=file,append=TRUE)
	HTML("<p>Available generic functions which does handle the class</p>",file=file,append=TRUE)
	HTML(info,file=file,append=TRUE,...)
	invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.libraryIQR" <- function (x,file=HTMLGetFile(), append=TRUE, ...)
{
    cat("\n",file=file,append=append,...)
    sQuote <- function(s) paste("`", s, "'", sep = "")
    db <- x$results
    out <- if (nrow(db) == 0)
        NULL
    else lapply(split(1:nrow(db), db[, "LibPath"]), function(ind) db[ind,
        c("Package", "Title"), drop = FALSE])
    first <- TRUE
    for (lib in names(out)) {
        HTML(paste(paste("<p>Packages in library ",
            sQuote(lib), ":</p>", sep = "")),file=file,append=TRUE,...)
        HTML(cbind(out[[lib]][, "Package"], out[[lib]][,
            "Title"]), file=file,append=TRUE,...)
        first <- FALSE
    }
    if (first) {
        HTML("<p>no packages found</p>",file=file, append=TRUE,...)    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.aov" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
      cat("\n",file=file,append=append,...)
    if (length(x) == 1)
        HTML(x[[1]], file=file)
    else NextMethod()
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.aovlist" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    nn <- names(x)
    for (i in nn) {
        HTMLli(paste(i, "\n<br>", sep = ""),file=file)
        HTML(x[[i]], file=file)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.glm.null" <- function (x, digits = max(3, getOption("digits") - 3), na.print = "",
    file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("\nCall:<font class=call> ",paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "</font>\n<br>\n", sep = ""),file=file)
    HTMLli("Deviance Residuals: \n<br>",file=file)
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q",
            "Max")
    }
    HTML.default(x$deviance.resid, digits = digits, na = "",file=file)
    HTMLli("No coefficients\n<br>")
    HTMLli(paste("\n(Dispersion parameter for ", x$family$family,
        " family taken to be ", x$dispersion, ")\n\n    Null deviance:<b> ",
        x$null.deviance, " </b>on <b>", x$df.null, " </b>degrees of freedom\n\n",
        "Residual deviance: <b>", x$deviance, " </b>on<b> ", x$df.residual,
        " </b>degrees of freedom\n\n", "Number of Fisher Scoring iterations<b>: ",
        x$iter, "</b>\n<br>\n", sep = ""),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.manova" <- function (x, digits = getOption("digits"),file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (length(stats <- x$stats)) {
        HTML.anova(stats,file=file)
    }
    else {
        HTML("<p>No error degrees of freedom</p>\n")
        HTML(data.frame(Df = x$Df, row.names = x$row.names),file=file)
    }
    invisible(x)
}



#----------------------------------------------------------------------------------------------------#

"HTML.summary.table" <- function (x, digits = max(1, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!inherits(x, "summary.table"))
        stop("x must inherit from class `summary.table'")
    if (!is.null(x$call)) {
        HTMLli(paste("Call:<font class='call'> ", x$call,"</font>"),file=file)
    }
    HTMLli(paste("Number of cases in table:<b>", x$n.cases, "</b>\n<br>"),file=file)
    HTMLli(paste("Number of factors:<b>", x$n.vars, "</b>\n<br>"),file=file)
    if (x$n.vars > 1) {
        HTMLli("Test for independence of all factors:\n<br>",file=file)
        ch <- x$statistic
        HTML(paste(" Chisq = <b>", format(round(ch, max(0, digits - log10(ch)))),
            "</b>, df = <b>", x$parameter, "</b>, p-value = <b>", format.pval(x$p.value,
                digits, eps = 0), "</b>\n<br>", sep = ""),file=file)
        if (!x$approx.ok)
            HTML("<p>Chi-squared approximation may be incorrect</p>\n",file=file)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#
"HTML.TukeyHSD" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<center><p><b>Tukey multiple comparisons of means</b></p>\n")
    HTML(paste("<p>", format(100 * attr(x, "conf.level"), 2), "% family-wise confidence level</p></center>\n",
        sep = ""),file=file)

    if (attr(x, "ordered"))
        HTML("<center><p>factor levels have been ordered</p></center>\n",file=file)
    HTMLli(paste("Fit: ", deparse(attr(x, "orig.call")), "\n<br>\n", sep = ""),file=file)
    attr(x, "orig.call") <- attr(x, "conf.level") <- attr(x, "ordered") <- NULL
	lapply(unclass(x),HTML,file=file,append=TRUE,...)
    #HTML.default(unclass(x), file=file,...)
    invisible(return(x))
}


#----------------------------------------------------------------------------------------------------#

"HTML.simple.list" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
	HTML(noquote(cbind("<-" = unlist(x))), file=file,append=TRUE,...)
}

#----------------------------------------------------------------------------------------------------#

"HTML.noquote" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- attr(x, "class"))) {
        cl <- cl[cl != "noquote"]
        attr(x, "class") <- (if (length(cl) > 0)
            cl
        else NULL)
    }
    HTML(x, file=file, append=TRUE,...)
}



###
### PACKAGES FUNCTIONS
###


### PACKAGE TS

#----------------------------------------------------------------------------------------------------#

"HTML.ar" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("Call:\n<font class='call'>", deparse(x$call), "</font>\n", sep = ""),file=file)
    nser <- NCOL(x$var.pred)
    if (nser > 1) {
        if (!is.null(x$x.intercept))
            res <- x[c("ar", "x.intercept", "var.pred")]
        else res <- x[c("ar", "var.pred")]
        res$ar <- aperm(res$ar, c(2, 3, 1))
        HTML(res, digits = digits,file=file)
    }
    else {
        if (x$order > 0) {
            HTMLli("Coefficients:\n",file=file)
            coef <- drop(round(x$ar, digits = digits))
            names(coef) <- seq(length = x$order)
            HTML.default(coef, file=file)
        }
        if (!is.null(xint <- x$x.intercept) && !is.na(xint))
            HTML(paste("<p>Intercept: <b>", format(xint, digits = digits),
                "</b> (", format(x$asy.se.coef$x.mean, digits = digits),
                ") ", "\n</p>", sep = ""),file=file)
        HTML(paste("<p>Order selected <b>", x$order, " </b>sigma^2 estimated as <b>",
            format(x$var.pred, digits = digits), "</b>\n<</p>"),file=file)
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.Arima" <- function (x, digits = max(3, getOption("digits") - 3), se = TRUE,
    file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("nCall:<font class='call'>", deparse(x$call, width.cutoff = 75), "</font>", sep = "\n"),file=file)
    HTMLli("Coefficients:\n<br>",file=file)
    coef <- round(x$coef, digits = digits)
    if (se && nrow(x$var.coef)) {
        ses <- rep(0, length(coef))
        ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
        coef <- matrix(coef, 1, dimnames = list(NULL, names(coef)))
        coef <- rbind(coef, s.e. = ses)
    }
    HTML.default(coef,file=file)
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS")
        HTML(paste("\n<p>sigma^2 estimated as <b>", format(x$sigma2, digits = digits),
            "</b>:  log likelihood = <b>", format(round(x$loglik, 2)),
            "</b>,  aic = <b>", format(round(x$aic, 2)), "</b>\n</p>", sep = ""),file=file)
    else HTML("<p>sigma^2 estimated as <b>", format(x$sigma2, digits = digits),
        "</b>:  part log likelihood =<b> ", format(round(x$loglik, 2)),
        "</b>\n</p>", sep = "")
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.arima0" <- function (x, digits = max(3, getOption("digits") - 3), se = TRUE,
    file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("\nCall:<font class='call'>", deparse(x$call, width.cutoff = 75), "</font>", sep = "\n"),file=file)
    HTMLli("Coefficients:\n<br>",file=file)
    coef <- round(x$coef, digits = digits)
    if (se && nrow(x$var.coef)) {
        ses <- rep(0, length(coef))
        ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
        coef <- matrix(coef, 1, dimnames = list(NULL, names(coef)))
        coef <- rbind(coef, s.e. = ses)
    }
    HTML.default(coef, file=file)
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS")
        HTML(paste("\n<p>sigma^2 estimated as <b>", format(x$sigma2, digits = digits),
            "</b>:  log likelihood = <b>", format(round(x$loglik, 2)),
            "</b>,  aic = <b>", format(round(x$aic, 2)), "</b>\n</p>", sep = ""),file=file)
    else HTML(paste("\n<p>sigma^2 estimated as <b>", format(x$sigma2, digits = digits),
        "</b>:  part log likelihood =<b> ", format(round(x$loglik, 2)),
        "</b>\n</p>", sep = ""),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.HoltWinters" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(paste("<p><b>Holt-Winters exponential smoothing", if (x$beta == 0)
        "without"
    else "with", "trend and", if (x$gamma == 0)
        "without"
    else paste(if (x$beta == 0)
        "with ", x$seasonal, sep = ""), "seasonal componenent.\n</b></p>"),file=file)

    HTMLli(paste("\nCall:\n", deparse(x$call), "\n<br>"),file=file)
    HTMLli("Smoothing parameters:\n<ul>",file=file)
    HTMLli(paste(" alpha: ", x$alpha, "\n"),file=file)
    HTMLli(paste(" beta: ", x$beta, "\n"),file=file)
    HTMLli(paste(" gamma: ", x$gamma, "\n<br>"),file=file)
    HTML("</ul>",file=file)
    HTMLli("Coefficients:\n",file=file)
    HTML(t(t(x$coefficients)),file=file)
}


#----------------------------------------------------------------------------------------------------#

"HTML.stl" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("Call:\n ",deparse(x$call),"\n<br>"),file=file)
    HTMLli("\nComponents\n",file=file)
    HTML(x$time.series, file=file,append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.StructTS" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("\nCall:", deparse(x$call, width.cutoff = 75), "\n", sep = " "),file=file)
    HTMLli("Variances:\n",file=file)
    HTML(x$coef,  digits=digits,file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.tskernel" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    y <- c(rev(x$coef[2:(x$m + 1)]), x$coef)
    i <- -x$m:x$m
    HTML(paste("<p>",attr(x, "name"), "</p>\n"),file=file)
    HTML(paste( paste("coef[", format(i), "] = ", format(y, digits = digits),sep = ""),collapse="<br>\n", sep = "\n<br>"),file=file)
}


### PACKAGE CTEST

#----------------------------------------------------------------------------------------------------#

"HTML.pairwise.htest" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("Pairwise comparisons using", x$method, "\n<br>\n<br>"),file=file)
    HTMLli(paste("data: <font class=dataname>", x$data.name,"</font>", "\n<br>\n<br>"),file=file)
    pp <- format.pval(x$p.value, 2, na.form = "-")
    attributes(pp) <- attributes(x$p.value)
    HTML(pp, file=file)
    HTMLli(paste("\nP value adjustment method:", x$p.adjust.method, "\n"),file=file)
}

#----------------------------------------------------------------------------------------------------#

"HTML.power.htest" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste(x$method,"<br>"), file=file)
    note <- x$note
    x[c("method", "note")] <- NULL
    HTML(paste(paste(formatC(names(x), width = 15, flag = "+"),
        format(x), sep = " = 	"), sep = "\n<br>",collapse="\n<br>"),file=file)
    if (!is.null(note))
        HTML(paste("\n<p>", "NOTE:", note, "\n</p>\n"),file=file)
    else HTMLbr(file=file)
}


#----------------------------------------------------------------------------------------------------#

"HTML.boot" <- function (x, digits = options()$digits, index = 1:ncol(boot.out$t), file=HTMLGetFile(),  append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    boot.out <- x
    sim <- boot.out$sim
    cl <- boot.out$call
    t <- matrix(boot.out$t[, index], nrow = nrow(boot.out$t))
    allNA <- apply(t, 2, function(t) all(is.na(t)))
    ind1 <- index[allNA]
    index <- index[!allNA]
    t <- matrix(t[, !allNA], nrow = nrow(t))
    rn <- paste("t", index, "*", sep = "")
    if (length(index) == 0)
        op <- NULL
    else if (is.null(t0 <- boot.out$t0)) {
        if (is.null(boot.out$call$weights))
            op <- cbind(apply(t, 2, mean, na.rm = TRUE), sqrt(apply(t,
                2, function(t.st) var(t.st[!is.na(t.st)]))))
        else {
            op <- NULL
            for (i in index) op <- rbind(op, boot::imp.moments(boot.out,
                index = i)$rat)
            op[, 2] <- sqrt(op[, 2])
        }
        dimnames(op) <- list(rn, c("mean", "std. error"))
    }
    else {
        t0 <- boot.out$t0[index]
        if (is.null(boot.out$call$weights)) {
            op <- cbind(t0, apply(t, 2, mean, na.rm = TRUE) -
                t0, sqrt(apply(t, 2, function(t.st) var(t.st[!is.na(t.st)]))))
            dimnames(op) <- list(rn, c("original", " bias  ",
                " std. error"))
        }
        else {
            op <- NULL
            for (i in index) op <- rbind(op, boot::imp.moments(boot.out,
                index = i)$rat)
            op <- cbind(t0, op[, 1] - t0, sqrt(op[, 2]), apply(t,
                2, mean, na.rm = TRUE))
            dimnames(op) <- list(rn, c("original", " bias  ",
                " std. error", " mean(t*)"))
        }
    }
    if (cl[[1]] == "boot") {
        if (sim == "parametric")
            HTML(as.title("PARAMETRIC BOOTSTRAP"),file=file)
        else if (sim == "antithetic") {
            if (is.null(cl$strata))
		HTML(as.title("ANTITHETIC BOOTSTRAP"),file=file)
            else
            HTML(as.title("STRATIFIED ANTITHETIC BOOTSTRAP"),file=file)

        }
        else if (sim == "permutation") {
            if (is.null(cl$strata))
		HTML(as.title("DATA PERMUTATION"),file=file)
           else HTML(as.title("STRATIFIED DATA PERMUTATION"),file=file)
        }
        else if (sim == "balanced") {
            if (is.null(cl$strata) && is.null(cl$weights))
                HTML(as.title("BALANCED BOOTSTRAP"),file=file)
            else if (is.null(cl$strata))
                HTML(as.title("BALANCED WEIGHTED BOOTSTRAP"),file=file)
            else if (is.null(cl$weights))
		HTML(as.title("STRATIFIED BALANCED BOOTSTRAP"),file=file)
            else HTML(as.title("STRATIFIED WEIGHTED BALANCED BOOTSTRAP"),file=file)
        }
        else {
            if (is.null(cl$strata) && is.null(cl$weights))
		HTML(as.title("ORDINARY NONPARAMETRIC BOOTSTRAP"),file=file)
            else if (is.null(cl$strata))
 		HTML(as.title("WEIGHTED BOOTSTRAP"),file=file)
             else if (is.null(cl$weights))
		HTML(as.title("STRATIFIED BOOTSTRAP"),file=file)
                else HTML(as.title("STRATIFIED WEIGHTED BOOTSTRAP"),file=file)
        }
    }
    else if (cl[[1]] == "tilt.boot") {
        R <- boot.out$R
        th <- boot.out$theta
        if (sim == "balanced")
		HTML(as.title("BALANCED TITLED BOOTSTRAP"),file=file)
        else HTML(as.title("TILTED BOOTSTRAP"),file=file)
        if ((R[1] == 0) || is.null(cl$tilt) || eval(cl$tilt))
            HTML("<p>Exponential tilting used\n</p>",file=file)
        else HTML("<p>Frequency Smoothing used\n</p>",file=file)
        i1 <- 1
        if (boot.out$R[1] > 0)
            HTML(paste("<p>First", R[1], "replicates untilted,\n</p>"),file=file)
        else {
            HTML(paste("<p>First ", R[2], " replicates tilted to ",
                signif(th[1], 4), ",\n</p>", sep = ""),file=file)
            i1 <- 2
        }
        if (i1 <= length(th)) {
            for (j in i1:length(th)) HTML(paste("<p>Next ", R[j +
                1], " replicates tilted to ", signif(th[j], 4),
                ifelse(j != length(th), ",\n", ".\n</p>"), sep = ""),file=file)
        }
        op <- op[, 1:3]
    }
    else if (cl[[1]] == "tsboot") {
        if (!is.null(cl$indices))
		HTML(as.title("TIME SERIES BOOTSTRAP USING SUPPLIED INDICES"),file=file)
            else if (sim == "model")
            HTML(as.title("MODEL BASED BOOTSTRAP FOR TIME SERIES"),file=file)
        else if (sim == "scramble") {
		HTML(as.title("PHASE SCRAMBLED BOOTSTRAP FOR TIME SERIES"),file=file)
            if (boot.out$norm)
                HTML("<p>Normal margins used.\n</p>",file=file)
            else HTML("<p>Observed margins used.\n</p>",file=file)
        }
        else if (sim == "geom") {
            if (is.null(cl$ran.gen))
                HTML(as.title("STATIONARY BOOTSTRAP FOR TIME SERIES"),file=file)
            else  HTML(as.title("POST-BLACKENED STATIONARY BOOTSTRAP FOR TIME SERIES"),file=file)
		HTML(paste("<p>Average Block Length of", boot.out$l,
                "\n</p>"),file=file)
        }
        else {
            if (is.null(cl$ran.gen))
		HTML("<p>BLOCK BOOTSTRAP FOR TIME SERIES</p>",file=file)
            else HTML("<p>POST-BLACKENED BLOCK BOOTSTRAP FOR TIME SERIES</p>",file=file)
            HTML(paste("<p>Fixed Block Length of", boot.out$l, "\n</p>"),file=file)
        }
    }
    else {
        cat("\n")
        if (sim == "weird") {
            if (!is.null(cl$strata))
                HTML(as.title("STRATIFIED BOOTSTRAP FOR CENSORED DATA"),file=file)
       }
        else if ((sim == "ordinary") || ((sim == "model") &&
            is.null(boot.out$cox))) {
            if (!is.null(cl$strata))
 		 HTML(as.title("STRATIFIED CASE RESAMPLING BOOTSTRAP FOR CENSORED DATA"),file=file)
        }
        else if (sim == "model") {
            if (!is.null(cl$strata))

		HTML(as.title("STRATIFIED MODEL BASED BOOTSTRAP FOR COX REGRESSION MODEL"),file=file)
        }
        else if (sim == "cond") {
            if (!is.null(cl$strata))
 	HTML(as.title("STRATIFIED CONDITIONAL BOOTSTRAP"),file=file)
            if (is.null(boot.out$cox))
                HTML("<p>FOR CENSORED DATA\n</p>\n",file=file)
            else HTML("<p>FOR COX REGRESSION MODEL\n</p>\n",file=file)
        }
    }
    HTMLli(paste("\nCall: ",deparse(cl)),file=file)

    HTMLli("Bootstrap Statistics :\n<br>",file=file)
    if (!is.null(op))
        HTML(op, digits = digits,file=file)
    if (length(ind1) > 0)
        for (j in ind1) HTML(paste("<p>WARNING: All values of t",
            j, "* are NA\n</p>", sep = ""),file=file)
    invisible(boot.out)
}

#----------------------------------------------------------------------------------------------------#

"HTML.simplex" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    simp.out <- x
    HTML("\n<p><b>Linear Programming Results\n</b></p>\n",file=file)
    cl <- simp.out$call
    HTMLli(paste("Call : ",deparse(cl)),file=file)
	HTML(paste("<p>", if (simp.out$maxi) "Maximization" else "Minimization", " Problem with Objective Function Coefficients\n</p>"),file=file)
    HTML(simp.out$obj,file=file)
    if (simp.out$solved == 1) {
        HTML("\n<p>\nOptimal solution has the following values\n</p>",file=file)
        HTML(simp.out$soln,file=file)
        HTML(paste("<p>The optimal value of the objective ", " function is ",
            simp.out$value, ".\n</p>", sep = ""),file=file)
    }
    else if (simp.out$solved == 0) {
        HTML("\n<p>\nIteration limit exceeded without finding solution\n</p>",file=file)
        HTML("<p>The coefficient values at termination were\n</p>",file=file)
        HTML(simp.out$soln,file=file)
        HTML(paste("<p>The objective function value was ", simp.out$value,
            ".\n</p>", sep = ""),file=file)
    }
    else HTML("\n<p>No feasible solution could be found\n</p>",file=file)
}

#----------------------------------------------------------------------------------------------------#

"HTML.saddle.distn" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    sad.d <- x
    cl <- sad.d$call
    rg <- range(sad.d$points[, 1])
    mid <- mean(rg)
    digs <- ceiling(log10(abs(mid)))
    if (digs <= 0)
        digs <- 4
    else if (digs >= 4)
        digs <- 0
    else digs <- 4 - digs
    rg <- round(rg, digs)
    level <- 100 * sad.d$quantiles[, 1]
    quans <- format(round(sad.d$quantiles, digs))
    quans[, 1] <- paste( format(level), "%     ", sep = "")
    HTML("\n<p><b>Saddlepoint Distribution Approximations\n</b></p>\n",file=file)
    HTMLli(paste("Call : ",paste(deparse(cl),collapse="")),file=file)
    HTML("\n<p>Quantiles of the Distribution\n</p>",file=file)
    HTML(t(t(quans)),file=file)
    HTML(paste("\n<p>\nSmoothing spline used ", nrow(sad.d$points),
        " points in the range ", rg[1], " to ", rg[2], ".</p>", sep = ""),file=file)
    if (sad.d$LR)
        HTMLli("Lugananni-Rice approximations used.",file=file)
       HTMLbr(file=file)
    invisible(sad.d)
}

#----------------------------------------------------------------------------------------------------#

"HTML.bootci" <- function (x, hinv = NULL, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    ci.out <- x
    cl <- ci.out$call
    ntypes <- length(ci.out) - 3
    nints <- nrow(ci.out[[4]])
    t0 <- ci.out$t0
    if (!is.null(hinv))
        t0 <- hinv(t0)
    digs <- ceiling(log10(abs(t0)))
    if (digs <= 0)
        digs <- 4
    else if (digs >= 4)
        digs <- 0
    else digs <- 4 - digs
    intlabs <- NULL
    basrg <- strg <- perg <- bcarg <- NULL
    if (!is.null(ci.out$normal))
        intlabs <- c(intlabs, "     Normal        ")
    if (!is.null(ci.out$basic)) {
        intlabs <- c(intlabs, "     Basic         ")
        basrg <- range(ci.out$basic[, 2:3])
    }
    if (!is.null(ci.out$student)) {
        intlabs <- c(intlabs, "   Studentized     ")
        strg <- range(ci.out$student[, 2:3])
    }
    if (!is.null(ci.out$percent)) {
        intlabs <- c(intlabs, "    Percentile     ")
        perg <- range(ci.out$percent[, 2:3])
    }
    if (!is.null(ci.out$bca)) {
        intlabs <- c(intlabs, "      BCa          ")
        bcarg <- range(ci.out$bca[, 2:3])
    }
    level <- 100 * ci.out[[4]][, 1]
    if (ntypes == 4)
        n1 <- n2 <- 2
    else if (ntypes == 5) {
        n1 <- 3
        n2 <- 2
    }
    else {
        n1 <- ntypes
        n2 <- 0
    }
    ints1 <- matrix(NA, nints, 2 * n1 + 1)
    ints1[, 1] <- level
    n0 <- 4
    for (i in n0:(n0 + n1 - 1)) {
        j <- c(2 * i - 6, 2 * i - 5)
        nc <- ncol(ci.out[[i]])
        nc <- c(nc - 1, nc)
        if (is.null(hinv))
            ints1[, j] <- ci.out[[i]][, nc]
        else ints1[, j] <- hinv(ci.out[[i]][, nc])
    }
    n0 <- 4 + n1
    ints1 <- format(round(ints1, digs))
    ints1[, 1] <- paste("\n<br>", level, "%  ", sep = "")
    ints1[, 2 * (1:n1)] <- paste("(", ints1[, 2 * (1:n1)], ",",
        sep = "")
    ints1[, 2 * (1:n1) + 1] <- paste(ints1[, 2 * (1:n1) + 1],
        ")  ")
    if (n2 > 0) {
        ints2 <- matrix(NA, nints, 2 * n2 + 1)
        ints2[, 1] <- level
        j <- c(2, 3)
        for (i in n0:(n0 + n2 - 1)) {
            if (is.null(hinv))
                ints2[, j] <- ci.out[[i]][, c(4, 5)]
            else ints2[, j] <- hinv(ci.out[[i]][, c(4, 5)])
            j <- j + 2
        }
        ints2 <- format(round(ints2, digs))
        ints2[, 1] <- paste("\n<br>", level, "%  ", sep = "")
        ints2[, 2 * (1:n2)] <- paste("(", ints2[, 2 * (1:n2)],
            ",", sep = "")
        ints2[, 2 * (1:n2) + 1] <- paste(ints2[, 2 * (1:n2) +
            1], ")  ")
    }
    R <- ci.out$R
    HTML(as.title("BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS"),file=file)
    HTML(paste("<p>Based on", R, "bootstrap replicates\n\n</p>"),file=file)
    HTMLli(paste("CALL : ",paste(deparse(cl),collapse=" ")),file=file)
    HTML("\n<p>Intervals : </p>",file=file)
    HTML(paste("\n<p>Level", intlabs[1:n1],"</p>"),file=file)
    HTML(t(ints1),file=file)
    if (n2 > 0) {
        HTML(paste("\n<p>\nLevel", intlabs[(n1 + 1):(n1 + n2)],"</p>"),file=file)
        HTML(t(ints2),file=file)
    }
    if (!is.null(cl$h)) {
        if (is.null(cl$hinv) && is.null(hinv))
            HTML("\n<p>Calculations and Intervals on Transformed Scale\n</p>",file=file)
        else HTML("\n<p>Calculations on Transformed Scale;  Intervals on Original Scale\n</p>",file=file)
    }
    else if (is.null(cl$hinv) && is.null(hinv))
        HTML("\n<p>Calculations and Intervals on Original Scale\n</p>",file=file)
    else HTML("\n<p>Calculations on Original Scale but Intervals Transformed\n</p>",file=file)
    if (!is.null(basrg)) {
        if ((basrg[1] <= 1) || (basrg[2] >= R))
            HTML("\n<p>Warning : Basic Intervals used Extreme Quantiles\n</p>",file=file)
        if ((basrg[1] <= 10) || (basrg[2] >= R - 9))
            HTML("\n<p>Some basic intervals may be unstable\n</p>",file=file)
    }
    if (!is.null(strg)) {
        if ((strg[1] <= 1) || (strg[2] >= R))
            HTML("\n<p>Warning : Studentized Intervals used Extreme Quantiles\n</p>",file=file)
        if ((strg[1] <= 10) || (strg[2] >= R - 9))
            HTML("\n<p>Some studentized intervals may be unstable\n</p>",file=file)
    }
    if (!is.null(perg)) {
        if ((perg[1] <= 1) || (perg[2] >= R))
            HTML("\n<p>Warning : Percentile Intervals used Extreme Quantiles\n</p>",file=file)
        if ((perg[1] <= 10) || (perg[2] >= R - 9))
            HTML("\n<p>Some percentile intervals may be unstable\n</p>",file=file)
    }
    if (!is.null(bcarg)) {
        if ((bcarg[1] <= 1) || (bcarg[2] >= R))
            HTML("\n<p>Warning : BCa Intervals used Extreme Quantiles\n</p>",file=file)
        if ((bcarg[1] <= 10) || (bcarg[2] >= R - 9))
            HTML("\n<p>Some BCa intervals may be unstable\n</p>",file=file)
    }
    invisible(ci.out)
}


#----------------------------------------------------------------------------------------------------#

### PACKAGE MVA (merged into stats)

#----------------------------------------------------------------------------------------------------#

"HTML.dist" <- function (x, diag = NULL, upper = NULL, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (is.null(diag))
        diag <- if (is.null(a <- attr(x, "Diag")))
            FALSE
        else a
    if (is.null(upper))
        upper <- if (is.null(a <- attr(x, "Upper")))
            FALSE
        else a
    size <- attr(x, "Size")
    df <- as.matrix(x)
    if (!upper)
        df[row(df) < col(df)] <- NA
    if (!diag)
        df[row(df) == col(df)] <- NA
    HTML(if (diag || upper)
        df
    else df[-1, -size], file=file, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.factanal" <- function (x, digits = 3, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("\nCall:\n", deparse(x$call), "\n<br>\n", sep = ""),file=file)
    HTMLli("Uniquenesses:\n<br>",file=file)
    HTML(round(x$uniquenesses, digits),file=file,append=TRUE,...)
    HTML(x$loadings, digits = digits,file=file,append=TRUE, ...)
    p <- nrow(x$loadings)
    factors <- x$factors
    if (!is.na(x$n.obs) && x$dof > 0) {
        dof <- x$dof
        stat <- (x$n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3) *
            x$criteria["objective"]
        HTMLli(paste("\n<p>Test of the hypothesis that", factors, if (factors ==
            1)
            "factor is"
        else "factors are", "sufficient.\n</p>"),file=file)
        HTML(paste("<p>The chi square statistic is <b>", round(stat, 2), " </b> on <b>",
            dof, if (dof == 1)
                " </b>degree"
            else "</b>degrees", "of freedom.\n<br>The p-value is <b>", signif(pchisq(stat,
                dof, lower.tail = FALSE), 3), "</b>\n</p>"),file=file)
    }
    else {
        HTML(paste("\n<p>The degrees of freedom for the model is <b>",
            x$dof, " </b>and the fit was <b>", round(x$criteria["objective"],
                4), "</b>\n</p>"),file=file)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.loadings" <- function (x, digits = 3, cutoff = 0.1, sort = FALSE, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    Lambda <- unclass(x)
    p <- nrow(Lambda)
    factors <- ncol(Lambda)
    if (sort) {
        mx <- max.col(abs(Lambda))
        ind <- cbind(1:p, mx)
        mx[abs(Lambda[ind]) < 0.5] <- factors + 1
        Lambda <- Lambda[order(mx, 1:p), ]
    }
    HTMLli("Loadings:\n<br>",file=file)
    fx <- format(round(Lambda, digits))
    names(fx) <- NULL
    nc <- nchar(fx[1])
    fx[abs(Lambda) < cutoff] <- paste(rep("&nbsp;", nc), collapse = "")
    HTML(fx, file=file, ...)
    vx <- colSums(x^2)
    varex <- rbind("SS loadings" = vx)
    if (is.null(attr(x, "covariance"))) {
        varex <- rbind(varex, "Proportion Var" = vx/p)
        if (factors > 1)
            varex <- rbind(varex, "Cumulative Var" = cumsum(vx/p))
    }
    HTMLbr(file=file)
    HTML(round(varex, digits),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.hclust" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(x$call))
        HTMLli(paste("Call : ", deparse(x$call), "\n<ul>\n", sep = ""),file=file)
    if (!is.null(x$method))
        HTMLli(paste("Cluster method :", x$method, "\n"),file=file)
    if (!is.null(x$dist.method))
        HTMLli(paste("Distance : ", x$dist.method, "\n"),file=file)
    HTMLli(paste("Number of objects: ", length(x$height) + 1, "\n"),file=file)
	HTML("</ul><br>&nbsp;<br>",file=file)
	invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.prcomp" <- function (x, print.x = FALSE, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Standard deviations:\n</p>",file=file,append=TRUE)
    HTML(x$sdev, file=file,append=TRUE,...)
    HTML("\n<p>Rotation:\n</p>")
    HTML(x$rotation, file=file,append=TRUE,...)
    if (print.x && length(x$x)) {
        HTML("\n<p>Rotated variables:\n</p>")
        HTML(x$x, file=file,append=TRUE,...)
    }
    invisible(x)
}



#----------------------------------------------------------------------------------------------------#

"HTML.princomp" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("Call: <font class=call>",deparse(x$call),"</font>"),file=file)
    HTML("\n<p>Standard deviations:\n</p>",file=file)
    HTML(t(as.matrix(x$sdev)), file=file,append=TRUE,...)
    HTML(paste("\n<p><b>", length(x$scale), " </b>variables and <b>", x$n.obs, " </b>observations.\n</p>"),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.prcomp" <- function (x, digits = min(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Importance of components:\n</p>",file=file)
    HTML(x$importance, digits = digits,file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.princomp" <- function (x, digits = 3, loadings = x$print.loadings, cutoff = x$cutoff, file=HTMLGetFile(), append=TRUE, ...)
{
    cat("\n",file=file,append=append,...)
    vars <- x$sdev^2
    vars <- vars/sum(vars)
    HTML("<p>Importance of components:\n</p>",file=file)
    HTML(rbind("Standard deviation" = x$sdev, "Proportion of Variance" = vars,
        "Cumulative Proportion" = cumsum(vars)),file=file)
    if (loadings) {
        HTMLli("Loadings:\n",file=file)
        cx <- format(round(x$loadings, digits = digits))
        cx[abs(x$loadings) < cutoff] <- substring("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
            1, nchar(cx[1, 1]))
        HTML(cx, quote = FALSE, file=file)
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
### PACKAGE EDA (merged into stats)
#----------------------------------------------------------------------------------------------------#

"HTML.medpolish" <- function (x, digits = getOption("digits"), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(paste("\n<p><b>Median Polish Results (Dataset: \"", x$name, "\")\n</b></p>",
        sep = ""),file=file)
    HTML(paste("\n<p>Overall:", x$overall, "\n</p>\n<p>Row Effects:\n</p>"),file=file)
    HTML(x$row, digits = digits, file=file,append=TRUE,...)
    HTML("\n<p>Column Effects:\n</p>",file=file)
    HTML(x$col, digits = digits, file=file)
    HTML("\n<p>Residuals:\n</p>",file=file)
    HTML(x$residuals, digits = max(2, digits - 2), file=file)
    HTMLbr(file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.tukeyline" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("Call:\n", deparse(x$call), "\n<br>\n", sep = ""),file=file)
    HTML("<p>Coefficients:\n</p>",file=file)
    print.default(format(coef(x), digits = digits),file=file)
    HTMLbr(file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.tukeysmooth" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(paste("<p><b>",attr(x, "kind"), " Tukey smoother resulting from ", deparse(attr(x,
        "call")), "\n",   if (twiced <- attr(x, "twiced")) " <-<-twiced<-<- ",
        if (!is.null(it <- attr(x, "iter"))) paste(" used", it, "iterations\n"),
        if (!is.null(ch <- attr(x, "changed"))) paste(if (!ch) " NOT ", " changed\n</b></p>")),file=file)
    if (length(class(x)) > 1)
        NextMethod()
    else {
        y <- x
        attributes(y) <- NULL
        HTML(y,file=file, append=TRUE)
        invisible(x)
    }
}


#----------------------------------------------------------------------------------------------------#
### PACKAGE EDA (merged into stats)
#----------------------------------------------------------------------------------------------------#

#
# 2008-05-23: Removed by Fernando H Rosa. Class appears to no longer exist on package stats
#
#"HTML.grob" <- function (x, file=HTMLGetFile(), append=TRUE,...)
#{
#    cat("\n",file=file,append=append,...)
#    cl <- class(get.value.grob(x))
#    HTML(paste(cl[1:(length(cl) - 1)], collapse = "&nbsp;"),file=file)
#    invisible(x)
#}

#----------------------------------------------------------------------------------------------------#

"HTML.unit" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(as.character(x), file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.viewport" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(class(x),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#
### PACKAGE LATTICE
#----------------------------------------------------------------------------------------------------#

"HTML.shingle" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("\n<p>Data:\n</p>",file=file)
    HTML(as.numeric(x),file=file)
    l <- levels(x)
    n <- nlevels(x)
    if (n < 1)
        HTML("\n<p>no intervals\n</p>",file=file)
    else {
        int <- data.frame(min = numeric(n), max = numeric(n),
            count = numeric(n))
        for (i in 1:n) {
            int$min[i] <- l[[i]][1]
            int$max[i] <- l[[i]][2]
            int$count[i] <- length(x[x >= l[[i]][1] & x <= l[[i]][2]])
        }
        HTML("\n<p>Intervals:\n</p>",file=file)
        HTML(int,file=file)
        olap <- numeric(n - 1)
        if (n > 2)
            for (i in 1:(n - 1)) olap[i] <- length(x[x >= l[[i]][1] &
                x <= l[[i]][2] & x >= l[[i + 1]][1] & x <= l[[i +
                1]][2]])
        HTML("\n<p>Overlap between adjacent intervals:\n</p>",file=file)
        HTML(olap,file=file)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

 "HTML.shingleLevel" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(do.call("rbind", x),file=file)
    invisible(x)
}



#----------------------------------------------------------------------------------------------------#
### PACKAGE MASS
#----------------------------------------------------------------------------------------------------#

"HTML.abbrev" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (is.list(x))
        x <- unlist(x)
    NextMethod("HTML")
}


#----------------------------------------------------------------------------------------------------#

"HTML.Anova" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    heading <- attr(x, "heading")
    if (!is.null(heading))
        HTML(paste("<p>",heading,"</p>", sep = " ",collapse="<br>"),file=file)
    attr(x, "heading") <- NULL
    HTML.data.frame(x,file=file)
}

#----------------------------------------------------------------------------------------------------#

"HTML.anova.loglm" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    y <- x
    y[, 5] <- round(y[, 5], 5)
    R <- array("", dim(x), dimnames(x))
    for (j in 1:5) {
        colj <- c(colnames(x)[j], format(y[, j]))
        R[, j] <- colj[-1]
        colnames(R)[j] <- colj[1]
    }
    R[1, 3:5] <- ""
    forms <- attr(x, "formulae")
    HTML("<p><b>LR tests for hierarchical log-linear models</b>\n</p>\n",file=file)
    for (i in seq(along = forms))
    HTML(paste(paste("<p>Model ", i, ":<br>", sep = ""), paste(deparse(forms[[i]]),collapse=""), "</p>"),file=file)
    HTMLbr(file=file)
    HTML(R,file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.correspondence" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(paste("<p>First canonical correlation(s):", format(x$cor, ...), "\n</p>"),file=file)
    rcn <- names(dimnames(x$Freq))
    HTML(paste("\n<p>", rcn[1], "scores:\n</p>"),file=file)
    HTML(x$rscore,file=file)
    HTML(paste("\n<p>", rcn[2], "scores:\n</p>"),file=file)
    HTML(x$cscore,file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.fitdistr" <- function (x, digits = getOption("digits"), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    ans <- format(rbind(x$estimate, x$sd), digits = digits)
    ans[1, ] <- sapply(ans[1, ], function(x) paste("", x))
    ans[2, ] <- sapply(ans[2, ], function(x) paste("(", x, ")",
        sep = ""))
    dn <- dimnames(ans)
    dn[[1]] <- rep("", 2)
    dn[[2]] <- paste(substring("  ", 1, (nchar(ans[2, ]) -
        nchar(dn[[2]]))%/%2), dn[[2]])
    dn[[2]] <- paste(dn[[2]], substring("  ", 1, (nchar(ans[2,
        ]) - nchar(dn[[2]]))%/%2))
    dimnames(ans) <- dn
    HTML(ans, file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.fractions" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    y <- attr(x, "fracs")
    att <- attributes(x)
    att$fracs <- att$class <- NULL
    x <- do.call("structure", c(list(y), att))
    NextMethod("HTML", file=file)
}


#----------------------------------------------------------------------------------------------------#

"HTML.gamma.shape" <- function (x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    y <- x
    x <- array(unlist(x), dim = 2:1, dimnames = list(c("Alpha ", "SE "), ""))
    NextMethod("HTML",file=file)
    invisible(y)
}

#----------------------------------------------------------------------------------------------------#

"HTML.glm.dose" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    M <- cbind(x, attr(x, "SE"))
    dimnames(M) <- list(names(x), c("Dose", "SE"))
    x <- M
    NextMethod("HTML",file=file)
}

#----------------------------------------------------------------------------------------------------#

"HTML.lda" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) {
        names(cl)[2] <- ""
        HTMLli(paste("Call: ",deparse(cl)),file=file)
    }
    HTML("\n<p>Prior probabilities of groups:\n</p>",file=file)
    HTML(x$prior, file=file,...)
    HTML("\n<p>Group means:\n</p>",file=file)
    HTML(x$means, file=file,...)
    HTML("\n<p>Coefficients of linear discriminants:\n</p>",file=file)
    HTML(x$scaling, file=file,...)
    svd <- x$svd
    names(svd) <- dimnames(x$scaling)[[2]]
    if (length(svd) > 1) {
        HTML("\n<p>Proportion of trace:\n</p>",file=file)
        HTML(round(svd^2/sum(svd^2), 4), file=file,append=TRUE,...)
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.loglm" <- function (x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("Call: <font class=call>",deparse(x$call),"</font>"),file=file)
    ts.array <- rbind(c(x$lrt, x$df, if (x$df > 0) 1 - pchisq(x$lrt,
        x$df) else 1), c(x$pearson, x$df, if (x$df > 0) 1 - pchisq(x$pearson,
        x$df) else 1))
    dimnames(ts.array) <- list(c("Likelihood Ratio", "Pearson"),
        c("X^2", "df", "P(> X^2)"))
    HTML("\n<p>Statistics:\n</p>",file=file)
    HTML(ts.array,file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.mca" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) HTMLli(paste("Call: ",deparse(cl)),file=file)

    HTML(paste("\n<p>Multiple correspondence analysis of <b>", nrow(x$rs),
        " </b>cases of <b> ", x$p, " </b>factors\n</p>"),file=file)

    p <- 100 * cumsum(x$d)/(x$p - 1)
    HTML(paste("\n<p>Correlations ",paste(round(x$d, 3),collapse=" "),"  cumulative % explained ", paste(round(p, 2),collapse=" "),"</p>" ),file=file)

    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.polr" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) HTMLli(paste("Call: ",deparse(cl)),file=file)
    if (length(coef(x))) {
        HTML("\n<p>Coefficients:\n</p>",file=file)
        HTML(coef(x), file=file,append=TRUE,...)
    }
    else {
        HTML("\n<p>No coefficients\n</p>",file=file)
    }
    HTML("\n<p>Intercepts:\n</p>",file=file)
    HTML(x$zeta, file=file,append=TRUE,...)
    HTML(paste("\n<p>Residual Deviance: <b>", format(x$deviance, nsmall = 2), "</b>\n</p>"),file=file)
    HTML(paste("<p>AIC:<b>", format(x$deviance + 2 * x$edf, nsmall = 2), "</b>\n</p>"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.qda" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) {
        names(cl)[2] <- ""
        HTMLli(paste("Call: ",deparse(cl)),file=file)
    }
    HTML("\n<p>Prior probabilities of groups:\n</p>",file=file)
    HTML(x$prior, file=file,...)
    HTML("\n<p>Group means:\n</p>",file=file)
    HTML(x$means, file=file,append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.ridgelm" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    scaledcoef <- t(as.matrix(x$coef/x$scales))
    if (x$Inter) {
        inter <- x$ym - scaledcoef %*% x$xm
        scaledcoef <- cbind(Intercept = inter, scaledcoef)
    }
    HTML(drop(scaledcoef), file=file,append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.rlm" <- function (x,file=HTMLGetFile(), append=TRUE, ...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) {
        HTMLli(paste("Call: ",paste(deparse(cl),collapse=" ")),file=file)
    }
    if (x$converged)
        HTML(paste("<p>Converged in <b>", length(x$conv), "</b> iterations\n</p>"),file=file)
    else HTML(paste("<p>Ran <b>", length(x$conv), " </b>iterations without convergence\n</p>"),file=file)
    coef <- x$coef
    HTML("\n<p>Coefficients:\n</p>",file=file)
    HTML(coef, file=file,append=TRUE,...)
    nobs <- length(x$resid)
    rdf <- nobs - length(coef)
    HTML(paste("\n<p>Degrees of freedom: <b>", nobs, " </b>total; <b>", rdf, " </b>residual\n</p>"),file=file)
    HTML(paste("<p>Scale estimate:<b>", paste(format(signif(x$s, 3)),collapse=" "), "</b>\n</p>"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.rms.curv" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(paste("<p><li>Parameter effects: c^theta x sqrt(FALSE) =<b>", round(x$pe,
        4), "</b>\n<br><li>", "Intrinsic: c^iota  x sqrt(FALSE) =<b>", round(x$ic,
        4), "\n</b></p>"),file=file, append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.loglm" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Formula:\n</p>",file=file)
    HTML(formula(x),file=file)
    HTML("\n<p>Statistics:\n</p>",file=file)
    HTML(x$tests,file=file)
    if (!is.null(x$oe)) {
        HTML("\n<p>Observed (Expected):\n</p>",file=file)
        HTML(x$oe, file=file)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.negbin" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
	    cat("\n",file=file,append=append,...)
	NextMethod(x,file=file)
	dp <- 2 - floor(log10(x$SE.theta))
    	HTML(paste("<p><li>Theta:<b> ", round(x$theta, dp), "</b>\n<li>Std. Err.:<b> ", round(x$SE.theta,  dp), "</b>\n</p>"),file=file)
    	if (!is.null(x$th.warn))
    	HTML(paste("<p>Warning while fitting theta:", x$th.warn, "\n</p>"),file=file)
	HTML(paste("\n<p><li> 2 x log-likelihood: ", format(round(x$twologlik, 3), nsmall = dp), "\n</p>"),file=file)
	invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.polr" <- function (x, digits = x$digits, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) {
        HTMLli(paste("Call: ",deparse(cl)),file=file)
    }
    coef <- format(round(x$coef, digits = digits))
    pc <- x$pc
    if (pc > 0) {
        HTML("\n<p>Coefficients:\n</p>",file=file)
        HTML(x$coef[seq(len = pc), ], file=file,append=TRUE, ...)
    }
    else {
        HTML("\n<p>No coefficients\n</p>",file=file)
    }
    HTML("\n<p>Intercepts:\n</p>",file=file)
    HTML(coef[(pc + 1):nrow(coef), ], file=file,append=TRUE, ...)
    HTML(paste("\n<p>Residual Deviance:<b>", format(x$deviance, nsmall = 2), "</b>\n</p>"),file=file)
    HTML(paste("\n<p>AIC:<b>", format(x$deviance + 2 * x$edf, nsmall = 2), "</b>\n</p>"),file=file)
    if (!is.null(correl <- x$correlation)) {
        cat("\n<p>Correlation of Coefficients:\n</p>",file=file)
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        HTML(correl[-1, -ncol(correl)], file=file, append=TRUE,...)
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.rlm" <- function (x, digits = max(3, .Options$digits - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTMLli(paste("\nCall: ",deparse(x$call)),file=file)
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
    if (rdf > 5) {
        HTML("<p>Residuals:\n</p>",file=file)
        if (length(dim(resid)) == 2) {
            rq <- apply(t(resid), 1, quantile)
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q",
                "Max"), colnames(resid))
        }
        else {
            rq <- quantile(resid)
            names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        }
        HTML(rq, file=file)
    }
    else if (rdf > 0) {
        HTML("<p>Residuals:\n</p>",file=file)
        HTML(resid,file=file)
    }
    if (nsingular <- df[3] - df[1])
        HTML(paste("\n<p>Coefficients: (", nsingular, " not defined because of singularities)\n</p>",sep = ""),file=file)
    else HTML("\n<p>Coefficients:\n</p>",file=file)
    HTML(format(round(x$coef, digits = digits)), file=file)
    HTML(paste("\n<p>Residual standard error:<b>", format(signif(x$sigma,
        digits)), " </b>on <b> ", rdf, " </b>degrees of freedom\n</p>"),file=file)
    if (!is.null(correl <- x$correlation)) {
        p <- dim(correl)[2]
        if (p > 1) {
            HTML("\n<p>Correlation of Coefficients:\n</p>",file=file)
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            HTML(correl[-1, -p, drop = FALSE], file=file)
        }
    }
    invisible(x)
}



#----------------------------------------------------------------------------------------------------#
### PACKAGE NNET
#----------------------------------------------------------------------------------------------------#

"HTML.multinom" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) {
        HTMLli(paste("Call: ",paste(deparse(cl),collapse="")),file=file)
    }
    HTML("\n<p>Coefficients:\n</p>",file=file)
    HTML(coef(x), file=file)
    HTML(paste("\n<p>Residual Deviance: <b>", format(x$deviance), "</b>\n</p>"),file=file)
    HTML(paste("<p>AIC:<b>", format(x$AIC), "</b>\n</p>"),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.nnet" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!inherits(x, "nnet"))
        stop("Not legitimate a neural net fit")
    HTML(paste("<p><b>a ", x$n[1], "-", x$n[2], "-", x$n[3], " network with ", length(x$wts), " weights.</b></p>", sep = ""),file=file)

    if (length(x$coefnames))
        HTML(paste("<p>inputs:", x$coefnames, "\noutput(s):", deparse(formula(x)[[2]]), "\n</p>"),file=file)
    HTML("<p>options were -</p>",file=file)
    tconn <- diff(x$nconn)
    if (tconn[length(tconn)] > x$n[2] + 1)
        HTMLli(" skip-layer connections ",file=file)
    if (x$nunits > x$nsunits && !x$softmax)
        HTMLli(" linear output units ",file=file)
    if (x$entropy)
        HTMLli(" entropy fitting ",file=file)
    if (x$softmax)
        HTMLli(" softmax modelling ",file=file)
    if (x$decay[1] > 0)
        HTMLli(paste(" decay=", x$decay[1], sep = ""),file=file)
    HTMLbr(file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.multinom" <- function (x, digits = x$digits, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    if (!is.null(cl <- x$call)) {
        HTMLli(paste("Call:",paste(deparse(cl),collapse=" ")),file=file)
    }
    HTML("\n<p>Coefficients:\n</p>",file=file)
    if (x$is.binomial) {
        HTML(cbind(Values = x$coefficients, "Std. Err." = x$standard.errors,
            "Value/SE" = x$Wald.ratios), file=file)
    }
    else {
        HTML(x$coefficients, file=file)
        HTML("\n<p>Std. Errors:\n</p>",file=file)
        HTML(x$standard.errors, file=file)
        if (!is.null(x$Wald.ratios)) {
            HTML("\n<O>Value/SE (Wald statistics):\n</p>",file=file)
            HTML(x$coefficients/x$standard.errors, file=file)
        }
    }
    HTML(paste("\n<p>Residual Deviance:<b>", format(x$deviance), "</b>\n</p>"),file=file)
    HTML(paste("\n<p>AIC:<b>", format(x$AIC), "</b>\n</p>"),file=file)
    if (!is.null(correl <- x$correlation)) {
        p <- dim(correl)[2]
        if (p > 1) {
            HTML("\n</p>Correlation of Coefficients:\n</p>",file=file)
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            HTML(correl[-1, -p], file= file)
        }
    }
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.nnet" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
     cat("\n",file=file,append=append,...)
     HTML(paste("<p><b>a ", x$n[1], "-", x$n[2], "-", x$n[3], " network with ", length(x$wts), " weights.</b></p>", sep = ""),file=file)

        HTML("<p>options were -</p>",file=file)
        tconn <- diff(x$nconn)
        if (tconn[length(tconn)] > x$n[2] + 1)
            HTMLli(" skip-layer connections ",file=file)
        if (x$nunits > x$nsunits && !x$softmax)
            HTMLli(" linear output units ",file=file)
        if (x$entropy)
            HTMLli(" entropy fitting ",file=file)
        if (x$softmax)
            HTMLli(" softmax modelling ",file=file)
        if (x$decay[1] > 0)
        HTMLli(paste(" decay=", x$decay[1], sep = ""),file=file)
    wts <- format(round(nnet::nnet(x), 2))
    lapply(split(wts, rep(1:x$nunits, tconn)), function(x) HTML(x,file=file))
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
### PACKAGE CLUSTER
#----------------------------------------------------------------------------------------------------#


"HTML.agnes" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Merge:\n</p>",file=file)
    HTML(x$merge, file=file,append=TRUE,...)
    HTML("<p>Order of objects:\n</p>",file=file)
    HTML(if (length(x$order.lab) != 0)
        x$order.lab
    else x$order, file=file,append=TRUE, ...)
    HTML("<p>Height:\n</p>",file=file)
    HTML(x$height, file=file,append=TRUE,...)
    HTML("<p>Agglomerative coefficient:\n</p>",file=file)
    HTML(x$ac, file=file,append=TRUE,...)
    HTML("\n<p>Available components:\n</p>",file=file)
    HTML(names(x), file=file,append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.clara" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Best sample:\n</p>",file=file)
    HTML(x$sample, file=file, append=TRUE,...)
    HTML("<p>Medoids:\n</p>",file=file)
    HTML(x$medoids, file=file,append=TRUE,...)
    HTML("<p>Clustering vector:\n</p>",file=file)
    HTML(x$clustering, file=file,append=TRUE,...)
    HTML("<p>Objective function:\n</p>",file=file)
    HTML(x$objective, file=file,append=TRUE,...)
    HTML("\n<p>Available components:\n</p>",file=file)
    HTML(names(x),file=file, append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.diana" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Merge:\n</p>",file=file)
    HTML(x$merge, file=file,append=TRUE,...)
    HTML("<p>Order of objects:\n</p>",file=file)
    HTML(if (length(x$order.lab) != 0)  x$order.lab    else x$order, file= file, append=TRUE,...)
    HTML("<p>Height:\n</p>",file=file)
    HTML(x$height, file=file,append=TRUE,...)
    HTML("<p>Divisive coefficient:\n</p>",file=file)
    HTML(x$dc,file=file, append=TRUE,...)
    HTML("\n<p>Available components:\n</p>",file=file)
    HTML(names(x),file=file,append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.dissimilarity" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Dissimilarities :\n</p>",file=file)
    HTML(as.vector(x),file=file,append=TRUE, ...)
    if (!is.null(attr(x, "na.message")))
        HTML(paste("<p>Warning : ", attr(x, "NA.message"), "\n</p>"),file=file)
    HTML(paste("<p>Metric : ", attr(x, "Metric"), "\n</p>"),file=file)
    HTML(paste("<p>Number of objects : ", attr(x, "Size"), "\n</p>"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.ellipsoid" <- function (x, digits = max(1, getOption("digits") - 2), file=HTMLGetFile(), append=TRUE,...)
{

    cat("\n",file=file,append=append,...)
    d <- length(x$loc)
    HTML(paste("<p>`ellipsoid' in <b> ", d, " </b>dimensions:<br> center = (<b>", paste(format(x$loc,
        digits = digits),collapse=" "), "</b>); squared ave.radius d^2 = <b>", format(x$d2,
        digits = digits), " </b>\n<br> and shape matrix =\n</p>"),file=file)
    HTML(x$cov, file=file, append=TRUE,...)
    HTML(paste("<p>&nbsp;  hence,", if (d == 2)
        " area "
    else " volume ", " = <b>", format(cluster::volume(x), digits = digits),
        "\n</b></p>"),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.fanny" <- function (x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(x$objective, file=file,append=TRUE,...)
    HTML("<p>Membership coefficients:\n</p>", file=file)
    HTML(x$membership, file=file,append=TRUE, ...)
    HTML("<p>Coefficients:\n</p>", file=file)
    HTML(x$coeff, file=file, append=TRUE,...)
    HTML("<p>Closest hard clustering:\n</p>", file=file)
    HTML(x$clustering, file=file,append=TRUE, ...)
    HTML("\n<p>Available components:\n</p>", file=file)
    HTML(names(x), file=file, append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.mona" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Revised data:\n</p>",file=file)
    HTML(x$data,file=file,  append=TRUE,...)
    HTML("<p>Order of objects:\n</p>",file=file)
    HTML(if (length(x$order.lab) != 0)  x$order.lab else x$order,file=file, append=TRUE,...)
    HTML("<p>Variable used:\n</p>",file=file)
    HTML(x$variable, file=file, append=TRUE,...)
    HTML("<p>Separation step:\n</p>",file=file)
    HTML(x$step,file=file, append=TRUE,...)
    HTML("\n<p>Available components:\n</p>",file=file)
    HTML(names(x),file=file,append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.pam" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Medoids:\n</p>",file=file)
    HTML(x$medoids,file=file, append=TRUE,...)
    HTML("<p>Clustering vector:\n</p>",file=file)
    HTML(x$clustering,file=file, append=TRUE,...)
    HTML("<p>Objective function:\n</p>",file=file)
    HTML(x$objective,file=file, append=TRUE,...)
    HTML("\n<p>Available components:\n</p>",file=file)
    HTML(names(x),file=file, append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.agnes" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Merge:\n</p>",file=file)
    HTML(x$merge, file=file, append=TRUE,...)
    HTML("<p>Order of objects:\n</p>",file=file)
    HTML(if (length(x$order.lab) != 0)
        x$order.lab
    else x$order, file=file, append=TRUE,...)
    HTML("<p>Height:\n</p>",file=file)
    HTML(x$height, file=file,append=TRUE, ...)
    HTML("<p>Agglomerative coefficient:\n</p>",file=file)
    HTML(x$ac, file=file, append=TRUE,...)
    HTML("<p>\n</p>",file=file)
    HTML(x$diss, file=file, append=TRUE,...)
    HTML("<p>\nAvailable components:\n</p>",file=file)
    HTML(names(x), file=file,append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.clara" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Best sample:\n</p>",file=file)
    HTML(x$sample, file=file, append=TRUE,...)
    HTML("<p>Medoids:\n</p>",file=file)
    HTML(x$medoids, file=file, append=TRUE,...)
    HTML("<p>Clustering vector:\n</p>",file=file)
    HTML(x$clustering, file=file,append=TRUE, ...)
    HTML("<p>Objective function:\n</p>",file=file)
    HTML(x$objective, file=file,append=TRUE, ...)
    HTML("<p>\nNumerical information per cluster:\n</p>",file=file)
    HTML(x$clusinfo, file=file, append=TRUE,...)
    if (length(x$silinfo) != 0) {
        HTML("<p>\nSilhouette plot information for best sample:\n</p>",file=file)
        HTML(x$silinfo[[1]], file=file,append=TRUE, ...)
        HTML("<p>Average silhouette width per cluster:\n</p>",file=file)
        HTML(x$silinfo[[2]], file=file,append=TRUE, ...)
        HTML("<p>Average silhouette width of best sample:\n</p>",file=file)
        HTML(x$silinfo[[3]], file=file,append=TRUE, ...)
    }
    HTML("<p>\n</p>",file=file)
    HTML(x$diss, file=file, append=TRUE,...)
    HTML("<p>\nAvailable components:\n</p>",file=file)
    HTML(names(x), file=file,append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.diana" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML("<p>Merge:\n</p>",file=file)
    HTML(x$merge, file=file, append=TRUE,...)
    HTML("<p>Order of objects:\n</p>",file=file)
    HTML(if (length(x$order.lab) != 0)
        x$order.lab
    else x$order, file=file, append=TRUE,...)
    HTML("<p>Height:\n</p>",file=file)
    HTML(x$height, file=file,append=TRUE, ...)
    HTML("<p>Divisive coefficient:\n</p>",file=file)
    HTML(x$dc, file=file, append=TRUE,...)
    HTML("<p>\n</p>",file=file)
    HTML(x$diss, file=file,append=TRUE, ...)
    HTML("<p>\nAvailable components:\n</p>",file=file)
    HTML(names(x), file=file,append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

 "HTML.summary.fanny" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n",file=file,append=append,...)
    HTML(x$objective, file=file, append=TRUE,...)
    HTML("<p>Membership coefficients:\n</p>",file=file)
    HTML(x$membership, file=file, append=TRUE, ...)
    HTML("<p>Coefficients:\n</p>",file=file)
    HTML(x$coeff, file=file, append=TRUE, ...)
    HTML("<p>Closest hard clustering:\n</p>",file=file)
    HTML(x$clustering, file=file, append=TRUE, ...)
    if (length(x$silinfo) != 0) {
        HTML("<p>\nSilhouette plot information:\n</p>",file=file)
        HTML(x$silinfo[[1]], file=file, append=TRUE, ...)
        HTML("<p>Average silhouette width per cluster:\n</p>",file=file)
        HTML(x$silinfo[[2]], file=file, append=TRUE, ...)
        HTML("<p>Average silhouette width of total data set:\n</p>",file=file)
        HTML(x$silinfo[[3]], file=file, append=TRUE, ...)
    }
    HTML("<p>\n</p>",file=file)
    HTML(x$diss, file=file, append=TRUE, ...)
    HTML("<p>\nAvailable components:\n</p>",file=file)
    HTML(names(x), file=file, append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.mona" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML.mona(x, file=file, append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.pam" <- function(x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file,append=append,...)
    HTML("<p>Medoids:\n</p>",file=file)
    HTML(x$medoids, file=file, append=TRUE, ...)
    HTML("<p>Clustering vector:\n</p>",file=file)
    HTML(x$clustering, file=file, append=TRUE, ...)
    HTML("<p>Objective function:\n</p>",file=file)
    HTML(x$objective, file=file, append=TRUE, ...)
    HTML("<p>\nNumerical information per cluster:\n</p>",file=file)
    HTML(x$clusinfo, file=file, append=TRUE, ...)
    HTML("<p>\nIsolated clusters:\n</p>",file=file)
    HTML("<p>L-clusters: ")
    HTML(names(x$isolation[x$isolation == "L"]),
        ...)
    HTML("<p>L*-clusters: ")
    HTML(names(x$isolation[x$isolation == "L*"]),
        ...)
    if (length(x$silinfo) != 0) {
        HTML("<p>\nSilhouette plot information:\n</p>",file=file)
        HTML(x$silinfo[[1]], file=file, append=TRUE, ...)
        HTML("<p>Average silhouette width per cluster:\n</p>",file=file)
        HTML(x$silinfo[[2]], file=file, append=TRUE, ...)
        HTML("<p>Average silhouette width of total data set:\n</p>",file=file)
        HTML(x$silinfo[[3]], file=file, append=TRUE, ...)
    }
    HTML("<p>\n</p>",file=file)
    HTML(x$diss, file=file, append=TRUE, ...)
    HTML("<p>\nAvailable components:\n</p>",file=file)
    HTML(names(x), file=file, append=TRUE, ...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
### PACKAGE MGCV
#----------------------------------------------------------------------------------------------------#

"HTML.gam" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML(x$family,file=file)
    HTML("<p>Formula:\n</p>",file=file)
    HTML(x$formula,file=file)
    if (x$dim == 0)
        HTML(paste("<p>Total model degrees of freedom <b>", x$nsdf, " </b>\n</p>"),file=file)
    else HTML(paste("\n<p>Estimated degrees of freedom:<b>", paste(x$edf,collapse=" , "), "</b>  total = <b>",
        paste(sum(x$edf) + x$nsdf,collapse=" , "), "</b>\n</p>"),file=file)
    gcv <- x$df.null * x$sig2/(x$df.null - sum(x$edf) - x$nsdf)
    HTML("\n<p>GCV score:</p> ",file=file)
    HTML(gcv,file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.gam" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML(x$family,file=file)
    HTML("<p>Formula:\n</p>",file=file)
    HTML(x$formula,file=file)
    if (length(x$p.coeff) > 0) {
        HTML("\n<p>Parametric coefficients:\n</p>",file=file)
        width <- max(nchar(names(x$p.coeff)))

        HTML("\n<p align=center><table cellspacing=0 border=1><td><table cellspacing=0> <tr class= firstline >        <th></th><th>Estimate</th><th>std.err.</th><th>t ratio</th><th>Pr(>|t[)</th></tr>\n",file=file)


        for (i in 1:length(x$p.coeff)) HTML(paste("<tr><td class=firstcolumn>",formatC(names(x$p.coeff)[i], width = width),"</td><td class=\"CellInside\">", formatC(x$p.coeff[i], width = 10,digits = 5),"</td><td class=\"CellInside\">", formatC(x$se[i], width = 10, digits = 4),"</td><td class=\"CellInside\">",formatC(x$p.t[i], width = 10, digits = 4), "</td><td class=\"CellInside\">",format.pval(x$p.pv[i]),"</td></tr>\n", sep = ""),file=file)

           HTML("\n</table></td></table></center>",file=file)

    }
	HTMLbr( file=file)
    if (x$m > 0) {
        HTML("<p>Approximate significance of smooth terms:\n</p>",file=file)
        width <- max(nchar(names(x$chi.sq)))

        HTML("\n<p align=center><table cellspacing=0 border=1><td><table cellspacing=0> <tr class= firstline > <th></th><th>edf</th><th>chi.sq</th><th>p-value</th></tr>\n",file=file)
        for (i in 1:x$m)

        HTML(paste("<tr><td class=firstcolumn>",formatC(names(x$chi.sq)[i], width = width),
        "</td><td class=CellInside>", formatC(x$edf[i], width = 10, digits = 4), "</td><td class=CellInside>",
            formatC(x$chi.sq[i], width = 10, digits = 5),"</td><td class=CellInside>",
            format.pval(x$s.pv[i]), "</td></tr>\n", sep = ""),file=file)

           HTML("\n</table></td></table></center>",file=file)

    }
    HTML(paste("\n<p>Adjusted r-sq. = <b>", formatC(x$r.sq, digits = 3, width = 5),
        " </b>   GCV score = <b>", formatC(x$gcv, digits = 5), "</b> \n<br>Scale estimate = <b>",
        formatC(x$scale, digits = 5, width = 8, flag = "-"),
        "    </b>     n = <b>", x$n, "</b>\n</p>", sep = ""),file=file)
        invisible(x)
}


#----------------------------------------------------------------------------------------------------#
### PACKAGE RPART
#----------------------------------------------------------------------------------------------------#


"HTML.rpart" <- function (x, minlength = 0, spaces = 2, cp, digits = getOption("digits"),
    file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    if (!inherits(x, "rpart"))
        stop("Not legitimate rpart object")
    if (!missing(cp))
        x <- rpart::prune.rpart(x, cp = cp)
    frame <- x$frame
    ylevel <- attr(x, "ylevels")
    node <- as.numeric(row.names(frame))
    # tree.depth is not exported by rpart anymore. Defining it locally:
    "Inttree.depth" <-
    function (nodes)
    {
        depth <- floor(log(nodes, base = 2) + 1e-7)
        as.vector(depth - min(depth))
    }
    depth <- Inttree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = " ")
    if (length(node) > 1) {
        indent <- substring(indent, 1, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")",
            sep = "")
    }
    else indent <- paste(format(node), ")", sep = "")
    tfun <- (x$functions)$print
    if (!is.null(tfun)) {
        if (is.null(frame$yval2))
            yval <- tfun(frame$yval, ylevel, digits)
        else yval <- tfun(frame$yval2, ylevel, digits)
    }
    else yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels(x, digits = digits, minlength = minlength, ...)
    n <- frame$n
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
        yval, term)
    omit <- x$na.action
    if (length(omit))
        HTML(paste("<p>n=<b>", n[1], "</b> (", naprint(omit), ")\n</p>\n", sep = ""),file=file)
    else HTML(paste("<p>n=<b>", n[1], "</b>\n</p>\n"),file=file)
    if (x$method == "class")
        HTML("<p>node), split, n, loss, yval, (yprob)\n</p>",file=file)
    else HTML("<p>node), split, n, deviance, yval\n</p>",file=file)
    HTML("<p>      * denotes terminal node\n\n</p>",file=file)
    HTML(paste("<xmp>", paste(z, sep = "\n",collapse="\n"),"</xmp>"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
### PACKAGE MODREG
#----------------------------------------------------------------------------------------------------#

"HTML.loess" <- function (x, digits = max(3, getOption("digits") - 3),file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    if (!is.null(cl <- x$call)) HTMLli(paste("Call: ",paste(deparse(cl),collapse=" ")),file=file)
    HTML(paste("\n<p>Number of Observations:<b>", x$n, "</b>\n</p>"),file=file)
    HTML(paste("<p>Equivalent Number of Parameters:<b>", format(round(x$enp,
        2)), "</b>\n</p>"),file=file)
    HTML(paste("<p>Residual", if (x$pars$family == "gaussian")
        " Standard Error: <b>"
    else " Scale Estimate:<b> ", format(signif(x$s, digits)), "</b>\n</p>"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.ppr" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    if (!is.null(cl <- x$call)) HTMLli(paste("Call:",paste(deparse(cl),collapse=" ")),file=file)
    mu <- x$mu
    ml <- x$ml
    HTML("\n<p>Goodness of fit:\n</p>",file=file)
    gof <- x$gofn
    names(gof) <- paste(1:ml, "terms")
    HTML(format(gof[mu:ml], ...), file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.smooth.spline" <- function (x, digits = getOption("digits"), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    if (!is.null(cl <- x$call)) HTMLli(paste("Call:",paste(deparse(cl),collapse=" ")),file=file)
    ip <- x$iparms
    cv <- cl$cv
    if (is.null(cv))
        cv <- FALSE
    else if (is.name(cv))
        cv <- eval(cv)
    HTML(paste("\n<p>Smoothing Parameter  spar=<b>", format(x$spar, digits = digits),
        "</b> lambda=<b>", format(x$lambda, digits = digits),"</b>", if (ip["ispar"] !=
            1) paste("(", ip["iter"], " iterations)", sep = ""), "\n</p>"),file=file)
    HTML(paste("<p>Equivalent Degrees of Freedom (Df):<b>", format(x$df, digits = digits),
        "</b>\n</p>"),file=file)
    HTML(paste("<p>Penalized Criterion:<b>", format(x$pen.crit, digits = digits),
        "</b>\n</p>"),file=file)
    HTML(paste ("<p>",if (cv) "PRESS:"
    else "GCV:", "<b>",format(x$cv.crit, digits = digits), "</b>\n</p>"),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.summary.loess" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
   if (!is.null(cl <- x$call)) HTMLli(paste("Call:",paste(deparse(cl),collapse=" ")),file=file)
	HTML(paste("\n<p>Number of Observations:<b>", x$n, "</b>\n</p>"),file=file)
    	HTML(paste("<p>Equivalent Number of Parameters:<b>", format(round(x$enp, 2)), "</b>\n</p>"),file=file)
    if (x$pars$family == "gaussian")
        HTML("<p>Residual Standard Error:</p>",file=file)
    else HTML("<p>Residual Scale Estimate:</p>",file=file)
        HTML(format(signif(x$s, digits)),file=file)
    HTML("<p>Trace of smoother matrix:</p>",file=file)
    HTML(format(round(x$trace.hat,2)), file=file)
    HTML("\n<p>Control settings:\n</p><ul>",file=file)
    HTMLli(paste("normalize: ", x$pars$normalize, "\n"),file=file)
    HTMLli(paste("  span     : ", format(x$pars$span), "\n"),file=file)
    HTMLli(paste("  degree   : ", x$pars$degree, "\n"),file=file)
    HTMLli(paste("  family   : ", x$pars$family),file=file)
    if (x$pars$family != "gaussian")
        HTMLli(paste("       iterations =", x$pars$iterations),file=file)
    	HTML("</ul>",file=file)
    HTML(paste("\n<p>  surface  : ", x$pars$surface, if (x$pars$surface == "interpolate")  paste("  cell =", format(x$pars$cell)),"</p>"),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.ppr" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML.ppr(x,file=file, ...)
    mu <- x$mu
    HTML("\n<p>Projection direction vectors:\n</p>",file=file)
    HTML(format(x$alpha, ...),file=file)
    HTML("\n<p>Coefficients of ridge terms:\n</p>",file=file)
    HTML(format(x$beta, ...), file=file)
    if (any(x$edf > 0)) {
        HTML("\n<p>Equivalent df for ridge terms:\n</p>")
        edf <- x$edf
        names(edf) <- paste("term", 1:mu)
        HTML(round(edf, 2),file=file, append=TRUE,...)
    }
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#
### PACKAGE SPLINES
#----------------------------------------------------------------------------------------------------#



"HTML.bSpline" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    value <- c(rep(NA, splines::splineOrder(x)), coef(x))
    names(value) <- format(splines::splineKnots(x), digits = 5)
    HTML(paste("<p> bSpline representation of spline", if (!is.null(form <- attr(x, "formula"))) paste (" for", paste(deparse(as.vector(form)),collapse=" "))  ,"</p>"),file=file)
    HTML(value, file=file,append=TRUE,...)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.polySpline" <- function (x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    coeff <- coef(x)
    dnames <- dimnames(coeff)
    if (is.null(dnames[[2]]))
        dimnames(coeff) <- list(format(splines::splineKnots(x)), c("constant",
            "linear", "quadratic", "cubic", paste(4:29, "th",
                sep = ""))[1:(dim(coeff)[2])])
    HTML(    paste(    "<p>Polynomial representation of spline ",    if (!is.null(form <- attr(x, "formula")))     	paste(" for ", paste(deparse(as.vector(form)),collapse=" ")  )    ,"</p>")    ,file=file    )
    HTML(coeff, file=file,append=TRUE,...)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#

"HTML.ppolySpline" <- function (x,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML("<p>periodic </p>",file=file)
    HTML(paste("\n<p>Period:<b>", format(x[["period"]]), "</b>\n</p>"),file=file)
    NextMethod("HTML",file=file)
    invisible(x)
}



#----------------------------------------------------------------------------------------------------#
### PACKAGE LSQ
#----------------------------------------------------------------------------------------------------#

"HTML.lqs" <- function (x, digits = max(3, getOption("digits") - 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
	if (!is.null(cl <- x$call)) HTMLli(paste("Call:",paste(deparse(cl),collapse=" ")),file=file)

	HTML("<p>Coefficients:\n</p>",file=file)
    HTML.default(format(coef(x), digits = digits), file=file)
    HTML(paste("\n<p>Scale estimates ", paste(format(x$scale, digits = digits),collapse=" "),
        "\n\n</p>"),file=file)
       invisible(x)
}


#----------------------------------------------------------------------------------------------------#
### PACKAGE NLS
#----------------------------------------------------------------------------------------------------#

"HTML.nls" <- function (x, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML("<p><b>Nonlinear regression model\n</b></p>",file=file)
    HTMLli(paste("Model: ", paste(deparse(formula(x)),collapse=" "), "\n"),file=file)
    HTMLli(paste("Data: ", as.character(x$data), "\n"),file=file)
    HTML(x$m$getAllPars(),file=file)
    HTMLli(paste("Residual sum-of-squares: ", format(x$m$deviance()),"\n"),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

"HTML.summary.nls" <- function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = p >
    4, signif.stars = getOption("show.signif.stars"), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    HTML(paste("<p>Formula:",paste(deparse(x$formula), collapse = " "),"</p>"),file=file)
    df <- x$df
    rdf <- df[2]
    HTML("\n<p>Parameters:\n</p>",file=file)
    HTML.coefmat(x$parameters, digits = digits, signif.stars = signif.stars,
        file=file,append=TRUE,...)
    HTML(paste("\n<p>Residual standard error:<b> ", format(signif(x$sigma,
        digits)), " </b>on <b>", rdf, " </b>degrees of freedom\n</p>"),file=file)
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- dim(correl)[2]
        if (p > 1) {
            HTML("\n<p>Correlation of Parameter Estimates:\n</p>",file=file)
            if (symbolic.cor)
                HTML(symnum(correl)[-1, -p],file=file)
            else {
                correl[!lower.tri(correl)] <- NA
                HTML(correl[-1, -p, drop = FALSE], file=file)
            }
        }
    }
    HTMLbr(file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#
### PACKAGE STEPFUN
#----------------------------------------------------------------------------------------------------#

"HTML.ecdf" <- function (x, digits = getOption("digits") - 2, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    numform <- function(x) paste(formatC(x, digits = digits), collapse = ", ")
    HTML(paste("<p><b>Empirical CDF</b></p> \n<li>Call:<font class='call'> ", paste(deparse(attr(x, "call")),collapse=" "),"</font>"), file=file)
    n <- length(xx <- eval(expression(x), envir = environment(x)))
    i1 <- 1:min(3, n)
    i2 <- if (n >= 4)
        max(4, n - 1):n
    else integer(0)
    HTML(paste(" x[1:", n, "] = ", numform(xx[i1]), if (n > 3)
        ", ", if (n > 5)
        " ..., ", numform(xx[i2]), "\n<br>", sep = ""),file=file)
    invisible(x)
}


#----------------------------------------------------------------------------------------------------#

 "HTML.stepfun" <- function (x, digits = getOption("digits") - 2, file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    numform <- function(x) paste(formatC(x, digits = digits), collapse = ", ")
    i1 <- function(n) 1:min(3, n)
    i2 <- function(n) if (n >= 4)
        max(4, n - 1):n
    else integer(0)
    HTML(paste("<p><b>Step function</b></p>\n<li>Call<font class='call'>: ",paste(deparse(attr(x, "call")) ,collapse=" "),"</font>"),file=file)
    env <- environment(x)
    n <- length(xx <- eval(expression(x), envir = env))
    HTML(paste(" x[1:", n, "] = ", numform(xx[i1(n)]), if (n > 3)
        ", ", if (n > 5)
        " ..., ", numform(xx[i2(n)]), "\n<br>", sep = ""),file=file)
    y <- eval(expression(c(yleft, y)), envir = env)
    HTML(paste(n + 1, " step heights = ", numform(y[i1(n + 1)]), if (n +
        1 > 3)
        ", ", if (n + 1 > 5)
        " ..., ", numform(y[i2(n + 1)]), "\n<br>", sep = ""),file=file)
    invisible(x)
}

#----------------------------------------------------------------------------------------------------#
### PACKAGE SURVIVAL
#----------------------------------------------------------------------------------------------------#

"HTML.cox.zph" <- function (x, digits = max(options()$digits - 4, 3), file=HTMLGetFile(), append=TRUE,...)
HTML(x$table, file=file,append=append,...)

#----------------------------------------------------------------------------------------------------#

"HTML.coxph.null" <- function (x, digits = max(options()$digits - 4, 3), file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
	if (!is.null(cl <- x$call)) HTMLli(paste("Call:",paste(deparse(cl),collapse=" ")),file=file)
    HTML(paste("<p>Null model  log likelihood=<b>", format(x$loglik), "</b>\n</p>"),file=file)
    omit <- x$na.action
    if (length(omit)) HTML(paste("<p>  n=<b>", x$n, " </b>(", naprint(omit), ")\n</p>", sep = ""),file=file)
    else HTML(paste("<p>  n=<b>", x$n, "</b>\n</p>"),file=file)
}

#----------------------------------------------------------------------------------------------------#
### XTABLE
#----------------------------------------------------------------------------------------------------#

"HTML.xtable" <- function(x,file=HTMLGetFile(), append=TRUE,...){
    cat("\n", file=file, append=append,...)
    cat(capture.output(print(x,type="html")),file=file,append=TRUE,sep="\n")
}

#----------------------------------------------------------------------------------------------------#
### UTILITARY FUNCTIONS
#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#

"HTML.title"<-
function(x, HR = 2,CSSclass=NULL,file=HTMLGetFile(), append=TRUE, ...)
{
	cat(paste("\n <h", HR, if(!is.null(CSSclass)) paste(" class=",CSSclass,sep="") ," > ", x, "</h", HR, ">\n", sep =
		""), file = file, append=append, sep = "")
}

#----------------------------------------------------------------------------------------------------#

"HTMLstem" <- function (x, file=HTMLGetFile(), append=TRUE,...)  HTML(paste("<pre>",paste(capture.output(stem(x)),collapse="<br>"),"</pre>"), file=file,append=append,...)



#----------------------------------------------------------------------------------------------------#

"HTMLbr"<- function(x=1,file=HTMLGetFile(), append=TRUE) { cat(paste("\n",rep("<br>&nbsp;",x),"\n",sep=""), append=append, file = file)}

#----------------------------------------------------------------------------------------------------#

"HTMLhr"<- function(file=HTMLGetFile(), Width = "100%", Size = "1",CSSclass=NULL,append=TRUE){ cat(paste("\n<hr ", ifelse(!is.null(CSSclass),paste("class=",CSSclass,sep=""),""), " width=", Width, " size=", Size, ">", sep = ""), file = file, append=append, sep = "")}

#----------------------------------------------------------------------------------------------------#

"HTMLli"<- function(txt="", file=HTMLGetFile(), append=TRUE) { cat(paste("\n<br><li>", txt, sep = ""), sep = "", append=append, file = file)}

#----------------------------------------------------------------------------------------------------#


"HTMLplot" <- function (Caption = "", file=HTMLGetFile(), append=TRUE, GraphDirectory = ".",
    GraphFileName = "", GraphSaveAs = "png", GraphBorder = 1, Align = "center",
    Width=500,Height=500,WidthHTML=NULL,HeightHTML=NULL,
    GraphPointSize=12,GraphBackGround="white",GraphRes=72,plotFunction=NULL,...)
{
## New version with code submitted by James Wettenhall <wettenhall@wehi.edu.au>
## Change  plotFunction by plotExpression...

    if (exists(".HTMLTmpEnv", where=.HTMLEnv))
    {
       GraphDirectory <- get(".HTML.outdir", envir=get(".HTMLTmpEnv", envir=.HTMLEnv))
    }

    cat("\n", file=file, append=append,...)
    if (GraphFileName == "") {
        nowd <- date()
        GraphFileName <- paste("GRAPH_", substring(nowd, 5, 7),
            substring(nowd, 9, 10), "_", substring(nowd, 12,
                13), substring(nowd, 15, 16), substring(nowd,
                18, 19), sep = "")
    }

     GraphFileName <- paste(GraphFileName, ".", GraphSaveAs, sep = "")
     AbsGraphFileName <- file.path(GraphDirectory, GraphFileName)

    if (GraphSaveAs=="png")
      {
        if (is.null(plotFunction))
          dev.print(device=png, file = AbsGraphFileName, width=Width,height=Height,pointsize=GraphPointSize,bg=GraphBackGround)
        else
        {
          if (exists("X11", envir =.GlobalEnv) && Sys.info()["sysname"] != "Windows" && Sys.info()["sysname"] != "Darwin")
            bitmap(file = AbsGraphFileName,bg=GraphBackGround,res=GraphRes)
          else
            png(filename = AbsGraphFileName, width=Width,height=Height,pointsize=GraphPointSize,bg=GraphBackGround)
          plotFunction()
          dev.off()
        }
      }
      else if (GraphSaveAs %in% c("jpg","jpeg"))
      {
        if (is.null(plotFunction))
          dev.print(device=jpeg, file = AbsGraphFileName, width=Width,height=Height,pointsize=GraphPointSize,bg=GraphBackGround)
        else
        {
          if (exists("X11", envir =.GlobalEnv) && Sys.info()["sysname"] != "Windows" && Sys.info()["sysname"] != "Darwin")
            bitmap(filename = AbsGraphFileName,bg=GraphBackGround,res=GraphRes,type="jpeg")
          else
            jpeg(filename = AbsGraphFileName, width=Width,height=Height,pointsize=GraphPointSize,bg=GraphBackGround)
          plotFunction()
          dev.off()
        }
      }
      else if (GraphSaveAs=="gif")
      {
        stop("Gif support was removed from base R because of patent restrictions. Use either jpg or png")
#
#        if (is.null(plotFunction))
#  Gif support was removed from base R because of patent restrictions.
#  see http://tolstoy.newcastle.edu.au/R/help/05/02/12809.html
#          dev.print(device=gif, file = AbsGraphFileName, width=Width,height=Height,pointsize=GraphPointSize,bg=GraphBackGround)
#
#        else
#        {
#          stop("When passing a plot function to HTMLplot, device must be jpg or png.")
#        }
      }
    else stop("GraphSaveAs must be either jpg, png or gif")

    cat(paste("<p align=", Align, "><img src='", GraphFileName,
        "' border=", GraphBorder, if (!is.null(Width)) paste(" width=",Width,sep="") else "",if (!is.null(HeightHTML)) paste(" height=",HeightHTML,sep=""), if(!is.null(WidthHTML)) paste(" width="),">", sep = "", collapse = ""),
        file = file, append=TRUE, sep = "")
    if (Caption != "") {
        cat(paste("<br><font class=caption>", Caption, "</font>"), file = file, append=TRUE, sep = "")
    }
    cat("</p>", file = file, append=TRUE, sep = "\n")
    if (substitute(file)=="HTMLGetFile()") try(assign(".HTML.graph", TRUE, envir = .HTMLEnv))
    invisible(return(TRUE))
}

#----------------------------------------------------------------------------------------------------#

"HTMLInsertGraph" <- function(GraphFileName="",Caption="",GraphBorder=1,Align="center",WidthHTML=500,HeightHTML=NULL,file=HTMLGetFile(), append=TRUE,...)
{
    cat("\n", file=file, append=append,...)
    cat(paste("<p align=", Align, "><img src='", GraphFileName, "' border=", GraphBorder, if (!is.null(WidthHTML)) paste(" width=",WidthHTML,sep="") else "",if (!is.null(HeightHTML)) paste(" height=",HeightHTML,sep="") else "",">", sep = "", collapse = ""),         file = file, append=TRUE, sep = "")
    if (Caption != "") cat(paste("<br><i class=caption>", Caption, "</i>"), file = file, append=TRUE, sep = "")
    invisible(return(TRUE))
}


#----------------------------------------------------------------------------------------------------#

"HTMLCSS" <- function(file=HTMLGetFile(), append=TRUE,CSSfile="R2HTML.css")
{

  cat(paste("\n<link rel=stylesheet type=text/css href=",CSSfile,">\n",sep=""),file=file,append=append)

}


#----------------------------------------------------------------------------------------------------#
"HTMLChangeCSS" <- function(newCSS="R2HTML",from=NULL){
	target=getwd()
	if(exists(".HTMLTmpEnv", .HTMLEnv))
        target=file.path(get(".HTML.outdir",envir=get(".HTMLTmpEnv", .HTMLEnv)))

	if (is.null(from)){
##		from=file.path(.find.package(package = "R2HTML"),"output")
                from=system.file(package = "R2HTML","output")
	}
	fromfile=file.path(from,paste(newCSS,"css",sep="."))
	if (!file.exists(fromfile)) stop(paste("Source CSS file",fromfile,"not found"))
	file.copy(fromfile,file.path(target,"R2HTML.css"),overwrite=TRUE)

}


"HTMLCommand" <- function(x,file=HTMLGetFile(),Num="",menu=FALSE,target="index<-main.html",append=TRUE,...)
	{
	cat("\n",file=file,append=append,...)
	if (menu==TRUE)
	cat(paste("<br><li><a class=command href='./",target,"#Num",Num,"' target=main> ",paste(x,collapse=""),"</a>",sep=""),file=file,append=TRUE,sep="")
	else {
	if (Num!="") cat(paste("<a name=Num",Num,">&nbsp;</a>",sep=""),file=file,append=TRUE,sep="")
	cat(paste("\n<p><xmp class=command>> ",x,"</xmp></p>\n",sep=""),file=file,append=TRUE,sep="")
	}
	}

#----------------------------------------------------------------------------------------------------#

"HTMLcode" <- function(x,...){
	tmpfic=tempfile()
	HTML(x,file=tmpfic,...)
	cat("\n",file=tmpfic,append=TRUE)
	tmptxt=readLines(tmpfic)
	unlink(tmpfic)
	return(paste(tmptxt,collapse="\n"))
}
#----------------------------------------------------------------------------------------------------#


"HTMLReplaceNA"<-
function(Vec, Replace = " ")
{
	Vec <- as.character(Vec)
	#Vec <- format( Vec, ... )
	for(i in 1:length(Vec))
	{
		try(if((Vec[i] == "NA") | (Vec[i] == "NaN") | is.na(Vec[i])){ Vec[i] <- Replace})
	}
	Vec
}


#----------------------------------------------------------------------------------------------------#
"HTML.cormat" <- function(x, file=HTMLGetFile(),  digits=2,append=TRUE,align="center",caption="",captionalign="bottom",classcaption="captiondataframe",classtable="cormat",useCSS=TRUE,...)
{
	cat("\n", file=file,append=append)
	x<-as.matrix(x)
	if (is.numeric(x)) x<-round(x,digits=digits)
	if (is.null(dimnames(x))) x <- as.data.frame(x)
	txt <- paste("<p align=",align,">")
	txtcaption <- ifelse(is.null(caption),"",paste("<caption align=",captionalign," class=",classcaption,">",caption,"</caption>",sep=""))
	cormat=x
	abscormat=abs(cormat)
	backcolors=matrix(grey(1-as.matrix(abscormat)),ncol=ncol(cormat))
	css = 10*round(abs(x),1)
	css=matrix(paste("cor",unlist(css),sep=""),ncol=ncol(x))
	diag(css)="cordiag"
	diag(backcolors)="#FFFFFF"
	forecolors=matrix("#000000",ncol=ncol(cormat),nrow=nrow(cormat))
	forecolors[abscormat>0.5]="#FFFFFF"
	forecolors[abscormat>0.8]="#F6FF6E"
	diag(forecolors)="#FFFFFF"
	forebold=matrix(FALSE,ncol=ncol(cormat),nrow=nrow(cormat))
	forebold[abscormat>0.9]=TRUE
	txt<- paste(txt,"<table cellspacing=0 cellpading=0 border=0 >",txtcaption,"<td valign=middle class=corbody><table cellspacing=0 border=0>")
	txt <- paste(txt,paste("\n<tr><td align=right class=corvarname>",dimnames(x)[[2]],"</td><td width=2>&nbsp;</td></tr>",collapse="\n"))
	txt <- paste(txt,"</table></td><td valign=top class=corsep>&nbsp;</td><td valign=top>")
	txt <- paste(txt, "<table cellspacing=0 cellpadding=0 border=1 ><td><table class=",classtable," cellspacing=0>", sep = "")
	for(i in 1:dim(x)[1]) {
		VecDebut <- c(rep(paste("\n\t<td align=right", sep = ""), dim(x)[2]))
		if (useCSS) VecAttrib=c(paste(" class= ",css[i,],">")) else  VecAttrib=c(paste("  bgcolor=",backcolors[i,],"><font color=",forecolors[i,],">",ifelse(forebold[i,],"<b>","")))
		VecMilieu <- HTMLReplaceNA(as.matrix(x[i,  ]))
		VecFin <-  rep("</td>", dim(x)[2] )
		txt <- paste(txt, "\n<tr>",paste(VecDebut,VecAttrib, VecMilieu, VecFin, sep = "", collapse = ""),"</tr>")
		}
	txt <- paste(txt, "</table></td></table></td></table>")
	cat(txt, "\n", file = file, sep = "", append=TRUE,...)
	invisible(return(x))

	}

#----------------------------------------------------------------------------------------------------#

"as.title"<-
function(x)
{
	if (!is.character(x)) {
		x <- try(as.character(x))
		if (!is.character(x)) stop("Input argument must be of character mode")
	}
	class(x) <- "title"
	return(x)
}


#----------------------------------------------------------------------------------------------------#
###   R2HTML CORE
#----------------------------------------------------------------------------------------------------#

"HTMLStart" <- function(outdir=tempdir(),filename="index",extension="html",echo=FALSE, autobrowse=FALSE, HTMLframe=TRUE, withprompt="HTML> ",CSSFile="R2HTML.css",BackGroundColor="FFFFFF",BackGroundImg="",Title="R output")
{
	if (outdir!=tempdir())
	{
	# Copy of CSS and logo, if outdir != tempdir
		file.copy(file.path(tempdir(),'R2HTML.css'), file.path(outdir,'R2HTML.css'))
		file.copy(file.path(tempdir(),'R2HTMLlogo.gif'), file.path(outdir,'R2HTMLlogo.gif'))
	}
	.HTMLTmpEnv <- new.env(parent=.GlobalEnv)
	assign(".HTMLTmpEnv",.HTMLTmpEnv,envir=.HTMLEnv)
	assign("oldprompt",getOption("prompt"),envir=.HTMLTmpEnv)
	assign("HTMLframe",HTMLframe,envir=.HTMLTmpEnv)
	assign(".HTML.outdir",outdir,envir=.HTMLTmpEnv)
	assign("HTMLtorefresh",file.path(outdir,paste(filename,extension,sep=".")),envir=.HTMLTmpEnv)
	options(prompt=withprompt)
	assign(".HTML.graph",FALSE,envir =.HTMLTmpEnv)

	# Creation of required HTML files

	try(.HTML.file <- HTMLInitFile(outdir = outdir,filename=filename,extension=extension,HTMLframe=HTMLframe, BackGroundColor = BackGroundColor, BackGroundImg = BackGroundImg, Title = Title,CSSFile=CSSFile,useLaTeX=TRUE))
	assign(".HTML.file", .HTML.file, .HTMLTmpEnv)


	ToHTML <- function(file,echo,HTMLframe,HTMLMenuFile,target,outdir)
	{
		NumCom<-0
		function(expr,value,ok,visible)
		{

		NumCom<<- NumCom+1

		if (NumCom>1){

			ToPrint<-TRUE

			if (get(".HTML.graph",envir=.HTMLTmpEnv)==TRUE)
				{
					ToPrint <- FALSE
					assign(".HTML.graph",FALSE,envir=.HTMLTmpEnv)
				}
			else
				{
					if (length(expr)>1) {if ((expr[[1]]=="=")||(expr[[1]]=="<-")) ToPrint<-FALSE}


					# Print the commands and/or their number
					if (echo) HTMLCommand(deparse(expr),file,NumCom) else cat(paste("<a name=Num",NumCom,">&nbsp</a>",sep=""),file=file,sep="",append=TRUE)
					if (HTMLframe) HTMLCommand(deparse(expr),HTMLMenuFile,NumCom,menu=TRUE,target=target)
					if (ToPrint) HTML(value,file=file)
				}
		}

		if (autobrowse) browseURL(url=get("HTMLtorefresh",envir=.HTMLTmpEnv))
		invisible(return(TRUE))
		}
	}
	on.exit(addTaskCallback(ToHTML(.HTML.file,echo=echo,HTMLframe=HTMLframe,HTMLMenuFile=file.path(outdir,paste(filename,"_menu.",extension,sep="")),target=paste(filename,"_main.",extension,sep=""),outdir=outdir),name="HTML"),add=TRUE)
	cat("\n *** Output redirected to directory: ", outdir)
	cat("\n *** Use HTMLStop() to end redirection.")
	invisible(return(TRUE))

}
#----------------------------------------------------------------------------------------------------#

"HTMLInitFile"<-function(outdir = tempdir(),filename="index",extension="html",
		HTMLframe=FALSE, BackGroundColor = "FFFFFF", BackGroundImg = "",
		Title = "R output",CSSFile="R2HTML.css",useLaTeX=TRUE,useGrid=TRUE)
{
if (HTMLframe==FALSE){
	file<-file.path(outdir,paste(filename,".",extension,sep=""))
	assign(".HTML.file",file,envir =.HTMLEnv)

  txt <- ifelse(useLaTeX,"<html xmlns:mml=\"http://www.w3.org/1998/Math/MathML\">","<html>")
  #<HEAD>
    txt <- c(txt, "<head>")
    txt <- c(txt, paste("<title>",Title,"</title>"))
    # css
    txt <- c(txt, paste("<link rel=stylesheet href=\"",CSSFile,"\" type=text/css>",sep=""))
    # LaTeX ?
    if (useLaTeX)   txt <- c(txt, "<object id=\"mathplayer\" classid=\"clsid:32F66A20-7614-11D4-BD11-00104BD3F987\"></object>\n<?import namespace=\"mml\" implementation=\"#mathplayer\"?>\n<script type=\"text/javascript\" src=\"ASCIIMathML.js\"></script>")
    # Grid?
    if (useGrid) {
      txt <- c(txt, HTMLgrid_references())
      txt <- c(txt, "<script>\n   nequations=0;\n</script>")
    }
  # </HEAD>
  txt <- c(txt, "</head>")
  # <BODY>
  body <- c("<body")
  if(useLaTeX) body=c(body," onload=\"translate()\"")
  body=c(body,paste(" bgcolor=",BackGroundColor))
   if (BackGroundImg!="") body = c(body, paste(" background=\"",BackGroundImg,"\"",sep=""))
   body <- c(body," >")
   body=paste(body,collapse="")
   txt <- c(txt, body)
   txt <- paste(txt, collapse="\n")
   cat(txt, file=file,append=FALSE)

	}
else	{
	filemenu<-paste(filename,"_menu.",extension,sep="")
	filemain<-paste(filename,"_main.",extension,sep="")
	absfilemenu<-file.path(outdir,filemenu)
	file<-absfilemain<-file.path(outdir,filemain)
	absfileindex<-file.path(outdir,paste(filename,".",extension,sep=""))
	assign(".HTML.file",absfilemain,envir =.HTMLEnv)

	cat(paste("<html><head>	\n <title>",Title,"</title>\n <meta http-equiv=content-type content=text/html;charset=iso-8859-1>\n <frameset cols=250,* border=1 frameborder=yes><frame src=",filemenu," name=menu scrolling=yes><frame src=",filemain," name=main scrolling=yes></frameset></body></html>"), append = FALSE, sep = "", file = absfileindex)

	cat("<html><head><link rel=stylesheet href=",CSSFile," type=text/css> </head><body bgcolor=\"#E5F5FF\">  <center> <img src=R2HTMLlogo.gif> <hr size=1></center><br>",sep="",append=FALSE,file=absfilemenu)

     txt <- ifelse(useLaTeX,"<html xmlns:mml=\"http://www.w3.org/1998/Math/MathML\">","<html>")
  #<HEAD>
    txt <- c(txt, "<head>")
    txt <- c(txt, paste("<title>",Title,"</title>"))
    # css
    txt <- c(txt, paste("<link rel=stylesheet href=\"",CSSFile,"\" type=text/css>",sep=""))
    # LaTeX ?
    if (useLaTeX)   txt <- c(txt, "<object id=\"mathplayer\" classid=\"clsid:32F66A20-7614-11D4-BD11-00104BD3F987\"></object>\n<?import namespace=\"mml\" implementation=\"#mathplayer\"?>\n<script type=\"text/javascript\" src=\"ASCIIMathML.js\"></script>")
   # Grid?
    if (useGrid) {
      txt <- c(txt, HTMLgrid_references())
      txt <- c(txt, "<script>\n   nequations=0;\n</script>")
    }  # </HEAD>
  txt <- c(txt, "</head>")
  # <BODY>
  body <- c("<body")
  if(useLaTeX) body=c(body," onload=\"translate()\"")
  body=c(body,paste(" bgcolor=",BackGroundColor))
   if (!is.null(BackGroundImg)) body = c(body, paste(" background=\"",BackGroundImg,"\"",sep=""))
   body <- c(body," >")
   body=paste(body,collapse="")
   txt <- c(txt, body)
   txt <- paste(txt, collapse="\n")
   cat(txt, file=absfilemain,append=FALSE)

}

	invisible(return(file))
}

#----------------------------------------------------------------------------------------------------#

"HTMLEndFile"<- function(file=HTMLGetFile())
{
	cat("\n<hr size=1>\n<font size=-1>\n\t Generated on: <i>", date(),
		"</i> - <b>R2HTML</b> \n<hr size=1>\n\t</body>\n</html>",
		sep = "", append=TRUE, file = file)
}


#----------------------------------------------------------------------------------------------------#

"HTMLStop"<-function()
{
	invisible(removeTaskCallback("HTML"))
	.HTMLTmpEnv <- get(".HTMLTmpEnv", envir=.HTMLEnv)
	options(prompt=get("oldprompt",envir=.HTMLTmpEnv))
	.tmp=get(".HTML.file",envir=.HTMLTmpEnv)
	HTMLEndFile(file=get(".HTML.file",envir=.HTMLTmpEnv))
	rm(".HTMLTmpEnv", envir=.HTMLEnv)
	invisible(return(.tmp))
}

#----------------------------------------------------------------------------------------------------#
# Function contributed by Gabor Grothendieck (ggrothendieck_at_gmail.com)

HTML2clip <- function(x, filename = file("clipboard", ifelse(.Platform$OS == "windows","w",stop("Writing to clipboard only supported on Windows"))), append = FALSE, ...) {
    HTML(x, file = filename, append = append, ...)
}

#----------------------------------------------------------------------------------------------------#


# "myunzip"   <-  function (zipname, dest)
# {
#     if (file.exists(zipname)) {
#       if (.Platform$OS.type=="unix")  system(paste(getOption("unzip"), "-oq", zipname, "-d", dest))
#       else .Internal(int.unzip(zipname, NULL, dest))
#     }
#     else stop(paste("zipfile", zipname, "not found"))
# }

".onLoad" <- function(lib,pkg)
{
	#cat("\nLoading R2HTML package...\n")
	#ps.options(bg="white")

  # Copy all the content of "output" directory to tempdir()
  # now we use a zip file as there are subdirectories...
   unzip(zipfile=file.path(lib,pkg,'output','R2HTMLstuff.zip'),exdir=tempdir())

  options(R2HTML.CSSdir=file.path(lib,pkg,"output"))
  options(R2HTML.sortableDF=FALSE)
  options(R2HTML.format.digits=2)
  options(R2HTML.format.nsmall=0)
  options(R2HTML.format.big.mark="")
  options(R2HTML.format.big.interval=3)
  options(R2HTML.format.decimal.mark=Sys.localeconv()[["decimal_point"]])
  options(R2HTML.grid.first=TRUE)
  options(R2HTML.grid.stuffbasepath="./")

}


options(R2HTML.sortableDF=FALSE)
options(R2HTML.format.digits=2)
options(R2HTML.format.nsmall=0)
options(R2HTML.format.big.mark="")
options(R2HTML.format.big.interval=3)
options(R2HTML.format.decimal.mark=Sys.localeconv()[["decimal_point"]])
options(R2HTML.grid.first=TRUE)
options(R2HTML.grid.stuffbasepath="./")

