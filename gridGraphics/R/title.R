
# C_title(main, sub, xlab, ylab, line, outer, ...)
C_title <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:7)])
    dev.set(playDev())
    mainArg <- GetTextArg(x[[2]],
                          par$cex.main*par$cex,
                          par$col.main, par$font.main, par$bg)
    main <- mainArg$txt
    subArg <- GetTextArg(x[[3]],
                          par$cex.sub*par$cex,
                          par$col.sub, par$font.sub, par$bg)
    sub <- subArg$txt
    xlabArg <- GetTextArg(x[[4]],
                          par$cex.lab*par$cex,
                          par$col.lab, par$font.lab, par$bg)
    xlab <- xlabArg$txt
    ylabArg <- GetTextArg(x[[5]], 
                          par$cex.lab*par$cex,
                          par$col.lab, par$font.lab, par$bg)
    ylab <- ylabArg$txt
    line <- x[[6]]
    outer <- x[[7]]
    if (outer) {
        depth <- gotovp(NA, "inner")
    } else {
        depth <- gotovp(if (is.na(par$xpd)) NA else TRUE, "plot")
    }        
    if (!is.null(main)) {
	if (outer) {
	    if (is.finite(line)) {
		vpos <- line
		adjy <- 0
	    } else {
		vpos <- 0.5*par$oma[3]
		adjy = 0.5
	    }
	} else {
	    if (is.finite(line)) {
		vpos = line
		adjy = 0
	    } else {
		vpos = 0.5*par$mar[3]
		adjy = 0.5
	    }
	}
        grid.text(main,
                  x=unit(par$adj, "npc"),
                  y=unit(1, "npc") +
                    unit(vpos*par$cex*par$cin[2], "in"),
                  vjust=adjy,
                  gp=gpar(cex=mainArg$pars$cex,
                          fontface=mainArg$pars$font, fontfamily=par$family,
                          col=mainArg$pars$col, lineheight=par$lheight),
                  name=grobname("main"))
    }
    if (!is.null(sub)) {
        cex <- if (is.null(subArg$pars$cex)) {
            par$cex
        } else {
            cex <- subArg$pars$cex
        }
        GMtext(sub, 1, line=par$mgp[1] + 1, at=0.5, las=0, xadj=0.5, yadj=0,
               mex=par$mex, cin=par$cin, cex=subArg$pars$cex,
               linecex=par$mex*par$cex, font=subArg$pars$font,
               family=par$family, col=subArg$pars$col,
               lheight=par$lheight, yLineBias=par$ylbias, label="sub")
    }
    if (!is.null(xlab)) {
        cex <- if (is.null(xlabArg$pars$cex)) {
            par$cex
        } else {
            cex <- xlabArg$pars$cex
        }
        GMtext(xlab, 1, line=par$mgp[1], at=0.5, las=0, xadj=0.5, yadj=0, 
               mex=par$mex, cin=par$cin, cex=xlabArg$pars$cex,
               linecex=par$mex*par$cex, font=xlabArg$pars$font,
               family=par$family, col=xlabArg$pars$col,
               lheight=par$lheight, yLineBias=par$ylbias,
               label="xlab")
    }
    if (!is.null(ylab)) {
        cex <- if (is.null(ylabArg$pars$cex)) {
            par$cex
        } else {
            cex <- ylabArg$pars$cex
        }
        GMtext(ylab, 2, line=par$mgp[1], at=0.5, las=0, xadj=0.5, yadj=0,
               mex=par$mex, cin=par$cin, cex=ylabArg$pars$cex,
               linecex=par$mex*par$cex, font=ylabArg$pars$font,
               family=par$family, col=ylabArg$pars$col,
               lheight=par$lheight, yLineBias=par$ylbias,
               label="ylab")
    }
    upViewport(depth)
}

GetTextArg <- function(x, cex, col, font, bg) {
    # Text may be specified as list(text, col=, cex=, font=)
    pars <- list(cex=cex, col=col, font=font)
    if (is.null(x) || is.na(x)) {
        txt <- NULL
    } else if (is.language(x) || is.character(x)) {
        txt <- x
    } else {
        if (is.list(x)) {
            names <- names(x)
            if (length(names)) {
                parnames <- names %in% c("col", "cex", "font")
                if (any(parnames)) {
                    if ("cex" %in% names) {
                        pars$cex = x$cex
                    }
                    if ("col" %in% names) {
                        pars$col <- FixupCol(x$col, NA, bg)
                    }
                    if ("font" %in% names) {
                        pars$font <- FixupFont(x$font, NA)
                    }
                } 
                nonpars <- x[!parnames]
                if (length(nonpars)) {
                    txt <- x[[length(nonpars)]]
                } else {
                    txt <- NULL
                }
            } else {
                txt <- x[[1]]
            }
        } else {
            stop("Unrecognised text argument type")
        }
    }
    list(txt=txt, pars=pars)
}
