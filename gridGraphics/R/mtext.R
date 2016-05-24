
# C_mtext(text, side, line, outer, at, adj, padj, cex, col, font, ...) */

C_mtext <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:11)])
    dev.set(playDev())
    text <- x[[2]]
    side <- x[[3]]
    line <- x[[4]]
    outer <- x[[5]]
    adj <- ComputeAdjValue(x[[7]], side, par$las)
    at <- ComputeAtValue(x[[6]], adj, side, par$las)
    padj <- ComputePAdjValue(x[[8]], side, par$las)
    # NOTE: default is not par$cex, but 1.0 
    # NOTE: yes, mtext() really does override 'col=NA', 'cex=NA', and 'font=NA'
    cex <- FixupCex(x[[9]], 1)
    # NOTE: deliberately reverse any auto scaling of 0.66 or 0.83 
    cex <- ifelse(is.finite(cex), cex, unadjustedCex(par))
    col <- FixupCol(x[[10]], NA, par$bg)
    col <- ifelse(is.na(col), par$col, col)
    font <- FixupFont(x[[11]], par$font)
    font <- ifelse(is.na(font), par$font, font)
    if (outer) {
        depth <- gotovp(NA, "inner")
    } else {
        # NOTE: there is a bug in C_mtext() in plot.c where it checks
        #       "if (outer)" when 'outer' is still an SEXP, so the
        #       result is ALWAYS TRUE, so xpd is ALWAYS set to 2
        depth <- gotovp(NA, "window")
        # depth <- gotovp(if (is.na(par$xpd)) NA else TRUE, "window")
    }
    name <- paste0("mtext-", switch(side, "bottom", "left", "top", "right"))
    if (outer)
        name <- paste(name, "outer", sep="-")
    GMtext(text, side, line, outer, at,
           las=par$las, xadj=adj, yadj=padj,
           mex=par$mex, cin=par$cin, cex=cex, linecex=par$cex,
           font=font, family=par$family,
           col=col, lheight=par$lheight, par$ylbias,
           label=name)
    upViewport(depth)
}

# Helpers for C_mtext()
unadjustedCex <- function(par) {
    # If par(mfrow/mfcol) is exactly 2x2 then there is a 0.83 scaling
    # If it is more than 3 in either dimension there is a 0.66 scaling
    # We want to reverse that scaling here
    nr <- par$mfrow[1]
    nc <- par$mfrow[2]
    if (nr == 2 && nc == 2) {
        par$cex/0.83
    } else if (nr > 2 || nc > 2) {
        par$cex/0.66
    } else {
        par$cex
    }
}

ComputeAdjValue <- function(adj, side, las) {
    if (is.finite(adj)) {
        adj
    } else {
	switch(las + 1,
               # las = 0
               0.5,
               # las = 1
               switch(side, 0.5, 1, 0.5, 0),
               # las = 2
               switch(side, 1, 1, 0, 0),
               # las = 3
               switch(side, 1, 0.5, 0, 0.5))
    }
}

ComputePAdjValue <- function(padj, side, las) {
    if (is.finite(padj)) {
        padj
    } else {
        switch(las+ 1,
               # las = 0
               0,
               # las = 1
               switch(side, 0, 0.5, 0, 0.5),
               # las = 2
               0.5,
               # las = 3
               switch(side, 0.5, 0, 0.5, 0))
    }
}

ComputeAtValue <- function(at, adj, side, las) {
    if (any(is.finite(at))) {
        unit(at, "native")
    } else {
	# If the text is parallel to the axis, use "adj" for "at"
	# Otherwise, centre the text
	unit(switch(las + 1,
	            # parallel to axis 
                    adj,
                    # horizontal 
                    switch(side, adj, 0.5, adj, 0.5),
                    # perpendicular to axis
                    0.5,
                    # vertical 
                    switch(side, 0.5, adj, 0.5, adj)), "npc")
    }
}

# Code to centralise the work that GMtext() does to mess around with
# the (x, y) locations that it is sent
# (so that those fiddly adjustments are not reproduced all over the place)

# NOTE that 'mgp' is in 'mex' units
# NOTE we use par("cin") rather than 'grid' "lines"
# NOTE that 'linecex' attempts to capture the fact that line height is based
#      on 'mex'*'cexbase' NOT 'cex'*'cexbase'

GMtext <- function(str, side, line, outer=FALSE, at, las, xadj, yadj,
                   mex, cin, cex, linecex, font, family,
                   col, lheight, yLineBias,
                   allowOverlap=TRUE, label) {
    if (side == 1) {
        if (las == 2 || las == 3) {
            angle <- 90
        } else {
            line <- line + 1/mex*(1 - yLineBias)
            angle <- 0
        }
        x <- at
        y <- unit(-line*cin[2]*linecex, "in")
    } else if (side == 2) {
        if(las == 1 || las == 2) {
	    angle <- 0
	} else {
	    line <- line + 1/mex*yLineBias
	    angle <- 90
	}
        x <- unit(-line*cin[2]*linecex, "in")
        y <- at
    } else if (side == 3) {
        if(las == 2 || las == 3) {
	    angle <- 90
	}
	else {
	    line <- line + 1/mex*yLineBias
	    angle <- 0
	}
        x <- at
        y <- unit(1, "npc") + unit(line*cin[2]*linecex, "in")
    } else if (side == 4) {
	if(las == 1 || las == 2) {
	    angle <- 0
	}
	else {
	    line <- line + 1/mex*(1 - yLineBias)
	    angle <- 90
	}
        x <- unit(1, "npc") + unit(line*cin[2]*linecex, "in")
        y <- at
    } else {
        stop("Invalid 'side'")
    }
    grid.text(str, x, y, hjust=xadj, vjust=yadj, rot=angle,
              gp=gpar(cex=cex, fontface=font, fontfamily=family,
                  col=col, lineheight=lheight),
              check.overlap=!allowOverlap,
              name=grobname(label))
}

