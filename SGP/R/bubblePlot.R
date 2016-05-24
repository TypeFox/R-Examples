`bubblePlot` <- function(
	bubble_plot_data.X,
	bubble_plot_data.Y,
	bubble_plot_data.SUBSET=NULL,
	bubble_plot_data.INDICATE=NULL,
    bubble_plot_data.BUBBLE_CENTER_LABEL=NULL,
	bubble_plot_data.SIZE,
	bubble_plot_data.LEVELS=NULL,
	bubble_plot_data.BUBBLE_TIPS_LINES,
	bubble_plot_labels.X=c("Growth", "Median Student Growth Percentile"),
	bubble_plot_labels.Y=c("Achievement", "Percent at/above Proficient"),
	bubble_plot_labels.SIZE=c(50, 100, 500, 1000),
	bubble_plot_labels.LEVELS=NULL,
	bubble_plot_labels.BUBBLE_TIPS_LINES=list("Median SGP (Count)", "Percent at/above Proficient"),
	bubble_plot_labels.BUBBLE_TITLES,
	bubble_plot_titles.MAIN="Growth and Achievement",
	bubble_plot_titles.SUB1="State School Performance",
	bubble_plot_titles.SUB2="Growth & Current Achievement",
	bubble_plot_titles.LEGEND1="School Size",
	bubble_plot_titles.LEGEND2_P1=NULL,
	bubble_plot_titles.LEGEND2_P2=NULL,
	bubble_plot_titles.NOTE = NULL,
	bubble_plot_configs.BUBBLE_MIN_MAX=c(0.03, 0.03),
	bubble_plot_configs.BUBBLE_X_TICKS=seq(0,100,10),
	bubble_plot_configs.BUBBLE_X_TICKS_SIZE=c(rep(0.6, 5), 1, rep(0.6, 5)),
	bubble_plot_configs.BUBBLE_X_BANDS=NULL,
	bubble_plot_configs.BUBBLE_X_BAND_LABELS=NULL,
	bubble_plot_configs.BUBBLE_Y_TICKS=seq(0,100,10),
	bubble_plot_configs.BUBBLE_Y_TICKS_SIZE=rep(0.6, 11),
	bubble_plot_configs.BUBBLE_Y_BANDS=NULL,
	bubble_plot_configs.BUBBLE_Y_BAND_LABELS=NULL,
	bubble_plot_configs.BUBBLE_SUBSET_INCREASE=0,
	bubble_plot_configs.BUBBLE_SUBSET_ALPHA=list(Transparent=0.3, Opaque=0.95),
	bubble_plot_configs.BUBBLE_COLOR="deeppink2",
    bubble_plot_configs.BUBBLE_COLOR_GRADIENT_REVERSE=FALSE,
	bubble_plot_configs.BUBBLE_TIPS=TRUE,
	bubble_plot_configs.BUBBLE_PLOT_DEVICE="PDF",
	bubble_plot_configs.BUBBLE_PLOT_FORMAT="print",
	bubble_plot_configs.BUBBLE_PLOT_LEGEND=FALSE,
	bubble_plot_configs.BUBBLE_PLOT_TITLE=TRUE,
	bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS=c("Growth", "Achievement"),
	bubble_plot_configs.BUBBLE_PLOT_EXTRAS=NULL,
    bubble_plot_configs.BUBBLE_PLOT_DIMENSION=NULL, ## List of WIDTH and HEIGHT
	bubble_plot_configs.BUBBLE_PLOT_NAME="bubblePlot.pdf",
	bubble_plot_configs.BUBBLE_PLOT_PATH=paste("Figures", sep=""),
	bubble_plot_pdftk.CREATE_CATALOG=FALSE) {


# Test for data to plot

if (length(bubble_plot_data.X)==0) {
	return("No data supplied for plotting. No plot produced")
}


# Test for installation of pdf2 package

if (bubble_plot_configs.BUBBLE_TIPS) {
	if (length(find.package("pdf2", quiet=TRUE)) > 0 & as.numeric(version$minor) < 14) {
		eval(parse(text="require(pdf2)"))
	} else {
		bubble_plot_configs.BUBBLE_TIPS <- FALSE
#		message("\tImplentation of BUBBLE_TIPS requires the installation of the package pdf2 from R-Forge: install.packages('pdf2',repos='http://R-Forge.R-project.org')")
	}
}


# Create directory for file

bubble_plot_configs.BUBBLE_PLOT_NAME <- gsub("/", "-", bubble_plot_configs.BUBBLE_PLOT_NAME)
if (!is.null(bubble_plot_configs.BUBBLE_PLOT_PATH)) {
        dir.create(bubble_plot_configs.BUBBLE_PLOT_PATH, recursive=TRUE, showWarnings=FALSE)
	file.path.and.name <- file.path(bubble_plot_configs.BUBBLE_PLOT_PATH, bubble_plot_configs.BUBBLE_PLOT_NAME)
} else {
	file.path.and.name <- bubble_plot_configs.BUBBLE_PLOT_NAME
}


# Calculate relevant quantities

if (!is.null(bubble_plot_labels.SIZE)) {
	numstud.range <- c(min(bubble_plot_labels.SIZE), max(bubble_plot_labels.SIZE))
} else {
	numstud.range <- c(25, 100)
}
num.sizes <- length(bubble_plot_labels.SIZE)
if (!missing(bubble_plot_data.BUBBLE_TIPS_LINES)) num.bubble.lines <- length(bubble_plot_data.BUBBLE_TIPS_LINES)
if (is.null(bubble_plot_data.LEVELS)) {
   num.levels <- 1; tmp.LEVELS <- rep(1, length(bubble_plot_data.X))
} else {
   num.levels <- length(unique(bubble_plot_labels.LEVELS)); tmp.LEVELS <- bubble_plot_data.LEVELS
}
if (!is.null(bubble_plot_configs.BUBBLE_COLOR)) {
    if (num.levels==1) {
	my.colors <- bubble_plot_configs.BUBBLE_COLOR
    } else {
       temp.colors <- rgb2hsv(col2rgb(bubble_plot_configs.BUBBLE_COLOR))
       my.colors <- hsv(h=temp.colors[1], s=1:num.levels/(num.levels+1), v=temp.colors[3])
       if (bubble_plot_configs.BUBBLE_COLOR_GRADIENT_REVERSE) my.colors <- rev(my.colors)
    }
} else {
     my.colors <- rev(rainbow_hcl(num.levels))
}

if (bubble_plot_configs.BUBBLE_PLOT_FORMAT=="print") {
     format.colors.background <- rgb(0.985, 0.985, 1.0)
     format.colors.border <- "grey20"
     format.colors.font <- c("grey20", rgb(0.985, 0.985, 1.0))
     format.colors.quadrant <- c(rgb(0.885, 0.885, 0.885), rgb(0.985, 0.985, 1.0))
} else {
     format.colors.background <- rgb(0.48, 0.48, 0.52)
     format.colors.border <- rgb(0.985, 0.985, 1.0)
     format.colors.font <- c(rgb(0.985, 0.985, 1.0), rgb(0.48, 0.48, 0.52))
     format.colors.quadrant <- c(rgb(0.885, 0.885, 0.885), rgb(0.985, 0.985, 1.0))
}


# Custom Color Function

bubblecolor <- function(x){
           temp <- character(length(x))
           for (i in 1:num.levels){
           temp[x == i] <- my.colors[i]
           }
           temp[is.na(x)] <- NA
           return(temp)
}


# Custom Bubble Size Function

bubblesize <- function(schoolsize, numstud.range) {
		slope <- (max.cex - min.cex)/(sqrt(numstud.range)[2] - sqrt(numstud.range)[1])
		temp <- slope*sqrt(schoolsize) - slope*sqrt(numstud.range)[2] + max.cex
		temp[temp < min.cex] <- min.cex; temp[temp > max.cex] <- max.cex
		return(temp)
}


# Custom Bubble Alpha Function

bubblealpha <- function(numbubbles, current.alpha) {
      if (numbubbles > 0 & numbubbles <= 100) return(1*current.alpha)
      if (numbubbles > 100 & numbubbles <= 250) return(0.85*current.alpha)
      if (numbubbles > 250 & numbubbles <= 500) return(0.7*current.alpha)
      if (numbubbles > 500 & numbubbles <= 1000) return(0.55*current.alpha)
      if (numbubbles > 1000) return(0.4*current.alpha)
}


# Indicator Tag Coordinate Function

indicate.tip <- function(x, y) {
          tmp.orientation <- character(2)
          if (y >= 0.8) {
                  tmp.orientation[2] <- "top"; tmp.y <- y-0.1
          } else {
                  tmp.orientation[2] <- "bottom"; tmp.y <- y+0.1
          }
          if (x <= 0.2 | (x >= 0.5 & x <= 0.8)) {
                  tmp.orientation[1] <- "left"; tmp.x <- x+0.1
          } else {
                  tmp.orientation[1] <- "right"; tmp.x <- x-0.1
          }
          list(x=tmp.x, y=tmp.y, orientation=tmp.orientation)
}


# Create viewports

if (bubble_plot_configs.BUBBLE_PLOT_LEGEND) {
    if (!is.null(bubble_plot_configs.BUBBLE_PLOT_DIMENSION)) {
        fig.width <- bubble_plot_configs.BUBBLE_PLOT_DIMENSION$WIDTH
        fig.height <- bubble_plot_configs.BUBBLE_PLOT_DIMENSION$HEIGHT
        text.buffer <- 0.1*fig.width/13
        text.start <- 0.7*fig.width/13
        if (!is.null(bubble_plot_configs.BUBBLE_MIN_MAX)) {
              min.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[1]*fig.width/13
              max.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[2]*fig.width/13
        } else {
              min.cex <- .01*fig.width/13
              max.cex <- .14*fig.width/13
        }
    } else {
        if (!is.null(bubble_plot_configs.BUBBLE_MIN_MAX)) {
              min.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[1]; max.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[2]
        } else {
              min.cex <- .01; max.cex <- .14
        }
        fig.width <- 13; fig.height <- 8.5; text.buffer <- 0.1; text.start <- 0.7
    }

    if (bubble_plot_configs.BUBBLE_PLOT_TITLE) {
        figure.vp <- viewport(layout = grid.layout(3, 3, widths = unit(c(0.8, 9.5, 2.7)*fig.width/13, rep("inches", 3)),
                              heights = unit(c(1.5, 6.2, 0.8)*fig.height/8.5, rep("inches", 3))),
                              gp=gpar(cex=fig.width/13))

        title.vp <- viewport(name="title.vp",
                    layout.pos.row=1, layout.pos.col=1:3,
                    xscale=c(0,1),
                    yscale=c(0,1),
                    gp=gpar(fill="transparent"))
    } else {
        figure.vp <- viewport(layout = grid.layout(3, 3, widths = unit(c(0.8, 9.5, 2.7)*fig.width/13, rep("inches", 3)),
                              heights = unit(c(0.2, 7.6, 0.8)*fig.height/8.5, rep("inches", 3))),
                              gp=gpar(cex=fig.width/13))
    }

        right.legend.vp <- viewport(name="right.top.legend.vp",
                  layout.pos.row=2, layout.pos.col=3,
                  xscale=c(0,1),
                  yscale=c(0,1),
                  gp=gpar(fill="transparent"))
} else {
     if (!is.null(bubble_plot_configs.BUBBLE_PLOT_DIMENSION)) {
        fig.width <- bubble_plot_configs.BUBBLE_PLOT_DIMENSION$WIDTH
        fig.height <- bubble_plot_configs.BUBBLE_PLOT_DIMENSION$HEIGHT
        text.buffer <- 0.1*fig.width/12
        text.start <- 0.7*fig.width/12
        if (!is.null(bubble_plot_configs.BUBBLE_MIN_MAX)) {
              min.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[1]*fig.width/12
              max.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[2]*fig.width/12
        } else {
              min.cex <- .01*fig.width/12
              max.cex <- .14*fig.width/12
        }
     } else {
        if (!is.null(bubble_plot_configs.BUBBLE_MIN_MAX)) {
              min.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[1]; max.cex <- bubble_plot_configs.BUBBLE_MIN_MAX[2]
        } else {
              min.cex <- .01; max.cex <- .14
        }
        fig.width <- 12; fig.height <- 8.5; text.buffer <- 0.1; text.start <- 0.7
     }

     if (bubble_plot_configs.BUBBLE_PLOT_TITLE) {
         figure.vp <- viewport(layout = grid.layout(3, 3, widths = unit(c(0.8, 10.7, 0.5)*fig.width/12, rep("inches", 3)),
                              heights = unit(c(1.5, 6.2, 0.8)*fig.height/8.5, rep("inches", 3))),
                              gp=gpar(cex=fig.width/12))

         title.vp <- viewport(name="title.vp",
                    layout.pos.row=1, layout.pos.col=1:3,
                    xscale=c(0,1),
                    yscale=c(0,1),
                    gp=gpar(fill="transparent"))
     } else {
         figure.vp <- viewport(layout = grid.layout(3, 3, widths = unit(c(0.8, 10.7, 0.5)*fig.width/12, rep("inches", 3)),
                              heights = unit(c(0.2, 7.6, 0.8)*fig.height/8.5, rep("inches", 3))),
                              gp=gpar(cex=fig.width/12))
    }
}

vaxis.vp <- viewport(name="vaxis.vp",
                       layout.pos.row=2, layout.pos.col=1,
                       xscale=c(0,1),
                       yscale=extendrange(bubble_plot_configs.BUBBLE_Y_TICKS, f=0.025),
                       gp=gpar(fill="transparent", cex=1.2))

chart.vp <- viewport(name="chart.vp",
                  layout.pos.row=2, layout.pos.col=2,
                  xscale=extendrange(bubble_plot_configs.BUBBLE_X_TICKS, f=0.025),
                  yscale=extendrange(bubble_plot_configs.BUBBLE_Y_TICKS, f=0.025),
                  gp=gpar(fill="transparent"))

haxis.vp <- viewport(name="haxis.vp",
                       layout.pos.row=3, layout.pos.col=2,
                       xscale=extendrange(bubble_plot_configs.BUBBLE_X_TICKS, f=0.025),
                       yscale=c(0,1),
                       gp=gpar(fill="transparent", cex=1.2))


# Set up device

if (bubble_plot_configs.BUBBLE_PLOT_DEVICE == "PDF") {
      pdf(file=file.path.and.name, width=fig.width, height=8.5, bg=format.colors.background, version="1.4")
}

if (bubble_plot_configs.BUBBLE_PLOT_DEVICE == "PNG") {
	Cairo(file=gsub(".pdf", ".png", file.path.and.name), width=fig.width, height=8.5, bg=format.colors.background, units="in", dpi=144, pointsize=10.5)
}

# Create plot (if bubble_plot_configs.BUBBLE_TIPS==TRUE)

if (bubble_plot_configs.BUBBLE_TIPS) {
    oldpar <- par(no.readonly = TRUE)
    plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), axes=FALSE, type = "n", xlab="", ylab="")
}


# Push figure.vp

pushViewport(figure.vp)


# Push chart.vp

pushViewport(chart.vp)

if (!is.null(bubble_plot_configs.BUBBLE_X_BANDS) | !is.null(bubble_plot_configs.BUBBLE_Y_BANDS)) {
	if (is.null(bubble_plot_configs.BUBBLE_X_BANDS) & !is.null(bubble_plot_configs.BUBBLE_Y_BANDS)) {
		for (i in seq_along(head(bubble_plot_configs.BUBBLE_Y_BANDS, -1))) {
			grid.roundrect(x=unit(0, "npc"), y=unit(bubble_plot_configs.BUBBLE_Y_BANDS[i], "native"),
				width=unit(1, "npc"), height=unit(diff(bubble_plot_configs.BUBBLE_Y_BANDS)[i], "native"),
				r=unit(0.1, "mm"), just=c("left", "bottom"), gp=gpar(fill=format.colors.quadrant[1], col=format.colors.background, lwd=0.4))
			grid.text(x=unit(0.05, "npc"), y=unit((bubble_plot_configs.BUBBLE_Y_BANDS[i] + diff(bubble_plot_configs.BUBBLE_Y_BANDS)[i]/2), "native"),
				bubble_plot_configs.BUBBLE_Y_BAND_LABELS[i], just="left",
				gp=gpar(cex=2, fontface=2, col=format.colors.quadrant[2]))
		}
	}
	if (!is.null(bubble_plot_configs.BUBBLE_X_BANDS) & is.null(bubble_plot_configs.BUBBLE_Y_BANDS)) {
		for (i in seq_along(head(bubble_plot_configs.BUBBLE_X_BANDS, -1))) {
			grid.roundrect(x=unit(bubble_plot_configs.BUBBLE_Y_BANDS[i], "native"), y=unit(0, "npc"),
				width=unit(diff(bubble_plot_configs.BUBBLE_Y_BANDS)[i], "native"), height=unit(1, "npc"),
				r=unit(0.1, "mm"), just=c("left", "bottom"), gp=gpar(fill=format.colors.quadrant[1], col=format.colors.background, lwd=0.4))
			grid.text(x=unit((bubble_plot_configs.BUBBLE_X_BANDS[i] + diff(bubble_plot_configs.BUBBLE_X_BANDS)[i]/2), "native"), y=unit(0.05, "npc"),
				bubble_plot_configs.BUBBLE_X_BAND_LABELS[i], just="left",
				gp=gpar(cex=2, fontface=2, col=format.colors.quadrant[2]))
		}
	}
	if (!is.null(bubble_plot_configs.BUBBLE_X_BANDS) & !is.null(bubble_plot_configs.BUBBLE_Y_BANDS)) {

	}
} else {
grid.rect(width=1, height=1, gp=gpar(fill=format.colors.quadrant[1], lwd=0.5, col=format.colors.border))
}

if (!is.null(bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS)) {
grid.text(x=0.05, y=0.15, paste("Lower", bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[1]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="left")
grid.text(x=0.05, y=0.08, paste("Lower",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[2]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="left")
grid.text(x=0.95, y=0.15, paste("Higher",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[1]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="right")
grid.text(x=0.95, y=0.08, paste("Lower",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[2]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="right")
grid.text(x=0.05, y=0.85, paste("Lower",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[1]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="left")
grid.text(x=0.05, y=0.92, paste("Higher",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[2]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="left")
grid.text(x=0.95, y=0.85, paste("Higher",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[1]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="right")
grid.text(x=0.95, y=0.92, paste("Higher",  bubble_plot_configs.BUBBLE_PLOT_BACKGROUND_LABELS[2]), gp=gpar(cex=1.8, fontface=2, col=format.colors.quadrant[2]), just="right")
}

# Add BUBBLE_PLOT_EXTRAS

base.line <- "grid.lines(x=unit(50, 'native'), y=c(0.03,0.97), gp=gpar(col='grey40', lwd=1.25, lty=2, alpha=0.5))"
if (!is.null(bubble_plot_configs.BUBBLE_PLOT_EXTRAS)) {
   for (i in c(base.line, bubble_plot_configs.BUBBLE_PLOT_EXTRAS)) {
      eval(parse(text=i))
   }
} else {
      eval(parse(text=base.line))
}

if (bubble_plot_configs.BUBBLE_TIPS) {
   if (!is.null(bubble_plot_data.SUBSET)) {
      grid.circle(x=bubble_plot_data.X, y=bubble_plot_data.Y, r=unit(bubblesize(bubble_plot_data.SIZE, numstud.range), rep("inches", length(bubble_plot_data.SIZE))),
               gp=gpar(col=rgb(0.4,0.4,0.4), lwd=0.05*bubble_plot_configs.BUBBLE_MIN_MAX[2]/0.12,
               fill=bubblecolor(unclass(tmp.LEVELS)), alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Transparent)), default.units="native")
      grid.circle(x=bubble_plot_data.X[bubble_plot_data.SUBSET], y=bubble_plot_data.Y[bubble_plot_data.SUBSET],
                  r=unit(bubble_plot_configs.BUBBLE_SUBSET_INCREASE+bubblesize(bubble_plot_data.SIZE[bubble_plot_data.SUBSET], numstud.range),
                        rep("inches", length(bubble_plot_data.SIZE[bubble_plot_data.SUBSET]))),
                  gp=gpar(lwd=0.750*bubble_plot_configs.BUBBLE_MIN_MAX[2]/0.12, fill=bubblecolor(unclass(tmp.LEVELS[bubble_plot_data.SUBSET])),
                  alpha=bubblealpha(length(bubble_plot_data.X[bubble_plot_data.SUBSET]), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque)), default.units="native")

     if (!is.null(bubble_plot_data.INDICATE)) {
          for (i in bubble_plot_data.INDICATE) {
              indicate.coordinates <- indicate.tip(as.numeric(convertX(unit(bubble_plot_data.X[i], "native"), "npc")),
                                               as.numeric(convertY(unit(bubble_plot_data.Y[i], "native"), "npc")))
              grid.segments(unit(indicate.coordinates$x, "npc"), unit(indicate.coordinates$y, "npc"),
                            unit(bubble_plot_data.X[i], "native"), unit(bubble_plot_data.Y[i], "native"),
                            gp=gpar(lwd=0.5))
              grid.circle(x=bubble_plot_data.X[i], y=bubble_plot_data.Y[i],
                      r=unit(c(1.0, 0.4)*bubblesize(bubble_plot_data.SIZE[i], numstud.range), rep("inches", length(bubble_plot_data.SIZE[i]))),
                      gp=gpar(lwd=c(0.5, 3.0), fill=bubblecolor(unclass(tmp.LEVELS[bubble_plot_data.INDICATE]))), default.units="native")
              grid.rect(x=unit(indicate.coordinates$x, "npc"), y=unit(indicate.coordinates$y, "npc"),
                             width=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.BUBBLE_TITLES[i]),
                             height=unit(1.5*text.buffer, "inches")+unit(1.0, "strheight", bubble_plot_labels.BUBBLE_TITLES[i]),
                             gp=gpar(col="grey20", lwd=0.7, fill=rgb(1.0, 0.94, 0.83, 0.6)), just=indicate.coordinates$orientation)
              if (indicate.coordinates$orientation[1]=="left") {
                        tmp.x <- unit(indicate.coordinates$x, "npc") + unit(text.buffer, "inches")
              } else {
                        tmp.x <- unit(indicate.coordinates$x, "npc") - unit(text.buffer, "inches")
              }
              if (indicate.coordinates$orientation[2]=="bottom") {
                        tmp.y <- unit(indicate.coordinates$y, "npc") + 0.75*unit(text.buffer, "inches")
              } else {
                        tmp.y <- unit(indicate.coordinates$y, "npc") - 0.75*unit(text.buffer, "inches")
              }
              grid.text(x=tmp.x, y=tmp.y, bubble_plot_labels.BUBBLE_TITLES[i], just=indicate.coordinates$orientation)
          }
     }

    if (!is.null(bubble_plot_data.BUBBLE_CENTER_LABEL)) {
            grid.text(x=bubble_plot_data.X, y=bubble_plot_data.Y, bubble_plot_data.BUBBLE_CENTER_LABEL,
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Transparent), font=2), default.units="native")
            grid.text(x=bubble_plot_data.X[bubble_plot_data.SUBSET], y=bubble_plot_data.Y[bubble_plot_data.SUBSET], bubble_plot_data.BUBBLE_CENTER_LABEL[bubble_plot_data.SUBSET],
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X[bubble_plot_data.SUBSET]), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque), font=2), default.units="native")
    }

    for (i in seq_along(bubble_plot_data.X[bubble_plot_data.SUBSET])) {
          pushViewport(viewport(x=unit(bubble_plot_data.X[bubble_plot_data.SUBSET][i], "native"),
                                y=unit(bubble_plot_data.Y[bubble_plot_data.SUBSET][i], "native"),
                                width=unit(1, "native"), height=unit(1, "native")))
          par(fig=gridFIG(), new = TRUE)
          tmp.bubble.txt <- character()
          for (j in seq(num.bubble.lines)) {
               tmp.bubble.txt <- c(tmp.bubble.txt, paste(bubble_plot_labels.BUBBLE_TIPS_LINES[[j]], ": ",
                                                         bubble_plot_data.BUBBLE_TIPS_LINES[[j]][bubble_plot_data.SUBSET][i], sep=""))
          }
          text(0.5, 0.5, "X", col=rgb(1,0,0,0.01), popup="PLACEHOLDER",
               cex=10*bubblesize(bubble_plot_data.SIZE[bubble_plot_data.SUBSET][i], numstud.range),
               annot.options=c(paste("/T (", bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.SUBSET][i], ")", sep=""),
                        paste("/Contents (", paste(tmp.bubble.txt, collapse="\n"), ")", sep="")))
          popViewport()
     } ## End for statement

   par(oldpar)
   } ## End SUBSET if statement

   else {
     grid.circle(x=bubble_plot_data.X, y=bubble_plot_data.Y, r=unit(bubblesize(bubble_plot_data.SIZE, numstud.range), rep("inches", length(bubble_plot_data.SIZE))),
               gp=gpar(col=rgb(0.2,0.2,0.2), lwd=0.05*bubble_plot_configs.BUBBLE_MIN_MAX[2]/0.12,
               fill=bubblecolor(unclass(tmp.LEVELS)), alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque)), default.units="native")

     if (!is.null(bubble_plot_data.INDICATE)) {
          for (i in bubble_plot_data.INDICATE) {
              indicate.coordinates <- indicate.tip(as.numeric(convertX(unit(bubble_plot_data.X[i], "native"), "npc")),
                                               as.numeric(convertY(unit(bubble_plot_data.Y[i], "native"), "npc")))
              grid.segments(unit(indicate.coordinates$x, "npc"), unit(indicate.coordinates$y, "npc"),
                            unit(bubble_plot_data.X[i], "native"), unit(bubble_plot_data.Y[i], "native"),
                            gp=gpar(lwd=0.5))
              grid.circle(x=bubble_plot_data.X[i], y=bubble_plot_data.Y[i],
                      r=unit(c(1.0, 0.4)*bubblesize(bubble_plot_data.SIZE[i], numstud.range), rep("inches", length(bubble_plot_data.SIZE[i]))),
                      gp=gpar(lwd=c(0.5, 3.0), fill=bubblecolor(unclass(tmp.LEVELS))), default.units="native")
              grid.rect(x=unit(indicate.coordinates$x, "npc"), y=unit(indicate.coordinates$y, "npc"),
                             width=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.BUBBLE_TITLES[i]),
                             height=unit(1.5*text.buffer, "inches")+unit(1.0, "strheight", bubble_plot_labels.BUBBLE_TITLES[i]),
                             gp=gpar(col="grey20", lwd=0.7, fill=rgb(1.0, 0.94, 0.83, 0.6)), just=indicate.coordinates$orientation)
              if (indicate.coordinates$orientation[1]=="left") {
                        tmp.x <- unit(indicate.coordinates$x, "npc") + unit(text.buffer, "inches")
              } else {
                        tmp.x <- unit(indicate.coordinates$x, "npc") - unit(text.buffer, "inches")
              }
              if (indicate.coordinates$orientation[2]=="bottom") {
                        tmp.y <- unit(indicate.coordinates$y, "npc") + 0.75*unit(text.buffer, "inches")
              } else {
                        tmp.y <- unit(indicate.coordinates$y, "npc") - 0.75*unit(text.buffer, "inches")
              }
              grid.text(x=tmp.x, y=tmp.y, bubble_plot_labels.BUBBLE_TITLES[i], just=indicate.coordinates$orientation)
          }
     }

     if (!is.null(bubble_plot_data.BUBBLE_CENTER_LABEL)) {
            grid.text(x=bubble_plot_data.X, y=bubble_plot_data.Y, bubble_plot_data.BUBBLE_CENTER_LABEL,
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque), font=2), default.units="native")
     }

     for (i in seq_along(bubble_plot_data.X)) {
         pushViewport(viewport(x=unit(bubble_plot_data.X[i], "native"),
                               y=unit(bubble_plot_data.Y[i], "native"),
                               width=unit(1, "native"), height=unit(1, "native")))
         par(fig=gridFIG(), new = TRUE)
         tmp.bubble.txt <- character()
         for (j in seq(num.bubble.lines)) {
              tmp.bubble.txt <- c(tmp.bubble.txt, paste(bubble_plot_labels.BUBBLE_TIPS_LINES[[j]], ": ",
                                                        bubble_plot_data.BUBBLE_TIPS_LINES[[j]][i], sep=""))
         }
         text(0.5, 0.5, "X", col=rgb(1,0,0,0.01), popup="PLACEHOLDER",
              cex=10*bubblesize(bubble_plot_data.SIZE[i], numstud.range),
              annot.options=c(paste("/T (", bubble_plot_labels.BUBBLE_TITLES[i], ")", sep=""),
                        paste("/Contents (", paste(tmp.bubble.txt, collapse="\n"), ")", sep="")))
         popViewport()
     } ## End for statement

   par(oldpar)
   } ## End SUBSET else statement
} ## End BUBBLE_TIPS if statement


else {
   if (!is.null(bubble_plot_data.SUBSET)){
      grid.circle(x=bubble_plot_data.X, y=bubble_plot_data.Y, r=unit(bubblesize(bubble_plot_data.SIZE, numstud.range), rep("inches", length(bubble_plot_data.SIZE))),
               gp=gpar(col=rgb(0.4,0.4,0.4), lwd=0.05*bubble_plot_configs.BUBBLE_MIN_MAX[2]/0.12,
               fill=bubblecolor(unclass(tmp.LEVELS)), alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Transparent)), default.units="native")
      grid.circle(x=bubble_plot_data.X[bubble_plot_data.SUBSET], y=bubble_plot_data.Y[bubble_plot_data.SUBSET],
                  r=unit(bubble_plot_configs.BUBBLE_SUBSET_INCREASE+bubblesize(bubble_plot_data.SIZE[bubble_plot_data.SUBSET], numstud.range),
                         rep("inches", length(bubble_plot_data.SIZE[bubble_plot_data.SUBSET]))),
                  gp=gpar(lwd=0.75*bubble_plot_configs.BUBBLE_MIN_MAX[2]/0.12, fill=bubblecolor(unclass(tmp.LEVELS[bubble_plot_data.SUBSET])),
                  alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque)), default.units="native")

   if (!is.null(bubble_plot_data.INDICATE)) {
          indicate.coordinates <- indicate.tip(as.numeric(convertX(unit(bubble_plot_data.X[bubble_plot_data.INDICATE], "native"), "npc")),
                                           as.numeric(convertY(unit(bubble_plot_data.Y[bubble_plot_data.INDICATE], "native"), "npc")))
          grid.segments(unit(indicate.coordinates$x, "npc"), unit(indicate.coordinates$y, "npc"),
                        unit(bubble_plot_data.X[bubble_plot_data.INDICATE], "native"), unit(bubble_plot_data.Y[bubble_plot_data.INDICATE], "native"),
                        gp=gpar(lwd=0.5))
          grid.circle(x=bubble_plot_data.X[bubble_plot_data.INDICATE], y=bubble_plot_data.Y[bubble_plot_data.INDICATE],
                  r=unit(c(1.0, 0.4)*bubblesize(bubble_plot_data.SIZE[bubble_plot_data.INDICATE], numstud.range), rep("inches", length(bubble_plot_data.SIZE[bubble_plot_data.INDICATE]))),
                  gp=gpar(lwd=c(0.5, 3.0), fill=bubblecolor(unclass(tmp.LEVELS[bubble_plot_data.INDICATE]))), default.units="native")
          grid.rect(x=unit(indicate.coordinates$x, "npc"), y=unit(indicate.coordinates$y, "npc"),
                         width=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.INDICATE]),
                         height=unit(1.5*text.buffer, "inches")+unit(1.0, "strheight", bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.INDICATE]),
                         gp=gpar(col="grey20", lwd=0.7, fill=rgb(1.0, 0.94, 0.83, 0.6)), just=indicate.coordinates$orientation)
          if (indicate.coordinates$orientation[1]=="left") {
                    tmp.x <- unit(indicate.coordinates$x, "npc") + unit(text.buffer, "inches")
          } else {
                    tmp.x <- unit(indicate.coordinates$x, "npc") - unit(text.buffer, "inches")
          }
          if (indicate.coordinates$orientation[2]=="bottom") {
                    tmp.y <- unit(indicate.coordinates$y, "npc") + 0.75*unit(text.buffer, "inches")
          } else {
                    tmp.y <- unit(indicate.coordinates$y, "npc") - 0.75*unit(text.buffer, "inches")
          }
          grid.text(x=tmp.x, y=tmp.y, bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.INDICATE], just=indicate.coordinates$orientation)
   }

   if (!is.null(bubble_plot_data.BUBBLE_CENTER_LABEL)) {
            grid.text(x=bubble_plot_data.X, y=bubble_plot_data.Y, bubble_plot_data.BUBBLE_CENTER_LABEL,
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque), font=2), default.units="native")
   }

   if (!is.null(bubble_plot_data.BUBBLE_CENTER_LABEL)) {
            grid.text(x=bubble_plot_data.X, y=bubble_plot_data.Y, bubble_plot_data.BUBBLE_CENTER_LABEL,
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Transparent), font=2), default.units="native")
            grid.text(x=bubble_plot_data.X[bubble_plot_data.SUBSET], y=bubble_plot_data.Y[bubble_plot_data.SUBSET], bubble_plot_data.BUBBLE_CENTER_LABEL[bubble_plot_data.SUBSET],
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X[bubble_plot_data.SUBSET]), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque), font=2), default.units="native")
   }
} else {
   grid.circle(x=bubble_plot_data.X, y=bubble_plot_data.Y, r=unit(bubblesize(bubble_plot_data.SIZE, numstud.range), rep("inches", length(bubble_plot_data.SIZE))),
               gp=gpar(col=rgb(0.2,0.2,0.2), lwd=0.05*bubble_plot_configs.BUBBLE_MIN_MAX[2]/0.12, fill=bubblecolor(unclass(tmp.LEVELS)),
               alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque)), default.units="native")

   if (!is.null(bubble_plot_data.INDICATE)) {
          indicate.coordinates <- indicate.tip(as.numeric(convertX(unit(bubble_plot_data.X[bubble_plot_data.INDICATE], "native"), "npc")),
                                           as.numeric(convertY(unit(bubble_plot_data.Y[bubble_plot_data.INDICATE], "native"), "npc")))
          grid.segments(unit(indicate.coordinates$x, "npc"), unit(indicate.coordinates$y, "npc"),
                        unit(bubble_plot_data.X[bubble_plot_data.INDICATE], "native"), unit(bubble_plot_data.Y[bubble_plot_data.INDICATE], "native"),
                        gp=gpar(lwd=0.5))
          grid.circle(x=bubble_plot_data.X[bubble_plot_data.INDICATE], y=bubble_plot_data.Y[bubble_plot_data.INDICATE],
                  r=unit(c(1.0, 0.4)*bubblesize(bubble_plot_data.SIZE[bubble_plot_data.INDICATE], numstud.range), rep("inches", length(bubble_plot_data.SIZE[bubble_plot_data.INDICATE]))),
                  gp=gpar(lwd=c(0.5, 3.0), fill=bubblecolor(unclass(tmp.LEVELS[bubble_plot_data.INDICATE]))), default.units="native")
          grid.rect(x=unit(indicate.coordinates$x, "npc"), y=unit(indicate.coordinates$y, "npc"),
                         width=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.INDICATE]),
                         height=unit(1.5*text.buffer, "inches")+unit(1.0, "strheight", bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.INDICATE]),
                         gp=gpar(col="grey20", lwd=0.7, fill=rgb(1.0, 0.94, 0.83, 0.6)), just=indicate.coordinates$orientation)
          if (indicate.coordinates$orientation[1]=="left") {
                    tmp.x <- unit(indicate.coordinates$x, "npc") + unit(text.buffer, "inches")
          } else {
                    tmp.x <- unit(indicate.coordinates$x, "npc") - unit(text.buffer, "inches")
          }
          if (indicate.coordinates$orientation[2]=="bottom") {
                    tmp.y <- unit(indicate.coordinates$y, "npc") + 0.75*unit(text.buffer, "inches")
          } else {
                    tmp.y <- unit(indicate.coordinates$y, "npc") - 0.75*unit(text.buffer, "inches")
          }
          grid.text(x=tmp.x, y=tmp.y, bubble_plot_labels.BUBBLE_TITLES[bubble_plot_data.INDICATE], just=indicate.coordinates$orientation)
   }

   if (!is.null(bubble_plot_data.BUBBLE_CENTER_LABEL)) {
            grid.text(x=bubble_plot_data.X, y=bubble_plot_data.Y, bubble_plot_data.BUBBLE_CENTER_LABEL,
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque), font=2), default.units="native")
   }

   if (!is.null(bubble_plot_data.BUBBLE_CENTER_LABEL)) {
            grid.text(x=bubble_plot_data.X, y=bubble_plot_data.Y, bubble_plot_data.BUBBLE_CENTER_LABEL,
                   gp=gpar(col=rgb(0.4,0.4,0.4), cex=bubblesize(bubble_plot_data.SIZE, numstud.range)/as.numeric(convertUnit(stringWidth(bubble_plot_data.BUBBLE_CENTER_LABEL), "inches")),
                   alpha=bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque), font=2), default.units="native")
   }
} ## End SUBSET else statement
} ## End BUBBLE_TIPS else statement

popViewport() ## Pop chart.vp


# Vertical Axis Viewport

pushViewport(vaxis.vp)

grid.rect(x=0.4, y=unit(text.start, "inches"), width=unit(2.0, "strheight", bubble_plot_labels.Y[1]),
          height=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.Y[1]), gp=gpar(fill=format.colors.border, lwd=0.5, col=format.colors.border),
          just=c("center", "bottom"))
grid.rect(x=0.4, y=unit(text.start, "inches")+unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.Y[1]),
          width=unit(2.0, "strheight", bubble_plot_labels.Y[1]),
          height=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.Y[2]), gp=gpar(fill=format.colors.background, lwd=0.5, col=format.colors.border),
          just=c("center", "bottom"))
grid.text(x=0.4, y=unit(text.start, "inches")+unit(text.buffer, "inches")+unit(0.5, "strwidth", bubble_plot_labels.Y[1]), bubble_plot_labels.Y[1],
          gp=gpar(col=format.colors.font[2]), rot=90, just="center")
grid.text(x=0.4, y=unit(text.start+3*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.Y[1])+unit(0.5, "strwidth", bubble_plot_labels.Y[2]),
          bubble_plot_labels.Y[2],
          gp=gpar(col=format.colors.font[1]), rot=90, just="center")
for (i in 2:(length(bubble_plot_configs.BUBBLE_Y_TICKS)-1)) {
     if (is.null(bubble_plot_configs.BUBBLE_Y_TICKS_SIZE)) {
          grid.text(x=0.925, y=bubble_plot_configs.BUBBLE_Y_TICKS[i], bubble_plot_configs.BUBBLE_Y_TICKS[i],
                    gp=gpar(col=format.colors.font[1], cex=0.65), just=c("right", "center"), default.units="native")
     } else {
          grid.text(x=0.925, y=bubble_plot_configs.BUBBLE_Y_TICKS[i], bubble_plot_configs.BUBBLE_Y_TICKS[i],
                    gp=gpar(col=format.colors.font[1], cex=bubble_plot_configs.BUBBLE_Y_TICKS_SIZE[i]), just=c("right", "center"), default.units="native")
}
}

popViewport() ## pop vaxis.vp


# Horizontal Axis Viewport

pushViewport(haxis.vp)

grid.rect(x=unit(text.start, "inches"), y=0.4, width=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.X[1]),
          height=unit(2.0, "strheight", bubble_plot_labels.X[1]), gp=gpar(fill=format.colors.border, lwd=0.5, col=format.colors.border),
          just=c("left", "center"))
grid.rect(x=unit(text.start, "inches")+unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.X[1]), y=0.4,
          width=unit(2*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.X[2]),
          height=unit(2.0, "strheight", bubble_plot_labels.X[1]), gp=gpar(fill=format.colors.background, lwd=0.5, col=format.colors.border),
          just=c("left", "center"))
grid.text(x=unit(text.start, "inches")+unit(text.buffer, "inches")+unit(0.5, "strwidth", bubble_plot_labels.X[1]), y=0.4, bubble_plot_labels.X[1],
          gp=gpar(col=format.colors.font[2]), just="center")
grid.text(x=unit(text.start+3*text.buffer, "inches")+unit(1.0, "strwidth", bubble_plot_labels.X[1])+unit(0.5, "strwidth", bubble_plot_labels.X[2]), y=0.4,
          bubble_plot_labels.X[2],
          gp=gpar(col=format.colors.font[1]), just="center")
for (i in seq_along(bubble_plot_configs.BUBBLE_X_TICKS)) {
     if (is.null(bubble_plot_configs.BUBBLE_X_TICKS_SIZE)) {
          grid.text(x=bubble_plot_configs.BUBBLE_X_TICKS[i], y=0.925, bubble_plot_configs.BUBBLE_X_TICKS[i],
                    gp=gpar(col=format.colors.font[1], cex=0.65), just=c("center", "top"), default.units="native")
     } else {
          grid.text(x=bubble_plot_configs.BUBBLE_X_TICKS[i], y=0.925, bubble_plot_configs.BUBBLE_X_TICKS[i],
                    gp=gpar(col=format.colors.font[1], cex=bubble_plot_configs.BUBBLE_X_TICKS_SIZE[i]), just=c("center", "top"), default.units="native")
}
}

popViewport() ## pop haxis.vp


# Right Legend Viewport

if (bubble_plot_configs.BUBBLE_PLOT_LEGEND) {
pushViewport(right.legend.vp)
grid.rect(width=0.9, height=1, gp=gpar(lwd=0.5, col=format.colors.border))


# Top legend (size)

y.coors <- (0.85+c(0, -0.0375, -0.075, -0.12, -0.175))[1:num.sizes]
grid.text(x=0.5, y=y.coors[1]+0.05, bubble_plot_titles.LEGEND1, gp=gpar(col=format.colors.font[1], fontface=2, cex=1.2))
if (!is.null(bubble_plot_data.SUBSET)) {
  bubble.legend.alpha <- 0.9
  bubble.legend.color <- rep(bubblecolor(unclass(tmp.LEVELS[bubble_plot_data.SUBSET]))[1], length=num.sizes)
} else {
  bubble.legend.alpha <- bubblealpha(length(bubble_plot_data.X), bubble_plot_configs.BUBBLE_SUBSET_ALPHA$Opaque)
  bubble.legend.color <- rep(sort(my.colors)[1], length=num.sizes)
}
for (i in 1:num.sizes){
grid.circle(x=0.25, y=y.coors[i], r=unit(bubblesize(bubble_plot_labels.SIZE[i], numstud.range), "inches"),
            gp=gpar(col="grey14", lwd=0.7, fill=bubble.legend.color[i], alpha=bubble.legend.alpha))
grid.text(x=0.35, y=y.coors[i], paste(bubble_plot_labels.SIZE[i], "Students"), gp=gpar(col=format.colors.font[1], cex=0.9), just="left")
}


# Bottom legend (color of bubbles)

if (!is.null(bubble_plot_data.LEVELS)){
num.levels <- length(unique(bubble_plot_labels.LEVELS))
y.coors <- seq(0.45, by=-.05, length=num.levels)
grid.text(x=0.5, y=y.coors[1]+0.1, bubble_plot_titles.LEGEND2_P1, gp=gpar(col=format.colors.font[1], fontface=2, cex=1.2))
grid.text(x=0.5, y=y.coors[1]+0.065, bubble_plot_titles.LEGEND2_P2, gp=gpar(col=format.colors.font[1], fontface=2, cex=1.2))
for (i in 1:num.levels){
grid.circle(x=0.15, y=y.coors[i], r=unit(0.1, "inches"), gp=gpar(col="grey14", lwd=0.7, fill=bubblecolor(i)))
grid.text(x=0.25, y=y.coors[i], bubble_plot_labels.LEVELS[i], gp=gpar(col=format.colors.font[1], cex=0.9), just="left")
}
}

if (!is.null(bubble_plot_data.LEVELS) & !is.null(bubble_plot_titles.NOTE)) {
	stop('\n\n\t\t Both NOTE and LEVELS can not be used simulateously.  Please choose one and proceed.\n')
}
if (is.null(bubble_plot_data.LEVELS) & !is.null(bubble_plot_titles.NOTE)){
	y.pos <- (nchar(bubble_plot_titles.NOTE)/300) * 0.35  # attempt to be adaptive with NOTE length...
	grid.text(x=0.5, y=y.pos, bubble_plot_titles.NOTE, gp=gpar(col=format.colors.font[1], fontface=3, cex=1.0))
}

popViewport() ## pop right.legend.vp
}


# Title Viewport

if (bubble_plot_configs.BUBBLE_PLOT_TITLE) {
pushViewport(title.vp)

grid.roundrect(width=unit(0.965, "npc"), height=unit(0.75, "npc"), r=unit(0.025, "snpc"), gp=gpar(col=format.colors.border, lwd=1.4))
grid.text(x=0.04, y=0.5, bubble_plot_titles.MAIN, gp=gpar(col=format.colors.font[1], fontface=2, fontfamily="Helvetica-Narrow", cex=3.4), just="left", default.units="native")
grid.text(x=0.96, y=0.65, bubble_plot_titles.SUB1, gp=gpar(col=format.colors.font[1], fontfamily="Helvetica-Narrow", cex=1.7), just="right", default.units="native")
grid.text(x=0.96, y=0.35, bubble_plot_titles.SUB2, gp=gpar(col=format.colors.font[1], fontfamily="Helvetica-Narrow", cex=1.7), just="right", default.units="native")

popViewport() ## pop title.vp
}


# End Viewport Creation

popViewport()


# Turn off device

if (bubble_plot_configs.BUBBLE_PLOT_DEVICE %in% c("PDF", "PNG")) {
    dev.off()
}


# Modify aspects of PDF bubbles

if (bubble_plot_configs.BUBBLE_TIPS) {
	temp_pdf <- readLines(file.path.and.name, encoding="UTF-8")
	temp_pdf <- temp_pdf[-which(temp_pdf=="/C [ 0 1 1 ]")]
	temp_pdf <- temp_pdf[-grep("PLACEHOLDER", temp_pdf)]
	writeLines(temp_pdf, file.path.and.name)
} ## End if BUBBLE_TIPS == TRUE


# Code for pdftk catalog creation

if (bubble_plot_pdftk.CREATE_CATALOG) {

	if (is.na(file.info(".pdftk_tmp")$isdir)){
                dir.create(".pdftk_tmp")
        }
	tmp.page.number <- length(list.files(".pdftk_tmp"))+1
	new.file.path.and.name <- file.path(".pdftk_tmp",
		paste(substr(paste("000000", as.character(tmp.page.number), sep=""), nchar(tmp.page.number), nchar(tmp.page.number)+7), ".pdf", sep=""))
	file.rename(file.path.and.name, new.file.path.and.name)
	if (tmp.page.number == 1) {
cat("InfoKey: Creator
InfoValue: R: A language and environment for statistical computing
InfoKey: Author
InfoValue: Rhode Island Department of Education/The National Center for the Improvement of Educational Assessment
InfoKey: Producer
InfoKey: Rhode Island Department of Education/The National Center for the Improvement of Educational Assessment
InfoKey: Title\n", file=file.path(".pdftk_tmp", ".meta_data.txt"))
cat(paste("InfoValue: ", bubble_plot_titles.SUB1, ": ", bubble_plot_titles.SUB2, "\n", sep=""), file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
cat(paste("BookmarkTitle: ", bubble_plot_configs.BUBBLE_PLOT_NAME, "\n", sep=""), file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
cat("BookmarkLevel: 1\n", file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
cat(paste("BookmarkPageNumber:", tmp.page.number, "\n"), file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
} else {
cat(paste("BookmarkTitle: ", bubble_plot_configs.BUBBLE_PLOT_NAME, "\n", sep=""), file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
cat("BookmarkLevel: 1\n", file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
cat(paste("BookmarkPageNumber:", tmp.page.number, "\n"), file=file.path(".pdftk_tmp", ".meta_data.txt"), append=TRUE)
	}
}
} ## END bubblePlot Function
