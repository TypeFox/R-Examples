##' A color-coded treatment time-line, with overlaid events
##' 
##' This individual-level graphic depicts horizontal time intervals of an
##' ongoing treatment course, color-coded by, e.g., agent.  Categorical events
##' which may occur during treatment, such as assessments of response, are
##' annotatated with color-coded arrows.
##' 
##' @param figlabel A string to be used as a LaTeX figure label
##' @param txs A data frame describing intervals of treatment
##' @param resp A data frame describing treatment response assessments
##' @param bsl A data frame of subject baseline characteristics
##' @param ptid The name of the unique patient identifier (default is "patnum")
##' @param condition An R expression giving a logical condition (a
##' predicate on subject baseline characteristics) used for selecting
##' the subjects to be plotted
##' @param formula A formula of the form \code{trt+resp~time|patnum},
##' interpreted as "plot treatment and response vs time, by patnum"
##' @param followed A list with difftime components \code{from} and
##' \code{to}, giving the minimum and maximum durations of follow-up
##' for patients to be plotted. This is necessary to prevent the
##' dwarfing of treatment courses for patients with short follow-up
##' when plotted alongside those of patients with extended follow-up
##' @param tx.start name of treatment start date column of \code{txs}
##' @param tx.end name of treatment end date column of \code{txs}
##' @param treatment name of treatment column of \code{txs}
##' @param response name of response column of \code{resp}
##' @param time name of time column of \code{resp}
##' @param timeunit The time unit desired for the horizontal axis
##' @param caption The figure caption may be provided explicitly or
##' (if NULL) constructed automatically
##' @param tx.key The plot legend for treatments
##' @param resp.key The plot legend for responses
##' @param cols.rows Trellis layout as \code{c(ncols, nrows)}
##' @param prefix.string Prefix string (including possibly a
##' directory) for cached plot output
##' @param xlim x-axis limits
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param filename Filename pattern for cached plot output
##' @note TODO: further notes
##' @author David C. Norris
##' @seealso TODO: List objects to See Also as \code{\link{help}}
##' @references TODO: Reference our white paper or pending publication
##' @keywords hplot
##' @examples
##' ## TODO: Provide an example
##' ## TODO: Document usage. If necessary, include sample data sets in package:VizOR.
##' @export timeline
## Note that there are 3 different column names used here: 'trt', 'resp' and 'time'.
## It may be worthwhile to require all of these to be provided in a formula, so that
## the data frame columns may be addressed more robustly.  For example, a formula
## such as resp ~ trt / time | patnum would permit all columns to be specified.
timeline <- function(figlabel, txs, resp, bsl, ptid='patnum', condition=TRUE,
           formula=trt+resp~time|patnum,
           followed=list(
             from=as.difftime(0, units='days'),
             to=as.difftime(Inf, units='days')),
           tx.start='start',
           tx.end='end',
           treatment=as.character(formula[[2]][[2]]),
           response=as.character(formula[[2]][[3]]),
           time=as.character(formula[[3]][[2]]),
           timeunit=units(resp[[time]]),
           caption=NULL,
           tx.key=key(txs[[treatment]]),
           resp.key=key(resp[[response]]),
           cols.rows=c(2,5), prefix.string="figs/plot",
           xlim=c(0, as.double(min(followed$to, max(resp[[time]])), units=timeunit)),
           xlab=paste("Time (", timeunit, ")", sep=""), ylab="",
           filename=paste(figlabel, "followed",
             sub("[.]", "_", followed$from), "-",
             sub("[.]", "_", followed$to), "y", sep="")){
    ## These assignments prove necessary to put tx.key and resp.key in local scope:
    tx.key <- tx.key
    resp.key <- resp.key
    ## Capture the condition as quote
    condition <- substitute(condition)
    if(is.null(caption)){
      caption <- paste("Patients with {\\tt ",
                       beNiceToLaTeX(paste(deparse(condition), collapse=" ")),
                       "}", sep="")
      if(!missing(followed)){
        followed.text <- ", and followed for"
        if(followed$from==0){
          followed.text <- paste(followed.text, "up to", followed$to, units(followed$to))
        } else if(followed$to==Inf) {
          followed.text <- paste(followed.text, "over", followed$from, units(followed$from))
        } else {
          unit.from <- units(followed$from)
          unit.to   <- units(followed$to)
          if(units(followed$from)==unit.to)
            unit.from <- ""
          followed.text <- paste(followed.text,
                                 followed$from, " ", unit.from, " to ",
                                 followed$to, " ", unit.to, sep="")
        }
        caption <- paste(caption, followed.text, sep="")
      }
      caption <- paste(caption, ".", sep="")
    }
    ## Calculate for each patient the duration of follow-up
    fu.duration <- aggregate(resp[[time]], by=list(pt=resp[,ptid]), FUN=max)
    fu.duration$x <- as.difftime(as.numeric(fu.duration$x), units=timeunit) # 'aggregate' disrespects difftime class
    ## Select patients with follow-up duration in given range followed.years
    #pts <- subset(fu.duration, followed$from < x & x <= followed$to)$pt
    pts <- fu.duration[followed$from < fu.duration$x & fu.duration$x <= followed$to,"pt"]
    ## Restrict to patients meeting given 'condition'
    pts.cond <- do.call("subset", args=list(x=bsl, subset=condition))[,ptid]
    pts <- intersect(pts, pts.cond)
    txs <- txs[txs[,ptid] %in% pts,]
    resp <- resp[resp[,ptid] %in% pts,]
    dataList <- list(resp=resp, txs=txs)
    trellis.device(pdf, file=paste(prefix.string, "-", filename, "%03d.pdf", sep=''),
                   width=6, height=9, onefile=FALSE, paper="special")
    formula <- NA ~ as.double(time, units=timeunit) | id
    formula[[3]][[3]] <- as.name(ptid)
    plots <-
      xyplot(formula,
             layout=cols.rows,
             as.table=TRUE,
             all.pages.skip=FALSE,
             data=dataList$resp,
             col='black',
             xlim=xlim,
             xlab=xlab, ylab=ylab,
             scales=list(y=list(draw=FALSE)),
             key=list(space="top", cols=2, rep=FALSE,
               points=list(
                 col=unname(unlist(tx.key)),
                 pch=rep(15, length(tx.key))),
               text=list(
                 lab=names(tx.key),
                 cex=0.8),
               points=list(
                 col=unname(unlist(resp.key)),
                 pch=rep(73, length(resp.key))),
               text=list(
                 lab=names(resp.key),
                 cex=0.8)
               ),
             panel=function(x, y, subscripts, panelData=dataList, ...) {
               ## Select the treatments received by this panel's patient
               this.pt <- as.character(panelData$resp[subscripts[1],ptid])
               tx <- panelData$txs[panelData$txs[,ptid] == this.pt,]
               for(i in seq_len(nrow(tx))){
                 ## Determine color for this treatment
                 tx.color <- tx.key[[as.character(tx[i,treatment])]]
                 grid.polygon(y=unit(c(0.45, 0.55, 0.55, 0.45), "npc"),
                              x=rep(
                                c(as.double(tx[i,tx.start], units=timeunit),
                                  as.double(tx[i,tx.end], units=timeunit)),
                                each=2),
                              gp=gpar(
                                fill=tx.color,
                                col=rgb(0,0,0,alpha=0.2)), # a subtle (quite transparent) border
                              default.units="native")
               }
               ## Plot responses last to put them on top layer
               this.resp <- panelData$resp[subscripts,c('time', response)]
               n <- nrow(this.resp)
               grid.polyline(x=unit(rep(as.double(this.resp[[time]], units=timeunit), each=2), "native"),
                             y=unit(rep(c(0.55,0.7), n), "npc"),
                             id=rep(1:n, each=2),
                             arrow=arrow(angle=20, length=unit(0.05, "npc"), ends="first", type="open"),
                             gp=gpar(
                               lwd=unit(1.5,"mm"),
                               col=unlist(resp.key[as.character(this.resp[[response]])]),
                               alpha=0.75) # sometimes, arrows overlap
                             )
             })
    print(plots, newpage=FALSE, prefix="plot")
    ## Close output device to write out all the files,
    ## diverting printed output to keep .tex file clean.
    ignore <- capture.output(dev.off())
    files <- list.files(path="figs",
                        pattern=glob2rx(paste("plot", "-", filename, "*.pdf", sep='')))
    cat("\\renewcommand{\\thefigure}{\\arabic{figure}-\\arabic{subfigure}}")
    cat("\\setcounter{subfigure}{1}")
    for(file in files){
      cat("\\begin{figure}\n\n")
      cat("\\includegraphics{", "figs/", file, "}\n\n", sep="") 
      cat("\\caption{\\label{fig:", figlabel, "}{\\bf ", caption, "}}", sep='')
      cat("\\addtocounter{figure}{-1}\n")
      cat("\\addtocounter{subfigure}{1}\n")
      cat("\\end{figure}\n\n")
      cat("\\clearpage\n\n")
    }
    cat("\\renewcommand{\\thefigure}{\\arabic{figure}}")
    cat("\\addtocounter{figure}{1}\n")
  }
