#### PostScript Goodies für R --- `a la /u/sfs/S/ps.goodies.S
####
####

## hidden in the name space -- FIXME? maybe more useful ?? ---
dev.latex <-
    function(file, DEV, height= 5+ main.space*1.25, width= 9.5,
             main.space = FALSE, lab.space = main.space,
             paper = "special", title = NULL,
             lab = c(10, 10, 7), mgp.lab = c(1.6, 0.7, 0),
             mar = c(4, 4, 0.9, 1.1), ...)
{
  ## Purpose: Setup for 1 LaTeX-includable picture SAVING on white space !
  ##	Calls  ps.do(.) ; par(.)  [ old par in global 'o.p']; USE  ps.end() !
  ## -------------------------------------------------------------------------
  ## Arguments: height & width in INCHES.   (5, 9.5) is 'horizontal look'
  ##		title: to be used in PostScript (-> for gv/ghostview !)
  ##		main.space & lab.space: if T, leave space for 'main' & 'x/ylab'
  ##		lab :  for  par(.);  (10,10,7): use more axis 'labels' ..
  ##		mgp.lab & mar :	 for par(.): these are values for 'lab.space=T'
  ## Note: FIRST fiddle with 'main.sp.', 'lab.sp.'  before 'mgp.lab' and 'mar'!
  ## -------------------------------------------------------------------------
  ## EXAMPLE:for(m in c(T,F)){str(ps.latex("q.ps",main=m));acf(hstart);ps.end()}
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: Sep 94; Sept. 95

  ## Cannot use  missing(.) here, as all arg.s *are* specified
  ## from the calling pdf.latex() etc ..
  frms <- formals()
  lab.def     <- identical(lab,     eval(frms[["lab"]]))
  mar.def     <- identical(mar,     eval(frms[["mar"]]))
  mgp.lab.def <- identical(mgp.lab, eval(frms[["mgp.lab"]]))

  if(!lab.def && !(length(lab)==3 && is.numeric(lab) && all(lab >=0)))
    stop("'lab' must be numeric vector >= 0, of length 3")
  if(!mgp.lab.def && !(length(mgp.lab)==3 && is.numeric(mgp.lab) &&
			    all(mgp.lab >=0) && all(diff(mgp.lab)<=0)))
    stop("'mgp.lab' must be non-increasing numeric vector >= 0, of length 3")
  if(!mar.def && !(length(mar)==4 && is.numeric(mar) && all(mar >=0)))
    stop("'mar' must be non-negative numeric vector of length 4")

  DEV(file=file, height=height, width=width, paper=paper, title = title, ...)
  ##=

  ## Now: just do the proper par(...) calls :
  mar.main.Extra  <- c(0,0, 3.2,0)
  mar.nolab.Minus <- c(1,1, 0.3,0)
  if(main.space && mar.def)
    mar <- mar + mar.main.Extra

  if(!lab.space) {
    mar <- mar - mar.nolab.Minus
    if(main.space)
      warning("'main.space' is TRUE, but 'lab.space' is FALSE ...")
  }
  o.p <- par(mar = mar, mgp= mgp.lab)
  o.p <- c(o.p, par(lab=lab)) # need 2 step for	 bug ?
  ## "frame 0 / GlobalEnv assignment deprecated: u.assign0("o.par.psl", o.p)
  invisible(list(old.par=o.p, new.par= par(c("mar","mgp","lab"))))
}

ps.latex <- function(file, height= 5+ main.space*1.25, width= 9.5,
		     main.space = FALSE, lab.space = main.space,
		     paper = "special", title = NULL,
		     lab = c(10, 10, 7), mgp.lab = c(1.6, 0.7, 0),
		     mar = c(4, 4, 0.9, 1.1), ...)
{
  dev.latex(DEV = ps.do, file=file, height=height, width=width,
            main.space=main.space, lab.space=lab.space, paper=paper,
            title=title, lab=lab, mgp.lab=mgp.lab, mar=mar, ...)
}

pdf.latex <- function(file, height= 5+ main.space*1.25, width= 9.5,
		     main.space = FALSE, lab.space = main.space,
		     paper = "special", title = NULL,
		     lab = c(10, 10, 7), mgp.lab = c(1.6, 0.7, 0),
		     mar = c(4, 4, 0.9, 1.1), ...)
{
  dev.latex(DEV = pdf.do, file=file, height=height, width=width,
            main.space=main.space, lab.space=lab.space, paper=paper,
            title=title, lab=lab, mgp.lab=mgp.lab, mar=mar, ...)
}


ps.do <- local({
    myfile <- NULL
    function(file, width = -1, height = -1,
		  onefile = FALSE, horizontal = FALSE, title = NULL, ...)
{
  ## Purpose: "Ghostview" device driver. --- to be "closed" by ps.end(..) --
  ## -------------------------------------------------------------------------
  ## Arguments: file, width, height : file name and dims in inch; 1 in:=2.54 cm
  ##		onefile = F  <==> Encapsulated PS  (Splus default: T, simple PS)

  ## -- new Args:  combining former   ps.do(.) and  ps.col(.) :

  ##	...  :	passed to ps.options
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, 1992-1995
  ##
  ## --->>>>>> CONSIDER	  'ps.latex'   instead	for pictures !

    myfile <<- file
##   if(length(l... <- list(...))) {
##     ## This does NOT work : pso are the *NEW*, not the *former* ones!
##     oldop <- ps.options()[names(l...)]
##     ps.options(...)
##     on.exit( do.call("ps.options", oldop) ) #- reset ps.options !
##   }

  if(is.null(title))
      title <- paste("R", paste(R.version[c("major", "minor")], collapse = "."),
		     "plot:", file)
  postscript(file = file, width = width, height = height, horizontal=horizontal,
	     onefile = onefile, title = title, print.it = FALSE, ...)
}## ps.do()
})## local(..)

ps.end <- function(call.gv = NULL, command = getOption("eps_view"),
		   debug = getOption("verbose"))
{
    ## Purpose:	 A "ghostview" device driver (almost).
    ## Author: Martin Maechler, Date:  May 26 1992, 15:32
    ## ----------------------------------------------------------------
    ## Arguments:   call.gv: If TRUE,  call ghostview.
    ##	  Default:	  Find out if ghostview already runs on this file,
    ##			  If yes, do not call it again.
    ## MUST be called after ps.do(..) or ps.latex() !
    ## Example:	 ps.end(com = "ghostview --media a4")
    ## ----------------------------------------------------------------
    ## Only if	postscript is running !! --
    if( names(dev.cur()) == "postscript")
	dev.off()
    if(.Platform $ OS.type == "unix") {
        .set.eps_view()
    } else { ## OS.type != "unix"  --- i.e. Windows :
	warning("using ps, ghostview,...is currently not implemented for non-Unix")
	return(FALSE)
    }
    ..ps.file <- environment(ps.do)$myfile
    if (is.null(call.gv)) {
	f <- u.sys(Sys.ps.cmd(), " | grep '", command, "' | grep -v grep")
	if(debug) { cat("ps.end(): f:\n");print(f) }
	call.gv <- length(f) == 0
	if(!call.gv) {
	    ##--- STILL does NOT work
	    ##--- if you work with two different pictures simultaneously
	    for(i in 1:length(f)) { #-- only NOT call if THIS ps.file .. --
		## find command in 'ps' output line (sub/gsub have no 'fixed=TRUE')
		ic <- regexpr(command, f[i], fixed=TRUE)
		## only keep the file name
		fil <- substr(f[i], ic + attr(ic,"match.length") + 1, 1e4)
		cat("ps.end(): fil:",fil,"\n")
		call.gv <- length(fil) < 1 || all(..ps.file != fil)
		if(!call.gv)
		    break #-- don't  call ghostview since it runs this file..
	    }
	}
    } else if(identical(call.gv, FALSE))
	fil <- "<unknown>"
    if (call.gv) {
	fil <- ..ps.file
	u.sys(command, " ", fil, "&", intern=FALSE)
    } else
    cat("\n >> switch to", sub(" .*", '', command),
        "(postscript viewer) window -- updated automagically!\n\n")
    invisible(fil)
}


###---  Using  pdf()  instead of postscript() --- otherwise "same" :

pdf.do <- local({
    myfile <- NULL
    function(file, paper = "default",
                   width = -1, height = -1, onefile = FALSE,
                   title = NULL, version = "1.4", quiet=FALSE, ...)
{
  ## Purpose: "PDF + view" device driver. --- to be "closed" by pdf.end(..) --
  ## -------------------------------------------------------------------------
  ## Arguments: file, width, height : file name and dims in inch; 1 in:=2.54 cm
  ##		onefile = FALSE <==> "Encapsulated"
  ##	...  :	passed to pdf.options
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, April 26, 2007 {built on much older ps.do()}

##   if(length(l... <- list(...))) {
##       ## ps.options also used for pdf -- in some way
##     oldop <- ps.options()[names(l...)]
##     ps.options(...)
##     on.exit( do.call("ps.options", oldop) ) #- reset ps.options !
##   }
    myfile <<- file

  if(missing(paper) && !missing(width) && !missing(height)) {
      if(!quiet)
	  message("as 'width' and 'height' are specified, setting 'paper = \"special\"")
      paper <- "special"
  }

  if(is.null(title))
      title <- paste("R", paste(R.version[c("major", "minor")], collapse = "."),
		     "plot:", file)
  ## default for 'paper' is now 'missing'
  pdf(file = file, version = version, paper = paper,
      width = width, height = height,
      onefile = onefile, title = title, ...)
}## pdf.do()
})## local(..)


pdf.end <- function(call.viewer = NULL, command = getOption("pdfviewer"),
		   debug = getOption("verbose"))
{
    ## Purpose:	 A "ghostview" device driver (almost).
    ## Author: Martin Maechler, Date:  April 26, 2007
    ## ----------------------------------------------------------------
    ## Arguments:   call.viewer: If TRUE,  call ghostview.
    ##	  Default:	  Find out if ghostview already runs on this file,
    ##			  If yes, do not call it again.
    ## MUST be called after pdf.do(..) or pdf.latex() !
    ## Example:	 pdf.end(com = "acroread")
    ## ----------------------------------------------------------------
    ## Only if	postscript is running !! --
    if( names(dev.cur()) == "pdf")
	dev.off()
    if(.Platform $ OS.type != "unix") {
	warning("using ps (process status) is currently not implemented for non-Unix")
	return(FALSE)
    }
    ..pdf.file <- environment(pdf.do)$myfile
    if (is.null(call.viewer)) {
        cmd <- basename(command)
	f <- u.sys(Sys.ps.cmd(), " | grep '", cmd, "' | grep -v grep")
	if(debug) { cat("pdf.end(): f:\n");print(f) }
	call.viewer <- length(f) == 0
	if(!call.viewer) {
	    ##--- STILL does NOT work
	    ##--- if you work with two different pictures simultaneously
	    for(i in 1:length(f)) { #-- only NOT call if THIS pdf.file .. --
		## find command in 'ps' output line (sub/gsub have no 'fixed=TRUE')
		ic <- regexpr(cmd, f[i], fixed=TRUE)
		## only keep the file name
		fil <- substr(f[i], ic + attr(ic,"match.length") + 1, 1e4)
		cat("pdf.end(): fil:",fil,"\n")
		call.viewer <- length(fil) < 1 || all(..pdf.file != fil)
		if(!call.viewer)
		    break #-- don't  call ghostview since it runs this file..
	    }
	}
    } else if(identical(call.viewer, FALSE))
	fil <- "<unknown>"
    if (call.viewer) {
	fil <- ..pdf.file
	u.sys(command, " ", fil, "&", intern=FALSE)
    } else {
	msg <- if(length(grep("acroread", command)))
	    " acroread -- and refresh via C-w M-f 1 !"
	else "	PDF viewer window and maybe refresh!"
	cat("\n >> switch to", msg,"\n\n")
    }
    invisible(fil)
}

## Alain Hauser <alain@huschhus.ch> --> ../man/cairoSwd.Rd
cairoSwd <- function(name, width, height, ...)
  cairo_pdf(filename = paste(name, "pdf", sep = "."),
            width = width, height = height)
