# ## Copyright (C) 2006  Antonio, Fabio Di Narzo
# ##
# ## This program is free software; you can redistribute it and/or modify
# ## it under the terms of the GNU General Public License as published by
# ## the Free Software Foundation; either version 2, or (at your option)
# ## any later version.
# ##
# ## This program is distributed in the hope that it will be useful, but
# ## WITHOUT ANY WARRANTY; without even the implied warranty of
# ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# ## General Public License for more details.
# ##
# ## A copy of the GNU General Public License is available via WWW at
# ## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# ## writing to the Free Software Foundation, Inc., 59 Temple Place,
# ## Suite 330, Boston, MA  02111-1307  USA.
# 
# #Custom environment S3 class. With methods for subclassing ('extend') and printing
# print.MyEnv <- function(x, ...) {
# 	cat("<",class(x)[1],">\n",sep="")
# }
# 
# Node <- function(parent=NULL, childs=NULL, opts=NULL, tkvar=NULL, ...) {
# 	this <- extend(MyEnv(), "Node", parent=NULL, childs=NULL, tkvar=tkvar, opts=opts)
# 	if(!is.null(parent))
# 		add(parent, this)
# 	if(!is.null(childs))
# 		for(nd in childs)
# 			add(this, nd)
# 	return(this)
# }
# 
# #Add a child to current node
# add.Node <- function(this, ...) {
# 	nodes <- list(...)
# 	if(length(nodes)>1) {
# 		for(node in nodes)
# 			add(this, node)
# 		return(this)
# 	}
# 	this$childs <- c(this$childs, nodes)
# 	nodes[[1]]$parent <- this
# }
# 
# #cmd: geometry manager command (as char. string)
# #opts: geometry manager options (as char. string)
# Frame <- function(cmd="pack", conf=NULL, ...) {
# 	extend(Node(...), "Frame", cmd=cmd, conf=conf)
# }
# 
# LabelFrame <- function(cmd="pack", conf=NULL, text="", ...) {
# 	extend(Frame(cmd="pack", conf=NULL, ...), "LabelFrame", text=text)
# }
# 
# #opts: widget type and options
# Widget <- function(...) {
# 	extend(Node(...), "Widget")
# }
# 
# #Creates all necessary (and properly linked) tcltk widgets, and shows them
# #'Frame' method
# build.Frame <- function(this, ...) {
# 	parent <- this$parent
# 	parent.frm <- parent$tkvar
# 	if(is.null(parent.frm))
# 		this$tkvar <- tktoplevel()
# 	else {
# 		args <- c(list(parent.frm), this$conf)
# 		this$tkvar <- do.call(tkframe, args)
# 	}
#   nds <- this$childs
#   if(!is.null(nds))
#   	for(nd in nds)
#     	build(nd)
# 	if(!is.null(parent.frm)) {
#   	cmd <- list(parent$cmd, this$tkvar)
#   	opts <- parent$opts
#   	do.call(tcl, c(cmd, opts))
# 	}
# }
# 
# build.LabelFrame <- function(this, ...) {
# 	parent <- this$parent
# 	parent.frm <- parent$tkvar
# 	if(is.null(parent.frm))
# 		this$tkvar <- tktoplevel()
# 	else {
# 		args <- c(list("labelframe", parent.frm), this$conf)
# 		this$tkvar <- do.call(tkwidget, args)
# 	}
#   nds <- this$childs
#   if(!is.null(nds))
#   	for(nd in nds)
#     	build(nd)
# 	if(!is.null(parent.frm)) {
#   	cmd <- list(parent$cmd, this$tkvar)
#   	opts <- parent$opts
#   	do.call(tcl, c(cmd, opts))
# 	}
# }
# 
# 
# #Creates all necessary (and properly linked) tcltk widgets, and shows them
# #'Widget' method
# build.Widget <- function(this, ...) {
# 	parent.frm <- this$parent$tkvar
# 	args <- c(list(parent.frm), this$opts)
#   tmp <- do.call(tkwidget, args)
#   this$tkvar <- tmp
#   cmd <- list(this$parent$cmd, this$tkvar)
#   opts <- this$parent$opts
#   do.call(tcl, c(cmd, opts))
# }
# 
# #From a named list of actions, instantiates a Frame with a button for each listed action
# makeButtonsFrame <- function(actions) {
# 	names <- names(actions)
# 	fr <- Frame(opts=list(side="left"))
# 	for(i in 1:length(actions))
# 		add(fr,Widget(opts=list(type="button", text=names[i], command=actions[[i]])))
# 	return(fr)
# }
# 
# #From a root frame, *builds* a complete dialog and shows it
# buildDialog <- function(title, rootNode) {
# 	build(rootNode)
# 	tkwm.title(rootNode$tkvar, title)
# }
# 
# 
# showDialog.aar <- function(x, ...) {
# 	frRoot <- Frame()
# 	vM <- tclVar(1)
# 	vD <- tclVar(1)
# 	vSteps <- tclVar(1)
# 	onFinish <- function() {
# 		res <- aar(x, m=as.numeric(tclObj(vM)), d=as.numeric(tclObj(vD)), steps=as.numeric(tclObj(vSteps)) )
# 		tkdestroy(frRoot$tkvar)
# 		assign("nlarModel", res, .GlobalEnv)
# 	}
# 	onCancel <- function()
# 		tkdestroy(frRoot$tkvar)
# 	frMain <- nlar.struct.Frame(vM, vD, vSteps)
# 	add(frRoot,
# 		frMain,
# 		makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
# 	)
# 	buildDialog(title="additive autoregressive model", frRoot)
# }
# 
# 
# showDialog.lstar <- function(x, ...) {
#   vML  <- tclVar(1)
#   vMH  <- tclVar(1)
#   vThDelay <- tclVar(0)
#   vTh <- tclVar(0)
#   vMaxit <- tclVar(3000)
#   
#   frTop <- Frame(opts=list(side="left"))
#   frLeft <- Frame()
#   add(frLeft,
#       namedSpinbox("low reg. AR order", vML, from=1, to=1000, increment=1, width=4),
#       namedSpinbox("high reg. AR order", vMH, from=1, to=1000, increment=1, width=4)
#       )
#   frRight <- Frame()
#   add(frRight,
#       namedEntry("treshold value", vTh, width=4),
#       namedSpinbox("treshold delay", vThDelay, from=0, to=1000)
#       )
#   frExtra <- Frame()
#   add(frExtra,
#       namedSpinbox("max iterations", vMaxit, from=100, to=1e5, increment=100, width=6)
#       )
#   add(frTop, frLeft, frRight, frExtra)
#   frRoot <- Frame()	
#   
#   onFinish <- function() {
#     mL <- as.numeric(tclObj(vML))
#     mH <- as.numeric(tclObj(vMH))
#     th <- as.numeric(tclObj(vTh))
#     thDelay <- as.numeric(tclObj(vThDelay))
#     maxit <- as.numeric(tclObj(vMaxit))
#     tkdestroy(frRoot$tkvar)
#     res <- lstar(x, mL=mL, mH=mH, th=th, thDelay=thDelay, control=list(maxit=maxit))
#     assign("nlarModel", res, .GlobalEnv)
#   }
#   onCancel <- function()
#     tkdestroy(frRoot$tkvar)
#   
#   bttnFrame <- makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
#   add(frRoot, frTop, bttnFrame)
#   buildDialog(title="LSTAR model", frRoot)
# 
# }
# 
# showDialog.lstar <- function(x, ...) {
#   vML  <- tclVar(1)
#   vMH  <- tclVar(1)
#   vThDelay <- tclVar(0)
#   vTh <- tclVar(0)
#   vMaxit <- tclVar(3000)
#   
#   frTop <- Frame(opts=list(side="left"))
#   frLeft <- Frame()
#   add(frLeft,
#       namedSpinbox("low reg. AR order", vML, from=1, to=1000, increment=1, width=4),
#       namedSpinbox("high reg. AR order", vMH, from=1, to=1000, increment=1, width=4)
#       )
#   frRight <- Frame()
#   add(frRight,
#       namedEntry("treshold value", vTh, width=4),
#       namedSpinbox("treshold delay", vThDelay, from=0, to=1000)
#       )
#   frExtra <- Frame()
#   add(frExtra,
#       namedSpinbox("max iterations", vMaxit, from=100, to=1e5, increment=100, width=6)
#       )
#   add(frTop, frLeft, frRight, frExtra)
#   frRoot <- Frame()	
#   
#   onFinish <- function() {
#     mL <- as.numeric(tclObj(vML))
#     mH <- as.numeric(tclObj(vMH))
#     th <- as.numeric(tclObj(vTh))
#     thDelay <- as.numeric(tclObj(vThDelay))
#     maxit <- as.numeric(tclObj(vMaxit))
#     tkdestroy(frRoot$tkvar)
#     res <- lstar(x, mL=mL, mH=mH, th=th, thDelay=thDelay, control=list(maxit=maxit))
#     assign("nlarModel", res, .GlobalEnv)
#   }
#   onCancel <- function()
#     tkdestroy(frRoot$tkvar)
#   
#   bttnFrame <- makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
#   add(frRoot, frTop, bttnFrame)
#   buildDialog(title="LSTAR model", frRoot)
# 
# }
# 
# 
# 
# showDialog.nnetTs <- function(x, ...) {
# 	frRoot <- Frame()
# 	vM <- tclVar(1)
# 	vD <- tclVar(1)
# 	vSteps <- tclVar(1)
# 	vSize <- tclVar(1)
# 	onFinish <- function() {
# 		res <- nnetTs(x, m=as.numeric(tclObj(vM)), d=as.numeric(tclObj(vD)), steps=as.numeric(tclObj(vSteps)) , size=as.numeric(tclObj(vSize)))
# 		tkdestroy(frRoot$tkvar)
# 		assign("nlarModel", res, .GlobalEnv)
# 	}
# 	onCancel <- function()
# 		tkdestroy(frRoot$tkvar)
# 	frMain <- nlar.struct.Frame(vM, vD, vSteps)
# 	add(frMain,
# 		Widget(opts=list(type="label", text="Hidden units")),
# 		Widget(opts=list(type="spinbox", from=1, to=100, increment=1, textvariable=vSize, width=4))
# 	)
# 	add(frRoot,
# 		frMain,
# 		makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
# 	)
# 	buildDialog(title="linear model", frRoot)
# }
# 
# showDialog.setar <- function(x, ...) {
#   vD <- tclVar(1)
#   vSteps <- tclVar(1)
#   vML  <- tclVar(1)
#   vMH  <- tclVar(1)
#   vThDelay <- tclVar(0)
#   vFixedTh <- tclVar(0)
#   vTh <- tclVar(0)
# 
#   frTop <- Frame(opts=list(side="left"))
#   frLeft <- Frame()
#   add(frLeft,
#       namedSpinbox("low reg. AR order", vML, from=1, to=1000, increment=1, width=4),
#       namedSpinbox("high reg. AR order", vMH, from=1, to=1000, increment=1, width=4))
#   frRight <- Frame()
#   add(frRight,
#       Widget(opts=list(type="checkbutton", text="fixed treshold", variable=vFixedTh)),
#       namedEntry("treshold value", vTh, width=4),
#       namedSpinbox("treshold delay", vThDelay, from=0, to=1000))
#   frExtra <- Frame()
#   add(frExtra,
#       namedSpinbox("time delay", vD, from=1, to=1000),
#       namedSpinbox("forecasting steps", vSteps, from=1, to=1000))
#   add(frTop, frExtra, frLeft, frRight)
#   frRoot <- Frame()
# 
#   onFinish <- function() {
#     d <- as.numeric(tclObj(vD))
#     steps <- as.numeric(tclObj(vSteps))
#     mL <- as.numeric(tclObj(vML))
#     mH <- as.numeric(tclObj(vMH))
#     fixedTh <- as.logical(tclObj(vFixedTh))
#     th <- as.numeric(tclObj(vTh))
#     thDelay <- as.numeric(tclObj(vThDelay))
#     tkdestroy(frRoot$tkvar)
#     if(fixedTh)
#       res <- setar(x, d=d, steps=steps, mL=mL, mH=mH, th=th, thDelay=thDelay)
#     else
#       res <- setar(x, d=d, steps=steps, mL=mL, mH=mH, thDelay=thDelay)
#     assign("nlarModel", res, .GlobalEnv)
#   }
#   onCancel <- function()
#     tkdestroy(frRoot$tkvar)
# 
#   bttnFrame <- makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
#   add(frRoot, frTop, bttnFrame)
#   buildDialog(title="SETAR model", frRoot)
# }
# 
# 
# showDialog.linear <- function(x, ...) {
# 	frRoot <- Frame()
# 	vM <- tclVar(1)
# 	vD <- tclVar(1)
# 	vSteps <- tclVar(1)
# 	onFinish <- function() {
# 		res <- linear(x, m=as.numeric(tclObj(vM)), d=as.numeric(tclObj(vD)), steps=as.numeric(tclObj(vSteps)) )
# 		tkdestroy(frRoot$tkvar)
# 		assign("nlarModel", res, .GlobalEnv)
# 	}
# 	onCancel <- function()
# 		tkdestroy(frRoot$tkvar)
# 	frMain <- nlar.struct.Frame(vM, vD, vSteps)
# 	add(frRoot,
# 		frMain,
# 		makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
# 	)
# 	buildDialog(title="linear model", frRoot)
# }
# 
# showDialog.star <- function(x, ...) {
#   vD <- tclVar(1)
#   vSteps <- tclVar(1)
#   vML  <- tclVar(1)
#   vMH  <- tclVar(1)
#   vThDelay <- tclVar(0)
#   vTh <- tclVar(0)
#   frTop <- Frame(opts=list(side="left"))
#   frLeft <- Frame()
#   add(frLeft,
#     namedSpinbox("low reg. AR order", vML, from=1, to=1000, increment=1, width=4),
#     namedSpinbox("high reg. AR order", vMH, from=1, to=1000, increment=1, width=4)
#   )
#   frRight <- Frame()
#   add(frRight,
#     namedSpinbox("treshold delay", vThDelay, from=0, to=1000)
#   )
#   frExtra <- Frame()
#   add(frExtra,
#           namedSpinbox("time delay", vD, from=1, to=1000),
#           namedSpinbox("forecasting steps", vSteps, from=1, to=1000)
#   )
#   add(frTop, frExtra, frLeft, frRight)
#   frRoot <- Frame()	
# 
#   onFinish <- function() {
#     d <- as.numeric(tclObj(vD))
#     steps <- as.numeric(tclObj(vSteps))
#     mL <- as.numeric(tclObj(vML))
#     mH <- as.numeric(tclObj(vMH))
#     thDelay <- as.numeric(tclObj(vThDelay))
#     tkdestroy(frRoot$tkvar)
#     res <- star(x, d=d, steps=steps, mL=mL, mH=mH, thDelay=thDelay)
#     assign("nlarModel", res, .GlobalEnv)
#   }
#   onCancel <- function()
#     tkdestroy(frRoot$tkvar)
#   bttnFrame <- makeButtonsFrame(list(Finish=onFinish, Cancel=onCancel))
#   add(frRoot, frTop, bttnFrame)
#   buildDialog(title="STAR model", frRoot)
# }


#   if(GUI) {
#     require(tcltk) || stop("tcltk package is required for displaying the GUI")
#     replot <- function(...) {
#       type <- as.character(tclObj(ttype))
#       h.new <- as.numeric(tclObj(th))
#       lag1.new <- as.numeric(tclObj(tlag1))
#       lag2.new <- as.numeric(tclObj(tlag2))
#       if(!is.na(h.new)) h <<- h.new
#       if(!is.na(h.new))
#         h <<- h.new
#       else
#         tclvalue(th) <- as.character(h)
#       if((!is.na(lag1.new))&(!is.na(lag2.new)))
#         lags <<- c(lag1.new, lag2.new)
#       else {
#         tclvalue(tlag1) <- as.character(lags[1])
#         tclvalue(tlag2) <- as.character(lags[2])
#       }
#       X <<- embedd(x, lags=c(-lags, 0))
#       xlab <<- paste("lag",lags[1])
#       ylab <<- paste("lag",lags[2])
#       mod <- sm.regression(X[,1:2], X[,3], h=rep(exp(h),2), display="none")
#       panel[[type]](mod$estimate)
#     }
#     types <- names(panel)
#     ttype <- tclVar(type)
#     tlag1 <- tclVar(lags[1])
#     tlag2 <- tclVar(lags[2])
#     th <- tclVar(h)
#     frLeft <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
#     frRight <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
#     frRoot <- Frame(opts=list(side="left"))
#     add(frRoot, frLeft, frRight)
# 
#     add(frLeft, Widget(opts=list(type="label", text="plot type:")))
#     for(i in 1:length(types))
#       add( frLeft, Widget(opts=list(type="radiobutton", command=replot, 
#                             text=types[i], value=types[i], variable=ttype)) )
# 
#     lagEntry1 <- Widget(opts=list(type="entry", textvariable=tlag1, width=2))
#     lagEntry2 <- Widget(opts=list(type="entry", textvariable=tlag2, width=2))
#     add(frRight, Widget(opts=list(type="label", text="first lag:")), lagEntry1,
#         Widget(opts=list(type="label", text="second lag:")), lagEntry2)
# 
#     h.start <- log(h) - 2
#     h.end <- log(h) + 2
#     hScale <- Widget(opts=list(type="scale", command=replot, from=h.start, to=h.end, 
#                        variable=th, resolution=(h.start-h.end)/100, orient="horiz"))
#     add(frRight, Widget(opts=list(type="label", text="kernel window:")), hScale)
# 
#     buildDialog(title="plot options", frRoot)
#     tkbind(lagEntry1$tkvar, "<Return>",replot)
#     tkbind(lagEntry2$tkvar, "<Return>",replot)
#     return(invisible(NULL))
#   }
# 
#   if(GUI) {
#     require(tcltk) || stop("tcltk package is required for displaying the GUI")
#     replot <- function(...) {
#       type <- as.character(tclObj(vType))
#       h.new <- exp(as.numeric(tclObj(vH)))
#       lag.new <- as.numeric(tclObj(vLag))
#       if(!is.na(h.new))
#         h <<- h.new
#       else
#         tclvalue(th) <- as.character(log(h))
#       if(!is.na(lag.new))
#         lag <<- lag.new
#       else
#         tclvalue(vLag) <- as.character(lag)
#       lags <- c(-lag, 0)
#       X <<- embedd(x, lags=lags)
#       xlab <<- paste("lag",lag)
#       ylab <<- paste("lag",0)
#       panel[[type]]()
#     }
#     types <- c("levels","persp","image","lines","points","regression")
#     vType <- tclVar(type)
#     vLag <- tclVar(lag)
#     vH <- tclVar(log(h))
# 	
#     frLeft <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
#     frRight <- Frame(opts=list(anchor="w"), conf=list(borderwidth=4, relief="raised"))
#     frRoot <- Frame(opts=list(side="left"))
#     add(frRoot, frLeft, frRight)
# 
#     add(frLeft, Widget(opts=list(type="label", text="plot type:")))
#     for(i in 1:length(types))
#       add( frLeft, Widget(opts=list(type="radiobutton", command=replot, 
#                             text=types[i], value=types[i], variable=vType)) )
# 
#     lagEntry <- Widget(opts=list(type="entry", textvariable=vLag, width=2))
#     add(frRight, Widget(opts=list(type="label", text="lag:")),lagEntry)
# 
#     h.start <- log(h) - 2
#     h.end <- log(h) + 2
#     F <- function (phi1, phi2, x_t, s_t) {
#       noRegimes <- dim(phi1)[1]
#       local <- array(0, c(noRegimes, dim(x_t)[1]))
#       local[1, ] <- x_t %*% phi1[1, ]
#       for (i in 2:noRegimes) local[i, ] <- (x_t %*% phi1[i, ]) * G(s_t, gamma = phi2[i - 1, 1], th = phi2[i - 1, 2])
#       return(apply(local, 2, sum))
#     }
#     hScale <- Widget(opts=list(type="scale", command=replot, from=h.start, to=h.end, 
#                        variable=vH, resolution=(h.start-h.end)/100, orient="horiz"))
#     add(frRight, Widget(opts=list(type="label", text="kernel window:")),hScale)
# 
#     buildDialog(title="plot options:", frRoot)
#     tkbind(lagEntry$tkvar, "<Return>", replot)
#     return(invisible(NULL))
#   }
