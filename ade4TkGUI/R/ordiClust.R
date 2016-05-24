################################
# tcl/Tk GUI for the ordiClust method
################################
"ordiClust" <- function(datatab=NULL, hscalef=1.2, vscalef=1.2, maxgr=20) {
  
  if (.Platform$OS.type == "windows") {
    hscale <- hscalef + 0.2
    vscale <- vscalef + 0.2
  } else {
    hscale <- hscalef
    vscale <- vscalef
  }
  
  tt <- tktoplevel()
  tkwm.title(tt,"OrdiClust")
  
  
  tkwm.title(tt,"OrdiClust")
  ## Scrollbar
  ##	
  #  scr <- tkscrollbar(tt, repeatinterval=5, command=function(...)tkyview(tt,...))
  ##  tkconfigure(tt, yscrollcommand=function(...)tkset(scr,...))
  #  tkpack(tt,side="left",fill="both",expand=TRUE)
  #  tkpack(scr,side="right",fill="y")
  #  
  
  #
  # Local variables
  #
  maxngr <- maxngr1 <- maxgr
  
  ordvar <- tclVar(2)
  distvar <- tclVar(1)
  clustvar <- tclVar(3)
  dfvar <- tclVar()
  dudivar <- tclVar()
  
  naxvar <- tclVar(2)
  xaxvar <- tclVar(1)
  yaxvar <- tclVar(2)
  
  nlevvar <- tclVar(1)
  hlevvar <- tclVar(0)
  colorvar <- tclVar(1)
  mgrvar <- tclVar(1)
  
  plotrow <- 0
  plotcol <- 0
  ploteig <- 0
  plotclust <- 0
  plotlev <- 0
  plotclass <- 0
  plotcurve <- 0
  
  distml <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  clustml <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
  
  hlev <- 0
  hlev1 <- 0
  nlev <- 0
  nax <- 2
  r1 <- 0
  h1 <- 0
  betrat <- 0
  
  ordiClust.dudi <- NULL
  ordiClust.factor <- NULL
  #
  # tkrplot plot function
  #
  "fplot" <- function() {
    #
    # Ellipses / convex hulls
    #
    if (plotclass) {
      if (!is.null(ordiClust.factor)) {
        xaxs <- tclvalue(xaxvar)
        yaxs <- tclvalue(yaxvar)
        nlevloc <- tclvalue(nlevvar)
        if (nlevloc == 1) fact1 <- as.factor(rep(1,nrow(ordiClust.dudi$li)))
        else fact1 <- ordiClust.factor
        if (xaxs != "") {
          xax <- as.numeric(xaxs)
        } else return()
        if (yaxs != "") {
          yax <- as.numeric(yaxs)
        } else  return()
        if ((xax >= 1) && (xax <= nax) && (yax >= 1) && (yax <= nax))
          if (tclvalue(mgrvar) == 1) {
            if (as.logical(tclObj(colorvar))) s.class(dfxy=ordiClust.dudi$li, fac=fact1, xax=xax, yax=yax, col=rainbow(nlevels(fact1)))
            else s.class(dfxy=ordiClust.dudi$li, fac=fact1, xax=xax, yax=yax)
          } else {
            if (as.logical(tclObj(colorvar))) s.class(dfxy=ordiClust.dudi$li, fac=fact1, chullSize=1, ellipseSize=0, starSize=0, xax=xax, yax=yax, col=rainbow(nlevels(fact1)))
            else s.class(dfxy=ordiClust.dudi$li, fac=fact1, chullSize=1, ellipseSize=0, starSize=0, xax=xax, yax=yax)
          }
      }
    } else if (ploteig) {
      if (!is.null(ordiClust.dudi)) {
        plotEig(ordiClust.dudi$eig)
      }
    } else if (plotrow) {
      if (!is.null(ordiClust.dudi)) {
        xaxs <- tclvalue(xaxvar)
        yaxs <- tclvalue(yaxvar)
        if (xaxs != "") {
          xax <- as.numeric(xaxs)
        } else return()
        if (yaxs != "") {
          yax <- as.numeric(yaxs)
        } else  return()
        if ((xax >= 1) && (xax <= nax) && (yax >= 1) && (yax <= nax))
          s.label(dfxy=ordiClust.dudi$li, xax=xax, yax=yax, plabels.boxes.draw=FALSE)
      }
    } else if (plotcol) {
      if (!is.null(ordiClust.dudi)) {
        xaxs <- tclvalue(xaxvar)
        yaxs <- tclvalue(yaxvar)
        if (xaxs != "") {
          xax <- as.numeric(xaxs)
        } else return()
        if (yaxs != "") {
          yax <- as.numeric(yaxs)
        } else  return()
        if ((xax >= 1) && (xax <= nax) && (yax >= 1) && (yax <= nax)) 
          s.label(dfxy=ordiClust.dudi$co, xax=xax, yax=yax, plabels.boxes.draw=FALSE)
      }
    } else if (plotclust) {
      if (class(h1) == "hclust") {
        plot(h1, hang=-1)
        if (plotlev) {
          if (as.logical(tclObj(colorvar))) abline(h=hlev, col="red", lwd=2)
          else abline(h=hlev, lwd=2)
        }
      }
    } else if (plotcurve) {
      if (!is.null(ordiClust.dudi))
        if (exists("betrat")) {
          plot(betrat, main="Between-groups / Total inertia ratio", xlab="Number of groups", ylab="Inertia ratio")
          lines(betrat)
          nlev <- as.numeric(tclvalue(nlevvar))
          if (as.logical(tclObj(colorvar))) {
            abline(h=betrat[nlev], col="red", lwd=2)
            abline(v=nlev, col="red", lwd=2)
          } else  {
            abline(h=betrat[nlev], lwd=2)
            abline(v=nlev, lwd=2)
          }
        }
    } else {
      plot(0, type="n")
    }
  }
  
  ################################
  # Function to save a graphic in a file
  ################################
  "outgraphOrdiClust" <- function()
  {
    #
    # Main dialog window with title and frames
    #
    tf <- tktoplevel()
    tkwm.title(tf,"Save graphic")
    #
    # Frames
    #
    frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
    frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
    frame3 <- tkframe(tf, relief="groove", borderwidth=2)	
    devframe <- tkframe(frame2, relief="groove", borderwidth=2)
    #
    # Tcl/Tk variables
    #
    done <- tclVar(0)
    formatvar <- tclVar(1)
    widthvar <- tclVar(6)
    heightvar <- tclVar(6)
    #
    # Save function
    #
    "savefic" <- function(formatvar, widthvar, heightvar)
    {
      outform <- tclvalue(formatvar)
      width <- as.numeric(tclvalue(widthvar))
      height <- as.numeric(tclvalue(heightvar))
      # odev <- dev.cur()
      if (outform == 1) { # postcript
        filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.ps", defaultextension=".ps",
                                           title="Save graph...", filetypes="{PostScript {.ps .eps}} {{All Files} {*.*}}"))
        if (filename != "") {
          postscript(file=filename, width=width, height=height)
        }
      } else if (outform == 2) { # pdf
        filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.pdf", defaultextension=".pdf",
                                           title="Save graph...", filetypes="{PDF {.pdf}} {{All Files} {*.*}}"))
        if (filename != "") {
          pdf(file=filename, width=width, height=height)
        }
      } else if (outform == 3) { # pictex
        filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.tex", defaultextension=".tex",
                                           title="Save graph...", filetypes="{PicTeX {.tex}} {{All Files} {*.*}}"))
        if (filename != "") {
          pictex(file=filename, width=width, height=height)
        }
      } else if (outform == 4) { # xfig
        filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.fig", defaultextension=".fig",
                                           title="Save graph...", filetypes="{XFig {.fig}} {{All Files} {*.*}}"))
        if (filename != "") {
          xfig(file=filename, width=width, height=height)
        }
      } else if (outform == 5) { # png
        filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.png", defaultextension=".png",
                                           title="Save graph...", filetypes="{PNG {.png}} {{All Files} {*.*}}"))
        if (filename != "") {
          png(filename=filename, width=width, height=height)
        }
      } else if (outform == 6) { # jpeg
        filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.jpeg", defaultextension=".jpeg",
                                           title="Save graph...", filetypes="{JPEG {.jpeg .jpg}} {{All Files} {*.*}}"))
        if (filename != "") {
          jpeg(filename=filename, width=width, height=height)
        }
      }
      # ndev <- dev.cur()
      # dev.set(odev)
      # dev.copy(which=ndev)
      # dev.off()
      fplot()
      dev.off()
      tkdestroy(tf)
    }
    #
    # Frames setup
    #
    tkgrid(tklabel(tf,text="Save current graphic", font="Times 18"), columnspan=2)
    
    tkgrid(tklabel(frame2,text="Output format : "), sticky="n")
    tkgrid(tkradiobutton(frame2, text="postscript", value=1, variable=formatvar), sticky="w")
    tkgrid(tkradiobutton(frame2, text="pdf", value=2, variable=formatvar), sticky="w")
    tkgrid(tkradiobutton(frame2, text="pictex", value=3, variable=formatvar), sticky="w")
    tkgrid(tkradiobutton(frame2, text="xfig", value=4, variable=formatvar), sticky="w")
    tkgrid(tkradiobutton(frame2, text="png", value=5, variable=formatvar), sticky="w")
    tkgrid(tkradiobutton(frame2, text="jpeg", value=6, variable=formatvar), sticky="w")
    tkgrid(frame2, rowspan=2, sticky="n")
    
    tkgrid(tklabel(frame3,text="Output size : "))
    width.entry <- tkentry(frame3, textvariable=widthvar, width=10)
    height.entry <- tkentry(frame3, textvariable=heightvar, width=10)
    tkgrid(tklabel(frame3,text="Width : "), width.entry)
    tkgrid(tklabel(frame3,text="Height : "), height.entry)
    tkgrid(frame3, column=1, row=1, sticky="n")
    
    save.but <- tkbutton(frame1, text="Save", command=function() savefic(formatvar, widthvar, heightvar))
    cancel.but <- tkbutton(frame1, text="Dismiss", command=function() tkdestroy(tf))
    tkgrid(save.but, cancel.but)
    tkgrid(frame1, column=1, row=2, sticky="n")
    
    tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tf, "<KeyPress-Return>", function() savefic(formatvar, widthvar, heightvar))
    tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
    tkwait.variable(done)
    if(tclvalue(done) == "2") return(0)
    tkdestroy(tf)
  }
  
  
  #
  # tkrplot init
  #
  plotFrame <- tkframe(tt, relief="groove", borderwidth=2)
  img <- tkrplot(plotFrame, fplot, hscale=hscale, vscale=vscale)
  
  #
  # Ordination method
  #
  "doOrd" <- function() {
    
    meth  <- as.numeric(tclvalue(ordvar))
    df1  <- tclvalue(dfvar)
    if (df1 == "") return()
    if (!is.null(ordiClust.dudi)) maxngr <<- min(maxngr1, nrow(ordiClust.dudi$li))
    tab <- eval(parse(text=df1), envir=env_ade4tkgui)
    #
    # cPCA
    #
    if (meth == 1) {
      nax <<- r1 <<- dudi.pca(as.data.frame(tab), scale=FALSE, scannf=FALSE, nf=2)$rank
      ordiClust.dudi <<- dudi.pca(as.data.frame(tab), scale=FALSE, scannf=FALSE, nf=r1)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      #
      # nPCA
      #
    } else if (meth == 2) {
      nax <<- r1 <<- dudi.pca(as.data.frame(tab), scannf=FALSE, nf=2)$rank
      ordiClust.dudi <<- dudi.pca(as.data.frame(tab), scannf=FALSE, nf=r1)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      #
      # COA
      #
    } else if (meth == 3) {
      nax <<- r1 <<- dudi.coa(as.data.frame(tab), scannf=FALSE, nf=2)$rank
      ordiClust.dudi <<- dudi.coa(as.data.frame(tab), scannf=FALSE, nf=r1)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      #
      # MCA
      #
    } else if (meth == 4) {
      nax <<- r1 <<- dudi.acm(as.data.frame(tab), scannf=FALSE, nf=2)$rank
      ordiClust.dudi <<- dudi.acm(as.data.frame(tab), scannf=FALSE, nf=r1)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
    }
    plotclass <<- 0
    plotrow <<- 0
    plotcol <<- 0
    ploteig <<- 1
    plotclust <<- 0
    plotcurve <<- 0
  }
  #
  # Set the number of axes on which cluster analysis will be computed
  # and plot row scores
  #
  "doChooseAxes" <- function() {

    meth  <- as.numeric(tclvalue(ordvar))
    df1  <- tclvalue(dfvar)
    if (df1 == "") return()
    naxloc  <- as.numeric(tclvalue(naxvar))
    tab <- eval(parse(text=df1), envir=env_ade4tkgui)
    #
    # cPCA
    #
    if (meth == 1) {
      nax <<- naxloc
      ordiClust.dudi <<- dudi.pca(as.data.frame(tab), scale=FALSE, scannf=FALSE, nf=nax)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      #
      # nPCA
      #
    } else if (meth == 2) {
      nax <<- naxloc
      ordiClust.dudi <<- dudi.pca(as.data.frame(tab), scannf=FALSE, nf=nax)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      #
      # COA
      #
    } else if (meth == 3) {
      nax <<- naxloc
      ordiClust.dudi <<- dudi.coa(as.data.frame(tab), scannf=FALSE, nf=nax)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      #
      # MCA
      #
    } else if (meth == 4) {
      nax <<- naxloc
      ordiClust.dudi <<- dudi.acm(as.data.frame(tab), scannf=FALSE, nf=nax)
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
    }
    plotclass <<- 0
    plotrow <<- 1
    plotcol <<- 0
    ploteig <<- 0
    plotclust <<- 0
    plotcurve <<- 0
  }
  #
  # Plot row scores
  #
  "plotr" <- function() {
    plotclass <<- 0
    plotrow <<- 1
    plotcol <<- 0
    ploteig <<- 0
    plotclust <<- 0
    plotcurve <<- 0
  }
  #
  # Plot column scores
  #
  "plotc" <- function() {
    plotclass <<- 0
    plotrow <<- 0
    plotcol <<- 1
    ploteig <<- 0
    plotclust <<- 0
    plotcurve <<- 0
  }
  
  #
  # Clustering method
  #
  "doClust" <- function() {
    nlev <<- as.numeric(tclvalue(nlevvar))
    distm  <- as.numeric(tclvalue(distvar))
    clustm  <- as.numeric(tclvalue(clustvar))
    df1  <- tclvalue(dfvar)
    if (df1 == "") if (is.null(ordiClust.dudi)) return()
    if (distm < 7) h1 <<- hclust(dist(ordiClust.dudi$li, method=distml[distm]), method=clustml[clustm])
    else h1 <<- hclust(dist.dudi(ordiClust.dudi), method=clustml[clustm])
    if (nlev > 1 && nlev <= nrow(ordiClust.dudi$li)) {
      ordiClust.factor <<- as.factor(cutree(h1, k=nlev))
      assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
      doCutk()
      doCuth()
    }
    if (plotcurve) doDraw()
    #		plotclass <<- 0
    plotrow <<- 0
    plotcol <<- 0
    ploteig <<- 0
    #		plotclust <<- 1
    #		plotcurve <<- 0
  }
  #
  # Cut tree (given a number of groups)
  #
  "doCutk" <- function() {
    nlev <<- as.numeric(tclvalue(nlevvar))
    if (class(h1) == "hclust")
      if (is.na(nlev)) {
        return()
      } else {
        if (nlev > 1 && nlev <= nrow(ordiClust.dudi$li)) {
          ordiClust.factor <<- as.factor(cutree(h1, k=nlev))
          assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
        }
        hh <- cutree(h1, h=h1$height)
        for (i in 1:ncol(hh))
          if (nlevels(ordiClust.factor) == nlevels(as.factor(hh[,i])))
            hlev1 <<- (as.numeric(dimnames(hh)[[2]][i])+as.numeric(dimnames(hh)[[2]][i+1]))/2
        hlev <<- hlev1
        tclvalue(hlevvar) <- hlev1
        if (nlev > 1 && nlev <= nrow(ordiClust.dudi$li)) {
          bet1 <- bca(ordiClust.dudi, as.factor(ordiClust.factor), scannf=FALSE)
          mc1 <- randtest(bet1)
          tkconfigure(bwg.label, text=paste(format(bet1$ratio*100, dig=2),"%",sep=""))
          tkconfigure(proba.label, text=format(mc1$pvalue, dig=3))
        }
      }
    if (plotclust) plotlev <<- 1
  }
  #
  # Cut tree (given a height in the tree)
  #
  "doCuth" <- function() {
    hlevloc <- as.numeric(tclvalue(hlevvar))
    if (class(h1) == "hclust")
      if (is.na(hlevloc)) {
        return()
      } else {
        ordiClust.factor <<- as.factor(cutree(h1, h=hlevloc))
        assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
        hlev <<- hlevloc
        tclvalue(nlevvar) <- nlevels(ordiClust.factor)
        nlev <<- as.numeric(tclvalue(nlevvar))
        if (nlev > 1 && nlev <= nrow(ordiClust.dudi$li)) {
          bet1 <- bca(ordiClust.dudi,as.factor(ordiClust.factor),scannf=FALSE)
          mc1 <- randtest(bet1)
          tkconfigure(bwg.label, text=paste(format(bet1$ratio*100, dig=2),"%",sep=""))
          tkconfigure(proba.label, text=format(mc1$pvalue, dig=3))
        }
      }
    if (plotclust) plotlev <<- 1
  }
  #
  # Do Monte-Carlo permutation test for between-groups analysis
  #
  "doTest" <- function() {
    nlev <- as.numeric(tclvalue(nlevvar))
    if (is.na(nlev)) {
      return()
    } else {
      if (nlev > 1 && nlev <= nrow(ordiClust.dudi$li)) {
        ordiClust.factor <<- as.factor(cutree(h1, k=nlev))
        assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
      }
      bet1 <- bca(ordiClust.dudi,as.factor(ordiClust.factor),scannf=FALSE)
      mc1 <- randtest(bet1)
      tkconfigure(bwg.label, text=paste(format(bet1$ratio*100, dig=2),"%",sep=""))
      tkconfigure(proba.label, text=format(mc1$pvalue, dig=3))
    }
    if (plotclust) plotlev <<- 1
  }
  #
  # Plot ellipses or convex hulls
  #
  "doClass" <- function() {
    plotrow <<- 0
    plotcol <<- 0
    ploteig <<- 0
    plotclust <<- 0
    plotlev <<- 0
    plotclass <<- 1
    plotcurve <<- 0
  }
  #
  # Draw the curve of Bet/Tot inertia ratio
  #
  "doDraw" <- function() {
    if (class(h1) == "hclust") {
      hhloc <- cutree(h1, h=h1$height[(length(h1$height)-(maxngr-2)):length(h1$height)])
      nchh <- ncol(hhloc)
      betrat <<- vector("numeric", nchh)
      for (i in 1:(nchh-1)) {
        f1loc <- as.factor(hhloc[,i])
        betrat[nchh-i+1] <<- bca(ordiClust.dudi, f1loc, scannf=FALSE)$ratio
      }
      plotrow <<- 0
      plotcol <<- 0
      ploteig <<- 0
      plotclust <<- 0
      plotlev <<- 0
      plotclass <<- 0
      plotcurve <<- 1
    }
  }
  #
  # Title
  #
  TFrame <- tkframe(tt, relief="groove")
  labh <- tklabel(TFrame, bitmap="questhead")
  tkgrid(tklabel(TFrame,text="OrdiClust", font="Times 18", foreground="red"), labh)
  tkgrid(TFrame, columnspan=2)
  tkbind(labh, "<Button-1>", function() print(help("ordiClust")))
  
  compFrame <- tkframe(tt, relief="groove")
  #
  # Input dataframe / dudi
  #	
  IOFrame <- tkframe(compFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(IOFrame,text="- Input (dataframe OR dudi) -", foreground="blue"), columnspan=5)
  df.entry <- tkentry(IOFrame, textvariable=dfvar)
  dfnr.label <- tklabel(IOFrame, width=4)
  dfnc.label <- tklabel(IOFrame, width=4)
  choosedf.but <- tkbutton(IOFrame, text="Set", command=function() {choosedf(df.entry, dfnr.label, dfnc.label); tclvalue(dudivar) <- ""})
  dudi.entry <- tkentry(IOFrame, textvariable=dudivar)
  choosedudi.but <- tkbutton(IOFrame, text="Set", command=function() {
    choosedudi(dudi.entry)
    if (tclvalue(dudivar) != "") {
      ordiClust.dudi <<- eval(parse(text=tclvalue(dudivar)))
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      tclvalue(naxvar) <<- as.character(ordiClust.dudi$nf)
      nax <<- ordiClust.dudi$nf
      plotrow <<- 0
      plotcol <<- 0
      ploteig <<- 1
      plotclust <<- 0
      plotlev <<- 0
      plotclass <<- 0
      plotcurve <<- 0
      tkrreplot(img)
    }
    tclvalue(dfvar) <- ""
  })
  
  tkgrid(tklabel(IOFrame,text="Input data frame : "), df.entry, choosedf.but, dfnr.label, dfnc.label, sticky="w")
  tkgrid(tklabel(IOFrame,text="Input dudi : "), dudi.entry, choosedudi.but, sticky="w")
  #
  # Ordination method
  #	
  ordFrame <- tkframe(compFrame, relief="groove", borderwidth=2)
  
  ordMethFrame <- tkframe(ordFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(ordMethFrame, text="- Ordination method -", foreground="blue"))
  tkgrid(tkradiobutton(ordMethFrame, text="cPCA", value=1, variable=ordvar))
  tkgrid(tkradiobutton(ordMethFrame, text="nPCA", value=2, variable=ordvar))
  tkgrid(tkradiobutton(ordMethFrame, text="COA", value=3, variable=ordvar))
  tkgrid(tkradiobutton(ordMethFrame, text="MCA", value=4, variable=ordvar))
  ord.but <- tkbutton(ordMethFrame, text="Submit", default="active", command=function() {doOrd(); tkrreplot(img)})
  tkgrid(ord.but)
  nax.but <- tkbutton(ordMethFrame, text="Set", default="active", command=function() {doChooseAxes(); tkrreplot(img)})
  nax.entry <- tkentry(ordMethFrame, textvariable=naxvar, width=4)
  tkgrid(tklabel(ordMethFrame,text="Number of axes : "), nax.entry, nax.but, sticky="w")
  #
  # Ordination graph
  #	
  ordGrFrame <- tkframe(ordFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(ordGrFrame, text="- Ordination graph -", foreground="blue"), columnspan=2)
  xax.entry <- tkentry(ordGrFrame, textvariable=xaxvar, width=4)
  tkgrid(tklabel(ordGrFrame,text="X-axis : "), xax.entry, sticky="w")
  yax.entry <- tkentry(ordGrFrame, textvariable=yaxvar, width=4)
  tkgrid(tklabel(ordGrFrame,text="Y-axis : "), yax.entry, sticky="w")
  plotrow.but <- tkbutton(ordGrFrame, text="Plot rows", default="active", command=function() {plotr(); tkrreplot(img)})
  plotcol.but <- tkbutton(ordGrFrame, text="Plot columns", default="active", command=function() {plotc(); tkrreplot(img)})
  tkgrid(plotrow.but, columnspan=2)
  tkgrid(plotcol.but, columnspan=2)
  
  tkgrid(tklabel(ordFrame,text="Ordination method :", font="Times 18", foreground="red"), columnspan=2)
  
  tkgrid(ordMethFrame, ordGrFrame)
  
  #
  # Distance matrix
  #	
  clustFrame <- tkframe(compFrame, relief="groove", borderwidth=2)
  
  clustDistFrame <- tkframe(clustFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(clustDistFrame, text="- Distance -", foreground="blue"))
  tkgrid(tkradiobutton(clustDistFrame, text="Euclidean", value=1, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustDistFrame, text="maximum", value=2, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustDistFrame, text="manhattan", value=3, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustDistFrame, text="canberra", value=4, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustDistFrame, text="binary", value=5, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustDistFrame, text="minkowski", value=6, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustDistFrame, text="dudi (ade4)", value=7, variable=distvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  #
  # Cluster method
  #	
  clustMethFrame <- tkframe(clustFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(clustMethFrame, text="- Cluster method -", foreground="blue"))
  tkgrid(tkradiobutton(clustMethFrame, text="ward", value=1, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustMethFrame, text="single", value=2, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustMethFrame, text="complete", value=3, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustMethFrame, text="average", value=4, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustMethFrame, text="mcquitty", value=5, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustMethFrame, text="median", value=6, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  tkgrid(tkradiobutton(clustMethFrame, text="centroid", value=7, variable=clustvar, command=function() {doClust(); tkrreplot(img)}), sticky="w")
  clust.but <- tkbutton(clustMethFrame, text="Submit", default="active", command=function() {
    plotclust <<- 1;
    plotcurve <<- 0;
    plotclass <<- 0;
    doClust();
    tkrreplot(img)})
  tkgrid(clust.but)
  #
  # Cut level
  #	
  cutFrame <- tkframe(compFrame, relief="groove", borderwidth=2)
  tkgrid(tklabel(cutFrame,text="Number of groups", font="Times 18", foreground="red"), columnspan=3)
  cutk.but <- tkbutton(cutFrame, text="Cut tree", default="active", command=function() {doCutk(); tkrreplot(img)})
  nlev.entry <- tkentry(cutFrame, textvariable=nlevvar, width=6)
  tkgrid(tklabel(cutFrame,text="Number of groups : "), nlev.entry, cutk.but, sticky="w")
  cuth.but <- tkbutton(cutFrame, text="Cut tree", default="active", command=function() {doCuth(); tkrreplot(img)})
  hlev.entry <- tkentry(cutFrame, textvariable=hlevvar, width=6)
  tkgrid(tklabel(cutFrame,text="Level height : "), hlev.entry, cuth.but, sticky="w")
  bwg.but <- tkbutton(cutFrame, text="Draw curve", default="active", command=function() {doDraw(); tkrreplot(img)})
  bwg.label <- tklabel(cutFrame, width=4)
  tkgrid(tklabel(cutFrame,text="Inertia ratio : "), bwg.label, bwg.but, sticky="w")
  proba.label <- tklabel(cutFrame, width=4)
  tkgrid(tklabel(cutFrame,text="BGA MCTest p-value : "), proba.label, sticky="w")
  
  tkgrid(tklabel(clustFrame,text="Cluster analysis :", font="Times 18", foreground="red"), columnspan=2)
  
  tkgrid(clustDistFrame, clustMethFrame)
  
  tkgrid(IOFrame)
  tkgrid(ordFrame)
  tkgrid(clustFrame)
  tkgrid(cutFrame)
  
  #
  # tkscale re-plot function
  #
  f <- function(...) {
    hlev <<- as.numeric(tclvalue(hlevvar))
    nlev <<- as.numeric(tclvalue(nlevvar))
    if (!is.null(ordiClust.dudi))
      if (nlev >= 1 && nlev <= nrow(ordiClust.dudi$li)) {
        if (class(h1) == "hclust") {
          ordiClust.factor <<- as.factor(cutree(h1, k=nlev))
          assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
        }
        else return(0)
        if (nlev >= 2 && nlev <= nrow(ordiClust.dudi$li)) {
          bet1 <- bca(ordiClust.dudi,as.factor(ordiClust.factor),scannf=FALSE)
          mc1 <- randtest(bet1)
          tkconfigure(bwg.label, text=paste(format(bet1$ratio*100, dig=2),"%",sep=""))
          tkconfigure(proba.label, text=format(mc1$pvalue, dig=3))
        } else return(0)
        if (plotclust) {
          plotlev <<- 1
          hh <- cutree(h1, h=h1$height)
          for (i in 1:ncol(hh))
            if (nlevels(ordiClust.factor) == nlevels(as.factor(hh[,i])))
              hlev1 <<- (as.numeric(dimnames(hh)[[2]][i])+as.numeric(dimnames(hh)[[2]][i+1]))/2
          hlev <<- hlev1
          tclvalue(hlevvar) <- hlev1				
        }
        tkrreplot(img)
        return(1)
      } else return(0)
  }
  
  #
  # scale
  #
  s <- tkscale(plotFrame, command=f, from=0, to=maxngr, variable=nlevvar, showvalue=TRUE, resolution=1, tickinterval=5, length=500, orient="horiz")
  tkgrid(img, columnspan=2)
  plotParFrame <- tkframe(plotFrame, relief="groove", borderwidth=2)
  color.but <- tkcheckbutton(plotParFrame,text="Colors", variable=colorvar, command=function() tkrreplot(img))
  ell.rbut <- tkradiobutton(plotParFrame, text="Ellipses", value=1, variable=mgrvar, command=function() tkrreplot(img))
  chul.rbut <- tkradiobutton(plotParFrame, text="Conv. hulls", value=2, variable=mgrvar, command=function() tkrreplot(img))
  tkgrid(color.but)
  tkgrid(ell.rbut)
  tkgrid(chul.rbut)
  tkgrid(plotParFrame, s)
  
  RCSFrame <- tkframe(plotFrame, relief="groove")
  cancel.but <- tkbutton(RCSFrame, text="Dismiss", command=function() tkdestroy(tt))
  submit.but <- tkbutton(RCSFrame, text="Submit", default="active", command=function() {doClass(); tkrreplot(img)})
  save.but <- tkbutton(RCSFrame, text="Save", default="active", command=function() {outgraphOrdiClust()})
  tkgrid(cancel.but, submit.but, save.but, ipadx=20)	
  tkgrid(RCSFrame, columnspan=2)
  
  tkgrid(compFrame, plotFrame)	
  
  
  if (!is.null(datatab)) {
    if (is.data.frame(datatab)) tclvalue(dfvar) <- as.character(substitute(datatab))
    if (is.dudi(datatab)) {
      ordiClust.dudi <- datatab
      assign("ordiClust.dudi", ordiClust.dudi, env_ade4tkgui)
      tclvalue(naxvar) <- as.character(ordiClust.dudi$nf)
      nax <- ordiClust.dudi$nf
      ploteig <- 1
      tkrreplot(img)
    }
  }
  
  #
  # Bind functions : up arrow
  #
  keyup <- function() {
    nlev <<- as.numeric(tclvalue(nlevvar))
    if (nlev < maxngr) {
      nlev <<- nlev + 1
      tclvalue(nlevvar) <- nlev
      if (!is.null(ordiClust.dudi)) if (class(h1) == "hclust") {
        if (nlev >= 2 && nlev <= nrow(ordiClust.dudi$li)) {
          ordiClust.factor <<- as.factor(cutree(h1, k=nlev))
          assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
          bet1 <- bca(ordiClust.dudi,as.factor(ordiClust.factor),scannf=FALSE)
          mc1 <- randtest(bet1)
          tkconfigure(bwg.label, text=paste(format(bet1$ratio*100, dig=2),"%",sep=""))
          tkconfigure(proba.label, text=format(mc1$pvalue, dig=3))
        }
        if (plotclust) {
          plotlev <<- 1
          hh <- cutree(h1, h=h1$height)
          for (i in 1:ncol(hh))
            if (nlevels(ordiClust.factor) == nlevels(as.factor(hh[,i])))
              hlev1 <<- (as.numeric(dimnames(hh)[[2]][i])+as.numeric(dimnames(hh)[[2]][i+1]))/2
          hlev <<- hlev1
          tclvalue(hlevvar) <- hlev1				
        }
      }
      tkrreplot(img)
    }
  }
  #
  # Bind functions : down arrow
  #
  keydown <- function() {
    nlev <<- as.numeric(tclvalue(nlevvar))
    if (nlev > 1) {
      nlev <<- nlev - 1
      tclvalue(nlevvar) <- nlev
      if (!is.null(ordiClust.dudi)) if (class(h1) == "hclust") {
        if (nlev >= 2 && nlev <= nrow(ordiClust.dudi$li)) {
          ordiClust.factor <<- as.factor(cutree(h1, k=nlev))
          assign("ordiClust.factor", ordiClust.factor, env_ade4tkgui)
          bet1 <- bca(ordiClust.dudi,as.factor(ordiClust.factor),scannf=FALSE)
          mc1 <- randtest(bet1)
          tkconfigure(bwg.label, text=paste(format(bet1$ratio*100, dig=2),"%",sep=""))
          tkconfigure(proba.label, text=format(mc1$pvalue, dig=3))
        } else {
          tkconfigure(bwg.label, text=paste(format(0, dig=2),"%",sep=""))
          tkconfigure(proba.label, text=format(NA, dig=3))
        }
        if (plotclust) {
          plotlev <<- 1
          hh <- cutree(h1, h=h1$height)
          for (i in 1:ncol(hh))
            if (nlevels(ordiClust.factor) == nlevels(as.factor(hh[,i])))
              hlev1 <<- (as.numeric(dimnames(hh)[[2]][i])+as.numeric(dimnames(hh)[[2]][i+1]))/2
          hlev <<- hlev1
          tclvalue(hlevvar) <- hlev1				
        }
      }
      tkrreplot(img)
    } else {
      tkconfigure(bwg.label, text=paste(format(0, dig=2),"%",sep=""))
      tkconfigure(proba.label, text=format(NA, dig=3))
    }
  }
  #
  # Bind functions : return
  #
  kret <- function() {
    doClass()
    tkrreplot(img)
  }
  
  tkbind(tt, "<KeyPress-Up>", keyup)
  tkbind(tt, "<KeyPress-Down>", keydown)
  tkbind(tt, "<KeyPress-Return>", kret)
  
  return(invisible())
}
