dotinfer <- function(dotdata, mtxdata) {
    require(tcltk) || stop("tcltk support is absent")
    require(tkrplot) || stop("tkrplot package is required")
    if (is.null(class(dotdata)) | class(dotdata) != "dotdata") {
        cat("The input object is not of class 'dotdata'\n")
        return(invisible())
    }
    if (!is.list(mtxdata)) {
        cat("Second argument must be a list\n")
        return(invisible())
    }
    if(!all(unlist(lapply(mtxdata, is.matrix)))) {
        cat("All list components of the second argument must be matrices\n")
        return(invisible())    
    } 
    nsp <- length(dotdata$Label)
    msg <- paste("Non-conformable arguments. \nAll matrices should have", 
           nsp, "X", nsp, "dimensions corresponding to the number of species", 
           "associated to the distributional data.")
    if(any(unlist(lapply(mtxdata, dim)) != nsp)) {
        cat(msg)
        return(invisible())    
    } 
    if(any(duplicated(names(mtxdata)))) {
        cat("Matrices in the second argument should be identified by unique names\n")
        return(invisible())
    }
    actsp <- 1 # Focus is put on the first species
    nghsp <- 1 # Neighbouring species to establish comparisons
    out <- list(sm = matrix(0, nsp, nsp), Label = dotdata$Label, occupancy = dotdata$occupancy, 
                coords = dotdata$coords, kind = "points")
    class(out) <- "dotinference"
    zoom <- FALSE
    wm <- mtxdata[[1]] 
    dendro <- hclust(as.dist(-1*wm), method = "single")
    vec <- wm[1,]
    vec[1] <- NA
    cbValue <- tclVar("0")
    seqsp <- order(vec, na.last = FALSE) 
    typmt <- tclVar("2")
    ruletext <- tclVar("Thresholding rule: >=")
    statusVar <- tclVar(paste("1", "/", nsp))
    mtxsel <- tclVar(names(mtxdata)[1])
    thresh <- tclVar()
    thr <- 0
    outputname <- tclVar("output_net")
    spsel <- tclVar(dotdata$Label[1])
    stblscore <- stablecouple(wm)$valref  
    cutoff  <- tclVar(80)
    typdn <- tclVar("1")
    kernel <- tclVar("gaussian")
    bw  <- tclVar("1")
    denscalc <- FALSE
    b <- c()
    ChangeMatrix <- function() {
      idmt <- which(names(mtxdata) == tclvalue(mtxsel))
      wm <<- mtxdata[[idmt]]
      if(tclvalue(typmt) == "1"){#Dissimilarity matrix
       stblaux <- stablecouple(wm, similarity = FALSE)
       dendro <<- hclust(as.dist(wm), method = "single")  
       if(tclvalue(cbValue) == "1") stblscore <<- stblaux$valref[stblaux$stpairs != 1:nsp]
       else stblscore <<- stblaux$valref
      }
      if(tclvalue(typmt) == "2"){#Similarity matrix
       stblaux <- stablecouple(wm)
       dendro <<- hclust(as.dist(-1*wm), method = "single")  
       if(tclvalue(cbValue) == "1") stblscore <<- stblaux$valref[stblaux$stpairs != 1:nsp]
       else stblscore <<- stblaux$valref
      }
      denscalc <<- FALSE
      ChangeSpecies()
      tkrreplot(densgr)
      tkrreplot(singlelink)
    }
    ChangeSpecies <- function() {
      actsp <<- which(dotdata$Label == tclvalue(spsel))
      nghsp <<- 1
      vec <- wm[actsp,]
      vec[actsp] <- NA
      if (tclvalue(typmt) == "1") seqsp <<- order(vec, na.last = FALSE)
      else seqsp <<- order(vec, na.last = FALSE, decreasing = TRUE)
      tkrreplot(dngr)
      tkrreplot(singlelink)
    }
    #Create the object dotinference
    Save <- function() {
      options(warn = -1)
      if(is.na(as.numeric(tclvalue(thresh)))) tclvalue(thresh) <- thr
      if(tclvalue(typmt) == "1") aa <- "<=" else aa <- ">="
      out$sm <- ifelse(eval(call(aa, wm, as.numeric(tclvalue(thresh)))), 1, 0) 
      if(tclvalue(outputname)=="") tclvalue(outputname) <- "output_net"
      assign(tclvalue(outputname), out, envir = .GlobalEnv)
      tkinsert(txt,"end","***********************\n")
      tkinsert(txt,"end",paste("Input Matrix:", tclvalue(mtxsel), "\n\n"))
      tkinsert(txt,"end",paste("Output dotinfer Object:", tclvalue(outputname),"\n\n"))
      tkinsert(txt,"end",paste("Thresholding:", aa, as.numeric(tclvalue(thresh))))
      tkinsert(txt,"end","\n**********************\n")
    }
    #####
    #Function that recovers the single linkage dendrogram    
    dendroplot <- function() {
      if(tclvalue(typmt) == "1") { #Dissimilarity matrix
            cutline <- thr
      }
      if(tclvalue(typmt) == "2"){#Similarity matrix
            cutline <- -1*thr
      }
      par(mar = rep(2, 4))
      plot(dendro, labels = rep("", nsp), axes = FALSE, main = "Single Linkage Dendrogram")    
      xred <- which(actsp == dendro$order)
      alter <- seqsp[nghsp]
      xblue <- which(alter == dendro$order)
      yred <-  row(dendro$merge)[dendro$merge == (-1*actsp)]
      yblue <- row(dendro$merge)[dendro$merge == (-1*alter)]
      segments(xblue, par()$usr[3], xblue, dendro$height[yblue], col = 4, lwd = 2)
      segments(xred, par()$usr[3], xred, dendro$height[yred], col = 2, lwd = 2) 
      abline(h = cutline, col = 5, lwd = 2) 
    }
    ##### End of function to obtain the single linkage dendrogram
    #Function to plot distributions
    dnplot <- function() {
        par(bg = "white")
        tclvalue(statusVar) <- paste(nghsp, "/", nsp) 
        ref <- array(FALSE, nrow(dotdata$coords))
        pts <- dotdata$occupancy[[actsp]] 
        if(zoom) ref[pts] <- TRUE else ref[] <- TRUE
        scdis1 <- 0.01*diff(range(dotdata$coords[,1]))
        scdis2 <- 0.01*diff(range(dotdata$coords[,2]))
        plot(dotdata$coords, bg = "white", col = gray(0.5), pch = 19, 
             xlim = (range(dotdata$coords[ref,1]) + c(-1*scdis1, scdis1)), 
             ylim = (range(dotdata$coords[ref,2]) + c(-1*scdis2, scdis2)))
        mtext(paste("Focus Species:", dotdata$Label[actsp]), line = 3, col = 2)
        if(tclvalue(typdn) == "1"){ 
          points(dotdata$coords[pts, ,drop = FALSE], col = 2, pch = 19, cex = 1.5)
        }  
        if(tclvalue(typdn) == "2") { 
          x0 <- dotdata$coords[dotdata$MSTsp[[actsp]][[3]], 1]
          y0 <- dotdata$coords[dotdata$MSTsp[[actsp]][[3]], 2]
          x1 <- dotdata$coords[dotdata$MSTsp[[actsp]][[4]], 1]
          y1 <- dotdata$coords[dotdata$MSTsp[[actsp]][[4]], 2]
          segments(x0, y0, x1, y1, lwd = 2.5, col = 2) 
        }
        if(tclvalue(typdn) == "3") { 
          x0 <- dotdata$coords[dotdata$MSTsp[[actsp]][[3]], 1]
          y0 <- dotdata$coords[dotdata$MSTsp[[actsp]][[3]], 2]
          x1 <- dotdata$coords[dotdata$MSTsp[[actsp]][[4]], 1]
          y1 <- dotdata$coords[dotdata$MSTsp[[actsp]][[4]], 2]
          segments(x0, y0, x1, y1, lwd = 2.5, col = 2) 
          points(dotdata$coords[pts, ,drop = FALSE], col = 2, pch = 19, cex = 1.5)
        }
        sp <- seqsp[nghsp]
        if(sp != actsp) {
          mtext(paste("Alter Species:", dotdata$Label[sp]), line = 2, col = 4)
          mtext(paste(tclvalue(mtxsel), "=", round(wm[actsp, sp], 5)), line = 1, cex = 0.8)
          pts <- dotdata$occupancy[[sp]]
          if(tclvalue(typdn) == "1"){ 
            points(dotdata$coords[pts, ,drop = FALSE], col = 4, pch = 19, cex = 1.2)
          }  
          if(tclvalue(typdn) == "2") { 
            x0 <- dotdata$coords[dotdata$MSTsp[[sp]][[3]], 1]
            y0 <- dotdata$coords[dotdata$MSTsp[[sp]][[3]], 2]
            x1 <- dotdata$coords[dotdata$MSTsp[[sp]][[4]], 1]
            y1 <- dotdata$coords[dotdata$MSTsp[[sp]][[4]], 2]
            segments(x0, y0, x1, y1, lwd = 2, col = 4) 
          }
          if(tclvalue(typdn) == "3") { 
            x0 <- dotdata$coords[dotdata$MSTsp[[sp]][[3]], 1]
            y0 <- dotdata$coords[dotdata$MSTsp[[sp]][[3]], 2]
            x1 <- dotdata$coords[dotdata$MSTsp[[sp]][[4]], 1]
            y1 <- dotdata$coords[dotdata$MSTsp[[sp]][[4]], 2]
            segments(x0, y0, x1, y1, lwd = 2, col = 4) 
            points(dotdata$coords[pts, ,drop = FALSE], col = 4, pch = 19, cex = 1.2)
          }
        }  
    }
    #####End of function to see the distributions
    #Interactive density plot, mainly based on the TCL version of Guido Masarotto
    replot <- function(...) {
        par(bg = "white", mar = c(0, 0, 1, 0))
        if (is.null(stblscore)) return() # too early...
        k <- as.character(tclObj(kernel))
        if (!denscalc) {x <- density(stblscore, kernel=k); b <<- x$bw; denscalc <<- TRUE}
        else x <- density(stblscore, bw=b*as.numeric(tclObj(bw)), kernel = k)
        eval(substitute(plot(x, main = "", axes = FALSE)))
        points(stblscore,rep(0,length(stblscore)))
        n <- length(x[[1]])
        areas <- diff(x[[1]][1:2])*(x[[2]][1:(n-1)] + x[[2]][2:n])*0.5 #trapezoid area
        cumareas <- c(0, cumsum(areas))
        cumareas[n] <- 1 #Force the sum to unity.
        if(tclvalue(typmt) == "1"){
            whref <- match(TRUE, cumareas >= as.numeric(tclvalue(cutoff))/100)
            thr <<- x[[1]][whref]
            tclvalue(thresh) <- thr
            pts <- x[[1]] <= thr
            mtext(paste("Cutoff <=", round(thr, 3)), line = 0)
            polygon(c(x[[1]][pts], thr), c(x[[2]][pts], tail(x[[2]], 1)), col = 2)
        }
        if(tclvalue(typmt) == "2"){
            whref <- match(TRUE, cumareas >= (1 - (as.numeric(tclvalue(cutoff))/100)))
            thr <<- x[[1]][whref]
            tclvalue(thresh) <- thr
            pts <- x[[1]] >= thr
            mtext(paste("Cutoff >=", round(thr, 3)), line = 0)
            polygon(c(x[[1]][pts], thr), c(x[[2]][pts], tail(x[[2]], 1)), col = 2)
        }
        abline(v = thr, col = 5)
    }
    #####End of function to display the density plot
    #Window to see the intensity matrix
    plotintensity <- function() {
        if(tclvalue(typmt) == "1") {fun1 <- "which.max"; fun2 <- "which.min"} 
        if(tclvalue(typmt) == "2") {fun1 <- "which.min"; fun2 <- "which.max"}
        vat <- function() {
          aux <- do.call(fun1, list(wm))
          nf <- aux %% nsp
          if(nf == 0) nf <- nsp
          I <- J <- c()
          P <- array(0, nsp)
          P[1] <- nf
          I <- nf
          J <- setdiff(1:nsp, nf)
          for(r in 2:nsp) {
            aux <- do.call(fun2, list(wm[I, J, drop = FALSE]))
            nc <-  ceiling(aux/length(I))
            P[r] <- J[nc]
            I <- c(I, J[nc])
            J <- setdiff(J, J[nc])
          }
          return(P) 
        }
        savepdf <- FALSE
        orvat <- vat()       
        orran <- c()
        orlab <- 1:nsp
        bound1 <- tclVar(min(wm))
        bound2 <- tclVar(max(wm))
        typgray <- tclVar("1")
        typor <- tclVar("2")
        sc <- tclVar("1")
        plotmtx <- function(){
          int <- range(as.numeric(tclvalue(bound1)), as.numeric(tclvalue(bound2)))
          plot.new()
          par(bg = "white", mar = rep(0, 4))
          if(length(int) != 2) return()
          plot.window(xlim = c(0, nsp), ylim = c(-nsp, 0))
          if(tclvalue(typor) == "1") {
             orran <<- sample(1:nsp)
             colcell <- wm[rep(orran, each = nsp) + rep(nsp*(orran - 1), nsp)]
          } 
          if(tclvalue(typor) == "2") colcell <- wm[rep(orvat, each = nsp) + rep(nsp*(orvat - 1), nsp)] 
          if(tclvalue(typor) == "3") colcell <- as.vector(wm) 
          if(diff(int) > 0) {
            graywm <- (colcell - int[1])/diff(int)
            graywm <- ifelse(graywm < 0, 0, graywm)
            graywm <- ifelse(graywm > 1, 1, graywm)
            if(fun1 == "which.min") graywm <- 1 - graywm #Similarity matrix. Zero means black
            if(tclvalue(typgray) == "2") {#Exponential decay
               cols <- gray(1 - exp(-1*as.numeric(tclvalue(sc))*graywm))
            } else cols <- gray(graywm)
          }
          if(diff(int)==0) { #Dichotomize
            if(fun1 == "which.min")
              cols <- ifelse(colcell >= int[1], gray(0), gray(1))
            if(fun1 == "which.max")
              cols <- ifelse(colcell <= int[1], gray(0), gray(1))
          }
          xleft <- rep(seq(0, nsp - 1), nsp)
          ybottom <- rep(seq(-1, -nsp), each = nsp)
          rect(xleft, ybottom, xleft + 1, ybottom + 1, col = cols, border = NA)  
          if(savepdf) {
            savepdf <<- FALSE
            filename <- tclvalue(tkgetSaveFile(initialfile = "VATmatrix.pdf",
                        defaultextension = ".pdf", title = "Save intensity matrix...",
                        filetypes = "{PDF {.pdf}} {{All Files} {*.*}}"))
            if (filename != "") pdf(file = filename) else return()
            plot.new()
            par(bg = "white", mar = rep(0, 4))
            plot.window(xlim = c(0, nsp), ylim = c(-nsp, 0))
            rect(xleft, ybottom, xleft + 1, ybottom + 1, col = cols, border = NA)  
            dev.off()
          }
        }
        intensity <- tktoplevel()
        tkwm.title(intensity, "Visual Assessment of Cluster Tendency")
        intens1 <- tkframe(intensity, relief = "groove", borderwidth = 2)
        intens2 <- tkframe(intensity, relief = "groove", borderwidth = 2)
        mtx <- tkrplot(intens2, hscale = 1.85, vscale = 1.85, fun = plotmtx)
        plotButton <- tkbutton(intens1, text = "Plot Intensity \nMatrix", command = function(...) tkrreplot(mtx)) 
        spButton <- tkbutton(intens1, text = "List Rearranged\nElements", command = function(...){
                             tt  <- tktoplevel()
                             xscr <- tkscrollbar(tt, repeatinterval=5,orient="horizontal",
                                                 command=function(...)tkxview(txt,...))
                             yscr <- tkscrollbar(tt, repeatinterval=5,
                                                 command=function(...)tkyview(txt,...))
                             txt <- tktext(tt,bg="white",font="courier",
                                     xscrollcommand=function(...)tkset(xscr,...),
                                     yscrollcommand=function(...)tkset(yscr,...), wrap="none")
                             tkgrid(txt,yscr)
                             tkgrid(xscr)
                             tkgrid.configure(yscr,sticky="ns")
                             tkgrid.configure(xscr,sticky="ew")
                             if(tclvalue(typor) == "1") {
                                splab <- dotdata$Label[orran]
                                tkinsert(txt,"end","**Items randomly arranged**\n")
                             }
                             if(tclvalue(typor) == "2") {
                                splab <- dotdata$Label[orvat]
                                tkinsert(txt,"end","**Items arranged by VAT algorithm**\n")
                             }
                             if(tclvalue(typor) == "3") {
                                splab <- dotdata$Label
                                tkinsert(txt,"end","**Items lexicographically arranged**\n")
                             }
                             for (i in (1:nsp)) tkinsert(txt,"end",paste(i,"=",splab[i],"\n"))
                             tkconfigure(txt, state="disabled")
                             tkfocus(txt)})
        copyPDF <- tkbutton(intens1, text = "Refresh and Save\nImage as PDF", 
                            command = function(...) {
                            savepdf <<- TRUE
                            tkrreplot(mtx)})
        bound1Entry <- tkentry(intens1, textvariable = bound1)
        bound2Entry <- tkentry(intens1, textvariable = bound2)
        scexp <- tkscale(intens1, from=1, to=100, showvalue=TRUE, variable=sc, resolution=1, orient="horiz")
        tkgrid(tklabel(intens1, text = "INPUT Adjacency Matrix\n", 
                       font = "Times 18",  foreground="blue"), columnspan = 2)
        tkgrid(tklabel(intens1, text = "Source", font = "Times 14"), columnspan = 2)
        tkgrid(tklabel(intens1, text = tclvalue(mtxsel)), columnspan = 2)
        tkgrid(tklabel(intens1, text = "\nRearrangement of Elements", font = "Times 14"), columnspan = 2)
        tkgrid(tkradiobutton(intens1, text="Random ordering", value=1, variable=typor), 
               columnspan = 2, sticky = "w")
        tkgrid(tkradiobutton(intens1, text="VAT algorithm", value=2, variable=typor),
               columnspan = 2, sticky = "w")
        tkgrid(tkradiobutton(intens1, text="Lexicographic ordering (by species label)", value=3, 
               variable=typor), columnspan = 2, sticky = "w")
        tkgrid(tklabel(intens1, text="\nRange of Brightness", font = "Times 14"), columnspan = 2)
        tkgrid(tklabel(intens1, text = "Bound 1 ="), bound1Entry, sticky = "e")
        tkgrid.configure(bound1Entry, sticky = "w")
        tkgrid(tklabel(intens1, text = "Bound 2 ="), bound2Entry, sticky = "e")
        tkgrid.configure(bound2Entry, sticky = "w")
        tkgrid(tklabel(intens1, text="\nBrightness Transfer Function", font = "Times 14"), columnspan = 2)
        tkgrid(tkradiobutton(intens1, text="Linear normalization",
               value=1, variable=typgray), sticky = "w")
        tkgrid(tkradiobutton(intens1, text="Natural exponential decay:\nf(x) = e^(-k*x), being k =",
                             value=2, variable=typgray), scexp, sticky = "w")
        tkgrid(tklabel(intens1, text="\n"), columnspan = 2)
        tkgrid(plotButton, spButton)
        tkgrid(copyPDF, tkbutton(intens1, text = "   Destroy all   \nthis", 
               command = function(...) tkdestroy(intensity)))
        tkgrid(tklabel(intens1, text="\nREMEMBER", font = "Times 14", foreground = "red"),
               columnspan = 2)
        txt1 <- "1) The stronger the association between items,\nthe darker the respective cell of the matrix" 
        txt2 <- paste("2) If both bounds are identical, the intensity\n",
                      "matrix will be a binary image at that threshold.", sep ="")
        txt3 <- paste("3) Select the appropriate kind of matrix\n",
                      "(i.e., similarity or dissimilarity)", sep ="")
        tkgrid(tklabel(intens1, text= txt1, font = "Times 12", foreground = "red"), 
               columnspan = 2, sticky = "w")
        tkgrid(tklabel(intens1, text= txt2, font = "Times 12", foreground = "red"), 
               columnspan = 2, sticky = "w")
        tkgrid(tklabel(intens1, text= txt3, font = "Times 12", foreground = "red"), 
               columnspan = 2, sticky = "w")
        tkgrid(tklabel(intens2, text = "OUTPUT Intensity Matrix", 
                       font = "Times 18",  foreground="blue"))
        tkgrid(mtx)
        tkgrid(intens1, intens2)
        tkgrid.configure(intens1, sticky = "n")
    }    
    #####End of window for tunning the intensity matrix
    #Window to see the network
    shownetwork <- function(sympmtx) {
      nettop <- tktoplevel()
      tkwm.title(nettop, "Sympatry Network")
      selsp <- tclVar("You have picked out no network node\n")
      or <- order(dendro$order) 
      X <- cos(or/nsp*2*pi)*250 + 300 
      Y <- sin(or/nsp*2*pi)*250 + 300 
      canvas <- tkcanvas(nettop, relief="raised", width=700, height=700, background = "white", 
                         borderwidth = 2)
      aux <- (row(sympmtx) < col(sympmtx)) & (sympmtx > 0)
      from <- row(sympmtx)[aux] 
      to <- col(sympmtx)[aux]
      chooseNode <- function(i){
        force(i)
        function(){
          tclvalue(selsp) <<- paste("You have picked out the node", i, "labelled\n", dotdata$Label[i])
          for(e in nodeEdges[[i]]) 
              tkitemconfigure(canvas, nodeItem[e$to], fill = "blue")
          tkitemconfigure(canvas, "text", fill = "black")
        }
      }
      moveNode <- function(i) {
        force(i)
        function(x, y) {
          x <- as.numeric(x)
          y <- as.numeric(y)
          width  <- as.numeric(tclvalue(tkwinfo("reqwidth",canvas)))
          height <- as.numeric(tclvalue(tkwinfo("reqheight",canvas)))
          if(x < 10) x <- 10
          if(x > width - 10) x <- width - 10
          if(y < 10) y <- 10
          if(y > height - 10) y <- height - 10
          for (e in nodeEdges[[i]]){
                  tkcoords(canvas,e$edgeItem,x,y,X[e$to],Y[e$to])
              }
          tkmove(canvas, nodeItem[i], x-X[i],y-Y[i])
          X[i] <<- x
          Y[i] <<- y
          tclvalue(selsp) <<- paste("You have picked out the network node",i,"labelled\n",dotdata$Label[i]) 
        }
      }
      nodeEdges <- vector("list",length(X))
      nodeItem <-  vector("character",length(X))
      for (i in seq(along=from)){
        fr <- from[i]
        t <- to[i]
        # add line to canvas
        e <- tkcreate(canvas, "line", X[fr],Y[fr],X[t],Y[t], width=2)
        nodeEdges[[fr]] <- c(nodeEdges[[fr]],list(list(to=t, edgeItem=e)))
        nodeEdges[[t]] <- c(nodeEdges[[t]],list(list(to=fr, edgeItem=e)))
      }
      for (i in seq(along=X)){
        # add the nodes
        p <- tkcreate(canvas,"oval",X[i]-10,Y[i]-10,X[i]+10,Y[i]+10, fill=gray(0.5))
        l <- tkcreate(canvas,"text", X[i] + 10, Y[i], text= i, anchor = "nw", font="10x20")
        tag <- paste("node",i,sep="")
        tkaddtag(canvas, tag, "withtag", p)
        tkaddtag(canvas, tag, "withtag", l)
        tkaddtag(canvas, "point", "withtag", p)
        tkaddtag(canvas, "text", "withtag", l)
        nodeItem[i] <- tag
        # animate them
        tkitembind(canvas, p, "<B1-Motion>", moveNode(i))
        tkitembind(canvas, p, "<ButtonPress-1>", chooseNode(i)) 
      }
      tkitembind(canvas, "point", "<B1-ButtonRelease>", function (){
                               tclvalue(selsp) <<- "You have picked out no network node\n" 
                               tkitemconfigure(canvas, "point", fill= gray(0.5))
                               tkitemconfigure(canvas, "current", fill="red")
                               })
      tkitembind(canvas, "point", "<Any-Enter>", function() tkitemconfigure(canvas, 
               "current", fill="red"))
      tkitembind(canvas, "point", "<Any-Leave>", function() tkitemconfigure(canvas, "current", fill= gray(0.5)))
      selectedsp <- tklabel(nettop, textvariable = selsp, font = "Times 16")
      tkgrid(selectedsp)
      tkgrid(canvas)
    }
    ##### End of window to see the network
    base <- tktoplevel()
    tkwm.title(base, "Sympatry Network Inference")
    l.frame <- tkframe(base)
    r.frame <- tkframe(base)
    #Frame with basic data
    setframe <-tkframe(r.frame, relief="groove", borderwidth=3)
    mtBox <- ttkcombobox(setframe, values = names(mtxdata), textvariable = mtxsel,
                         state="readonly")
    spBox <- ttkcombobox(setframe, values = dotdata$Label, textvariable = spsel,
                         state="readonly")
    #Binding function when user changes selection
    tkbind(mtBox, "<<ComboboxSelected>>", ChangeMatrix)
    tkbind(spBox, "<<ComboboxSelected>>", ChangeSpecies)
    tkgrid(tklabel (setframe, text="Choose Affinity Matrix", font = "Times 14", 
           foreground = "blue"))
    tkgrid(mtBox)
    disButton <- tkradiobutton(setframe, command= function(...) {
                         ChangeMatrix()
                         tkrreplot(densgr)
                         tclvalue(ruletext) <<- "Thresholding rule: <="}, text="Dissimilarity",
                         value=1, variable=typmt)
    simButton <- tkradiobutton(setframe, command= function(...) {
                         ChangeMatrix()
                         tkrreplot(densgr)
                         tclvalue(ruletext) <<- "Thresholding rule: >="}, text="Similarity",
                         value=2, variable=typmt)
    tkgrid(disButton, sticky = "w")
    tkgrid(simButton, sticky = "w")
    tkgrid(tklabel (setframe, font = "Times 14", text="Choose Focus Species",
           foreground = "blue"))
    tkgrid(spBox)
    tkgrid(tklabel(setframe, text="Distribution Plot", font = "Times 14", 
           foreground = "blue"))
    tkgrid(tkradiobutton(setframe, command = function() tkrreplot(dngr), text="Only dots",
           value=1, variable=typdn), sticky = "w")
    tkgrid(tkradiobutton(setframe, command= function() tkrreplot(dngr), text="Only MST",
           value=2, variable=typdn), sticky = "w")
    tkgrid(tkradiobutton(setframe, command= function() tkrreplot(dngr), text="MST and Dots",
           value=3, variable=typdn), sticky = "w")
    tkgrid(tklabel(setframe, text="Threshold Manager", font = "Times 14", 
           foreground = "blue"))
    checkButton <- tkcheckbutton(setframe, text = "Discard selfmatchings", variable=cbValue,
                                 command = function(...) {
                                           if(tclvalue(typmt) == "1"){#Dissimilarity matrix
                                              stblaux <- stablecouple(wm, similarity = FALSE)
                                              if(tclvalue(cbValue) == "1") stblscore <<- stblaux$valref[stblaux$stpairs != 1:nsp]
                                              else stblscore <<- stblaux$valref
                                           }
                                           if(tclvalue(typmt) == "2"){#Similarity matrix
                                             stblaux <- stablecouple(wm)
                                             if(tclvalue(cbValue) == "1") stblscore <<- stblaux$valref[stblaux$stpairs != 1:nsp]
                                             else stblscore <<- stblaux$valref
                                           }
                                           tkrreplot(densgr)
                                           tkrreplot(singlelink)
                                           })
    tkgrid(checkButton)
    tkgrid(tklabel(setframe, text="Bandwith"))
    tkgrid(tkscale(setframe, command= function(...) tkrreplot(densgr), from=0.05, to=2.00,
                   showvalue=FALSE, variable=bw, resolution=0.05, orient="horiz"))
    tkgrid(tkscale(setframe, command= function(...) {tkrreplot(densgr); tkrreplot(singlelink)}, from=0, to=100,
                   showvalue=TRUE, variable= cutoff, resolution=1, orient="horiz"))
    tkgrid(tklabel(setframe, text="Strongest relationships (% AUC)"))
    spec.frm <- tkframe(setframe, borderwidth=2)
    densgr <- tkrplot(spec.frm, fun = replot, hscale = 0.5, vscale = 0.5)
    tkgrid(densgr)
    #Frame for dendrogram
    proxframe <-tkframe(l.frame) 
    singlelink <- tkrplot(proxframe, fun = dendroplot, hscale = 1.2, vscale = 1.2)
    dendroButton <- tkbutton(proxframe, text = "Plot dendrogram on the current graphic device", command = function() {
                                             if(tclvalue(typmt) == "1") 
                                               plot(dendro, labels = dotdata$Label, main = "Single Linkage Clustering",
                                                    ylab = paste("Affinity Measure:", tclvalue(mtxsel)), xlab = "")
                                             if(tclvalue(typmt) == "2") 
                                               plot(dendro, labels = dotdata$Label, main = "Single Linkage Clustering",
                                                    ylab = paste("Association Measures: -1*", tclvalue(mtxsel), sep =""),
                                                    xlab = "")
                                             })
    outframe <- tkframe(proxframe, relief="groove", borderwidth=3) 
    createButton <- tkbutton(outframe, text = "Do", command = Save, padx = 20) 
    saveEntry <- tkentry(outframe, textvariable = outputname, width = 20)
    logdata <- tklabel(outframe, text = "") 
    tkgrid(singlelink)
    tkgrid(dendroButton)
    tkgrid(tklabel(outframe, text = "Save Binary Sympatry Network", font = "Times 14"), 
           logdata, columnspan = 3)
    tkgrid.configure(logdata, columnspan = 3)
    cutoffEntry <- tkentry(outframe, textvariable = thresh, width = 10)
    rule <- tklabel(outframe, textvariable = ruletext)
    tkgrid(rule, cutoffEntry)
    tkgrid.configure(cutoffEntry, sticky = "w", columnspan = 2)
    tkgrid.configure(rule, sticky = "e")
    tkgrid(tklabel(outframe, text = "Output Name (dotinfer object) ="), saveEntry, createButton, sticky = "e")
    tkgrid.configure(saveEntry, sticky = "w")
    tkgrid.configure(createButton, sticky = "w")
    tkgrid(outframe)
    #Frame for plotting species distributions
    dnframe <-tkframe(l.frame, relief="groove", borderwidth=3)
    dngr <- tkrplot(dnframe, fun = dnplot, hscale = 1.2, vscale = 1.2)
    fB <- tkbutton(dnframe, text = "<<", command = function(...) { 
                                         nghsp <<- 1
                                         tkrreplot(dngr)
                                         tkrreplot(singlelink)}, padx = 20) 
    lB <- tkbutton(dnframe, text = ">>", command = function(...) {
                                         nghsp <<- nsp
                                         tkrreplot(dngr)
                                         tkrreplot(singlelink)}, padx = 20)
    prevB <- tkbutton(dnframe, text = "<", command = function(...) {
                                         nghsp <<- pmax(1, nghsp - 1)
                                         tkrreplot(dngr)
                                         tkrreplot(singlelink)}, padx = 20) 
    nextB <- tkbutton(dnframe, text = ">", command = function(...) {
                                         nghsp <<- pmin(nghsp + 1, nsp)
                                         tkrreplot(dngr)
                                         tkrreplot(singlelink)}, padx = 20) 
    statusLabel <- tklabel(dnframe, textvariable = statusVar)
    tkgrid(dngr, columnspan = 5)
    tkgrid(fB, prevB, statusLabel, nextB, lB)
    subframe <- tkframe(dnframe)
    zoomout <- tkbutton(subframe, text = "Zoom Out", command = function(...) {
                                         zoom <<- FALSE  
                                         tkrreplot(dngr)}, padx = 20)  
    zoomin  <- tkbutton(subframe, text = "Zoom In", command = function(...) {
                                         zoom <<- TRUE  
                                         tkrreplot(dngr)}, padx = 20)  
    tkgrid(zoomout, tkbutton(subframe, text = "Copy to Clipboard", command = function(...) { 
                    tkrreplot(dngr)}, padx = 20), zoomin)
    tkgrid(subframe, columnspan = 5)
    #Frame to report activities
    textframe <- tkframe(l.frame, width = 480, height = 100)
    xscr <- tkscrollbar(textframe, repeatinterval=5,orient="horizontal",
                       command=function(...)tkxview(txt,...))
    yscr <- tkscrollbar(textframe, repeatinterval=5, 
                        command=function(...)tkyview(txt,...))
    txt <- tktext(textframe,bg="white",font="courier", xscrollcommand=function(...)tkset(xscr,...),
                  yscrollcommand=function(...)tkset(yscr,...), wrap="none")
    tkgrid(txt,yscr)
    tkgrid(xscr)
    tkgrid.configure(yscr,sticky="ns")
    tkgrid.configure(xscr,sticky="ew")
    tkconfigure(txt, width = 90, height = 10)
    tkgrid(dnframe, proxframe, sticky = "n")
    tkgrid(textframe, columnspan = 2)
    tkgrid(spec.frm)
    tkgrid(setframe, sticky = "n")
    tkgrid(l.frame, r.frame)
    #Menu
    topmenu <- tkmenu(base)
    tkconfigure(base, menu = topmenu)
    displaymenu <- tkmenu(topmenu, tearoff = FALSE)
    quitmenu <- tkmenu(topmenu)
    kernelmenu <- tkmenu(displaymenu, tearoff = FALSE)
    matrixmenu <- tkmenu(displaymenu, tearoff = FALSE)
    for (i in c("gaussian", "epanechnikov", "rectangular",
                 "triangular", "cosine")){
        tkadd(kernelmenu, "radiobutton", label = i, value=i, variable = kernel,
              command = function() tkrreplot(densgr))
    }
    tkadd(matrixmenu, "command", label = "Visual Assessment of Cluster Tendency", command = plotintensity)
    tkadd(displaymenu, "cascade", label = "Select kernel", menu = kernelmenu)
    tkadd(displaymenu, "cascade", label = "Intensity Matrix", menu = matrixmenu)
    tkadd(displaymenu,"command",label="Inferred Network", command = function(...) {
                      options(warn = -1)
                      if(is.na(as.numeric(tclvalue(thresh)))) tclvalue(thresh) <- thr
                      if(tclvalue(typmt) == "1") 
                         sympnetw <- ifelse(wm <= as.numeric(tclvalue(thresh)), 1, 0)
                      if(tclvalue(typmt) == "2") 
                         sympnetw <- ifelse(wm >= as.numeric(tclvalue(thresh)), 1, 0)
                      shownetwork(sympnetw)})
    tkadd(topmenu,"cascade",label="Display", menu=displaymenu)
    tkadd(topmenu,"command",label="Quit", command = function() tkdestroy(base))
}    

