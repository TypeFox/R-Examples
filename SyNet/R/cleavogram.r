cleavogram <- function () {
    require(tcltk) || stop("tcltk support is absent")
    require(tkrplot) || stop("tcltk support is absent")
    truelimits <- c()
    gr <- 0 #By default, no branch is considered selected
    usrCoords <- NULL
    targetseg <- NULL
    rbValue <- tclVar("10%")
    dens <- tclVar(0.5)
    miniwindow <- tclVar(1)
    transratio <- tclVar(0.3)
    mineccen <- tclVar(1)
    maxeccen <- tclVar(3)
    meaneccen <- tclVar(2)
    cvobj <- tclVar()
    done <- tclVar()
    color <- "yellow" #To select a single species
    linewidth <- tclVar(1.5)
    focussp <- tclVar("")
    xl <- yl <-c()
    searchgr <- function(dataref) {
        options(warn = -1)
        p1 <- dataref[,1] >= as.numeric(tclvalue(dens))
        p2 <- dataref[,2] >= as.numeric(tclvalue(transratio))
        p3 <- dataref[,3] <= ifelse(tclvalue(mineccen) %in% c(1,2,3), as.integer(tclvalue(mineccen)), Inf)
        p4 <- dataref[,4] <= ifelse(tclvalue(maxeccen) %in% c(1,2,3), as.integer(tclvalue(maxeccen)), Inf)
        aux <- ifelse(is.na(as.numeric(tclvalue(meaneccen))), Inf, as.numeric(tclvalue(meaneccen)))
        p5 <- dataref[,5] <= aux
        tclvalue(meaneccen) <<- aux
        whichr <- ((p1 & p2) & (p3 & p4)) & p5
        whichr <- replace(whichr, is.na(whichr), FALSE)
        return(whichr)
    }
    locategr <- function(x, y) {
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
        if (!is(cleavo, "cleavogram")) return()
        coords <- matrix(unlist(cleavo$components[, 6:9]), ncol = 4)
        Fracx <- as.numeric(x)/as.numeric(tkwinfo("reqwidth",cleavogram))
        Fracy <- 1 - as.numeric(y)/as.numeric(tkwinfo("reqheight",cleavogram))
        Newx <- truelimits[1] + Fracx*(truelimits[2] - truelimits[1])
        Newy <- truelimits[3] + Fracy*(truelimits[4] - truelimits[3])
        inx <- Newx >= coords[,1] & Newx <= (coords[,4] + 1)
        iny <- Newy >= (coords[,2] - 0.5) & Newy <= (coords[,3] + 0.5)
        targetseg <<- NULL
        gr <<- match(TRUE, inx & iny, nomatch = 0)
        if(gr > 0) targetseg <<- c(coords[gr,1] , mean(coords[gr,2:3]), coords[gr,4] + 1, mean(coords[gr,2:3]))
        tkrreplot(cleavogram)
        if(tclvalue(miniwindow) == "1") tkrreplot(miniplot)
    }
    locatesp <- function() {
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
        if (!is(cleavo, "cleavogram")) return()
        whsp <- match(tclvalue(focussp), cleavo$leaves, nomatch = 0)
        if(whsp == 0) return()
        brcol <- which((as.numeric(cleavo$components[,7]) <= whsp) & (as.numeric(cleavo$components[,8]) >= whsp))
        for(i in brcol){
          xl <<- c(xl, as.numeric(cleavo$components[i,c(6, 9)]) + c(0,1))
          yl <<- c(yl, rep(mean(as.numeric(cleavo$components[i,7:8])), 2))
        }
        xl <<- c(xl, tail(xl, 1), NA)
        yl <<- c(yl, whsp, NA)
        tkrreplot(cleavogram)
        if(tclvalue(miniwindow) == "2") tkrreplot(miniplot)
    }
    ya <- function (grsel = NULL) {
        par(bg = "white", mar = rep(0, 4))
        plot.new()
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
        if (!is(cleavo, "cleavogram")) return()
        if(is.null(usrCoords)) usrCoords <<- c(1, cleavo$nsub + 3.5, 0, cleavo$nsp + 0.5)
        plot.window(xlim = usrCoords[1:2], ylim = usrCoords[3:4])
        matrix(as.numeric(apply(cleavo$components, 1, tail, n = 4)), ncol= 4, byrow = TRUE) -> coords
        #Vertical lines
        segments(cleavo$verticals[,1], cleavo$verticals[,2], cleavo$verticals[,3], cleavo$verticals[,4])
        tips <- coords[,4] == cleavo$nsub
        if(any(tips)) segments(cleavo$nsub + 1, coords[tips,2], cleavo$nsub + 1, coords[tips,3])
        #Horizontal lines
        elong <- array(1, nrow(cleavo$components))
        elong[grep("Rem", rownames(cleavo$components))] <- 0.5
        yc <- apply(coords[,2:3], 1, mean)
        segments(coords[,1] , yc, coords[,4] + elong, yc)
        text(cleavo$namlast + 0.5*(cleavo$namlast <= cleavo$nsub), 1:cleavo$nsp, cleavo$leaves, col = 2, cex = 0.75, pos = 4)
        selbr <- searchgr(cleavo$components)
        coords <- matrix(unlist(cleavo$components[selbr, 6:9]), ncol = 4)
        points(0.5*(coords[,1] + coords[,4] + 1), 0.5*(coords[,2] + coords[,3]), col = 2, pch = 19)
        lines(xl, yl, col=color, lwd = as.numeric(tclvalue(linewidth)))
        if(!is.null(targetseg)) segments(targetseg[1],targetseg[2],targetseg[3],targetseg[4], col = "blue", lwd = 3)
        truelimits <<- par()$usr
    } 
    ya.small <- function () {
        par(bg = "white", mar = rep(0, 4))
        plot.new()
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
        if (!is(cleavo, "cleavogram")) return()
        if(tclvalue(miniwindow) == "1") {
          ptn <- cleavo 
          if(gr > 0) listofsp <- cleavo$leaves[cleavo$components[gr,7]:cleavo$components[gr, 8]]
          else listofsp <- NULL
          if(ptn$kind == "grids") {
              plot.window(xlim = range(ptn$coords[,1]), ylim = range(-1*ptn$coords[,2])) #chequear exactitud
              x <- as.numeric(ptn$coords[,1])
              y <- as.numeric(ptn$coords[,2])                           
              v <- seq(0.5, max(x) + 0.5, by = 1)
              h <- seq(-0.5, -1*max(y) - 0.5, by = -1)
              abline(v = v)
              abline(h = h)
              rect(x - 0.5, -1*y - 0.5, x + 0.5, -1*y + 0.5, col = gray(0.7))
              dnes <- unlist(ptn$occupancy[listofsp])
              whforgr <- unique(dnes)
              #freq <- c()
              #for(i in whforgr) freq <- c(freq, sum(dnes == i)/length(listofsp))
              #rect(x[whforgr] - 0.5, -1*y[whforgr] - 0.5, x[whforgr] + 0.5, -1*y[whforgr] + 0.5, col = gray(1 - freq))
              rect(x[whforgr] - 0.5, -1*y[whforgr] - 0.5, x[whforgr] + 0.5, -1*y[whforgr] + 0.5, col = "blue")
          }
          else if(ptn$kind == "points") {
              plot.window(xlim = range(ptn$coords[,1]), ylim = range(ptn$coords[,2]))
              points (ptn$coords, col = gray(0.9), pch = 19, asp = 1)
              grsel <- unique(unlist(ptn$occupancy[listofsp]))
              points(ptn$coords[grsel, ,drop = FALSE], pch = 19)
          }
        }
        else if (tclvalue(miniwindow) == "2") {
          plot.window(xlim = c(1, cleavo$nsub + 0.5), ylim = c(0,cleavo$nsp + 0.5))
          segments(cleavo$verticals[,1], cleavo$verticals[,2],cleavo$verticals[,3], cleavo$verticals[,4])
          matrix(as.numeric(apply(cleavo$components, 1, tail, n = 4)),ncol = 4, byrow = TRUE) -> coords
          segments(coords[,1] , apply(coords[,2:3], 1, mean), coords[,4] + 1, apply(coords[,2:3], 1, mean) )
          rect(usrCoords[1], usrCoords[3], usrCoords[2], usrCoords[4], border = "red")
          lines(xl, yl, col=color, lwd = as.numeric(tclvalue(linewidth)))
        }
    }
    savecleavogram <- function() {
        tf <- tktoplevel()
        tkwm.title(tf, "Save Cleavogram")
        frame1 <- tkframe(tf, relief = "groove", borderwidth = 2)
        frame2 <- tkframe(tf, relief = "groove", borderwidth = 2)
        done <- tclVar(0)
        formatvar <- tclVar(1)
        savefic <- function(formatvar) {
            outform <- tclvalue(formatvar)
            if (outform == "1") {
                filename <- tclvalue(tkgetSaveFile(initialfile = "Cleavo.ps",
                  defaultextension = ".ps", title = "Save cleavogram...",
                  filetypes = "{PostScript {.ps .eps}} {{All Files} {*.*}}"))
                if (filename != "") postscript(file = filename)
            }
            else if (outform == "2") {
                filename <- tclvalue(tkgetSaveFile(initialfile = "Cleavo.pdf",
                  defaultextension = ".pdf", title = "Save cleavogram...",
                  filetypes = "{PDF {.pdf}} {{All Files} {*.*}}"))
                if (filename != "") pdf(file = filename)
            }
            else if (outform == "3") {
                filename <- tclvalue(tkgetSaveFile(initialfile = "Cleavo.jpeg",
                  defaultextension = ".jpeg", title = "Save cleavogram...",
                  filetypes = "{JPEG {.jpeg .jpg}} {{All Files} {*.*}}"))
                if (filename != "") jpeg(file = filename)
            }
            ya() #Print the image
            dev.off()
            tkdestroy(tf)
        }
        tkgrid(tklabel(tf, text = "Save Current Graphic", font = "Times 18"),
            columnspan = 2)
        tkgrid(tklabel(frame2, text = "Output format : "), sticky = "n")
        tkgrid(tkradiobutton(frame2, text = "postscript", value = 1,
            variable = formatvar), sticky = "w")
        tkgrid(tkradiobutton(frame2, text = "pdf", value = 2,
            variable = formatvar), sticky = "w")
        tkgrid(tkradiobutton(frame2, text = "jpeg", value = 3,
            variable = formatvar), sticky = "w")
        tkgrid(frame2, rowspan = 2, sticky = "n")
        save.but <- tkbutton(frame1, text = "Save", command = function() savefic(formatvar))
        cancel.but <- tkbutton(frame1, text = "Cancel", command = function() tkdestroy(tf))
        tkgrid(save.but, cancel.but)
        tkgrid(frame1, column = 1, row = 2, sticky = "n")
        tkbind(tf, "<KeyPress-Return>", function() savefic(formatvar))
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
        tkwait.variable(done)
        if (tclvalue(done) == "2")
            return(0)
        tkdestroy(tf)
    }
    flatpartition <- function() {
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux)
        else {
          tkmessageBox(title = "NO DATA", message = "Please, provide some valid input", icon = "info", type = "ok")
          return()
        }
        if (!is(cleavo, "cleavogram")) {
          tkmessageBox(title = "NO DATA", message = "Please, provide some valid input", icon = "info", type = "ok")
          return()
        }
        done <- tclVar(0)
        outpartition <- tclVar("out_partition")
        formatvar <- tclVar(3)
        out <- c()
        batch <- c()
        batch <- grep("Rem", rownames(cleavo$components))
        batch <- union(batch, which(cleavo$components[,9] == max(cleavo$components[,9])))
        flpt <- tktoplevel()
        tkwm.title(flpt, "Generate Flat Partition")
        frame1 <- tkframe(flpt, relief = "groove", borderwidth = 2)
        frame2 <- tkframe(flpt)
        xscr <- tkscrollbar(flpt, repeatinterval=5,orient="horizontal",
                              command=function(...)tkxview(txt,...))
        yscr <- tkscrollbar(flpt, repeatinterval=5,
                              command=function(...)tkyview(txt,...))
        txt <- tktext(flpt,bg="white",font="courier",
            xscrollcommand=function(...)tkset(xscr,...),yscrollcommand=function(...)tkset(yscr,...),
            wrap="none")
        savefic <- function(formatvar) {
            out <<- c()
            textshow <- ""
            selbr <- searchgr(cleavo$components)
            fromto <- cleavo$components[selbr, 7:8, drop = FALSE]
            whichgrn <- which(selbr)
            fixbr <- c()
            whichgrf <- array(TRUE, sum(selbr))
            outform <- tclvalue(formatvar)
            textshow <- switch(outform, 
                              "1" = "Forward, first occurrence",
                              "2" = "Backward, first occurrence",
                              "3" = "Forward, conditioned advance")
            if (outform == "1" & any(selbr)) {
              while(any(whichgrf)) {
                first <- match(TRUE, whichgrf)
                fixbr <- c(fixbr, whichgrn[first])
                #FALSE for the largest group that matches the criteria and also
                #for the subordinated groups
                whichgrf[(fromto[first,1] <= fromto[,1]) & (fromto[first,2] >= fromto[,2])] <- FALSE
              }
            }
            else if (outform == "2" & any(selbr)) {
              while(any(whichgrf)) {
                first <- match(TRUE, whichgrf)
                whichgrf[first] <- FALSE
                if (sum((fromto[first,1] <= fromto[,1]) & (fromto[first,2] >= fromto[,2])) > 1) next
                fixbr <- c(fixbr, whichgrn[first])
                }
            }
            else if (outform == "3" & any(selbr)) {
                while(any(whichgrf)) {
                  first <- match(TRUE, whichgrf)
                  whichgrf[first] <- FALSE
                  desc <- ((fromto[first,1] <= fromto[,1]) & (fromto[first,2] >= fromto[,2])) & whichgrf
                  if (sum(desc) < 2) {fixbr <- c(fixbr, whichgrn[first]); whichgrf[desc] <- FALSE; next}
                  nosubd <- TRUE
                  for (i in fromto[desc,1])  if(any(i > fromto[desc,2])) {nosubd <- FALSE; break}
                  #I need to check if there are disjoint entities through a traversing in depth.
                  #Groups are characterized by the consecutive indices of their members.
                  #Two groups are disjoint if their intervals do not overlap.
                  if (nosubd) {whichgrf[desc] <- FALSE; fixbr <- c(fixbr, whichgrn[first])}
                }
            }
            #Next, descendant branches from those holded by fixbr are discarded
            nomore <- c()
            for (i in fixbr) {
              fcond <- which(as.integer(cleavo$components[,7]) >= as.integer(cleavo$components[i,7]))
              scond <- which(as.integer(cleavo$components[,8]) <= as.integer(cleavo$components[i,8]))
              nomore <- c(nomore, intersect(fcond, scond))
            }
            brend <- union(setdiff(batch, nomore), fixbr)
            UC <- diad <- iso <- inter <- 1
            options(warn = -1)
            tkinsert(txt, "end", " *********  FLAT  PARTITION  *********\n\n")
            tkinsert(txt, "end", "****  CRITERIA FOR COHESIVENESS  ****\n\n")
            tkinsert(txt, "end", " [Density/Transitivity/(Min;Max;Mean Eccentricities)]\n\n")
            tkinsert(txt, "end", paste("    ", round(as.numeric(tclvalue(dens)), 2),"/",
                                  round(as.numeric(tclvalue(transratio)), 2),"/(", tclvalue(mineccen), ";",
                                  tclvalue(maxeccen), ";", round(as.numeric(tclvalue(meaneccen)), 2), ")"))
            tkinsert(txt, "end", "\n*****************************************\n\n")
            tkinsert(txt, "end", " *********  SEARCH DIRECTION  *********\n\n")
            tkinsert(txt, "end", paste("      ", textshow, "\n\n"))
            tkinsert(txt, "end", "\n*****************************************\n\n")
            tkinsert(txt, "end", "************  FINAL STATUS  ************\n\n")
            for (i in brend) {
              dif <- diff(as.integer(cleavo$components[i,7:8]))
              if (dif > 1) {
                tkinsert(txt, "end", paste("\nUC ", UC, ":\n  "))
                spp <- seq(as.integer(cleavo$components[i,7]), as.integer(cleavo$components[i,8]), by = 1)
                for (j in spp) {
                  tkinsert(txt, "end", paste(cleavo$leaves[j], "\n  "))
                  out <<- rbind(out, c(cleavo$leaves[j], paste("UC ", UC)))
                }
                UC <- UC + 1
              }
              else if (dif == 1) {
                tkinsert(txt, "end", paste("\nDiad ", diad, ":\n  "))
                sp1 <- cleavo$leaves[as.integer(cleavo$components[i,7])]
                sp2 <- cleavo$leaves[as.integer(cleavo$components[i,8])]
                tkinsert(txt, "end", paste(sp1, "--" , sp2))
                out <<- rbind(out, c(sp1, paste("Diad ", diad)))
                out <<- rbind(out, c(sp2, paste("Diad ", diad)))
                diad <- diad + 1
              }
              else if (dif == 0 & cleavo$namlast[cleavo$components[i,8]] > cleavo$nsub) {
                tkinsert(txt, "end", paste("\nIsoalted ", iso, ":\n  "))
                sp1 <- cleavo$leaves[as.integer(cleavo$components[i,7])]
                tkinsert(txt, "end", sp1)
                out <<- rbind(out, c(sp1, "Isolated"))
                iso <- iso + 1
              }
              else if (dif == 0 & cleavo$namlast[cleavo$components[i,8]] <= cleavo$nsub) {
                tkinsert(txt, "end", paste("\nRemoved at Sub-Network", as.numeric(cleavo$components[i,6]) - 2,": \n  "))
                sp1 <- cleavo$leaves[as.integer(cleavo$components[i,7])]
                tkinsert(txt, "end", sp1)
                out <<- rbind(out, c(sp1, "Intermediary"))
                inter <- inter + 1
              }
            } #End for
            tkinsert(txt, "end", "\n*****************************************\n\n")
            tkfocus(txt)
            rslt <- c()
            rslt$kind <- cleavo$kind
            rslt$status <- cbind(Taxa = out[,1], Status = out[,2])
            rslt$occupancy <- cleavo$occupancy
            rslt$coords <- cleavo$coords
            class(rslt) <- "nampartition"
            if(tclvalue(outpartition)=="") tclvalue(outpartition) <- "out_partition"
            assign(tclvalue(outpartition), rslt, pos = 1)
        }
        outEntry <- tkentry(frame1, textvariable = outpartition, width = 20)
        save.but <- tkbutton(frame2, text = " Run! ", command = function() savefic(formatvar), padx = 25)
        cancel.but <- tkbutton(frame2, text = "Cancel", command = function() tkdestroy(flpt), padx = 25)
        tkgrid(tklabel(flpt, text = "Flat Partition from Cleavogram in", font = "Times 16",
                      foreground = "blue"))
        tkgrid(tklabel(flpt, text = aux, font = "Times 16", foreground = "blue"))
        tkgrid(tklabel(frame1, text = "Output Name:", font = "Times 14"))
        tkgrid(outEntry)
        tkgrid(tklabel(frame1, text = "Search Options:", font = "Times 14"))
        tkgrid(tkradiobutton(frame1, text = "Forward, first occurrence", value = 1,
            variable = formatvar), sticky = "w")
        tkgrid(tkradiobutton(frame1, text = "Backward, first occurrence", value = 2,
            variable = formatvar), sticky = "w")
        tkgrid(tkradiobutton(frame1, text = "Forward, conditioned advance", value = 3,
            variable = formatvar), sticky = "w")
        tkgrid(save.but, cancel.but)
        tkgrid(frame1)
        tkgrid(frame2)
        tkgrid(txt,yscr)
        tkgrid(xscr)
        tkgrid.configure(yscr,sticky="ns")
        tkgrid.configure(xscr,sticky="ew")
        tkbind(flpt, "<Destroy>", function() tkdestroy(flpt))
        tkbind(flpt, "<KeyPress-Return>", function() savefic(formatvar))
        tkbind(flpt, "<KeyPress-Escape>", function() tkdestroy(flpt))
    } #End of flat partition function
    zoom <- function(f) {
        fac <- f*as.numeric(unlist(strsplit(tclvalue(rbValue), "%")))/100
        cx <- fac*(usrCoords[2] - usrCoords[1])/2
        cy <- fac*(usrCoords[4] - usrCoords[3])/2
        usrCoords <<- usrCoords + c(cx, -cx, cy, -cy)
        tkrreplot(cleavogram)
        tkrreplot(miniplot)
    }
    advance <- function(n) {
        if (n == 1) usrCoords[1:2] <<- usrCoords[1:2] + 0.10*(usrCoords[2] - usrCoords[1]) #The visual effect is 10% Rightward
        else if (n == 2) usrCoords[1:2] <<- usrCoords[1:2] - 0.10*(usrCoords[2] - usrCoords[1]) #Leftward
        else if (n == 3) usrCoords[3:4] <<- usrCoords[3:4] + 0.10*(usrCoords[4] - usrCoords[3])
        else if (n == 4) usrCoords[3:4] <<- usrCoords[3:4] - 0.10*(usrCoords[4] - usrCoords[3])
        tkrreplot(cleavogram)
        tkrreplot(miniplot)
    }
    menurightbutton <- function(x, y) {
        menuright <- tkmenu(cleavogram, tearoff = FALSE)
        if (.Platform$OS.type == "windows") {
              tkadd(menuright, "command", label = "Copy", command = function() tkrreplot(cleavogram))
              tkadd(menuright, "separator")
          }
        tkadd(menuright, "command", label = "Save", command = savecleavogram)
        Newx <- as.integer(tkwinfo("rootx",cleavogram)) + as.integer(x)
        Newy <- as.integer(tkwinfo("rooty",cleavogram)) + as.integer(y)
        tcl("tk_popup", menuright, Newx, Newy)
    }
    loadbranch <- function(x, y){
        locategr(x, y)
        if(gr == 0) return()
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
        if (!is(cleavo, "cleavogram")) return()
        tf <- tktoplevel()
        tkwm.title(tf, "Spatial Expression of Selected Branch")
        #tkwm.geometry(tf, "770x605")
        #tkwm.resizable(tf, 0, 0)
        listofsp <- cleavo$leaves[cleavo$components[gr, 7]:cleavo$components[gr, 8]]
        whforgr <- unique(unlist(cleavo$occupancy[listofsp]))
        global <- TRUE
        ref <- rep(TRUE, nrow(cleavo$coords))
        targetsp <- NULL
        coldata <- gray(0.7)
        colgr <- "black"
        colsp <- "red"
        sizedata <- tclVar(1)
        sizegr <- tclVar(1.5)
        sizesp <- tclVar(1)
        frameEdit <- tkframe(tf, relief = "groove", borderwidth = 2)
        frameList <- tkframe(frameEdit)
        tlbsp <- tklistbox(frameList)
        scrsp <- tkscrollbar(frameList, repeatinterval = 5, command = function(...) tkyview(tlbsp,
            ...))
        tkconfigure(tlbsp, yscrollcommand = function(...) tkset(scrsp,
            ...))
        for(i in listofsp) tkinsert(tlbsp, "end", i)
        tkselection.set(tlbsp,0) #Default selection is the first item
        drawmap <- function () {
          par(bg = "white")
          if(cleavo$kind == "none") {plot.new(); return()}
          if(global) ref[] <- TRUE else {ref[] <- FALSE; ref[whforgr] <- TRUE}
          numc <- tclvalue(tkcurselection(tlbsp))
          targetsp <- tclvalue(tkget(tlbsp, numc))
          spsel <- unlist(cleavo$occupancy[targetsp])
          if(cleavo$kind == "grids") {
              x <- as.numeric(cleavo$coords[,1])
              y <- -1*as.numeric(cleavo$coords[,2])
              plot(x, y, xlab = "X_coord", ylab = "-(Y_coord)", main = paste("Item selected:",targetsp),
                  xlim = range(x[ref]) + c(-1, 1), ylim = range(y[ref]) + c(-1, 1))
              abline(v = seq(0.5, max(x[ref]) + 0.5, by = 1))
              abline(h = seq(-0.5, min(y[ref]) - 0.5, by = -1))
              rect(x - 0.5, y - 0.5, x + 0.5, y + 0.5, col = coldata, border = coldata, 
                   lwd = as.numeric(tclvalue(sizedata)))
              rect(x[whforgr] - 0.5, y[whforgr] - 0.5, x[whforgr] + 0.5, y[whforgr] + 0.5, col = colgr, 
                   border = colgr, lwd = as.numeric(tclvalue(sizegr)))
              rect(x[spsel] - 0.5, y[spsel] - 0.5, x[spsel] + 0.5, y[spsel] + 0.5, col = colsp,
                   border = colsp, lwd = as.numeric(tclvalue(sizesp)))
          }
          if(cleavo$kind == "points") {
              plot(cleavo$coords, xlab = "LONGITUDE", ylab = "LATITUDE", main = paste("Item selected:",targetsp),
                  col = coldata, pch = 19, asp = 1, cex = as.numeric(tclvalue(sizedata)),
                  xlim = range(cleavo$coords[ref,1]), ylim = range(cleavo$coords[ref,2]))
              points(cleavo$coords[whforgr, ,drop = FALSE], col = colgr, pch = 19, cex = as.numeric(tclvalue(sizegr)))
              points(cleavo$coords[spsel, ,drop = FALSE], col = colsp, pch = 19, cex = as.numeric(tclvalue(sizesp)))
          }
        } #End drawmap function
        framePlot <- tkframe(tf, relief = "groove", borderwidth = 2)
        cancel.but <- tkbutton(framePlot, text = "Exit", command = function() tkdestroy(tf), padx = 30)
        alldata.but <- tkbutton(framePlot, text = "Full Extent", command = function() {
                                                                global <<- TRUE
                                                                tkrreplot(mapview)},padx = 15)
        onlygroup.but <- tkbutton(framePlot, text = "Focus on Active Group", command = function() {
                                                                  global <<- FALSE
                                                                  tkrreplot(mapview)}, padx = 15)
        copyimg.but <- tkbutton(framePlot, text = "Copy to Clipboard", command = function() tkrreplot(mapview))
        mapview <- tkrplot(framePlot, fun = drawmap, 1.5, 1.5)
        canvasdata <- tkcanvas(frameEdit, width="80", height="25", bg=coldata)
        colorButtondata <- tkbutton(frameEdit, text="Change",command= function() {
                                auxcol <- coldata
                                tkconfigure(colorButtondata, state = "disabled")
                                coldata <<- tclvalue(tcl("tk_chooseColor", title="Choose a color"))
                                if (nchar(coldata)==0) coldata <<- auxcol
                                tkconfigure(canvasdata,bg=coldata)
                                tkconfigure(colorButtondata, state = "active")
                                tkrreplot(mapview)
                                })
        canvasgr <- tkcanvas(frameEdit, width="80", height="25", bg=colgr)
        colorButtongr <- tkbutton(frameEdit, text="Change",command= function() {
                                auxcol <- colgr
                                tkconfigure(colorButtongr, state = "disabled")
                                colgr <<- tclvalue(tcl("tk_chooseColor", title="Choose a color"))
                                if (nchar(colgr)==0) colgr <<- auxcol
                                tkconfigure(colorButtongr, state = "active")
                                tkconfigure(canvasgr,bg=colgr)
                                tkrreplot(mapview)
                                })
        canvassp <- tkcanvas(frameEdit, width="80", height="25", bg=colsp)
        colorButtonsp <- tkbutton(frameEdit, text="Change",command= function() {
                                auxcol <- colsp
                                tkconfigure(colorButtonsp, state = "disabled")
                                colsp <<- tclvalue(tcl("tk_chooseColor", title="Choose a color"))
                                if (nchar(colsp)==0) colsp <<- auxcol
                                tkconfigure(colorButtonsp, state = "active")
                                tkconfigure(canvassp,bg=colsp)
                                tkrreplot(mapview)
                                })
        tkgrid(tklabel(frameList, text = "Items", font = "Times 17", foreground = "blue"), columnspan=2)
        tkgrid(tlbsp, scrsp, columnspan = 2)
        tkgrid.configure(scrsp, sticky = "wns")
        tkgrid.configure(tlbsp, sticky = "e")
        tkgrid(frameList, columnspan = 2)
        tkgrid(tklabel(frameEdit, text = "\nColor and Dot Size\nSetting",
               font = "Times 17", foreground = "blue"), columnspan = 2)
        tkgrid(tklabel(frameEdit, text = "Universe of records\n(global dataset)",
               font = "Times 12", foreground = "blue"), columnspan = 2)
        tkgrid(canvasdata, colorButtondata)
        tkgrid.configure(colorButtondata, sticky = "w")
        tkgrid.configure(canvasdata, sticky = "e")
        tkgrid(tkscale(frameEdit, from = 0.5, to = 4, resolution = 0.1,
               orient = "horiz", variable = sizedata, command = function(...) tkrreplot(mapview),
               showvalue = FALSE), columnspan = 2)
        tkgrid(tklabel(frameEdit, text = "\nPooled set of records\n(selected group)",
               font = "Times 12", foreground = "blue"), columnspan = 2)
        tkgrid(canvasgr, colorButtongr)
        tkgrid.configure(colorButtongr, sticky = "w")
        tkgrid.configure(canvasgr, sticky = "e")
        tkgrid(tkscale(frameEdit, from = 0.5, to = 4, resolution = 0.1, orient = "horiz",
               variable = sizegr, command = function(...) tkrreplot(mapview),
               showvalue = FALSE), columnspan = 2)
        tkgrid(tklabel(frameEdit, text = "\nSpecies set of records\n(selected species)",
               font = "Times 12", foreground = "blue"), columnspan = 2)
        tkgrid(canvassp, colorButtonsp)
        tkgrid.configure(colorButtonsp, sticky = "w")
        tkgrid.configure(canvassp, sticky = "e")
        tkgrid(tkscale(frameEdit, from = 0.5, to = 4, resolution = 0.1, orient = "horiz",
               variable = sizesp, command = function(...) tkrreplot(mapview),
               showvalue = FALSE), columnspan = 2)
        tkgrid(tklabel(frameEdit, text =""), columnspan = 2)
        tkgrid(mapview, columnspan = 4, sticky = "n")
        tkgrid(alldata.but, onlygroup.but, copyimg.but, cancel.but)
        tkgrid.configure(onlygroup.but, sticky = "ew")
        tkgrid.configure(alldata.but, sticky = "ew")
        tkgrid.configure(copyimg.but, sticky = "ew")
        tkgrid.configure(cancel.but, sticky = "ew")
        tkgrid(frameEdit, framePlot)
        tkbind(tlbsp, "<ButtonRelease-1>", function() tkrreplot(mapview))
        tkbind(tlbsp, "<KeyRelease>", function() tkrreplot(mapview))
        tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
        tkbind(tf, "<Destroy>", function() tkdestroy(tf))
    }
    showresult <- function(){
        aux <- as.character(tclvalue(cvobj))
        if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
        if (!is(cleavo, "cleavogram")) return()
        tt  <- tktoplevel()
        tkwm.title(tt, "Full Report...")
        xscr <- tkscrollbar(tt, repeatinterval=5,orient="horizontal",
                              command=function(...)tkxview(txt,...))
        yscr <- tkscrollbar(tt, repeatinterval=5,
                              command=function(...)tkyview(txt,...))
        txt <- tktext(tt,bg="white",font="courier",
            xscrollcommand=function(...)tkset(xscr,...),yscrollcommand=function(...)tkset(yscr,...),
            wrap="none")
        tkgrid(txt,yscr)
        tkgrid(xscr)
        tkgrid.configure(yscr,sticky="ns")
        tkgrid.configure(xscr,sticky="ew")
        options(warn = -1)
        tkinsert(txt, "end", "****  CRITERIA FOR COHESIVENESS  ****\n\n")
        tkinsert(txt, "end", " [Density/Transitivity/(Min;Max;Mean Eccentricities)]\n\n")
        aux <- ifelse(is.na(as.numeric(tclvalue(meaneccen))), Inf, as.numeric(tclvalue(meaneccen)))
        tclvalue(meaneccen) <<- aux
        tkinsert(txt, "end", paste("    ", round(as.numeric(tclvalue(dens)), 2),"/",
                            round(as.numeric(tclvalue(transratio)), 2),"/(", tclvalue(mineccen), ";",
                            tclvalue(maxeccen), ";", round(as.numeric(tclvalue(meaneccen)), 2), ")"))
        tkinsert(txt, "end", "\n*****************************************\n\n")
        tkinsert(txt, "end", "*********  CODES FOR TAXA NAMES  *********\n\n")
        for (i in 1:cleavo$nsp) tkinsert(txt, "end", paste(i, ":" , cleavo$leaves[i], "\n"))
        tkinsert(txt, "end", "\n*****************************************\n\n")
        tkinsert(txt, "end", "*********  LIST OF BRANCHES  *********\n\n")
        tkinsert(txt, "end", "Branch ID [Cohesiveness parameters] --> [From Taxa ID] -- [To Taxa ID]\n\n")
        brnam <- rownames(cleavo$components)
        i <- 1
        apply(cleavo$components, 1, function(x) {
                                    aux <- as.character(round(as.numeric(x[1:5]), 2))
                                    tkinsert(txt,"end", paste("\nBranch", brnam[i],"[", aux[5], "/", aux[4],
                                              "/(", aux[1], ";", aux[2], ";", aux[3], ")] --> "))
                                    tkinsert(txt, "end", paste(x[7], "--", x[8]))
                                    i <<- i + 1  })
        tkinsert(txt, "end", "\n*****************************************\n\n")
        tkfocus(txt)
      }  #End showresult
      loadcombo <- function(combosp, aux) {
            if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
            if (!is(cleavo, "cleavogram")) return()
            tclvalue(focussp) <<- ""
            tkconfigure(combosp, values = sort(cleavo$leaves))
      }
      #Main frames
      m <- tktoplevel()
      tkwm.title(m, "Cleavogram Management")
      left <- tkframe(m, relief = "groove", borderwidth = 4)
      right <- tkframe(m, relief = "groove", borderwidth = 4)
      #Menu
      topMenu <- tkmenu(m)
      tkconfigure(m,menu=topMenu)
      fileMenu <- tkmenu(topMenu,tearoff=FALSE)
      analysisMenu <- tkmenu(topMenu, tearoff = FALSE)
      zoomMenu <- tkmenu(topMenu, tearoff = FALSE)
      zoomplus <- tkmenu(topMenu, tearoff = FALSE)
      zoomminus <- tkmenu(topMenu, tearoff = FALSE)
      tkadd(fileMenu,"command",label="Choose cleavogram...",command=function() {
                                                          auxen <- tclvalue(cvobj)
                                                          selectcv(entryCleav)
                                                          if(auxen != tclvalue(cvobj)){
                                                              usrCoords <<- NULL
                                                              xl <<- NULL
                                                              yl <<- NULL
                                                              gr <<- 0
                                                              targetseg <<- NULL
                                                              loadcombo(combosp, tclvalue(cvobj))
                                                              tkrreplot(cleavogram)
                                                              tkrreplot(miniplot)}
                                                          })
      tkadd(fileMenu,"command",label="Quit",command=function() tkdestroy(m))
      tkadd(topMenu,"cascade",label="Data",menu=fileMenu)
      tkadd(analysisMenu, "command", label = "Filter by criteria", command = function() tkrreplot(cleavogram))
      tkadd(analysisMenu, "command", label = "Full report", command = showresult)
      tkadd(analysisMenu, "command", label = "Flat partitions", command = flatpartition)
      tkadd(analysisMenu, "command", label = "View NAM partitions", command = viewsp)
      tkadd(topMenu, "cascade", label="Analysis",menu=analysisMenu)
      zoomMenu <- tkmenu(topMenu,tearoff=FALSE)
      tkadd(zoomMenu,"checkbutton",label="10%",variable=rbValue,onvalue="10%")
      tkadd(zoomMenu,"checkbutton",label="25%",variable=rbValue,onvalue="25%")
      tkadd(zoomMenu,"checkbutton",label="50%",variable=rbValue,onvalue="50%")
      tkadd(zoomMenu,"checkbutton",label="75%",variable=rbValue,onvalue="75%")
      tkadd(zoomMenu,"checkbutton",label="100%",variable=rbValue,onvalue="100%")
      tkadd(topMenu, "cascade", label="Zoom", menu=zoomMenu)
      #left frame
      IOFrame <- tkframe(left, relief = "groove", borderwidth = 2)
      entryCleav <- tkentry(IOFrame, textvariable = cvobj)
      tkgrid(tklabel(left, text = " Input ",  font = "Times 14", foreground = "blue"), columnspan = 2)
      tkgrid(tklabel(IOFrame, text = "Input Cleavogram Data : ", foreground = "red"), entryCleav)
      tkgrid.configure(entryCleav, sticky = "w")
      tkgrid(IOFrame, columnspan = 2)
      eccenFrame <- tkframe(left, relief = "groove", borderwidth = 2)
      tkgrid(tklabel(eccenFrame, text = "Eccentricity\n  (Upper Bounds)  ",
             font = "Times 12", foreground = "blue"))
      eccenMin <- ttkcombobox(eccenFrame, values = c("1", "2", "3", "NONE"),
                  textvariable = mineccen, state = "readonly", width = 10)
      eccenMax <- ttkcombobox(eccenFrame, values = c("1", "2", "3", "NONE"),
                  textvariable = maxeccen, state = "readonly", width = 10)
      eccenMean <- tkentry(eccenFrame, textvariable = meaneccen, width = 10)
      tkgrid(tklabel(eccenFrame, text = "Minimum",
             font = "Times 11", foreground = "red"))
      tkgrid(eccenMin)
      tkgrid(tklabel(eccenFrame, text = "Maximum\n(Graph Diameter)",
             font = "Times 11", foreground = "red"))
      tkgrid(eccenMax)
      tkgrid(tklabel(eccenFrame, text = "Mean",
             font = "Times 11", foreground = "red"))
      tkgrid(eccenMean)
      ratioFrame <- tkframe(left, relief = "groove", borderwidth = 2)
      transSlider <- tkscale(ratioFrame, from = 0.00, to = 1.00, resolution = 0.01, orient = "horiz",
                     variable = transratio)
      densSlider <- tkscale(ratioFrame, from = 0.00, to = 1.00, resolution = 0.01, orient = "horiz",
                    variable = dens)
      tkgrid(tklabel(ratioFrame, text = "Linking Ratios\n(Lower Bounds)",
             font = "Times 12", foreground = "blue"))
      tkgrid(tklabel(ratioFrame, text = "Global Clustering\nCoeffic.(Transitivity)",
             font = "Times 11", foreground = "red"))
      tkgrid(transSlider)
      tkgrid(tklabel(ratioFrame, text = "Graph Density", font = "Times 11", foreground = "red"))
      tkgrid(densSlider)
      tkgrid(tklabel(left, text = "Cohesiveness Criteria",  font = "Times 14", foreground = "blue"),
             columnspan = 2)
      tkgrid(eccenFrame, ratioFrame)
      tkgrid(tklabel(left, text = "Cleavogram Panel",  font = "Times 14",
             foreground = "blue"), columnspan = 2)
      controlFrame <- tkframe(left, relief = "groove", borderwidth = 2)
      rightButton <- tkbutton(controlFrame, text = " > ", font = "Times 14", foreground = "red",
                     pady = 3, command = function() advance(1))
      leftButton <- tkbutton(controlFrame, text = " < ", font = "Times 14", foreground = "red",
                    pady = 3, command = function() advance(2))
      upButton <- tkbutton(controlFrame, text = "/\\", font = "Times 14", foreground = "red",
                  padx = 6, pady = 3, command = function() advance(3))
      downButton <- tkbutton(controlFrame, text = "\\/", font = "Times 14", foreground = "red",
                    padx = 6, pady = 3, command = function() advance(4))
      redrawButton <- tkbutton(controlFrame, text = "Re-\nDraw", font = "Times 10", foreground = "red", command = function () {usrCoords <<- NULL; tkrreplot(cleavogram); tkrreplot(miniplot)})
      zoominButton <- tkbutton(controlFrame, text = "Zoom\nIn", font = "Times 10", foreground = "red",
                      command = function() zoom(1))
      zoomoutButton <- tkbutton(controlFrame, text = "Zoom\nOut", font = "Times 10", foreground = "red",
                       command = function() zoom(-1))
      tkgrid(upButton, downButton, leftButton, rightButton, zoomoutButton, zoominButton, redrawButton)
      tkgrid(controlFrame, columnspan = 2)
      miniplFrame <- tkframe(left, relief = "groove", borderwidth = 2)
      miniplot <- tkrplot(miniplFrame, fun = ya.small, 0.7, 0.7)
      tkgrid(miniplot)
      tkgrid(miniplFrame, columnspan = 2, sticky = "w")
      spradioButton <- tkradiobutton(left, text = "See Geographical Map", value = 1,
                       variable = miniwindow, command = function() tkrreplot(miniplot))
      clradioButton <- tkradiobutton(left, text = "See Cleavogram Map", value = 2,
                       variable = miniwindow, command = function() tkrreplot(miniplot))
      tkgrid(spradioButton, clradioButton)
      #right frame
      cleavoFrame <- tkframe(right)
      cleavogram <- tkrplot(cleavoFrame, fun = ya, 1.85, 1.65)
      tkgrid(cleavogram)
      tkconfigure(cleavogram, cursor = "hand2")
      tkbind(cleavogram, "<Button-3>", menurightbutton)
      tkbind(cleavogram, "<B1-Motion>", locategr)
      tkbind(cleavogram, "<Button-1>", locategr)
      tkbind(cleavogram, "<Double-ButtonPress-1>", loadbranch)
      tkgrid(cleavoFrame)
      focusFrame <- tkframe(right)
      combosp <- ttkcombobox(focusFrame, width = 20,
                textvariable = focussp, state = "readonly")
      canvas <- tkcanvas(focusFrame, width="80", height="25", bg=color)
      colorButton <- tkbutton(focusFrame, text="Set Line\nColor",command = function() {
                              auxcol <- color
                              tkconfigure(colorButton, state = "disabled")
                              color <<- tclvalue(tcl("tk_chooseColor", title="Choose a color"))
                              if (nchar(color)==0) color <<- auxcol
                              tkconfigure(colorButton, state = "active")
                              tkconfigure(canvas,bg=color)
                              tkrreplot(cleavogram)
                              if(tclvalue(miniwindow) == "2")
                                  tkrreplot(miniplot)
                              })
      lwidth <- tkscale(focusFrame, from = 1.00, to = 5.00, resolution = 0.1, orient = "horiz",
                        variable = linewidth, label = "Set Line Width:", command = function(...){
                                                            tkrreplot(cleavogram)
                                                            if(tclvalue(miniwindow) == "2")
                                                                tkrreplot(miniplot)
                        })
      setButton <- tkbutton(focusFrame, text = "  Add  \n Line ", command = locatesp)
      refButton <- tkbutton(focusFrame, text = " Clear \nLines", command = function() {
                                                          color <<- "red"
                                                          tkconfigure(canvas,bg=color)
                                                          xl <<- NULL
                                                          yl <<- NULL
                                                          tclvalue(linewidth) <- "1.5"
                                                          tkrreplot(cleavogram)
                                                          if(tclvalue(miniwindow) == "2")
                                                            tkrreplot(miniplot)
                                                          })
      histButton <- tkbutton(focusFrame, text = "Current\nLines", command = function(){
                            aux <- as.character(tclvalue(cvobj))
                            if (aux %in% ls(globalenv())) cleavo <- get(aux) else return()
                            if (!is(cleavo, "cleavogram")) return()
                            if(is.null(yl)) return()
                            actlines <- unique(yl[which(is.na(yl)) - 1])
                            linehist <- tktoplevel()
                            tkwm.title(linehist, "List of added lines")
                            xscr <- tkscrollbar(linehist, repeatinterval=5,orient="horizontal",
                                                  command=function(...)tkxview(txt,...))
                            yscr <- tkscrollbar(linehist, repeatinterval=5,
                                                  command=function(...)tkyview(txt,...))
                            txt <- tktext(linehist,bg="white",font="courier", width = 30, height = 10,
                                xscrollcommand=function(...)tkset(xscr,...),
                                yscrollcommand=function(...)tkset(yscr,...),
                                wrap="none")
                            tkgrid(txt, yscr)
                            tkgrid.configure(yscr, sticky = "wns")
                            tkgrid(xscr)
                            tkgrid.configure(xscr, sticky = "ewn")
                            for(i in actlines) tkinsert(txt, "end", paste(cleavo$leaves[i], "\n"))
                            })
      tkgrid(tklabel(focusFrame, foreground = "blue",
            text = "Select a Traversing Line\n(Track a single node on the cleavogram)"),
            combosp, setButton, refButton, histButton, canvas, colorButton, lwidth)
      tkgrid(focusFrame)
      tkgrid(left, right)
      tkgrid.configure(right, sticky = "n")
      tkbind(m, "<Destroy>", function() tkdestroy(m))
      tkbind(m, "<KeyPress-Escape>", function() tkdestroy(m))
}
