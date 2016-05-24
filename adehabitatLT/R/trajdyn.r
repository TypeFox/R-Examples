## inspired from functions of the package sudoku

"trajdyn" <- function(x, burst = attr(x[[1]],"burst"), hscale=1,
                      vscale=1, recycle = TRUE,
                      display = c("guess", "windows", "tk"), ...)
{
    ## checks that x is ltraj
    if (!inherits(x, "ltraj"))
        stop("x should be of class 'ltraj'")

    ## creates the environment where we will store the tk variables
    e1 <- new.env(parent = baseenv())

    ## type of trajectories
    typeII <- attr(x,"typeII")

    ## supprimer les NA
    x <- lapply(x, function(i) {
        jj <- i[!is.na(i$x),]
        attr(jj, "id") <- attr(i,"id")
        attr(jj, "burst") <- attr(i,"burst")
        return(jj)
    })
    class(x) <- c("ltraj","list")
    attr(x, "typeII") <- typeII
    attr(x, "regular") <- is.regular(x)

    ## now, assign all the required variables in environment e1
    u <- x
    assign("x", x[burst = burst], envir=e1)
    assign("v", x[burst = burst], envir=e1)
    assign("ajouli", FALSE, envir=e1)
    assign("ajoupo", FALSE, envir=e1)
    assign("ajoubu", FALSE, envir=e1)
    assign("addpoints", TRUE, envir=e1)
    assign("addlines", TRUE, envir=e1)
    assign("lim", TRUE, envir=e1)
    assign("buadd", burst, envir=e1)
    assign("K",1, envir=e1)
    assign("N",nrow(get("x", envir=e1)[[1]]), envir=e1)
    assign("cusr", rep(0 + NA, 4), envir=e1)
    assign("cplt", rep(0 + NA, 4), envir=e1)

    ## remove all unnecessary warnings
    opt <- options(warn=-1)
    on.exit(options(opt))

    ## which kind of display: windows or tcl/tk?
    dsp <- substring(match.arg(display), 1, 1)
    if (dsp == "g")
        dsp <- switch(.Platform$OS.type, windows = "w", "t")
    if (dsp == "t" && !requireNamespace("tkrplot", quietly = TRUE))
        stop("'tkrplot' package needed\n")

    if (dsp == "t")
        assign("hoho", 1, envir=e1)


    ## fonction replot de base (used for plotting the trajectories)
    replot <- function() {
        ## no margins
        opar <- par(mar=c(0,0,0,0), bg="white")

        ## get the burst
        tmptmp <- get("x", envir=e1)
        attr(tmptmp[[1]], "id") <- " "
        assign("x", tmptmp, envir=e1)

        ## and the limits
        if (get("lim", envir=e1)) {
            assign("xlim", range(get("x", envir=e1)[[1]]$x), envir=e1)
            assign("ylim", range(get("x", envir=e1)[[1]]$y), envir=e1)
        }

        ## plot.ltraj
        plot(get("x", envir=e1), id = attr(get("x", envir=e1)[[1]],"id"),
             addlines=FALSE, addp=FALSE, final=FALSE,
             xlim = get("xlim", envir=e1), ylim = get("ylim", envir=e1), ...)

        ## and assign the par values
        assign("cusr", par("usr"), envir=e1)
        assign("cplt", par("plt"), envir=e1)

        ## show the date
        scatterutil.sub(as.character(get("x", envir=e1)[[1]]$date[get("K", envir=e1)]),
                        1, "topleft")

        ## should we add other bursts?
        if (get("ajoubu", envir=e1)) {
            lapply(u[burst=get("buadd", envir=e1)], function(zz) {
                if (get("addpoints", envir=e1))
                    points(zz[,c("x","y")], pch=16, col="grey")
                if (get("addlines", envir=e1))
                    lines(zz[,c("x","y")], pch=16, col="grey")})
        }

        ## add points and lines?
        if (get("addpoints", envir=e1))
            points(get("x", envir=e1)[[1]][1:get("K", envir=e1),c("x","y")], pch=16)
        if (get("addlines", envir=e1))
            if (get("K", envir=e1)>1)
                lines(get("x", envir=e1)[[1]][1:get("K", envir=e1),c("x","y")], lwd=2)

        ## If distance is measured (left click)
        if (get("ajouli", envir=e1))
            lines(c(get("a1", envir=e1)[1],
                    get("a2", envir=e1)[1]),c(get("a1", envir=e1)[2],
                              get("a2", envir=e1)[2]), lwd=2, col="red")

        ## if relocation is identified (right click)
        if (get("ajoupo", envir=e1))
            points(get("a5", envir=e1)[1],
                   get("a5", envir=e1)[2], pch=16, col="red", cex=1.7)

        ## current relocation (end of the displayed trajectory
        iti <- unlist(get("x", envir=e1)[[1]][get("K", envir=e1),c("x","y")])
        points(iti[1],iti[2], col="blue", pch=16, cex=1.4)
        par(opar)
    }

    ## text to display to help the user
    help.txt <- paste("\n-------- to obtain this help, type 'h' ------------------",
                      "n/p            -- Next/Previous relocation",
                      "a              -- show all relocations",
                      "g              -- Go to...",
                      "0-9            -- show a given part of the path",
                      "b              -- change Burst",
                      "i              -- add/remove other bursts on the graph",
                      "z/o            -- Zoom in/Out",
                      "Left-Click     -- measure the distance between two points",
                      "Right-Click    -- identify a relocation",
                      "r/l            -- add or remove points/Lines",
                      "q              -- Quit",
                      "---------------------------------------------------------",
                      "\n", sep = "\n")
    assign("D",0, envir=e1)
    assign("a1", 0, envir=e1)
    assign("a2", 0, envir=e1)

    ## in case the display is tcltk
    if (dsp == "t") {

        ## initialize the tcltk environment
        tt <- tcltk::tktoplevel()
        tcltk::tkwm.title(tt, "Exploration of Animal Movements")
        img <- tkrplot::tkrplot(tt, replot, hscale = hscale, vscale = vscale)
        txt <- tcltk::tktext(tt, bg = "white", font = "courier 10")
        scr <- tcltk::tkscrollbar(tt, repeatinterval = 5,
                                  command = function(...) tcltk::tkyview(txt, ...))
        tcltk::tkconfigure(txt, yscrollcommand = function(...) tcltk::tkset(scr, ...))
        tcltk::tkpack(img, side = "top")
        tcltk::tkpack(txt, side = "left", fill = "both", expand = TRUE)
        tcltk::tkpack(scr, side = "right", fill = "y")

        iw <- as.numeric(tcltk::tcl("image", "width", tcltk::tkcget(img, "-image")))
        ih <- as.numeric(tcltk::tcl("image", "height", tcltk::tkcget(img, "-image")))
    }

    ## afunction for display of the results
    showz <- function() switch(dsp, w = replot(),
                               t = {tkrplot::tkrreplot(img)})

    ## show the help
    type <- function(s) switch(dsp, w = cat(s), t = {
      tcltk::tkinsert(txt, "end", s)
      tcltk::tksee(txt, "end")
    })
    type(help.txt)

    ## to allow the conversion user coordinates/trajectory coordinates
    cc <- function(x, y) {
      if (dsp == "t") {
        x <- (as.double(x) - 1)/iw
        y <- 1 - (as.double(y) - 1)/ih
      }
      px <- (x - get("cplt", envir=e1)[1])/(get("cplt", envir=e1)[2] - get("cplt", envir=e1)[1])
      py <- (y - get("cplt", envir=e1)[3])/(get("cplt", envir=e1)[4] - get("cplt", envir=e1)[3])
      ux <- px * (get("cusr", envir=e1)[2] - get("cusr", envir=e1)[1]) + get("cusr", envir=e1)[1]
      uy <- py * (get("cusr", envir=e1)[4] - get("cusr", envir=e1)[3]) + get("cusr", envir=e1)[3]
      c(ux,uy)
    }

    ## what to do with the mouse on windows
    mm.w <- function(buttons, x, y) {

        ## if a distance is measured
        if (buttons == 0) {

            ## if first click, store the coordinates of the point
            i<-get("D", envir=e1)
            if (i == 0) {
                assign("a1",  cc(x,y), envir=e1)
                assign("D", 1, envir=e1)
            }

            ## if second click, calculate the distance between the two points and display the line
            if (i == 1) {
                assign("a2", cc(x,y), envir=e1)
                assign("D", 0, envir=e1)
                di <- sqrt(sum((get("a2", envir=e1)-get("a1", envir=e1))^2))
                cat(paste("distance:",round(di,6),"\n"))
                lines(c(get("a1", envir=e1)[1],get("a2", envir=e1)[1]),c(get("a1", envir=e1)[2],get("a2", envir=e1)[2]), lwd=2, col="red")
            }
            return()
        }

        ## if a relocation is identified, print info and display the relocation
        if (buttons == 2) {
            w <- get("v",envir=e1)[[1]][1:get("K", envir=e1),]
            assign("a3", cc(x,y), envir=e1)
            di <- sqrt((w$x-get("a3", envir=e1)[1])^2 + (w$y-get("a3", envir=e1)[2])^2)
            print(w[which.min(di),])
            cat("\n")
            points(w[which.min(di),c("x","y")], pch=16, col="red", cex=1.7)
        return()
        }
    }

    ### mouse on tcl/tk for distance calculation
    mm.t <- function(x, y) {

        ## first click, store the coordinates
        i<-get("D", envir=e1)
        if (i == 0) {
            assign("a1", cc(x,y), envir=e1)
            assign("D", 1, envir=e1)
        }

        ## second click, show the line and compute the distance
        if (i == 1) {
            assign("a2", cc(x,y), envir=e1)
            assign("D", 0, envir=e1)
            di <- sqrt(sum((get("a2", envir=e1)-get("a1", envir=e1))^2))
            type(paste("distance:",di,"\n"))
            assign("ajouli", TRUE, envir=e1)
            showz()
            assign("ajouli", FALSE, envir=e1)
        }
        return()
    }


    ### mouse on tcl/tk for relocation identification
    mm.t2 <- function(x, y) {
        w <- get("v",envir=e1)[[1]][1:get("K", envir=e1),]
        assign("a3", cc(x,y), envir=e1)
        di <- sqrt((w$x-get("a3", envir=e1)[1])^2 + (w$y-get("a3", envir=e1)[2])^2)
        assign("a5", unlist(w[which.min(di),c("x","y")]), envir=e1)
        assign("ajoupo", TRUE, envir=e1)
        showz()
        assign("ajoupo", FALSE, envir=e1)
        tmp <- w[which.min(di),]
        se <-unlist(lapply((max(nchar(names(tmp))+
                                nchar(sapply(tmp,as.character))+1) -
                            nchar(names(tmp))-nchar(sapply(tmp, as.character))),
                           function(zz) paste(rep(" ",zz),
                                              collapse="")))
        so<-unlist(lapply(1:length(tmp),
                          function(i) paste(paste(names(tmp)[i],
                                                  as.character(tmp[1,i]),
                                                  sep = se[i]),"\n")))
        type(paste("Relocation",row.names(w)[which.min(di)],":\n"))
        sapply(so,type)
        type("\n")
        return()
    }

    ## current location of the mouse
    mm.mouse <- function(buttons, x, y) {
        assign("a8", cc(x,y), envir=e1)
        return()
    }
    mm.mouset <- function(x, y) {
        assign("a8", cc(x,y), envir=e1)
        return()
    }


    ## keyboard interaction
    kb <- function(A) {
        key <- tolower(A)

        ## quit
        if (key == "q") {
            if (dsp=="t")
                tcltk::tkdestroy(tt)
            return("OK - Finished")
        }

        ## parts of the trajectories
        if (key %in% c(0:9)) {
            if (key > 0)
                assign("K", round(seq(1,get("N", envir=e1),length=11))[as.numeric(key)+1], envir=e1)
            if (key == 0)
                assign("K", 1, envir=e1)
            showz()
        }

        ## zoom
        if (key == "z") {
            assign("tmppx", (get("cusr", envir=e1)[1:2]-
                             get("cusr", envir=e1)[1])/2, envir=e1)

            assign("xlim",  c((get("a8", envir=e1)[1] -
                               (get("tmppx", envir=e1)[2] -
                                get("tmppx", envir=e1)[1])/2),
                              (get("a8", envir=e1)[1] +
                               (get("tmppx", envir=e1)[2] -
                                get("tmppx", envir=e1)[1])/2)), envir=e1)

            assign("tmppy", (get("cusr", envir=e1)[3:4]-get("cusr", envir=e1)[3])/2, envir=e1)
            assign("ylim", c((get("a8", envir=e1)[2] - (get("tmppy",envir=e1)[2] - get("tmppy",envir=e1)[1])/2),
                             (get("a8", envir=e1)[2] + (get("tmppy",envir=e1)[2] - get("tmppy",envir=e1)[1])/2)), envir=e1)

            assign("lim", FALSE, envir=e1)
            showz()
        }

        ## unzoom
        if (key == "o") {
            assign("lim", TRUE, envir=e1)
            showz()
        }

        ## next relocation
        if (key == "n") {
            if (get("K", envir=e1)<=get("N", envir=e1))
                assign("K", get("K", envir=e1)+1, envir=e1)
            if (get("K", envir=e1)>get("N", envir=e1)) {
                if (recycle) assign("K", 1, envir=e1)
                if (!recycle) {
                    assign("K", get("N", envir=e1), envir=e1)
                    cat("End of burst !\n")
                }
            }
            showz()
        }

        ## add lines
        if (key == "l") {
            assign("addlines", !get("addlines", envir=e1), envir=e1)
            showz()
        }

        ## go to...
        if (key == "g") {
            if (dsp == "w") {
                recom <- TRUE
                while (recom) {
                    rr <- readline("Enter a relocation number: ")
                    recom <- FALSE
                    if (!(rr%in%row.names(get("x", envir=e1)[[1]]))) {
                        cat("invalid number\n")
                        recom <- TRUE
                    }
                }
                assign("K", which(row.names(get("x", envir=e1)[[1]])==as.numeric(rr)), envir=e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(row.names(get("x", envir=e1)[[1]])[1])
                tu <- tcltk::tktoplevel(tt, width=500, height=50)
                tcltk::tkwm.title(tu, "Enter a relocation number")
                tcltk::tkwm.resizable(tu, 0, 0)
                en <- tcltk::tkentry(tu, textvariable=lv, width=50)
                submit.but <- tcltk::tkbutton(tu, text="    OK     ",
                                              command=function() {
                                                  rr <- tcltk::tclvalue(lv)
                                                  if (!(rr%in%row.names(get("x", envir=e1)[[1]]))) {
                                                      tcltk::tkmessageBox(message="invalid number",
                                                                          type="ok")
                                                  } else {
                                                      assign("K", which(row.names(get("x", envir=e1)[[1]])==as.numeric(rr)), envir=e1)
                                                      showz()
                                                      tcltk::tkdestroy(tu)}})
                tcltk::tkpack(en, side = "top", fill = "both")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
            }
        }

        ## add relocations
        if (key == "r") {
            assign("addpoints", !get("addpoints", envir=e1), envir=e1)
            showz()
        }

        ## change burst burst
        if (key == "b") {
            assign("K", 1, envir=e1)
            if (dsp == "w") {
                assign("hoho", select.list(unlist(lapply(u, function(y) attr(y, "burst")))), envir=e1)
                type(paste("Choice of the burst:", get("hoho", envir=e1),"\n\n"))
                assign("x",u[burst=get("hoho", envir=e1)], envir=e1)
                assign("v",u[burst=get("hoho", envir=e1)], envir=e1)
                assign("N", nrow(get("x", envir=e1)[[1]]), envir=e1)
                showz()
            }
            if (dsp == "t") {
                lv <- tcltk::tclVar(unlist(lapply(u, function(y) attr(y, "burst"))))
                bubu <- unlist(lapply(u, function(y) attr(y, "burst")))
                tu <- tcltk::tktoplevel(tt)
                tcltk::tkwm.title(tu, "Choose a burst of relocations")
                tcltk::tkwm.resizable(tu, 0, 0)
                tfr <- tcltk::tkframe(tu)
                tli <- tcltk::tklistbox(tfr, bg = "white", font = "courier 12",
                                        listvariable = lv)
                scr2 <- tcltk::tkscrollbar(tfr, repeatinterval = 5,
                                           command = function(...) tcltk::tkyview(tli, ...))
                tcltk::tkconfigure(tli, yscrollcommand = function(...) tcltk::tkset(scr2, ...))
                submit.but <- tcltk::tkbutton(tu, text="    OK     ",
                                              command=function() {
                                                  assign("hoho", ifelse(nchar(tcltk::tclvalue(tcltk::tkcurselection(tli)))==0, 1,
                                                                        as.numeric(tcltk::tclvalue(tcltk::tkcurselection(tli)))+1), envir=e1)
                                                  type(paste("Choice of the burst:", bubu[get("hoho", envir=e1)],"\n\n"))
                                                  tcltk::tkdestroy(tu)})
                tcltk::tkpack(tli, side = "left", fill = "both", expand = TRUE)
                tcltk::tkpack(scr2, side = "right", fill = "y")
                tcltk::tkpack(tfr, side = "right", fill = "y")
                tcltk::tkpack(submit.but, side = "bottom")
                tcltk::tkwait.window(tu)
                assign("x",u[burst=bubu[get("hoho", envir=e1)]], envir=e1)
                assign("v",u[burst=bubu[get("hoho", envir=e1)]], envir=e1)
                assign("N", nrow(get("x", envir=e1)[[1]]), envir=e1)
                showz()
            }
        }

        ## add bursts
        if (key == "i") {
            if (get("ajoubu", envir=e1)) {
                assign("ajoubu", FALSE, envir=e1)
                showz()
            } else {
                if (dsp == "w") {
                    assign("buadd", select.list(unlist(lapply(u,
                                                              function(y) attr(y, "burst"))),
                                                multiple=TRUE), envir=e1)
                    if (length(get("buadd", envir=e1)>0)) {
                        type(paste("show bursts:", paste(get("buadd", envir=e1), collapse=" "),"\n\n"))
                        assign("ajoubu", TRUE, envir=e1)
                        showz()
                    }
                }
                if (dsp == "t") {
                    lv <- tcltk::tclVar(unlist(lapply(u, function(y) attr(y, "burst"))))
                    bubu <- unlist(lapply(u, function(y) attr(y, "burst")))
                    tu <- tcltk::tktoplevel(tt)
                    tcltk::tkwm.title(tu, "Choose one or several bursts")
                    tcltk::tkwm.resizable(tu, 0, 0)
                    tfr <- tcltk::tkframe(tu)
                    tli <- tcltk::tklistbox(tfr, bg = "white", font = "courier 12",
                                            listvariable = lv, selectmode="multiple")
                    scr2 <- tcltk::tkscrollbar(tfr, repeatinterval = 5,
                                               command = function(...) tcltk::tkyview(tli, ...))
                    tcltk::tkconfigure(tli, yscrollcommand = function(...) tcltk::tkset(scr2, ...))
                    submit.but <- tcltk::tkbutton(tu, text="    OK     ",
                                                  command=function() {
                                                      argg <- ifelse(nchar(tcltk::tclvalue(tcltk::tkcurselection(tli)))==0,
                                                                     1,0)
                                                      if (argg==0) {
                                                          assign("ajoubu", TRUE, envir=e1)
                                                          assign("buadd", bubu[as.numeric(unlist(strsplit(tcltk::tclvalue(tcltk::tkcurselection(tli)), " ")))+1], envir=e1)
                                                          type(paste("show bursts:", paste(get("buadd", envir=e1), collapse=" "),"\n\n"))
                                                          showz()
                                                          tcltk::tkdestroy(tu)}})
                    tcltk::tkpack(tli, side = "left", fill = "both", expand = TRUE)
                    tcltk::tkpack(scr2, side = "right", fill = "y")
                    tcltk::tkpack(tfr, side = "right", fill = "y")
                    tcltk::tkpack(submit.but, side = "bottom")
                    tcltk::tkwait.window(tu)
                    assign("x", u[burst=bubu[get("hoho", envir=e1)]], envir=e1)
                    assign("v", u[burst=bubu[get("hoho", envir=e1)]], envir=e1)
                    assign("N", nrow(get("x", envir=e1)[[1]]), envir=e1)
                    showz()
                }
            }
        }

        ## previous relocation
        if (key == "p") {
            if (get("K", envir=e1)>1)
                assign("K", get("K", envir=e1)-1, envir=e1)
            if (get("K", envir=e1)==1) {
                if (recycle)
                    assign("K", get("N", envir=e1), envir=e1)
                if (!recycle) {
                    assign("K", 1, envir=e1)
                    cat("Beginning of burst!\n")
                }
            }
            showz()
        }

        ## all relocations
        if (key == "a") {
            assign("K", get("N", envir=e1), envir=e1)
            showz()
        }

        ## help
        if (key == "h")
            type(help.txt)
        return()
    }

    showz()
    toto <- switch(dsp, w = getGraphicsEvent("", onKeybd = kb, onMouseDown = mm.w,
                        onMouseMove = mm.mouse),
                   t ={tcltk::tkbind(tt, "<Key>", kb)
                       tcltk::tkbind(img, "<Button-1>", mm.t)
                       tcltk::tkbind(img, "<Motion>", mm.mouset)
                       tcltk::tkbind(img, "<Button-3>", mm.t2)
                       tcltk::tkwait.window(tt)})
}
