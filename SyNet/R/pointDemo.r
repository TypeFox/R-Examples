pointDemo <- function ()
{
    #This function can be easily optimized in terms of the number of involved code lines.
    #I just keep this suboptimnal code to make more transparent the different operations. 
    require(tcltk) || stop("tcltk support is absent")
    color <- tclVar(0)
    species <- tclVar(0)
    d <- data.frame()
    runsp1 <- NULL
    runsp2 <- NULL
    l1 <- list(x = NULL, y = NULL)
    l2 <- list(x = NULL, y = NULL)
    plot(0, xlim = c(0, 10), ylim = c(0, 10), axes = F, xlab = "",
        ylab = "", main = "", type = "n")
    mst <- function(x) {
        numpoints <- nrow(x)
        caso <- order(x[lower.tri(x)])
        anc <- vector("integer", numpoints)
        path <- 0
        nanc <- 0
        inmst <- vector("logical", numpoints)
        cl <- ceiling(-0.5 + numpoints - sqrt(0.25 + numpoints*(numpoints-1) - 2*caso))
        rw <- numpoints*(1 - cl) + 0.5*cl*(cl + 1) + caso
        flag <- 0
        i <- 0
        while (flag < (numpoints - 1))
        {
          i <- i + 1
          aux <- 2*inmst[cl[i]] + inmst[rw[i]]
          if(aux == 0) {inmst[c(cl[i], rw[i])] <- TRUE; nanc <- nanc + 1; anc[c(cl[i], rw[i])] <- nanc;  path <- path + x[rw[i], cl[i]]; flag <- flag + 1; next}
          if(aux == 1) {inmst[cl[i]] <- TRUE; anc[cl[i]] <- anc[rw[i]];  path <- path + x[rw[i], cl[i]]; flag <- flag + 1; next }
          if(aux == 2) {inmst[rw[i]] <- TRUE; anc[rw[i]] <- anc[cl[i]];  path <- path + x[rw[i], cl[i]]; flag <- flag + 1; next }
          if(anc[cl[i]] != anc[rw[i]]) {anc[anc == anc[rw[i]]] <- anc[cl[i]]; path <- path + x[rw[i], cl[i]]; flag <- flag + 1 }
        }
        return (path = path)
    }
    drawing <- function(runsp1, runsp2) {
        dat <- data.frame(matrix(rep(0, 3*(runsp1 + runsp2)), ncol = 3))
        names(dat) <- c("IDsp", "x", "y")
        rect(0,0,10,1, col = "gray")
        text(5,0.5, "Adding Sp1 records...", col = "red")
        l1 <<- locator(runsp1, type = "p", col = "red", pch = 19)
        rect(0,0,10,1, col = "gray")
        text(5,0.5, "Adding Sp2 records...", bg = "white", col = "green")
        l2 <<- locator(runsp2, type = "p", col = "green", pch = 19)
        dat[,1] <- c(rep("sp1", runsp1), rep("sp2", runsp2))
        dat[,2] <- c(l1$x, l2$x)
        dat[,3] <- c(l1$y, l2$y)
        dat
    }
    movept <- function(move = NULL) {
        aux <- as.integer(tclvalue(species))
        if (move == 0) ifelse (aux, l2$y <<- l2$y + 0.2, l1$y <<- l1$y + 0.2)
        if (move == 1) ifelse (aux, l2$x <<- l2$x + 0.2, l1$x <<- l1$x + 0.2)
        if (move == 2) ifelse (aux, l2$y <<- l2$y - 0.2, l1$y <<- l1$y - 0.2)
        if (move == 3) ifelse (aux, l2$x <<- l2$x - 0.2, l1$x <<- l1$x - 0.2)
        opar <- par(mar = c(0,0,3,0))
        on.exit(opar)
        plot(0, xlim = c(0, 10), ylim = c(0, 10), axes = F, xlab = "",
                ylab = "", main = "", type = "n")
        rect(0,0,10,1, col = gray(0.7))
        if (aux)
            text(5,0.5, "Moving Sp2 records...", bg = "white", col = "green")
        else 
            text(5,0.5, "Moving Sp1 records...", bg = "white", col = "red")
        points(l1, col = "red", pch = 19)
        points(l2, col = "green", pch = 19)        
        d[,2] <<- c(l1$x, l2$x)
        d[,3] <<- c(l1$y, l2$y)    
        sp1l <- mst(as.matrix(dist(d[(1:runsp1),2:3]))) #It is redundant.
        sp2l <- mst(as.matrix(dist(d[(runsp1+1):(runsp1 + runsp2),2:3]))) #It is redundant. 
        sp12l <- mst(as.matrix(dist(d[,2:3])))
        mtext(paste("MST(Sp1): ", round(sp1l, 4)), line = 2, col = "red")
        mtext(paste("MST(Sp2): ", round(sp2l, 4)), line = 1, col = "green")
        mtext(paste("MST(Sp1 U Sp2): ", round(sp12l, 4)), line = 0)
    }
    rotate <- function(ckw) {
        aux <- as.integer(tclvalue(species))
        opar <- par(mar = c(0,0,3,0))
        on.exit(opar)
        plot(0, xlim = c(0, 10), ylim = c(0, 10), axes = F, xlab = "",
                ylab = "", main = "", type = "n")
        rect(0,0,10,1, col = gray(0.7))
        centrl1 <- sapply(l1, mean)
        centrl2 <- sapply(l2, mean)
        auxx1 <- l1$x - centrl1[1]
        auxy1 <- l1$y - centrl1[2]
        auxx2 <- l2$x - centrl2[1]
        auxy2 <- l2$y - centrl2[2]
        r1 <- sqrt((l1$x - centrl1[1])^2 + (l1$y - centrl1[2])^2)
        r2 <- sqrt((l2$x - centrl2[1])^2 + (l2$y - centrl2[2])^2)
        if (aux) {
            text(5,0.5, "Rotating Sp2 records...", bg = "white", col = "green")
            l2$x <<- auxx2*cos(ckw*1/6*pi) + auxy2*sin(ckw*1/6*pi) + centrl2[1]
            l2$y <<- -1*auxx2*sin(ckw*1/6*pi) + auxy2*cos(ckw*1/6*pi) + centrl2[2]
        }
        else { 
            text(5,0.5, "Rotating Sp1 records...", bg = "white", col = "red")
            l1$x <<- auxx1*cos(ckw*1/6*pi) + auxy1*sin(ckw*1/6*pi) + centrl1[1]
            l1$y <<- -1*auxx1*sin(ckw*1/6*pi) + auxy1*cos(ckw*1/6*pi) + centrl1[2]
        }
        points(l1, col = "red", pch = 19)
        points(l2, col = "green", pch = 19)        
        d[,2] <<- c(l1$x, l2$x)
        d[,3] <<- c(l1$y, l2$y)    
        sp1l <- mst(as.matrix(dist(d[(1:runsp1),2:3]))) #It is redundant.
        sp2l <- mst(as.matrix(dist(d[(runsp1+1):(runsp1 + runsp2),2:3]))) #It is redundant.
        sp12l <- mst(as.matrix(dist(d[,2:3])))
        mtext(paste("MST(Sp1): ", round(sp1l, 4)), line = 2, col = "red")
        mtext(paste("MST(Sp2): ", round(sp2l, 4)), line = 1, col = "green")
        mtext(paste("MST(Sp1 U Sp2): ", round(sp12l, 4)), line = 0)
    }
    refresh <- function() {
        opar <- par(mar = c(0,0,3,0))
        on.exit(opar)
        plot(0, xlim = c(0, 10), ylim = c(0, 10), axes = F, xlab = "",
            ylab = "", main = "", type = "n")
        runsp1 <<- as.integer(tclvalue(tkget(EntryRun1)))
        if(is.na(runsp1)) return()
        if(runsp1 < 1) return()
        runsp2 <<- as.integer(tclvalue(tkget(EntryRun2)))
        if(is.na(runsp2)) return() 
        if(runsp2 < 1) return()
        d <<- drawing(runsp1 = runsp1, runsp2 = runsp2)
        print(d)
        mtext(paste("MST(Sp1): ", round(mst(as.matrix(dist(d[(1:runsp1),2:3]))), 4)), line = 2, col = "red")
        mtext(paste("MST(Sp2): ", round(mst(as.matrix(dist(d[(runsp1+1):(runsp1 + runsp2),2:3]))), 4)),
              line = 1, col = "green")
        mtext(paste("MST(Sp2 U Sp1): ", round(mst(as.matrix(dist(d[,2:3]))), 4)), line = 0)
    }
    m <- tktoplevel()
    tkwm.title(m, "Comparing species ranges")
    runFrame <- tkframe(m, relief = "groove", borderwidth = 4)
    tkgrid(tklabel(runFrame, text = "Drawing Points", font = "Times 14", foreground = "blue"), 
           columnspan = 3)
    EntryRun1 <- tkentry(runFrame, width = 4, bg = "white", textvariable = tclVar("10"))
    EntryRun2 <- tkentry(runFrame, width = 4, bg = "white", textvariable = tclVar("10"))
    tkgrid(tklabel(runFrame, text = "Size Sp1", foreground = "red", background = "white"), EntryRun1)
    tkgrid(tklabel(runFrame, text = "Size Sp2", foreground = "green", background = "white"), EntryRun2)
    runButton <- tkbutton(runFrame, command = refresh, text = "Draw!")
    tkgrid(runButton, columnspan = 3, ipadx = 20)
    tkgrid(runFrame)
    moveFrame <- tkframe(m, relief = "groove", borderwidth = 4)
    tkgrid(tklabel(moveFrame, text = "Moving Points", font = "Times 14", foreground = "blue"),
           columnspan = 3)
    sp1Button <- tkradiobutton(moveFrame, variable = species, value = 0, text = "Species 1")
    sp2Button <- tkradiobutton(moveFrame, variable = species, value = 1, text = "Species 2")
    tkgrid(tklabel(moveFrame, text= "Acitve set of points", font = "Times 12", foreground = "blue"),
           columnspan = 3, sticky = "w")
    tkgrid(sp1Button)
    tkgrid(sp2Button)
    tkgrid(tklabel(moveFrame, text= "Spatial manager", font = "Times 12", foreground = "blue"),
           columnspan = 3, sticky = "w")
    upButton <- tkbutton(moveFrame, command = function() movept(0), text = " Up ", padx = 20)
    rightButton <- tkbutton(moveFrame, command = function () movept(1), text = "Right", padx = 20)
    downButton <- tkbutton(moveFrame, command = function () movept(2), text = " Down ", padx = 20)
    leftButton <- tkbutton(moveFrame, command = function () movept(3), text = " Left ", padx = 20)
    rotButton1 <-  tkbutton(moveFrame, command = function() rotate(1), text = " 30 CW ", padx = 20)
    rotButton2 <-  tkbutton(moveFrame, command = function() rotate(-1), text = " 30 CCW", padx = 20)
    tkgrid(tklabel(moveFrame,text="<    >"), leftButton, rightButton)
    tkgrid(tklabel(moveFrame,text="Rotate"), rotButton1, rotButton2)
    tkgrid(tklabel(moveFrame,text="\\/  /\\"), upButton, downButton)
    tkgrid(moveFrame)
}

