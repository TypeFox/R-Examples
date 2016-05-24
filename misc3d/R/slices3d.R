if (! exists("gray.colors"))
    gray.colors <- function (n, start = 0.3, end = 0.9, gamma = 2.2)
        gray(seq(from = start^gamma, to = end^gamma, length = n)^(1/gamma))

vslice <- function(vol, which, k, tpt = 1) {
    if (length(dim(vol)) == 4)
        switch(which,
               x = vol[k,,,tpt],
               y = vol[,k,,tpt],
               z = vol[,,k,tpt])
    else
        switch(which,
               x = vol[k,,],
               y = vol[,k,],
               z = vol[,,k])
}

slices3d <- function(vol1, vol2=NULL, rlim1=c(-Inf, Inf), rlim2=NULL,
                     col1=gray.colors(512), col2=NULL,
                     main="Three Planes View", scale = 0.8,
                     alpha=1, cross = TRUE,
                     layout=c("counterclockwise", "clockwise")){

    mkimg <- function(which) {
        switch(which,
               x = { i <- 1; j <- 2; k <- 3 },
               y = { i <- 2; j <- 1; k <- 3 },
               z = { i <- 3; j <- 1; k <- 2 })
        f <- function() {
            opar = par(mar=c(0,0,0,0))
            on.exit(par(opar))
            if(!(is.array(col)))
                image(vslice(vol, which, bb[i],bb[4]), col=col, zlim = rlim1)
            else{
                v <- switch(which,
                           x = matrix(1:(d[2]*d[3]), nrow=d[2]),
                           y = matrix(1:(d[1]*d[3]), nrow=d[1]),
                           z = matrix(1:(d[1]*d[2]), nrow=d[1]))
                image(v, col=vslice(col, which, bb[i],bb[4]))
            }
            lines(rep(bb[j]/d[j],100), seq(0,1,len=100))
            lines(seq(0,1,len=100), rep(bb[k]/d[k],100))
        }
        tkrplot(tt, f, hscale = 0.8, vscale = 0.8)
    }
    mkscale <- function(i) {
        f <- function(...) {
             b <- as.numeric(tclvalue(bbv[[i]]))
             if (b != bb[i]) {
                bb[i] <<- b
                if (cross || i == 4)
                    for (j in 1:3) tkrreplot(img[[j]])
                else  tkrreplot(img[[i]])
                tkconfigure(l2, text=bb[i])
            }
        }
        fr <- tkframe(tt)
        s <- tkscale(fr, command=f, from=1, to=d[i], resolution=1,
                variable=bbv[[i]], showvalue=FALSE, orient="horiz")
        l1 <- tklabel(fr, text = dn[i])
        l2 <- tklabel(fr, textvariable = bbv[[i]])
        tkgrid(l1, s, l2)
        fr
    }
    move <- function(which){
        if(lay=="clockwise"){
            switch(which,
                   x = { i <- 1; j <- 2; k <- 3 },
                   y = { i <- 2; j <- 1; k <- 3 },
                   z = { i <- 3; j <- 1; k <- 2 })
        }
        else{
            switch(which,
                   y = { i <- 1; j <- 2; k <- 3 },
                   x = { i <- 2; j <- 1; k <- 3 },
                   z = { i <- 3; j <- 1; k <- 2 })
        }
        tkbind(img[[i]],"<Button-1>", function(x,y){
            wid <- as.integer(tkwinfo("width",img[[i]]))
            hei <- as.integer(tkwinfo("height",img[[i]]))
            if(lay=="clockwise" || which=="z")
                bb[j] <<- as.numeric(x)/wid*d[j]
            else
                bb[i] <<- as.numeric(x)/wid*d[i]
            bb[k] <<- d[k] - as.numeric(y)/hei*d[k]
                
            
            for (j in 1:3){
                tkrreplot(img[[j]])
                tclvalue(bbv[[j]]) <<- as.character(round(bb[j]))
            }
        })
    }
    
    overlay <- function(vol1, vol2, rlim1, rlim2, col1, col2, alpha){
        choose1 <- vol1 <= rlim1[2] & vol1 >= rlim1[1]
        vol1 <- floor((length(col1) - .01) *
                    (vol1 - min(vol1))/(max(vol1) - min(vol1)) + 1)
        vol1c <- col1[vol1]
        vol1c[!choose1] <- "white"

        choose2 <- vol2 <= rlim2[2] & vol2 >= rlim2[1]
        vol2 <- floor((length(col2) - .01) *
                    (vol2 - min(vol2))/(max(vol2) - min(vol2)) + 1)
        vol2c <- col2[vol2]
        vol2c[!choose2] <- "transparent"
        alpha <- as.vector(ifelse(choose2, alpha, 0))
        col <- t(col2rgb(vol1c)) * (1 - alpha) + t(col2rgb(vol2c)) * alpha
        array(rgb(col, maxColorValue=255), dim=dim(vol1))
    }

    if (! require(tkrplot)) stop("tkrplot is required.");

    if(missing(rlim1))
        rlim1 <- range(vol1,na.rm = TRUE)
    if(is.null(vol2)){
        vol <- vol1
        col <- col1
    }
    else{
        if(!all(dim(vol1 == vol2)))
            stop("two layers have to have the same dimensions")
        if(missing(rlim2))
            rlim2 <- range(vol2,na.rm = TRUE)
        col <- overlay(vol1, vol2, rlim1, rlim2, col1, col2, alpha)
        vol <- array(0, dim=dim(vol1))
    }

    lay <- match.arg(layout)
    layout <- switch(lay, counterclockwise = c(2,1,3), clockwise = c(1,2,3))
    direct <- c("x", "y", "z")
    d <- dim(vol)
    #dn <- c("x", "y", "z", "t")
    dn <- c(direct, "t")
    tt <- tktoplevel()
    tktitle(tt) <- main
    bb <- c(round(d[1:3]) / 2, 1)
    bbv <- lapply(bb, tclVar)
    s <- lapply(layout, mkscale)
    img <- lapply(direct[layout], mkimg)
    tkgrid(img[[1]], img[[2]])
    tkgrid(s[[1]],s[[2]])
    tkgrid(img[[3]])
    if (length(d) == 4 && d[4] > 1)
        tkgrid(s[[3]], mkscale(4))
    else tkgrid(s[[3]])
    lapply(direct[layout], move)

    environment()
}
