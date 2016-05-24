shardsplot <- function(object, plot.type = c("eight", "four", "points", "n"),
    expand = 1, stck = TRUE, grd = FALSE, standardize = FALSE, data.or = NA,
    label = FALSE, plot = TRUE, classes = 0, vertices = TRUE,
    classcolors = "rainbow", wghts = 0, xlab = "Dimension 1", 
    ylab = "Dimension 2", xaxs = "i", yaxs = "i", plot.data.column = NA,
    log.classes = FALSE, revert.colors = FALSE, ...)
{
    diss <- FALSE
    plot.type <- match.arg(plot.type) # "delaunay" not yet implemented
    
    if (plot.type == "delaunay"){
        stop("plot type delaunay not yet implemented") 
        grd <- FALSE
    }
    
    if (class(object)=="som"){
        preimages <- object$code
        nzx <- object$x
        nzy <- object$y
        Z0 <- t(matrix(1:(nzx*nzy),nzx,nzy))
        nobs <- object$code.sum$nobs
        nobs.col <- nobs+1
        maxn <- max(nobs.col)
        cl.ord <- nobs.col
        if (!is.na(plot.data.column)){
            if(length(plot.data.column)!=1 || plot.data.column < 1 || plot.data.column > NCOL(preimages)){
                stop(paste("plot.data.column must be an integer between 1 and ", NCOL(preimages)))
            }
            cl.ord <- preimages[,plot.data.column]
            if(log.classes){
                cl.ord <- log(cl.ord)
            }
            cl.ord <- 99*(cl.ord - min(cl.ord))/diff(range(cl.ord))+1
        }
    }
    else if(class(object) == "EDAM"){
        preimages <- object$preimages
        Z0 <- object$Z
        if (classes[1]==0) cl.ord <- as.numeric(object$cl.ord)
        if (classes[1]!=0) cl.ord <- as.numeric(classes[as.vector(t(object$Z.old.terms))])
        if (ncol(preimages)==nrow(preimages) && preimages==t(preimages)){
            diss <- TRUE
            EV.dist <- preimages
            preimages <- cbind(c(1:nrow(preimages)),c(1:nrow(preimages)))
            rownames(preimages) <- rownames(object$preimages)
        }
    }
    if (standardize){
        EV.sd <- sqrt(apply(preimages,2,var))
        EV.sd[!EV.sd] <- 1
        preimages <- t(t(preimages)/EV.sd)
    }
    if (wghts[1]!=0) preimages <- preimages*kronecker(t(wghts),rep(1,nrow(preimages)))
    nzx <- ncol(Z0)
    nzy <- nrow(Z0)
    nzx.ex <- round(nzx*expand)
    nzy.ex <- round(nzy*expand)
    Cells0 <- cbind(kronecker(c(1:nzy),rep(1,nzx)),rep(c(1:nzx),nzy))
    Cells.ex <- Cells0*expand
    Cells.exldru <- Cells.ex
    Cells.exlurd <- Cells.ex
    EV <- preimages

    if (!diss) EV.dist <- distmirr(dist(EV, method = "euclidean"))
    dist.max <- 0
    dist.min <- max(EV.dist)
    for (i in 1:(nzy*nzx)){
        sub.mat <-  t(matrix(rep(Cells0[i,], nzx*nzy), 2, nzx*nzy))
        rel.dists <- EV.dist[which(diag((Cells0-sub.mat)%*%t(Cells0-sub.mat)) <= 2),i]
        rel.dists <- rel.dists[which(rel.dists > 0)]
        if (max(rel.dists) > dist.max) dist.max <- max(rel.dists)
        if (min(rel.dists) < dist.min) dist.min <- min(rel.dists)
    }
    Cells.ny <- (nzy*expand)
    Cells.nx <- (nzx*expand)
    if (any(expand!=1) || stck){
        for(j in 1:nzx){
            nb.dists <- rep(0, nzy-1)
            for (i in (2:nzy)){
              nb.dists[i-1] <- EV.dist[Z0[i,j],Z0[i-1,j]]
            }
            nb.dc <- sum(nb.dists)
            nb.pos <- c(expand,rep(0,nzy-1))
            for (l in 1:(nzy-1)){
                nb.pos[l+1] <- expand+(sum(nb.dists[1:l])*(Cells.ny-expand))/nb.dc
            }
            Cells.ex[Z0[,j],1] <- nb.pos
        }
        for (i in 1:nzy){
            nb.dists <- rep(0, nzx-1)
            for (j in (2:nzx)){
                nb.dists[j-1] <- EV.dist[Z0[i,j],Z0[i,j-1]]
            }
            nb.dc <- sum(nb.dists)
            nb.pos <- c(expand,rep(0, nzx-1))
            for (l in 1:(nzx-1)){
                nb.pos[l+1] <- expand+(sum(nb.dists[1:l])*(Cells.nx-expand))/nb.dc
            }
            Cells.ex[Z0[i,],2] <- nb.pos
        }
        if (min(nzx,nzy) > 2){
            for (j in 1:(nzy-2)){
                len <- min(nzx,nzy-j+1)
                nb.dists <- rep(0, len-1)
                rel.numbers <- rep(0,len)
                rel.numbers[1] <- Z0[j,1]
                for (i in (2:len)){
                    rel.numbers[i] <- Z0[i+j-1,i]
                    nb.dists[i-1] <- EV.dist[rel.numbers[i],rel.numbers[i-1]]
                }
                nb.dc <- sum(nb.dists)
                nb.pos <- cbind(c(j*expand,rep(0,len-1)),c(expand,rep(0,len-1)))
                for (l in 1 :(len-1)){
                    nb.pos[l+1,1] <- j*expand+(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                    nb.pos[l+1,2] <- expand+(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                }
                Cells.exldru[rel.numbers,] <- nb.pos
            }
            for (i in 2:(nzx-2)){
                len <- min(nzx-i+1,nzy)
                nb.dists <- rep(0, len-1)
                rel.numbers <- rep(0,len)
                rel.numbers[1] <- Z0[1,i]
                for (j in (2:len)){
                    rel.numbers[j] <- Z0[j,j+i-1]
                    nb.dists[j-1] <- EV.dist[rel.numbers[j],rel.numbers[j-1]]
                }
                nb.dc <- sum(nb.dists)
                nb.pos <- cbind(c(expand,rep(0,len-1)),c(i*expand,rep(0,len-1)))
                for (l in 1 :(len-1)){
                    nb.pos[l+1,1] <- expand+(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                    nb.pos[l+1,2] <- i*expand+(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                }
                Cells.exldru[rel.numbers,] <- nb.pos
            }
            for (j in 1:(nzy-2)){
                len <- min(nzx,nzy-j+1)
                nb.dists <- rep(0, len-1)
                rel.numbers <- rep(0,len)
                rel.numbers[1] <- Z0[j,nzx]
                for (i in (2:len)){
                    rel.numbers[i] <- Z0[i+j-1,nzx-i+1]
                    nb.dists[i-1] <- EV.dist[rel.numbers[i],rel.numbers[i-1]]
                }
                nb.dc <- sum(nb.dists)
                nb.pos <- cbind(c(j*expand,rep(0,len-1)),c(Cells.nx,rep(0,len-1)))
                for (l in 1 :(len-1)){
                    nb.pos[l+1,1] <- j*expand+(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                    nb.pos[l+1,2] <- Cells.nx-(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                }
                Cells.exlurd[rel.numbers,] <- nb.pos
            }
            for (i in 2:(nzx-1)){
                len <- min(i,nzy)
                nb.dists <- rep(0, len-1)
                rel.numbers <- rep(0,len)
                rel.numbers[1] <- Z0[1,i]
                for (j in (2:len)){
                    rel.numbers[j] <- Z0[j,i-j+1]
                    nb.dists[j-1] <- EV.dist[rel.numbers[j],rel.numbers[j-1]]
                }
                nb.dc <- sum(nb.dists)
                nb.pos <- cbind(c(expand,rep(0,len-1)),c(i*expand,rep(0,len-1)))
                for (l in 1:(len-1)){
                    nb.pos[l+1,1] <- expand+(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                    nb.pos[l+1,2] <- i*expand-(sum(nb.dists[1:l])*(len*expand-expand))/nb.dc
                }
                Cells.exlurd[rel.numbers,] <- nb.pos
            }
        }
        if (plot.type!="four") Cells.ex <- (Cells.ex+Cells.exldru+Cells.exlurd)/3
        if (grd) Cells.ex <- round(Cells.ex)
    }
    if (plot){
        plot(c(expand,expand*nzy), c(expand,expand*nzx), type = "n", 
            xlab = xlab, ylab = ylab, xaxs = xaxs, yaxs = yaxs, ...)
        if (vertices && plot.type!="delaunay"){
            for (i in (1:nzx)){
                for (j in (1:nzy)){
                    if (j>1) lines(Cells.ex[Z0[c(j,j-1),c(i,i)],], col="gray")
                    if (i>1) lines(Cells.ex[Z0[c(j,j),c(i,i-1)],], col="gray")
                }
            }
        }
    }
    if (class(object)!="som" | (class(object)=="som" & !is.na(plot.data.column))){
        vec.col.ord <- if(is.character(classcolors) && length(classcolors)==1){
            temp.col <- switch(classcolors,
                "rainbow" = rainbow(max(cl.ord))[cl.ord],
                "topo"    = topo.colors(max(cl.ord))[cl.ord],
                "gray"    = gray(1:max(cl.ord)/max(cl.ord))[cl.ord],
                stop("argument classcolors only support 'rainbow', 'topo', and 'gray'")
            )
            } else classcolors[cl.ord]
    } else {
        if (any(is.na(data.or))){
            vec.col.ord <- c(rgb(1,1,1),
                hsv(1,1, ((maxn*1.5):1) / (maxn*1.5)))[nobs.col]
        } else {
            code.classes <- rep(1,nzx*nzy)
            i <- 0
            nob <- nrow(data.or)
            dimen <- ncol(data.or)-1
            obj.codes <- object$visual$x+object$visual$y*nzx+1
            class.vec <- as.factor(data.or[,dimen+1])
            for (i in 1:(nzx*nzy)){
                class.table <- table(class.vec[(obj.codes==i)])
                if (length(class.table)>0) code.classes[i] <- which(class.table==max(class.table))[1]
            }
            vec.col.ord1 <-
                if(is.character(classcolors) && length(classcolors)==1){
                    temp.col<-  switch(classcolors,
                        "rainbow" = rainbow(max(code.classes))[code.classes],
                        "topo"    = topo.colors(max(code.classes))[code.classes],
                        "gray"    = gray(1:max(code.classes)/max(code.classes))[code.classes],
                        stop("argument classcolors only support 'rainbow', 'topo', and 'gray'")
                    )
                    if(revert.colors) rev(temp.col) else temp.col
                } else classcolors[cl.ord]
            vec.col.ord <- rep(0,nzx*nzy)
            rgb.mat <- col2rgb(vec.col.ord1)/255
            for (i in 1:(nzx*nzy)){
                rgb.act <- rgb.mat[,i]
                V.act <- max(rgb.act)
                rgb.min <- min(rgb.act)
                S.act <- 0
                if (V.act>0) S.act <- (V.act-rgb.min)/V.act
                H.act <- 0
                if (S.act>0){
                    if (V.act==rgb.act[1]) H.act <- (rgb.act[2]-rgb.act[3])/(V.act-rgb.min)
                    else if (V.act==rgb.act[2]) H.act <- 2 + (rgb.act[3]-rgb.act[1])/(V.act-rgb.min)
                    else if (V.act==rgb.act[3]) H.act <- 4 + (rgb.act[1]-rgb.act[2])/(V.act-rgb.min)
                    H.act <- H.act*60
                    if (H.act<0) H.act <- H.act+360
                    H.act <- H.act/360
                }
                S.act <- (maxn*1.5-nobs[i])/(maxn*1.5)
                vec.col.ord[i] <- hsv(H.act, V.act, S.act)
            }
            vec.col.ord[!nobs] <- "#FFFFFF"
        }
    }
    for (i in 1:(nzy*nzx)){
        sub.mat <-  t(matrix(rep(Cells0[i,], nzx*nzy), 2, nzx*nzy))
        Cells0msubmat <- Cells0 - sub.mat
        rel.numbers <- which(diag(Cells0msubmat %*% t(Cells0msubmat)) <= 2)
        rel.dists <- EV.dist[rel.numbers, i]
        rel.numbers <- rel.numbers[rel.dists > 0]
        rel.dists <- rel.dists[rel.dists > 0]
        rel.points <- matrix(rep(Cells.ex[i,], length(rel.dists)),
            length(rel.dists), 2, byrow = TRUE)
        rel.points <- rel.points + Cells0msubmat[rel.numbers,] *
            (0.5-((rel.dists-dist.min) / (2*dist.max)))
        rel.coords <- Cells0msubmat[rel.numbers,]
        if(length(rel.dists) == 3){
            rel.points <- rbind(rel.points, Cells0[i,])
            rel.coords <- rbind(rel.coords, c(0,0))
            rel.numbers <- c(rel.numbers, i)
        }
        rel.count <- 0
        circle.mat <- matrix(c(-1, -1, -1, 0, -1, 1, 0, 1, 1, 1, 1, 0, 1, -1, 0, -1),
            8, 2, byrow = TRUE)
        draw.points <- matrix(rep(Cells.ex[i,], 8), 8, 2, byrow = TRUE)
        for(l in 1:8){
            for(j in 1:nrow(rel.coords)){
                if (t(rel.coords[j,]-circle.mat[l,]) %*% (rel.coords[j,]-circle.mat[l,]) == 0){
                    rel.count <- rel.count+1
                    draw.points[l,] <- rel.points[j,]
                }
            }
        }
        if (cl.ord[1] && class(object)!="som"){
            if (Cells0[i,1] > 1 && Cells0[i,2] > 1){
                if (cl.ord[Z0[Cells0[i,1]-1, Cells0[i,2]]] == cl.ord[i] &&
                    cl.ord[Z0[Cells0[i,1], Cells0[i,2]-1]] == cl.ord[i]){
                    if (cl.ord[Z0[Cells0[i,1]-1,Cells0[i,2]-1]] != cl.ord[i]){
                        draw.points[1,] <- c(draw.points[2,1], draw.points[8,2])
                    }
                }
            }
            if (Cells0[i,1] > 1 && Cells0[i,2] < nzx){
                if (cl.ord[Z0[Cells0[i,1]-1,Cells0[i,2]]] == cl.ord[i] &&
                    cl.ord[Z0[Cells0[i,1],Cells0[i,2]+1]] == cl.ord[i]){
                    if (cl.ord[Z0[Cells0[i,1]-1,Cells0[i,2]+1]] != cl.ord[i]){
                        draw.points[3,] <- c(draw.points[2,1], draw.points[4,2])
                    }
                }
            }
            if (Cells0[i,1] < nzy && Cells0[i,2] > 1){
                if (cl.ord[Z0[Cells0[i,1]+1,Cells0[i,2]]] == cl.ord[i] &&
                    cl.ord[Z0[Cells0[i,1],Cells0[i,2]-1]] == cl.ord[i]){
                    if (cl.ord[Z0[Cells0[i,1]+1,Cells0[i,2]-1]] != cl.ord[i]){
                        draw.points[7,] <- c(draw.points[6,1], draw.points[8,2])
                    }
                }
            }
            if (Cells0[i,1] < nzy && Cells0[i,2] < nzx){
                if (cl.ord[Z0[Cells0[i,1]+1, Cells0[i,2]]] == cl.ord[i] &&
                    cl.ord[Z0[Cells0[i,1],Cells0[i,2]+1]] == cl.ord[i]){
                        if (cl.ord[Z0[Cells0[i,1]+1, Cells0[i,2]+1]] != cl.ord[i]){
                            draw.points[5,] <- c(draw.points[6,1], draw.points[4,2])
                    }
                }
            }
        }
        if (plot.type == "four"){
            draw.points[1,] <- c(draw.points[2,1], draw.points[8,2])
            draw.points[3,] <- c(draw.points[2,1], draw.points[4,2])
            draw.points[7,] <- c(draw.points[6,1], draw.points[8,2])
            draw.points[5,] <- c(draw.points[6,1], draw.points[4,2])
        }
        draw.points[(draw.points[,2] < 1),2] <- expand
        draw.points[(draw.points[,1] < 1),1] <- expand
        if (plot && plot.type!="delaunay" && plot.type!="n"){
            if (plot.type!="points") {
                stopifnot(length(vec.col.ord) == length(cl.ord))            
                polygon(draw.points, col = vec.col.ord[i], ...)
                if (label)
                    text(mean(draw.points[,1]), mean(draw.points[,2]),
                        rownames(preimages)[i], ...)
            }
            if (plot.type=="points")
                points(Cells.ex[i,1], Cells.ex[i,2], col=vec.col.ord[i], pch = 19, ...)
        }
        if (label && plot.type!="delaunay")
            text(mean(draw.points[,1]), mean(draw.points[,2]), rownames(preimages)[i], ...)
    }
#    if (plot && plot.type=="delaunay"){
#        lnames <- 0
#        if (label) lnames <- rownames(preimages)
#        cont.shardsplot(coords=Cells.ex, classes=cl.ord, ydistmat=EV.dist,
#            radius=radius, smallest=smallest, percentage=percentage, convexify=convexify,
#            label=label, vertices=vertices, classcolors=classcolors, pnew=FALSE,
#            vec.col.ord=vec.col.ord, lnames=lnames, ...)
#    }
    Cells.ex.dist <- distmirr(dist(Cells.ex, method = "euclidean"))
    results <- list(Cells.ex = Cells.ex, S = as.numeric(TopoS(EV.dist,Cells.ex.dist)))
    class(results) <- "EDAM.ex"
    invisible(results)
}
