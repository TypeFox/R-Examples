.arcp <- function(xy)
{
    if (nrow(xy) < 3)
        return(0);
    x.segmat <- cbind(xy, rbind(xy[2:nrow(xy), ],
                                xy[1, ]))
    abs(sum(x.segmat[,1] * x.segmat[,4] - x.segmat[,3]
            * x.segmat[,2])) / 2
}


clusthr <-function(xy, id=NULL)
{
    ## Verifications
    if (ncol(xy)!=2)
        stop("xy should have two columns")
    if (is.null(id))
        id<-rep(1, nrow(xy))
    id<-factor(id)
    if (length(id)!=nrow(xy))
        stop("id should have the same length as xy")

    ## splits the coordinates into one component per animal
    lixy<-split(xy, id)

    ## the output list
    res<-list()

    ## The function clubase is used to compute the length
    ## of the output vectors. It relies on the C function "longfacclust".
    ## It is needed to reserve memory for the main C function "clusterhrr"
    ## called below
    clubase <- function(xy)
    {
        nr <- as.integer(nrow(xy))
        xy <- as.double(t(as.matrix(xy)))
        len2 <- integer(1)
        toto <- .C("longfacclustr", as.double(xy),
                   as.integer(nr), integer(1),
                   PACKAGE = "adehabitat")[[3]]
        return(toto)
    }

    ## for each animal, an home range is desired
    for (i in names(lixy)){
        x<-lixy[[i]]

        ## computation of the output vectors
        len <- clubase(x)

        ## Computation of the tree: call to the C function "clusterhrr"
        toto <- .C("clusterhrr",
                   as.double(t(as.matrix(x))), as.integer(nrow(x)),
                   integer(len), integer(len), integer(len),
                   as.integer(len), PACKAGE = "adehabitat")

        facso <- toto[[3]]              # contains indices of the step of the
                                        # algorithm (First step =1, second = 2)

        nolocso <- toto[[4]]            # contains the indices of the
                                        # relocations clustered at the step
                                        # indicated by the factor step

        cluso <- toto[[5]]              # contains the indices of the cluster
                                        # in which the relocations clustered
                                        # at the step i are clustered

        ## output
        re <- data.frame(step = facso, clust = cluso, reloc = nolocso)
        res[[i]] <- list(xy = x, results = re)
    }

    ## Output
    class(res) <- "clusthr"
    return(res)
}



print.clusthr <- function(x, ...)
{
    if (!inherits(x, "clusthr"))
        stop("x should be of class \"ichr\"")
    cat("**********************************************\n*\n")
    cat("*  Home range estimation by cluster analysis\n")
    cat("*       (Kenward et al. 2001)\n*\n\n")
    cat("Home range estimated for the individuals:\n")
    print(names(x), quote = FALSE)
    cat("\nEach individual is a component of the list.\n")
    cat("For each individual, the following components are available:\n")
    cat("$xy: the relocations\n")
    cat("$results: a data.frame with the following columns:\n")
    cat("         $step: The step number of the algorithm\n")
    cat("         $clust: The cluster number assigned to some relocations\n")
    cat("         $reloc: The relocation(s) which is (are) assigned\n")
    cat("                 to the cluster \"clust\" at step \"step\"\n\n")

}


plot.clusthr <- function(x, whi = names(x), pch = 21,
                         bgpts = "white", colpts="black", cex=0.7,
                         plotit = TRUE, colpol = "grey",...)
{
    ## Verifications
    if (!inherits(x, "clusthr"))
        stop("x should be of class \"clusthr\"")
    x <- x[whi]
    class(x) <- "clusthr"

    ## Graphical settings
    if (plotit) {
        if (length(whi)>1) {
            opar <- par(mfrow = n2mfrow(length(whi)), mar=c(0,0,2,0))
            on.exit(par(opar))
        }
    }

    ## For each animal
    restep <- lapply(whi, function(i) {

        ## The main graph, with relocations
        if (plotit) {
            plot(x[[i]]$xy, asp=1, ty="n", main=names(x[i]),
                 axes=(length(whi)==1),...)
            box()
            points(x[[i]]$xy, pch= pch, bg = bgpts, col = colpts, cex=cex)
        }

        ## local variables
        step <- x[[i]]$results$step
        clust <- x[[i]]$results$clust
        reloc <- x[[i]]$results$reloc

        ##################
        ## Computes the home ranges for each step of the algorithm
        ## as an object of class "area"

        ## The relocations clustered at step one
        liclu <- list()
        liclu[[1]] <- reloc[step==1]

        ## The home range at step one
        poltot <- list()
        pc <- data.frame(id=factor(rep(1,3)),x[[i]]$xy[reloc[step==1],])
        pc <- as.area(pc)
        attr(pc,"nlocs") <- 3

        ## poltot contains the home range
        poltot[1] <- list(pc)



        ## For each step:
        for (j in 2:max(step)) {

            ## The relocations clustered at step j
            relocj <- reloc[step==j]
            r1 <- relocj[1]

            ## are the relocations already clustered (in which case the step is
            ## a merging of two clusters)
            oussa <- unlist(lapply(1:length(liclu), function(o)
                                   r1%in%liclu[[o]]))

            ## If it is a merging of two clusters
            if (any(oussa)) {
                ## we merge the two sets of relocations
                liclu[[which(oussa)]] <- liclu[[which(oussa)]][-c(which(liclu[[which(oussa)]]%in%relocj))]
            }

            ## update liclu according to the new cluster
            liclu[clust[step==j][1]] <- list(c(unlist(liclu[clust[step==j][1]]), relocj))

            ## computes the convex polygons around all current clusters
            kkk <- lapply(1:length(liclu), function(m) {
                k <- liclu[[m]]
                xy2 <- x[[i]]$xy[k,]
                pol <- xy2[chull(xy2[,1], xy2[,2]),]
                id <- rep(m, nrow(pol))
                pol <- data.frame(id, pol)
                attr(pol,"nlocs") <- nrow(xy2)
                return(pol)
            })

            ## Computes the number of relocations already clustered
            nlocc <- sum(unlist(lapply(kkk, function(w) attr(w,"nlocs"))))

            ## creates an object of class area with the polygons
            ## described above
            pol2 <- do.call("rbind", kkk)
            pol2[,1] <- factor(pol2[,1])
            class(pol2) <- c("area", "data.frame")
            attr(pol2, "nlocs") <- nlocc

            ## the home range is stored in poltot
            poltot[j] <- list(pol2)
        }

        ## In case of relocations located at the same place,
        ## the home range does not change.
        ## We delete the repetitions of the same home range
        ## and store the results in poltot2
        poltot2 <- list()
        poltot2[[1]] <- poltot[[1]]
        k <- 2
        for (p in 2:length(poltot)) {
            if (!identical(poltot[[p]], poltot[[p-1]])) {
                poltot2[[k]] <- poltot[[p]]
                k <- k+1
            }
        }

        ## And plots the result
        if (plotit) {

            ## The color of the polygons
            if (!is.na(colpol)) {
                foncol <- get(colpol, pos=".GlobalEnv")
                if (colpol=="grey") {
                    colp <- grey((1:length(poltot2))/(length(poltot2)+1))
                } else {
                    colp <- foncol(length(poltot2))
                }
            } else {
                colp <- NA
            }

            ## plots the polygons and the points
            lapply(length(poltot2):1, function(h) {
                ii <- poltot2[[h]]
                lapply(split(ii[,2:3], ii[,1]), function(u)
                       polygon(u, col=colp[h],border="black"))            })
            points(x[[i]]$xy, pch= pch, bg = bgpts, col = colpts, cex=cex)
        }

        ## returns the home ranges
        return(poltot2)
    })

    ## The function returns the home ranges
    invisible(restep)
}





clusthr.area <- function(x, percent = seq(20, 100, by = 5),
                         unin = c("m", "km"), unout = c("ha", "km2", "m2"),
                         plotit=TRUE)
{
    ## Verifications
    if (!inherits(x, "clusthr"))
        stop("x should be of class \"clusthr\"")
    unin <- match.arg(unin)
    unout <- match.arg(unout)

    ## graphical settings
    if (plotit) {
        opar <- par(mfrow=n2mfrow(length(x)))
        on.exit(par(opar))
    }

    ## Computes the home range
    u <- plot(x, plotit=FALSE)

    ## for each animal
    li <- lapply(1:length(u), function(d) {

        ## gets the home range
        o <- u[[d]]
        ## number of relocs
        nlo <- unlist(lapply(o, function(y) attr(y, "nlocs")))

        ## computes the area of the polygons
        ou <- unlist(lapply(o, function(y) {
            lib <- split(y[,2:3], y[,1])
            ji <- sum(unlist(lapply(lib, function(r) {
                class(r) <- c("data.frame")
                names(r) <- c("X","Y")
                r <- rbind(r,r[1,])
                pol <- Polygon(as.matrix(r))
                spdf <- SpatialPolygons(list(Polygons(list(pol), 1)))
                lar <- unlist(lapply(polygons(spdf)@polygons,
                                     function(x) unlist(lapply(x@Polygons, function(y)
                                                               .arcp(y@coords)))))
                lhol <- unlist(lapply(polygons(spdf)@polygons,
                                      function(x) unlist(lapply(x@Polygons, function(y)
                                                                y@hole))))
                sum(lar[!lhol])-sum(lar[lhol])
            })))
            return(ji)
        }))

        ## Depending on the ooutput units, change
        if (unin == "m") {
            if (unout == "ha")
                ou <- ou/10000
            if (unout == "km2")
                ou <- ou/1e+06
        }
        if (unin == "km") {
            if (unout == "ha")
                ou <- ou * 100
            if (unout == "m2")
                ou <- ou * 1e+06
        }

        ## percentage of relocations included in the home range
        nlo <- 100 * nlo/nrow(x[[d]]$xy)
        rere <- data.frame(nlo,ou)

        ## finds the home range containing the percentage
        ## of relocations the closest to the specified percent
        if (!is.null(percent)) {
            rere <- unlist(lapply(percent, function(e) {
                if (any(nlo<e)) {
                    res <- ou[max(which(nlo<=e))]
                } else {
                    warning(paste(e,
                                  "% contour could not be created.\n More data points are probably needed."))
                    res <- NA
                }
                return(res)
            }))
        }
        return(rere)
    })

    ## plots the results
    if (plotit) {
        lapply(1:length(x), function(i) {
            if (is.null(percent)) {
                mm <- li[[i]]
            } else {
                mm <- data.frame(percent, li[[i]])
            }
            plot(mm,
                 ty="l", xlab="Home range level",
                 ylab=paste("Home range size (",unout,")", sep=""),
                 main=names(x)[i])
            points(mm, pch=16, cex=0.5)
        })
    }

    ## in the case where percent is not specified,
    ## all the home range size are returned
    if (!is.null(percent)) {
        li <- as.data.frame(do.call("cbind", li))
        names(li) <- names(x)
        row.names(li) <- as.character(percent)
        class(li) <- c("hrsize", "data.frame")
        attr(li, "units") <- unout
    }
    return(li)
}






getverticesclusthr <- function(x, whi=names(x), lev=95)
{
    ## Verifications
    if (!inherits(x, "clusthr"))
        stop("x should be of class \"clusthr\"")
    x <- x[whi]
    class(x) <- "clusthr"

    ## Computes the home ranges
    uu <- plot(x, plotit=FALSE)

    ## gets the home range for a given level (percent)
    res2 <- lapply(1:length(uu), function(r) {
        y <- uu[[r]]
        nlo <- unlist(lapply(y, function(z) attr(z, "nlocs")))
        nlo <- 100 * nlo/nrow(x[[r]]$xy)
        if (any(nlo<lev)) {
            res <- max(which(nlo<=lev))
        } else {
            stop(paste(lev,"% contour could not be created.\n More data points are probably needed."))
        }
        return(y[[res]])
    })

    ## output of class "kver" (as with the function getverticeshr)
    names(res2) <- names(x)
    class(res2) <- "kver"
    return(res2)
}




kver.rast <- function(kv, asc)
{
    ## Verification
    if (!inherits(kv,"kver"))
        stop("kv should be of class \"kver\"")
    if (!inherits(asc,"asc"))
        stop("asc should be of class \"asc\"")

    ## For each animal, use of the function hr.rast to rasterize the
    ## objects of class "area", and sum the results over the columns
    ## of the resulting object of class "kasc"
    li <- lapply(kv, function(z) {

        ## rasterization: one map per polygon of the object of class
        ## area
        ka <- hr.rast(z,asc)
        class(ka) <- "data.frame"
        ka <- as.matrix(ka)

        ## if the pixel is contained in at least one polygon,
        ## 1 and 0 otherwise
        a <- apply(ka,1,function(a) {
            if (!all(is.na(a))) {
                return(1)
            } else {
                return(NA)
            }})
        ## output
        as2 <- matrix(a, ncol = ncol(asc))
        as2 <- getascattr(asc,as2)
        return(as2)
    })

    ## The resulting list is converted to an object of class "kasc"
    names(li) <- names(kv)
    return(as.kasc(li))
}


kver2shapefile <- function(kv, which=names(kv))
{
    ## Verifications
    if (!inherits(kv, "kver"))
        stop("x should be of class \"kver\"")
    whi <- which
    kv <- kv[whi]

    ## how many parts in the home range, for each home-range level
    ## in the object kver?
    nlev <- sum(unlist(lapply(kv, function(x) length(unique(x[,1])))))


    nlo <- 0
    Idatt <- numeric(0)
    Names <- character(0)

    for (i in 1:length(kv)) {

        ## gets the object area for the ith animal and converts it into a df
        class(kv[[i]]) <- "data.frame"

        ## The ID are re-defined (polygons belonging to the
        ## home range of different animals will be merged in the same object
        ## and a variable "Names" will be used to separate the home range
        ## in the shapefile
        kv[[i]][,1] <- as.numeric(factor(kv[[i]][,1]))+nlo

        ## list of polygons from the object area
        lkv <- split(kv[[i]], kv[[i]][,1])


        ## verifies that the polygons are closed (same first and last vertex)
        kv[[i]] <- do.call("rbind", lapply(lkv, function(x) {
            if (abs(sum(unlist(x[1,2:3]-x[nrow(x),2:3])))>1e-16) {
                return(rbind(x, x[1,]))
            } else {
                return(rbind(x))
            }}))

        ## The attributes of the shapefile (ID of the polygons and
        ## names of the animals to identify the home ranges
        Idatt <- c(Idatt, (nlo+1):(nlo+length(unique(kv[[i]][,1]))))
        Names <- c(Names, rep(names(kv)[i], length(unique(kv[[i]][,1]))))

        ## update nlo
        nlo <- nlo+length(unique(kv[[i]][,1]))
    }

    ## Output
    shp <- do.call("rbind", kv)
    names(shp) <- c("Id","X","Y")
    att <- data.frame(Id=Idatt, Names=Names)
    shp.file = convert.to.shapefile(shp, att, "Id", 5)
    return(shp.file)
}

