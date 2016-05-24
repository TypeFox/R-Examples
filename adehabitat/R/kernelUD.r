"kernelUD" <- function(xy, id=NULL, h="href", grid=40, same4all=FALSE,
                       hlim=c(0.1, 1.5), kern = c("bivnorm", "epa"),
                       extent=0.5)
{
    ## Verifications
    kern <- match.arg(kern)
    if (ncol(xy)!=2)
      stop("xy should have 2 columns")
    if ((!is.null(id))&(length(id)!=nrow(xy)))
      stop("id should have the same length as xy")
    if ((!is.numeric(h))&(h!="href")&(h!="LSCV"))
      stop("h should be numeric or equal to either \"href\" or \"LSCV\"")
    if ((h == "LSCV")&(kern == "epa"))
      stop("LSCV is not implemented with an Epanechnikov kernel")
    if (is.null(id))
      id<-rep(1, nrow(xy))
    id<-factor(id)
    if (min(table(id))<5)
    stop("At least 5 relocations are required to fit an home range")
    if (is.list(grid)) {
        if (is.null(names(grid)))
            stop("when grid is a list, it should have named elements")
        nn <- names(grid)
        lev <- levels(id)
        if (length(lev) != length(nn))
            stop("the length of the grid list should be equal to the number of levels of id")
        if (!all(lev%in%nn))
            stop("some ID levels do not have corresponding grids")
    }

    ## split xy into a list where each animal is an element
    lixy<-split(xy, id)
    sorties<-list()
    typh<-h
    htmp<-h
    gr<-grid

    ## If the same grid is wanted for all animals:
    ## First computes this grid
    if (same4all) {
        if (!is.list(gr)) {
            ## if the grid is not given
            if (length(as.vector(gr))==1) {

                ## the "core" grid
                if (!is.numeric(gr))
                    stop("grid should be an object of class asc or a number")
                xli<-range(xy[,1])
                yli<-range(xy[,2])
                xli<-c(xli[1]-extent*abs(xli[2]-xli[1]),
                       xli[2]+extent*abs(xli[2]-xli[1]))
                yli<-c(yli[1]-extent*abs(yli[2]-yli[1]),
                       yli[2]+extent*abs(yli[2]-yli[1]))
                xygg<-data.frame(x=xli, y=yli)
                grid<-ascgen(xygg, nrcol=grid)
                cellsize<-attr(grid, "cellsize")
                lx<-nrow(grid)*cellsize
                ly<-ncol(grid)*cellsize
                ref<-lx
                if (ly>lx)
                    ref<-ly
                xll<-attr(grid, "xll")
                yll<-attr(grid, "yll")

                ## One adds empty rows and columns to the "core" grid
                xll<-xll-lx/2
                yll<-yll-ly/2
                arajlig<-ceiling((lx/2)/cellsize)
                arajcol<-ceiling((ly/2)/cellsize)
                mrajlig<-matrix(0, ncol=ncol(grid), nrow=arajlig)
                grid<-rbind(mrajlig, grid, mrajlig)
                mrajcol<-matrix(0, ncol=arajcol, nrow=nrow(grid))
                grid<-cbind(mrajcol, grid, mrajcol)

                ## We add the attributes
                attr(grid, "xll")<-xll
                attr(grid, "yll")<-yll
                attr(grid, "cellsize")<-cellsize
                attr(grid, "type")<-"numeric"
                class(grid)<-"asc"
            }
        } else {
            stop("when same4all=TRUE, a list of grid cannot be passed as \"grid\"")
        }
    }


    ## UD estimation for each animal
    for (i in 1:nlevels(id)) {

        df<-lixy[[i]]

        ## 1. Computation of h
        varx<-var(df[,1])
        vary<-var(df[,2])
        sdxy<-sqrt(0.5*(varx+vary))
        n<-nrow(df)
        ex<-(-1/6)
        href<-sdxy*(n^ex)
        if (kern=="epa")
            href <- href*1.77
        if (h=="href") {
            htmp<-href
        }
        if (h=="LSCV") {
            hvec<-seq(hlim[1]*href, hlim[2]*href, length=100)
            CV<-.C("CVmise", as.integer(nrow(df)), as.double(df[,1]),
                   as.double(df[,2]),
                   as.double(hvec), double(length(hvec)),
                   as.integer(length(hvec)), PACKAGE="adehabitat")[[5]]
            htmp<-hvec[CV==min(CV)]
            if ((CV[CV==min(CV)]==CV[1])|(CV[CV==min(CV)]==CV[length(CV)]))
                warning("The algorithm did not converge \nwithin the specified range of hlim: try to increase it")
        }


        ## 2. The grid if not the same for all
        if (!is.list(gr)) {
            if (length(as.vector(gr))==1) {
                if (!is.numeric(gr))
                    stop("grid should be an object of class asc or a number")

                if (!same4all) {

                    ## the "core" grid
                    grid<-matrix(0, ncol=gr, nrow=gr)
                    rgx<-range(df[,1])
                    rgy<-range(df[,2])
                    lx<-rgx[2]-rgx[1]
                    ly<-rgy[2]-rgy[1]
                    ref<-lx
                    if (ly>lx)
                        ref<-ly

                    xll<-rgx[1]
                    yll<-rgy[1]
                    cellsize<-ref/ncol(grid)

                    ## One adds empty rows and columns to the "core" grid
                    xll<-xll-lx*extent
                    yll<-yll-ly*extent
                    arajlig<-ceiling((lx*extent)/cellsize)
                    arajcol<-ceiling((ly*extent)/cellsize)
                    mrajlig<-matrix(0, ncol=ncol(grid), nrow=arajlig)
                    grid<-rbind(mrajlig, grid, mrajlig)
                    mrajcol<-matrix(0, ncol=arajcol, nrow=nrow(grid))
                    grid<-cbind(mrajcol, grid, mrajcol)

                    ## We add the attributes of the grid
                    attr(grid, "xll")<-xll
                    attr(grid, "yll")<-yll
                    attr(grid, "cellsize")<-cellsize
                    attr(grid, "type")<-"numeric"
                    class(grid)<-"asc"
                }
            }
        } else {
            grid <- gr[[names(lixy)[i]]]
        }

        grille<-grid
        xylo<-getXYcoords(grid)
        xg<-xylo$x
        yg<-xylo$y


        ## Kernel estimation in itself (the C function called
        ## depends on the choosed kernel)
        if (kern=="bivnorm") {
            toto<-.C("kernelhr", double(nrow(grid)*ncol(grid)),as.double(xg),
                     as.double(yg),
                     as.integer(ncol(grid)), as.integer(nrow(grid)),
                     as.integer(nrow(df)), as.double(htmp),
                     as.double(df[,1]), as.double(df[,2]),
                     PACKAGE="adehabitat")
        }
        if (kern=="epa") {
            toto<-.C("kernepan", double(nrow(grid)*ncol(grid)),as.double(xg),
                     as.double(yg),
                     as.integer(ncol(grid)), as.integer(nrow(grid)),
                     as.integer(nrow(df)), as.double(htmp),
                     as.double(df[,1]), as.double(df[,2]),
                     PACKAGE="adehabitat")
        }

        ## output
        UD<-matrix(toto[[1]], nrow=nrow(grid), byrow=TRUE)
        UD<-getascattr(grid, UD)
        if (typh=="LSCV") {
            CV<-data.frame(h=hvec, CV=CV)
            convergence<-min(CV[,2])!=CV[1,2]
            htmp<-list(CV=CV, convergence=convergence, h=htmp)
        }
        attr(UD, "UD") <- "simple"
        sorties[[names(lixy)[i]]]<-list(UD=UD, h=htmp, locs=df, hmeth=typh)
    }
    ## general output
    class(sorties)<-c("khrud", "khr")
    return(sorties)
  }

