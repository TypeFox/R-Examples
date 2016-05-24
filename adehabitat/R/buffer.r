## buffer.point.unic is used to compute a "mask" for buffers:
## It returns a matrix containing a rasterized circle of radius
## equal to the buffer size, with the same resolution as the map

".buffer.point.unic" <- function(x, md)
{
    ## gets the cell size
    res<-attr(x, "cellsize")
    ## the number of pixels of the "mask" matrix
    nmax<-ceiling(md/res)

    ## creates the mini-matrix of NA
    calc<-matrix(NA,ncol = 2*nmax+1, nrow=2*nmax+1)
    cour<-nmax+1
    ## cour contains the index of the row (and of the column)
    ## of the pixel located at the center of the map

    ## each pixel of the map takes the value 1 if its center
    ## is located at a distance < md of the center of the map
    for (i in 1:nrow(calc)) {
        for (j in 1:ncol(calc)) {
            d1<-(((cour-i-0.5)*res)^2+((cour-j-0.5)*res)^2)
            d2<-(((cour-i-0.5)*res)^2+((cour-j+0.5)*res)^2)
            d3<-(((cour-i+0.5)*res)^2+((cour-j-0.5)*res)^2)
            d4<-(((cour-i+0.5)*res)^2+((cour-j+0.5)*res)^2)

            if (min(c(d1,d2,d3,d4))<= md^2)
                calc[i,j]<-1
        }
    }

    ## output
    return(calc)
}


"buffer" <- function(pts, x, dist)
{
    ## Verifications
    if (inherits(x, "asc"))
        x<-as.kasc(list(toto=x))
    if (inherits(x, "kasc"))
        x<-storemapattr(x)
    if (!inherits(x, "mapattr"))
        stop("non convenient format for x")
    res<-attr(x, "cellsize")
    nmax<-ceiling(dist/res)

    ## removes the missing values
    pts <- pts[!is.na(pts[,1]),]
    pts <- pts[!is.na(pts[,2]),]

    ## Computation of the "mask" with .buffer.point.unic (see above)
    ## and replace the NA by 0
    calc<-.buffer.point.unic(x, dist)
    calc0<-calc
    calc0[is.na(calc0)]<-0

    ## Count the points in each pixel of the map: allows to
    ## identify the pixels containing at least one point
    asc<-count.points(pts, x)
    vasc<-as.vector(asc)

    ## Computation of an idlig and idcol of the same length as vasc.
    ## idasc will allow to identify the values of idlig et idcol
    ## for which the cells are >0
    idasc<-1:length(vasc)
    idcons<-idasc[vasc>0]
    idlig<-as.vector(row(asc))
    idcol<-as.vector(col(asc))

    ## keeps only the indices of the rows and columns for the
    ## pixels containing at least one point
    ligcons<-idlig[idcons]
    colcons<-idcol[idcons]

    ## output matrix
    ## temporarily adds rows and columns to the map, so that it is
    ## sufficiently large to place the "mask" matrix (so that we
    ## add +2*nmax)
    sorties<-matrix(0, nrow=(attr(x, "ncol")+2*nmax),
                    ncol=(attr(x, "nrow")+2*nmax))

    ## for each non-empty pixel, adds the mask on the map
    for (i in 1:length(idcons)) {

        ## places the "mask" matrix on the area, depending on
        ## the non-empty pixel location
        car<-matrix(0, nrow=(attr(x, "ncol")+2*nmax),
                    ncol=(attr(x, "nrow")+2*nmax))
        car[c(ligcons[i]:(ligcons[i]+2*nmax)),
            c(colcons[i]:(colcons[i]+2*nmax))]<-
                car[c(ligcons[i]:(ligcons[i]+2*nmax)),
                    c(colcons[i]:(colcons[i]+2*nmax))]+
                        calc0
        ## ...and adds this temporary matrix to the output map
        sorties<-sorties+car
    }

    ## remove the additional rows and columns added previously
    sorties<-sorties[c((nmax+1):(nrow(sorties)-nmax)),
                     c((nmax+1):(ncol(sorties)-nmax))]

    ## The matrix is either 0 (outside the buffer) or 1 (inside)
    sorties<-matrix(as.numeric(sorties!=0), ncol=attr(x, "nrow"))

    ## Output
    sorties[sorties==0]<-NA
    attr(sorties, "cellsize")<-attr(x, "cellsize")
    attr(sorties, "xll")<-attr(x, "xll")
    attr(sorties, "yll")<-attr(x, "yll")
    attr(sorties, "type")<-"numeric"
    class(sorties)<-"asc"
    return(sorties)
  }

