## distance of point to linear interpolation of DVH curve
dvhDistance <-
function(x, DV) {
    ## relevant dose + volume data as matrix
    oneDistance <- function(x, oneDV) {
        if(any(is.na(oneDV))) {
            return(list(dstMin=NA_real_, dstMinRel=NA_real_, dstCnstr=NA_real_,
                        ptMinD=NA_real_, ptMinV=NA_real_))
        }

        DVmat <- if(oneDV$volRel &&  oneDV$doseRel) {
            data.matrix(x$dvh[ , c("doseRel", "volumeRel")])
        } else if(  oneDV$volRel && !oneDV$doseRel) {
            data.matrix(x$dvh[ , c("dose",    "volumeRel")])
        } else if( !oneDV$volRel &&  oneDV$doseRel) {
            data.matrix(x$dvh[ , c("doseRel", "volume")])
        } else if( !oneDV$volRel && !oneDV$doseRel) {
            data.matrix(x$dvh[ , c("dose",    "volume")])
        }

        ## there may be duplicate points due to precision in DVH file
        DVmat <- unique(DVmat)

        ## constraint point
        DVcpt <- c(oneDV$D, oneDV$V)

        #sp::spDistsN1(DVmat, DVcpt, longlat=FALSE)
        #DVlines   <- sp::Line(DVmat)
        #DVlinesL  <- sp::Lines(list(DVlines), ID="A")
        #DVlinesSL <- sp::SpatialLines(list(DVlinesL))
        #DVpt      <- sp::SpatialPoints(matrix(DVcpt, nrow=1))
        #rgeos::gDistance(DVlinesSL, DVpt)
        ## orthogonal projections of DV onto subspaces a'x
        ## defined by segments of x (DVcpt is treated as column vector)
        ## subspaces are affine -> shift constraint to origin first - for each subspace
        cptMat <- -sweep(DVmat, 2, DVcpt, "-")
        DVbase <- diff(DVmat)              # basis vectors of segments
        lens   <- sqrt(rowSums(DVbase^2))  # length of basis vectors
        
        ## basis vectors scaled to unit length
        DVbu <- diag(1/lens, ncol=length(lens)) %*% DVbase
        
        ## orthogonal projections of all shifted constraint points
        ## onto all subspaces in in subspace coords
        opAll <- DVbu %*% t(cptMat[-nrow(cptMat), , drop=FALSE])
        op    <- diag(opAll)  # just the relevant projections (corresponding subscpace)

        ## check if projection is outside of DVH segment -> not in [0, segmentLength]
        outside <- (op < 0) | (op > lens)

        ## projections in DV coords aa'x (corresponding pairs of row vectors)
        ## general case (>= 1D)
        # opOCL <- Map(tcrossprod, split(DVbu, f=seq_len(nrow(DVbu))),
        #                          split(op,   f=seq_len(length(op))))
        # opOC0 <- t(do.call(cbind, opOCL))
        
        ## special case 1D
        opOC0 <- DVbu * op

        ## add original coords back to move away from origin
        opOC <- opOC0 + DVmat[-nrow(DVmat), ]
        
        ## keep only those projections that are on the DVH segment
        opOCin <- opOC[!outside, , drop=FALSE]

        ## distance of each point in DV to
        ## all points in DVmat as well as all non-outside orthogonal projections
        dstSupport <- sqrt(rowSums(sweep(DVmat,  2, DVcpt, "-")^2))
        dstProj    <- sqrt(rowSums(sweep(opOCin, 2, DVcpt, "-")^2))
        
        ## catch pathological cases
        if(!is.finite(min(dstSupport)) && !is.finite(min(dstProj))) {
            dstMin <- NA_real_
            ptMin  <- c(NA_real_, NA_real_)
        } else if(min(dstSupport) < suppressWarnings(min(dstProj))) {
            dstMin <- min(c(dstSupport, dstProj))
            ptMin  <- DVmat[which.min(dstSupport), ]
        } else {
            dstMin <- min(c(dstSupport, dstProj))
            ptMin  <- opOCin[which.min(dstProj), ]
        }

        ## get x for y=0 and y for x=0 (two-point form)
        yAxisX <- 0
        yAxisY <- -oneDV$D * ((ptMin[2]-oneDV$V)/(ptMin[1]-oneDV$D)) + oneDV$V
        xAxisX <- -oneDV$V * ((ptMin[1]-oneDV$D)/(ptMin[2]-oneDV$V)) + oneDV$D
        xAxisY <- 0

        ## get distance from constraint to axes along line with ptMin
        dstCnstrXaxis <- sqrt((oneDV$D-xAxisX)^2 + (oneDV$V-xAxisY)^2)
        dstCnstrYaxis <- sqrt((oneDV$D-yAxisX)^2 + (oneDV$V-yAxisY)^2)
        
        ## choose axis with smaller distance
        dstCnstr  <- min(c(dstCnstrXaxis, dstCnstrYaxis))
        dstMinRel <- dstMin / dstCnstr

        return(list(dstMin=dstMin, dstMinRel=dstMinRel, dstCnstr=dstCnstr,
                    ptMinD=ptMin[1], ptMinV=ptMin[2]))
    }

    ## get distance for each DV constraint point
    Map(oneDistance, list(x), split(DV, f=seq_len(nrow(DV))))
}
