##' This stitches together the incisions and tears by inserting new
##' points in the tears and creating correspondences between new
##' points.
##'
##' @title Stitch together incisions and tears in an AnnotatedOutline
##' @param a \code{\link{AnnotatedOutline}} object 
##' @return
##' \item{\code{Rset}}{the set of points on the rim}
##' \item{\code{i0}}{the index of the landmark}
##' \item{\code{P}}{a new set of meshpoints}
##' \item{\code{V0}}{indicies of the apex of each tear}
##' \item{\code{VF}}{indicies of the forward vertex of each tear}
##' \item{\code{VB}}{indicies of the backward vertex of each tear}
##' \item{\code{TFset}}{list containing indicies of points in each foward tear}
##' \item{\code{TBset}}{list containing indicies of points in each backward tear}
##' \item{\code{gf}}{new forward pointer list}
##' \item{\code{gb}}{new backward pointer list}
##' \item{\code{h}}{correspondence mapping}
##' \item{\code{hf}}{correspondence mapping in forward direction for
##' points on boundary}
##' \item{\code{hb}}{correspondence mapping in backward direction for
##' points on boundary}
##' @author David Sterratt
##' @export
StitchedOutline <- function(a) {
  
  r <- computeTearRelationships(a, a$V0, a$VB, a$VF)

  ## If not set, set the landmark marker index. Otherwise
  ## check it
  Rset <- r$Rset
  if (!(a$i0 %in% Rset)) {
    print(a$i0)
    print(Rset)
    stop("Fixed Point is not in rim")
  }

  P <- a$P
  V0 <- a$V0
  VF <- a$VF
  VB <- a$VB
  TFset <- r$TFset
  TBset <- r$TBset
  gf <- a$gf
  gb <- a$gb
  hf <- r$hf
  hb <- r$hb
  h <- r$h
  
  ## Insert points on the backward tears corresponding to points on
  ## the forward tears
  sF <-      stitch.insert.points(P, V0, VF, VB, TFset, TBset,
                                               gf, gb, hf, hb, h,
                                               "Forwards")

  ## Insert points on the forward tears corresponding to points on
  ## the backward tears
  sB <- with(sF,
             stitch.insert.points(P, V0, VB, VF, TBset, TFset,
                                  gb, gf, hb, hf, h,
                                  "Backwards"))
  ## Extract data from object
  P <- sB$P
  gf <- sB$gb
  gb <- sB$gf
  hf <- sB$hb
  hb <- sB$hf
  h <- sB$h

  ## Link up points on rim
  h[Rset] <- hf[Rset]
  
  ## Make sure that there are no chains of correspondences
  while (!all(h==h[h])) {
   h <- h[h]
  }

  s <- merge(list(Rset=Rset, i0=a$i0,
                  VF=VF, VB=VB, V0=V0,
                  TFset=TFset, TBset=TBset,
                  P=P, h=h, hf=hf, hb=hb,
                  gf=gf, gb=gb), a)
  class(s) <- addClass("stitchedOutline", a)
  return(s)
}

## Inner function responsible for inserting the points
stitch.insert.points <- function(P, V0, VF, VB, TFset, TBset, gf, gb, hf, hb, h,
                                 dir) {
  M <- length(V0)                       # Number of tears
  ## Iterate through tears to insert new points
  for (j in 1:M) {
    ## Compute the total path length along each side of the tear
    Sf <- path.length(V0[j], VF[j], gf, hf, P)
    Sb <- path.length(V0[j], VB[j], gb, hb, P)
    message(paste("Tear", j, ": Sf =", Sf, "; Sb =", Sb))

    ## For each point in the forward path, create one in the backwards
    ## path at the same fractional location
    message(paste("  ", dir, " path", sep=""))
    for (i in setdiff(TFset[[j]], c(V0[j], VF[j]))) {
      sf <- path.length(V0[j], i, gf, hf, P)
      ## If the point isn't at the apex, insert a point
      if (sf > 0) {
        message(paste("    i =", i,
                                "; sf/Sf =", sf/Sf,
                                "; sf =", sf))
        for (k in TBset[[j]]) {
          sb <- path.length(V0[j], k, gb, hb, P)
          message(paste("      k =", format(k, width=4),
                                  "; sb/Sb =", sb/Sb,
                                  "; sb =", sb))
          if (sb/Sb > sf/Sf) {
            break;
          }
          k0 <- k
          sb0 <- sb
        }

        ## If this point does not point to another, create a new point
        if ((hf[i] == i)) {
          f <- (sf/Sf*Sb-sb0)/(sb-sb0)
          message(paste("      Creating new point: f =", f))
          p <- (1-f) * P[k0,] + f * P[k,]

          ## Find the index of any row of P that matches p
          n <- anyDuplicated(rbind(P, p), fromLast=TRUE) 
          if (n == 0) {
            ## If the point p doesn't exist
            P <- rbind(P, p)
            ## Update forward and backward pointers
            n <- nrow(P)                    # Index of new point
            gb[n]     <- k
            gf[n]     <- gf[k]
            gb[gf[k]] <- n
            gf[k]     <- n

            ## Update correspondences
            hf[n] <- n
            hb[n] <- n
            h[i] <- n
            h[n] <- n
          } else {
            message(paste("      Point", n, "already exists"))
            h[i] <- n
            h[n] <- n
          }
        } else {
          ## If not creating a point, set the point to point to the forward pointer 
          h[i] <- hf[i]
        }
      }
    } 
  }
  return(list(P=P, hf=hf, hb=hb, gf=gf, gb=gb, h=h))
}

##' Plot flat \code{\link{StitchedOutline}}. If the optional argument
##' \code{stitch} is \code{TRUE} the user markup is displayed.
##'
##' @title Flat plot of AnnotatedOutline
##' @param x \code{\link{AnnotatedOutline}} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param stitch If \code{TRUE}, plot stitch
##' @param lwd Line width
##' @param ... Other parameters
##' @method flatplot stitchedOutline
##' @author David Sterratt
##' @export
flatplot.stitchedOutline <- function(x, axt="n", ylim=NULL,
                                     stitch=TRUE, lwd=1,
                                     ...) {
  NextMethod()

  if (stitch) {
    with(x, {
      for (TF in TFset) {
        lines(P[TF,], col=getOption("TF.col"), lwd=lwd)
      }
      for (TB in TBset) {
        lines(P[TB,], col=getOption("TB.col"), lwd=lwd)
      }
      for (j in 1:length(h)) {
        if (h[j] != j) {
          lines(P[c(j, h[j]),], col=getOption("stitch.col"), lwd=lwd)
        }
      }
      
      for (j in 1:length(hf)) {
        if (hf[j] != j) {
          lines(P[c(j, hf[j]),], col=getOption("stitch.col"), lwd=lwd)
        }
      }
      for (j in 1:length(hb)) {
        if (hb[j] != j) {
          lines(P[c(j, hb[j]),], col=getOption("V.stitch.col"), lwd=lwd)
        }
      }
    })
  }
}
