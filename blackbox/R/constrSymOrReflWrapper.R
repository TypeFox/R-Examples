constrSymOrReflWrapper <- function(points, edge, logLcheck=-Inf, refpoint=c(), lower, upper) {
  ## a function to compute the symmetric of a point refpoints or its orthogonal reflexion through an edge
  ## points must be of class matrix
  ## boundcheck gives a minimal predicted log L value
  ## the edge is described by d points; it is used only is refpoint is empty
  if (length(refpoint)>0) { ## symmetric through refpoint
    resu <- t(apply(points, 1, constrSymOrRefl, logLcheck=logLcheck, refpoint=refpoint, lower=lower, upper=upper, edge=edge))
  } else { ##symmetric through edge. An old try with obvious deficiencies as an exploration of parameter space
    ## example: symmetricPts(matrix(c(0, 0, 0), ncol=3), diag(3)) is reflexion of origin through unit simplex
    ## weakness: say edge is ((Nm=40, g=0.55), (Nm=50, g=0.45)); the mirror image of (Nm=10, g=0.5) is
    ##symmetricPts(matrix(c(10, 0.5), ncol=2), matrix(c(50, 0.45, 40, 0.55), ncol=2)) => c(2.54..., 9.86...) huge g
    ## FR->FR should use sweep on next line !!!
    vv <- t(apply(edge[-1, , drop=FALSE], 1, function(v) {v-edge[1, ]})) ## edge eq y=edge[1, ]+sum_j a_j vv_j
    orthog <- qr.Q(qr(t(vv))) ## column vectors = orthonormal basis for space spanned by vv: t(orthog) %*% orthog=I
    resu <- t(apply(points, 1, constrSymOrRefl, logLcheck=logLcheck, orthog=orthog, lower=lower, upper=upper, edge=edge))
  }
  return(resu)
} ## end def constrSymOrReflWrapper
