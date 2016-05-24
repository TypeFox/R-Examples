# screeningStep.r


ScreeningStep <- function(y.tilde, gram, cg.all, nm, v, r, q0=0.1, scale=1) {

  #############################################################################
  # This function conduct the GS-step of the graphet screening
  #
  # Args:
  #   y.tilde=X'%*%y
  #   gram: the regularized gram matrix
  #   cg.all: a list whose kth element is a matrix of k columns. Its rows contains all the connected
  #       subgraph with k nodes.
  #   nm: the maximal subgraph invesgated in the screening step
  #   v, r: the essential tuning parameter of graphlet screening, 
  #       tied to the sparse level and the signal strength, respectively.
  #   q0: the minimal screening parameter.
  #   scale: q(D,F)=q^{max}(D,F)*scale, default is scale=1
  #
  # Returns:                                                       
  #   survivor: a logical vector, where TRUE means retained as a protential signal.
  #   
  ###############################################################################

  p <- length(y.tilde)
  # univeriate screening
  q1 <- (v+r)^2/r^2/4
  tau <- sqrt(2 * r * log(p))
  tau1 <- sqrt(q1) * tau * sqrt(scale)
  survivor <- (abs(y.tilde) > tau1)
  indices<-which(survivor)
  if (nm>1){
    # multivariate screening
    for (ii in 2:nm) {
      if (ii > 2/v) {
        break
      }
      #find all the CG with ii nodes in the graph
      cg.ii <- cg.all[[ii]]
      #if there is no such cg.all, break
      if (length(cg.ii) == 0){
        break
      }
      nii <- dim(cg.ii)[1]
      indicatorii <- rep(FALSE, p)
      for (jj in 1:nii){
        #find the indices to be tested.
        ds <- setdiff(cg.ii[jj, ], indices)
        d <- length(ds)  
        if ((d > 0) & (d <= 1/v)) {
          #define the test statistic.
          es <- intersect(cg.ii[jj, ], indices)
          if (length(es) > 0) {
            ncmatrx <- gram[ds, ds] - gram[ds, es]%*%solve(gram[es, es]) %*% gram[es, ds]
          } else {
            ncmatrx <- gram[ds,ds]
          }
          if (d == 1) {
            ww <- ncmatrx
          } else if (d == 2) {
            ww <- d * min(eigen(ncmatrx, only.values =TRUE)$values)
          } else {
              www <- rep(0,2^(d-1));             
              # for a length d vector, let the sign of the first coordinate to be 1, there are 2^(d-1) combinations.
              # by symmetry of Lfun, the the above assumption is ok.
              for (kk in 0:(2^(d-1)-1)) {               
                idx <-  VectorizeBase(kk, 2, d-1) # express k on base 2 as a row vector
                signs  <- c(1,2*idx - 1) # convert to -1, 1 to indicate sign                   
                dvec <- rep(0,d)
                amat <- diag(signs)
                bvec <- rep(1,d)
                result.qp <- solve.QP(ncmatrx, dvec, amat, bvec)  
                www[kk+1] <- result.qp$value
              }
            ww<-min(www)          
          }           
          ww.condition <- (ww > (v/r)*(2*d-(d-sqrt(d^2-1))*((d/2) != floor(d/2))))
          if (ww.condition == 1) {
            coor <- ww*r
            deltadelta <- v^2/coor/4-sqrt((coor^2-2*coor*d*v+v^2)*(coor^2-2*coor*d*v+v^2>0))+sqrt((coor^2-2*coor*d*v)*(coor^2-2*coor*d*v>0))
            teb <- scale*2*log(p)*(coor*5/4-v*d/2-sqrt((coor^2-2*coor*d*v)*(coor^2-2*coor*d*v > 0))+deltadelta*(d/2 != floor(d/2)))
          } else {
            teb <- q0
          }
          if (length(es) == 0) {
            test.stat <- t(y.tilde[cg.ii[jj, ]])%*%ginv(gram[cg.ii[jj, ], cg.ii[jj, ]])%*%y.tilde[cg.ii[jj, ]]
          } else {
            test.stat <- t(y.tilde[cg.ii[jj, ]])%*%ginv(gram[cg.ii[jj, ], cg.ii[jj, ]])%*%y.tilde[cg.ii[jj, ]]-t(y.tilde[es])%*%ginv(gram[es, es])%*%y.tilde[es]
          }
          if (test.stat > teb) {
            # as long as a node pass in one test, it will be included in the survival indices.
            #survivor[ds] <- TRUE
            indicatorii[ds] <- TRUE
          }
        }      
      }
    survivor <- as.logical(survivor + indicatorii)
    indices <- which(survivor) 
    }
  }
  return(survivor)
}



