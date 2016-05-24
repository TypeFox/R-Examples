#########################################################################
# Categorical Network Class Methods
# Joint Probability Calculations

probMatAddNode <- function(object, node, pmat, imat) {
  pars <- object@pars[[node]]
  ncats <- length(object@cats[[node]])
  id <- which(colnames(pmat)==node)
  if(!is.null(pmat) && nrow(pmat) > exp(log(2)*14)) {
    warning("The table exceeds ", nrow(pmat), " rows.")
    return(NULL)
  }
  if(length(id) > 0)
    return(list(pmat,imat))    
  if(is.null(pars)) {
    pm <- pmat
    im <- imat
    pmat <- NULL
    imat <- NULL
    for(c in 1:ncats) {
      if(is.null(pm)) {
        pmc <- cbind(pm, object@probs[[node]][c])
        imc <- cbind(im, c)
      }
      else {
        pmc <- cbind(pm, rep(object@probs[[node]][c], nrow(pm)))
        imc <- cbind(im, rep(c, nrow(pm)))
      }
      pmat <- rbind(pmat, pmc)
      imat <- rbind(imat, imc)
    }
    colnames(pmat) <- c(colnames(pm), node)
    colnames(imat) <- c(colnames(im), node)
    return(list(pmat,imat))
  }
  ipar <- NULL
  res <- TRUE
  for(par in pars) {
    id <- which(colnames(pmat)==par)
    if(length(id) < 1) {
      res <- probMatAddNode(object, par, pmat, imat)
      if(is.null(res))
        break
      pmat <- res[[1]]
      imat <- res[[2]]
      id <- ncol(pmat)
    }
    ipar <- c(ipar, id[1])
  }
  if(is.null(res))
    return(NULL)
  pm <- pmat
  im <- imat
  pmat <- NULL
  imat <- NULL
  for(j in 1:nrow(pm)) {
    prow <- pm[j,]
    irow <- im[j,]
    pl <- object@probs[[node]]
    for(k in 1:length(pars)) {
      pl <- pl[[im[j,ipar[k]]]]
    }
    for(nc in 1:ncats) {
      pmat <- rbind(pmat, c(prow, pl[nc]))
      imat <- rbind(imat, c(irow, nc))
    }
  }
  colnames(pmat) <- c(colnames(pm), node)
  colnames(imat) <- c(colnames(im), node)
  return(list(pmat,imat))
}

setMethod("cnJointProb", "catNetwork",
          function(object, nodes) {
            if(!is(object, "catNetwork"))
              stop("catNetwork object is required.")

            if(is.character(nodes))
              nodes <- sapply(nodes, function(cc) {
                id <- which(object@nodes == cc)
                if(length(id)>0)
                  return(id[1])
                return(-1)
              })
            nodes <- as.integer(nodes)
            for(node in nodes)
              if(node < 0 || nodes > object@numnodes)
                stop("Invalid nodes")
              
            pmat <- NULL
            imat <- NULL
            for(node in nodes) {
              res <- probMatAddNode(object, node, pmat, imat)
              if(is.null(res))
                break
              pmat <- res[[1]]
              imat <- res[[2]]
            }
            if(is.null(res)) 
              return(NULL)
            pjoint <- cbind(imat, p=apply(pmat, 1, "prod"))
            return(pjoint)
          })

## calculates P(x|y)
## x and y should be named
setMethod("cnCondProb", "catNetwork",
          function(object, x, y) {
            if(!is(object, "catNetwork"))
              stop("catNetwork object is required.")

            xnodes <- names(x)
            ynodes <- names(y)
            if(is.character(xnodes))
              xnodes <- sapply(xnodes, function(cc) {
                id <- which(object@nodes == cc)
                if(length(id)>0)
                  return(id[1])
                return(-1)
              })
            xnodes <- as.integer(xnodes)
            if(length(xnodes) < 1)
              stop("x cannot be empty")
            if(is.character(ynodes))
              ynodes <- sapply(ynodes, function(cc) {
                id <- which(object@nodes == cc)
                if(length(id)>0)
                  return(id[1])
                return(-1)
              })
            ynodes <- as.integer(ynodes)
            for(i in 1:length(ynodes)) {
              id <- which(xnodes == ynodes[i])
              if(length(id) >= 1) {
                if(x[id[1]] != y[i])
                  stop("Wrong expression")
                ## P(X, z=a|Y, z=a) = P(X|Y), remove z
                x <- x[-id[1]]
                y <- y[-i]
                xnodes <- xnodes[-id[1]]
                ynodes <- ynodes[-i]
              }
            }
            nodes <- c(xnodes, ynodes)
            vals <- c(x, y)
            for(i in 1:length(nodes)) {
              node <- nodes[i]
              if(node < 0 || nodes > object@numnodes)
                stop("Invalid nodes")
              if(is.character(vals[i])) {
                id <- which(object@cats[[node]] == vals[i])
                if(length(id) < 0)
                  stop("Invalid value for node ", node)
                vals[i] <- id[1]
              }
              vals[i] <- as.integer(vals[i])
              if(vals[i] < 0 || vals[i] > length(object@cats[[node]]))
                stop("Invalid value for node ", node)
            }
            xvals <- vals[1:length(xnodes)]
            yvals <- NULL
            if(length(xnodes) < length(nodes))
              yvals <- vals[(length(xnodes)+1):length(nodes)]
            
            pmat <- NULL
            imat <- NULL
            for(node in nodes) {
              res <- probMatAddNode(object, node, pmat, imat)
              if(is.null(res))
                break
              pmat <- res[[1]]
              imat <- res[[2]]
            }

            if(!is.null(res)) {
              pjoint <- apply(pmat, 1, "prod")
              ixnodes <- sapply(xnodes, function(nn) which(colnames(imat)==nn)[1])
              iynodes <- sapply(ynodes, function(nn) which(colnames(imat)==nn)[1])
              
              py <- 1
              if(length(iynodes) > 0) {
                jy <- NULL
                for(j in 1:nrow(imat))
                  if(prod(imat[j,iynodes] == yvals))
                    jy <- c(jy, j)
                py <- sum(pjoint[jy])
                jx <- NULL
                for(j in 1:nrow(imat))
                  if(prod(imat[j,ixnodes] == xvals)*prod(imat[j,iynodes] == yvals))
                    jx <- c(jx, j)
              }
              else {
                jx <- NULL
                for(j in 1:nrow(imat))
                  if(prod(imat[j,ixnodes] == xvals))
                    jx <- c(jx, j)
              }
              px <- sum(pjoint[jx])
              return(px/py)
            }

            ## approximate from a sample
            warning("The required probability will be approximated")
            ps <- cnSamples(object, floor(1024*exp(log(object@maxcats)*(1+object@maxpars))), as.index = TRUE)
            px <- 0
            py <- 0
            for(j in 1:nrow(ps)) {
              if(prod(ps[j, ynodes] == yvals)) {
                py <- py+1
                if(prod(ps[j,xnodes] == xvals))
                  px <- px+1
              }
            }
            if(py < 1) {
              warning("Can't calculate the probability")
              return(-1)
            }
            return(px/py)
          })

setMethod("cnJointKLdist", "catNetwork",
          function(object1, object2,...) {
            if(!is(object1, "catNetwork") || !is(object2, "catNetwork"))
              stop("catNetwork object is required.")
            if(object1@numnodes != object2@numnodes)
              stop("Number of nodes should be equal.")
            pm1 <- cnJointProb(object1)
            pm2 <- cnJointProb(object2)
            p1 <- pm1[,ncol(pm1)]
            p2 <- pm2[,ncol(pm2)]
            probs <- p1
            probs[p2==0] <- 0
            p2[p2==0] <- 1
            p1[p1==0] <- 1
            return(sum(probs*log(p1/p2)))
          })
