data.prep <- function(obj, Y=NULL, group.size=3)
  {
    if(is.null(Y)) {
        ctrl <- obj$mdata$Y[obj$mdata$Tr==0]
        trt <- obj$mdata$Y[obj$mdata$Tr==1]
      }
    else {
        ctrl <- Y[obj$index.control]
        trt <- Y[obj$index.treated]
      }
    
    if(group.size==3){
      n.i <- length(trt)
      idx <- seq(2,n.i,2)
      trt <- trt[-idx]
      n.t <- length(trt)
      grp.dta <- c()
      j <- 1
      for(i in 1:n.t) {
          tmp <- c(trt[i], ctrl[j:(j+1)])	
          grp.dta <- c(grp.dta, tmp)
          j <- j + 2
          rm(tmp)	
        }
    }
    else {
      n.i <- length(trt)
      idx <- seq(2,n.i,2)
      trt <- trt[-idx]
      trt <- trt[-idx]
      n.t <- length(trt)
      grp.dta <- c()
      j <- 1

      for(i in 1:n.t) {
          tmp <- c(trt[i], ctrl[j:(j+2)])	
          grp.dta <- c(grp.dta, tmp)
          j <- j + 3
          rm(tmp)	
        }
    }
    
    id <- rep(1:n.t, each=group.size)
    zeros <- rep(0,group.size-1)
    trt.ind <- c(1,zeros)
    treat <- rep(trt.ind, times=n.t)
    mctrl <- list()
    mctrl$Y <- grp.dta
    mctrl$id <- id
    mctrl$treat <- treat
    mctrl
  }
