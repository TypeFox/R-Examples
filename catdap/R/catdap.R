.packageName <- "catdap"

#===============================================

catdap1 <- function(cdata, response.names=NULL, plot=1)  {

#===============================================

    if (plot != 1 && plot != 2) plot <- 0

#   cdata                    # original data
    n <- dim(cdata)[2]       # total number of variable
    nsamp <- dim(cdata)[1]   # sample size
    catdata <- array(0, dim=c(n, nsamp))
    for (i in 1:n) catdata[i, ] <- cdata[, i]

#   title                   # variable names
    title <- names(cdata)
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

    dname <- list()
    for (i in 1:n)
      if (is.numeric(cdata[, i]) == FALSE) 
        dname[[i]] <- levels(ordered(cdata[, i]))

    catdap1.out <- catdap01(catdata, title, dname, response.names, plot) 

    class(catdap1.out) <- "catdap1"
    return(catdap1.out)
}

#================================================

catdap1c <- function(ctable, response.names=NULL, plot=1) {

#================================================

    if (plot != 1 && plot != 2) plot <- 0

#  ctable                  # original contingency table data
    id <- dim(ctable)
    n <- length(id)
    y <- as.array(ctable)
    yn <- length(y)
    nsamp <- 0
    for (i in 1:yn)  nsamp <- nsamp + y[[i]]
    cdata <- array(0, dim=c(n, nsamp))

    nn <-  yn
    x <- c(0:(yn-1))
    for (i in n:1) {
      nn <-  nn / id[i]
      nd <- 0
      for (j in 1:yn) {
        m <-  as.integer(x[j] / nn)
        if (y[[j]] != 0) for (k in 1:y[[j]]) cdata[i, nd+k] <- m + 1
        nd <- nd + y[[j]]
        x[j] <- x[j] %% nn
      }
    }

#   title                   # variable names
    title <- names(dimnames(ctable))
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

    dname <- dimnames(ctable)
    if (is.null(dname)) {
      dname <- rep(list(), n)
      for (i in 1:n) dname[[i]] <- as.character(1:id[i])
    }

  catdap1c.out <- catdap01(cdata, title, dname, response.names, plot) 

  class(catdap1c.out) <- "catdap1"
  return(catdap1c.out)
}

#==================================================

catdap01 <- function(cdata, title, dname, response.names, plot) {   # PROGRAM CATDAP 01

#==================================================

    n <- dim(cdata)[1]       # total number of variable
    nsamp <- dim(cdata)[2]   # sample size

#   response variable numbers
    ier <- 0
    if (is.null(response.names) ) {
      ires <- c(1:n)
      l <- n
    } else {
      l <- length(response.names)
      ires <- NULL
      m <- 0
      for (i in 1:l)
        for (j in 1:n)
          if (response.names[i] == title[j]) {
            m <- m + 1
            ires <- c(ires, j)
          }
      if (m == 0) stop(" Error : response variable name is wrong.")
      if (m < l)  ier <- 3001
      l <- m
    }

#   conv                    # variable number and recode numbers
    recode <- 0
    iconv <- array(0, dim=c(20, recode))

#   ex                      # variable numbers of explanatory variables
    n1 <- 0
    iex <- 0

#   minmax
    item1 <- rep(0, n)
    item2 <- rep(0, n)
    for (i in 1:n) item1[i] <- as.integer(min(cdata[i, ]))
    for (i in 1:n) item2[i] <- as.integer(max(cdata[i, ]))
    nc <- rep(0, n)
    for (i in 1:n) nc[i] <- item2[i] - item1[i] + 1
    n4 <- max(nc)

    if (length(dname) == 0) {
      for (i in 1:n)
        if (is.numeric(cdata[i, ])) {
          j1 <- item1[i]
          j2 <- item2[i]
          dname[[i]] <- c(j1:j2)
        }
    } else {
      for (i in 1:n)
        if (is.null(dname[[i]]) == TRUE) {
          j1 <- item1[i]
          j2 <- item2[i]
          dname[[i]] <- c(j1:j2)
        }
    }

#   skip                     # variable number and code number to delete
    iskip1 <- 0
    iskip <- 0
    
    z <- .Call("catdap1",
                 as.integer(nsamp),
                 as.integer(n),
                 as.integer(l),
                 as.integer(recode),
                 as.integer(n1),
                 as.integer(iskip1),
                 as.integer(n4),
                 as.integer(item1),
                 as.integer(item2),
                 as.integer(iconv),
                 as.integer(ires),
                 as.integer(iex),
                 as.integer(iskip),
                 as.integer(cdata))

    nc <- z[[1L]]
    ia <- array(z[[2L]], dim=c(n4, n4, l, n))
    p <- array(z[[3L]], dim=c(n4, n4, l, n))
    total <- array(z[[4L]], dim=c(n4, n))
    aic <- array(z[[5L]], dim=c(l, n))
    for (i in 1:l)  aic[i, ires[i]] <- NA
    ord <- array(z[[6L]], dim=c(l, n-1))
    ia <- aperm(ia, c(2, 1, 3, 4))
    p <- aperm(p, c(2, 1, 3, 4))

    cross.table <- list()
    for (i in 1:l) {
      i1 <- 1
      i2 <- nc[ires[i]]
      cross.table[[i]] <- list()
      for (j in 1:n) {
        if (j == ires[i]){
          cross.table[[i]][[j]] <- NULL
        } else {
          j1 <- 1
          j2 <- nc[j]
          cross.table[[i]][[j]] <- list()
          cross.table[[i]][[j]]$res <- title[ires[i]]
          cross.table[[i]][[j]]$n <- array(ia[j1:j2, i1:i2, i, j], dim=c(j2,i2), dimnames = list(dname[[j]], dname[[ires[i]]]))
          cross.table[[i]][[j]]$p <- array(p[j1:j2, i1:i2, i, j], dim=c(j2,i2))
        }
      }  # for (j in 1:n) end
    }  # for (i in 1:l) end

    etotal <- matrix(list(), n)
    for (j in 1:n) {
      j1 <- 1
      j2 <- nc[j]
      etotal[j] <- list(total[j1:j2, j])
    }

    catdap01.out <- list(tway.table=cross.table, title=title, total=etotal, aic=aic, aic.order=ord, plot=plot, ier=ier)

    return(catdap01.out)
}

#========================

print.catdap1 <- function(x,...) {

#========================

    l <- dim(x$aic)[1]
    n <- dim(x$aic)[2]
    m <- n
    iexp <- 0
    ires <- rep(0, l)

#--------------------------------------------------------------#
    cat("\n < Summary of AIC's for the two-way tables >\n\n")
#--------------------------------------------------------------#
    cat("         Explanatory \n           variables ")
    for (i in 1:l) cat(sprintf("\t%20s", x$tway.table[[i]][[1]]$res))
    cat("\n")
    for (i in 1:(l+2)) cat("---------------------")
    cat("\n")
    for (j in 1:m) {
      jj <- j
      cat(sprintf("%20s", x$title[jj]))
      for (j in 1:l) {
        if (is.na(x$aic[j,jj]) == TRUE) {
          cat("\t                 -  ")
          ires[j] <- jj
        } else {
          if (is.null(x$aic[j,jj]) == FALSE) cat(sprintf("\t            %8.2f", x$aic[j,jj]))
        }
      }
      cat("\n")
    }
    nc <- min(m, n-1)
    if (nc < 6) {
      nc1 <- nc
      nr <- 1
    } else if (nc < 11) {
      nc1 <- as.integer((nc+1) / 2)
      nr <- 2
    } else {
      nc1 <- 5
      nr <- 2
    }

    dname <- list()
    item <- rep(0, n)
    for (i in 1:n)
      if (i != ires[1]) {
        dname[[i]] <- dimnames(x$tway.table[[1]][[i]]$n)[[1]]
        dname[[ires[1]]] <- dimnames(x$tway.table[[1]][[i]]$n)[[2]]
      }
    for (i in 1:n)
      item[i] <- length(dname[[i]])

    for (i in 1:l) {
      i1 <- 1
      i3 <- ires[i]
      i2 <- item[i3]
#---------------------------------------------------------------------------------------#
      cat("\n\n < List of explanatory variables arranged in ascending order of AIC >")
#---------------------------------------------------------------------------------------#
      cat(sprintf("\n\n Response variable  : %s\n\n", x$title[i3]))
      cat("     No.\t  Explanatory\tNumber of categories\t  A I C\t\tDifference of AIC\n")
      cat("        \t  variable   \tof exp. var.\n")
      for (j in 1:9) cat("----------")
      cat("\n")
      as <- 0.0
      for (j in 1:nc) {
        nv <- x$aic.order[i, j]
        aic <- x$aic[i, nv]
        if (j != 1) as <- aic - as
        cat(sprintf("%8i %20s\t%8i\t\t%8.2f\t%8.2f\n", j,x$title[nv], item[nv], aic, as))
        as <- aic
      }  # for (j in 1:nc) end
#--------------------------------------------------------------------------#
      cat("\n < Two-way tables arranged in ascending order of AIC >\n\n")
#--------------------------------------------------------------------------#

      cat(sprintf("\t\t\t\t( %s )\n", x$title[i3]))
      nsamp <- 0
      ptt <- 0.0
      pt <- rep(0, i2)
      cat("        ")
      for (ii in i1:i2) {
        cat(sprintf("\t%s\t", dname[[i3]][ii]))
        nsamp <- nsamp + x$total[[i3]][ii]
      }
      for (ii in i1:i2) {
        pt[ii] <- x$total[[i3]][ii] / nsamp * 100.0
        ptt <- ptt + pt[ii]
      }

      cat("\t\tTotal\n")
      for (j in 1:nc) {
        nv <- x$aic.order[i, j]
        j1 <- 1
        j2 <- item[nv]
        cat(sprintf("( %s )\n", x$title[nv]))
        for( k in j1:j2) {
          cat(sprintf("%s\t", dname[[nv]][k]))
          tc <- 0
          for (ii in i1:i2)
            cat(sprintf("%8i ( %5.1f )", x$tway.table[[i]][[nv]]$n[k,ii], x$tway.table[[i]][[nv]]$p[k,ii]))
          cat(sprintf("%8i ( %5.1f )\n", x$total[[nv]][k], ptt))
        }
        cat("Total\t")
        for (ii in i1:i2) cat(sprintf("%8i ( %5.1f )", x$total[[i3]][ii], pt[ii]))
        cat(sprintf("%8i ( %5.1f )\n\n", nsamp, ptt))
      }  # for (j in 1:nc) end
    }  # for (i in 1:l) end

      if (x$ier == 3001)  cat(" caution : some response variable name is wrong.\n")

      if (x$plot > 0) {
        old.par <- par(no.readonly=TRUE)
        nc <- 0
        for (i in 1:l) {
          i1 <- 1
          if (i > 1) par(ask=TRUE)
          if (x$plot == 1)
            plot.single1(x$title, dname, ires[i], item, x$aic.order[i, ], x$aic[i, ], x$tway.table[[i]], x$total, old.par, nc)
          if (x$plot == 2)
            plot.single2(x$title, dname, ires[i], item, x$aic.order[i, ], x$aic[i, ], x$tway.table[[i]], old.par, nc)
        } 
        par(old.par)

        par(ask=TRUE)
        plot.grshade( x$title, x$aic, ires, nsamp, old.par)
        par(old.par)
      }
}


#==============================================================

catdap2 <-
 function(data, pool=NULL, response.name, accuracy=NULL, nvar=NULL,
 explanatory.names=NULL, pa1=1, pa2=4, pa3=10, print.level=0, plot=1) {

#==============================================================

    if (print.level != 1) print.level <- 0
    if (plot != 1 && plot != 2) plot <- 0

#   data			# original data
    nsamp <- dim(data)[1]	# total number of variable
    n <- dim(data)[2]	# sample size

#   title                   # variable names
    title <- names(data)
    if (is.null(title)) {
      title <- rep(1:n)
      title <- as.character(title)
    }

#   response variable number
    ier <- 0
    ires <- 0
    for (i in 1:n) if (response.name == title[i]) ires <- i
    if (ires == 0) stop(" Error : response variable name is wrong.")

#   accuracy                 # accuracy of measurement
    xx <- accuracy
    if (is.null(accuracy)) xx <- rep(0, n)

#   pooling                    # the way of pooling categories of the variable
    if (is.null(pool))  pool <- rep(1, n)

    cdata <- array(0, dim=c(n, nsamp))
    dname <- list()
    n1 <- 0
    for (i in 1:n)
      if (is.numeric(data[, i])) {
        if (pool[i] != 2) {
          n1 <- n1 + 1
          cdata[i, ] <- data[, i]
          dname[[i]] <- NULL
        } else if (pool[i] == 2) {
          cv <- convi(data[,i])
          cdata[i, ] <- cv$cdata
          dname[[i]] <- cv$cname
        }
      } else {                         # character
        cdata[i, ] <- as.numeric(data[, i])
        pool[i] <- 2
        xx[i] <- 0.0
        dname[[i]] <- levels(ordered(data[, i]))
      }
    if (n1 == n) dname <- NULL

#   minmax                   # minimum and maximum of the code values
    item1 <- rep(0, n)
    item2 <- rep(0, n)

    for (i in 1:n) {
      if (i == ires) {
        if (pool[i] < 0) {   # continuous response variable
          item1[i] <- 1
          item2[i] <- -pool[ires]
          if (xx[i] == 0) {
            xx[i] <- (max(cdata[i, ]) - min(cdata[i, ])) / 50
            cat(sprintf("\n Warning : %d-th accuracy is corrected to %8.2f\n", i, xx[i]))
          }
          pool[i] <- 2
        } else if (pool[i] < 2) {
          stop(" Error : 'pool' is less than 0 for the response variable.")

        } else {   # categorical response variable
          item1[i] <- as.integer(min(cdata[i, ]))
          item2[i] <- as.integer(max(cdata[i, ]))
          pool[i] <- 2
          xx[i] <- 0
        }
      } else { # if (i != ires)
        if (pool[i] == 2) {
          xx[i] <- 0
          item1[i] <- as.integer(min(cdata[i, ]))
          item2[i] <- as.integer(max(cdata[i, ]))
        } else {
          if (pool[i] < 0) pool[i] <- 1
          if (xx[i] == 0) {   # continuous explanatory variable
            xx[i] <- (max(cdata[i, ]) - min(cdata[i, ])) / 50
            cat(sprintf("\n Warning : %d-th accuracy is corrected to %8.2f\n", i, xx[i]))
          }
          item1[i] <- 0
          item2[i] <- 0
        } # if(pool[i] == 2) end
      }  # if (i == ires) end
    } # for (i in 1:n) end

    l <- 1
    recode <- 0
    iskip1 <- 0

    it <- 0
    for (i in 1:n)
      if (xx[i] != 0) {
        it <- 1
        break
      }
    xx.in <- rep(0, n)
    for (i in 1:n) xx.in[i] <- xx[i]

    ipart <- 0
    iconv <- rep(0, 20)
    face <- ires
    iskip <- rep(0, 20)
    isk <- rep(0, 2)
    sk <- rep(0, 20)
    ida <- array(0, dim=c(n, nsamp))
    da <- array(0, dim=c(n, nsamp))
    if (it == 0) ida <- cdata
    if (it == 1) da <- cdata

    if (is.null(nvar)) nvar <- n
    if (nvar <2) stop(" Error : 'nvar' is greater than or equal to 2.")
    if (nvar > n) stop(" Error : 'nvar' is less than or equal to the number of variable.")
    novv <- 0
    for (i in 1:n)
      if (pool[i] == 0 || pool[i] > 0) novv <- novv + 1
    nov <- min(nvar, novv)

#   additional contingency table 
    icls <- rep(0, max(10, nov+1))
    iwa <- 0
    if (is.null(explanatory.names)) {
      icl <- 0
    } else {
      icl <- 1
      lex <- length(explanatory.names)
      if (lex > (n-1)) {
        iwa <- 3001
        icl <- 0
      } else {
        icls[1] <- lex + 1
        icls[2] <- ires
        m <- 0
        for (i in 1:lex)
          for (j in 1:n)
            if (explanatory.names[i] == title[j]) {
              icls[2 + i] = j
              m <- m + 1
            }
        if (m != lex) {
          iwa <- 3001
          icl <- 0
        } else if (m > max((nov-1),8)) {
          iwa <- 3002
          icl <- 0
        }
      }
    }

    n11 <- item2[ires] - item1[ires] + 1
    n2 <- 0
    for (i in 1:n) {
      n1 <- item2[i] - item1[i] + 1
      n2 <- max(n2, n1)
    }
    n33 <- 100
    if (nsamp / 2 + 1 < n33) n33 <- nsamp / 2 + 2
    n33 <- max(n33, n2)
    if (it == 0) n33=n2

    ikr <- n * pa2
    jkr <- pa3
    icl1 <- icl + 1

### border check
    ikkk <- n33 ** 3 * pa1 * pa1
    n31 <- n * ikkk * icl1
    ddmax <- 2 ** 24 - 1
    ikkkm <- as.integer(ddmax / n / icl1)
    if (ikkk < 2**7) ikkk = 2 ** 7
    if (ikkk > ikkkm) ikkk <- ikkkm
    
    eps <- 1.0e-10   # an error tolerance in the difference of AIC's

    z <- .Call("catdap2m",
	     as.integer(nsamp),
	     as.integer(n),
	     as.integer(l),
	     as.integer(recode),
	     as.integer(iskip1),
	     as.integer(it),
	     as.integer(ipart),
	     as.integer(nov),
	     as.integer(icl),
	     as.integer(item1),
	     as.integer(item2),
	     as.integer(pool),
	     as.integer(iconv),
	     as.integer(ires),
	     as.integer(iskip),
	     as.integer(isk),
	     as.double(sk),
	     as.integer(icls),
	     as.double(xx.in),
	     as.integer(ida),
	     as.double(da),
	     as.integer(n11),
	     as.integer(n33),
	     as.integer(ikr),
	     as.integer(jkr),
	     as.integer(ikkk),
               as.double(eps))

    ier <- z[[22L]]
 
    if (ier[1] == 0) {
      lk77 <- z[[14L]]
      iab <- array(z[[1L]], dim=c(n, n11, n33))

      totalc <- z[[2L]]
      ttrr <- array(z[[3L]], dim=c(n, n33))
      ab <- array(z[[4L]], dim=c(n, n33))
      iaa <- array(z[[5L]], dim=c(n, n33))
      pa <- array(z[[6L]], dim=c(n, n11, n33))
      idata <- z[[7L]]
      ite <- z[[8L]]
      aic <- z[[9L]]
      aaam <- z[[10L]][2:lk77]
      caa <- array(z[[11L]], dim=c(ikr, jkr))
      caa <- caa[1:lk77, 1:jkr]
      icaa <- z[[12L]][1:lk77]

      nnia <- z[[13L]][2:lk77]
      morder <- z[[15L]][1:lk77]

      iby <- array(z[[16L]], dim=c(nov ,ikkk, icl1))
      ibc <- array(z[[17L]], dim=c(n11, ikkk, icl1))
      pbc <- array(z[[18L]], dim=c(n11,ikkk,icl1))
      aic1 <- z[[19L]]
      iabse <- array(z[[20L]], dim=c(nov, n33, icl1))
      baic <- z[[21L]]

      j1 <- 1
      if (is.null(dname)) {
        for (j in 1:n) {
          j2 <- ite[j]
          dname[[j]] <- c(j1:j2)
        }
      } else {
        if (length(dname) == n) {
          for (j in 1:n)  if (is.null(dname[[j]])) {
            j2 <- ite[j]
            dname[[j]] <- c(j1:j2)
          }
        } else {
          for (j in 1:(n-1))  if (is.null(dname[[j]])) {
            j2 <- ite[j]
            dname[[j]] <- c(j1:j2)
          }
          j2 <- ite[n]
          dname[[n]] <- c(j1:j2)
        }
      }

      iab <- aperm(iab, c(1, 3, 2))
      pa <- aperm(pa, c(1, 3, 2))
      cross.table <- list()
      i1 <- 1
      i2 <- ite[ires]
      for (j in 1:n) {
        if (j == ires){
          cross.table[[j]] <- NA
        } else {
          j1 <- 1
          j2 <- ite[j]
          cross.table[[j]] <- list()
          cross.table[[j]]$res <- title[ires]
          cross.table[[j]]$n <- array(iab[j, j1:j2, i1:i2], dim=c(j2,i2), dimnames = list(dname[[j]], dname[[ires]]))
          cross.table[[j]]$p <- array(pa[j, j1:j2, i1:i2], dim=c(j2,i2))
        }
      }

      total <- matrix(list(), n)
      for (j in 1:n) {
        j1 <- 1
        j2 <- ite[j]
        if (j == ires) total[j] <- list(totalc[j1:j2])
        if (j != ires) total[j] <- list(ttrr[j, j1:j2])
      }

      nint <- max(ite) + 1
      interval <- list()
      for (j in 1:n) {
        j3 <- idata[j]
        j1 <- 1
        j2 <- ite[j3]
        if (xx[j3] != 0.0) interval[[j3]] <- ab[j3, j1:(j2+1)]
        if (xx[j3]==0.0) interval[[j3]] <- as.integer(iaa[j3, j1:j2])
      }
      aic.order <- idata[2:n]

      aorder <- NULL
      for (i in 1:(lk77-1))
        if (morder[i] > 1) aorder <- c(aorder, morder[i] - 1)
      nsub <- length(aorder)

      nv <- rep(0, nsub)
      ncc <- rep(0, nsub)
      aaic <- rep(0, nsub)
      icaa <- icaa[2:lk77]

      for (i in 1:nsub) {
        j <- aorder[i]
        nv[i] <- icaa[j] - 1
        ncc[i] <- nnia[j]
        aaic[i] <- aaam[j]
      }

      cexp <- list()
      for (i in 1:nsub) {
        j <- aorder[i] + 1
        k <- nv[i] + 1
        cexp[[i]] <- aic.order[caa[j, 2:k] - 1]
      }

      iby <- aperm(iby, c(2, 1, 3))
      ibc <- aperm(ibc, c(2, 1, 3))
      pbc <- aperm(pbc, c(2, 1, 3))

      np1 <- ncc[1]
      nn <- icaa[aorder[1]]

      ctbl <- list()
      ctbl.cnum <- list()
      ctbl.p <- list()
      caic <- rep(0, icl1)
      cinterval <- list()
      idx <- 0
	  
      ctbl.cnum[[icl1]] <- iby[1:np1, 2:nn, 1]
      ctbl[[icl1]] <- ibc[1:np1, i1:i2, 1]
      ctbl.p[[icl1]] <- pbc[1:np1, i1:i2, 1]
      caic[icl1] <- aic1[1]

      if (icl > 0) {   ## The output of the additional analysis 
        np1 <- 1
        idx <- idx + 1
        cint <- NULL
        intval <- list()
        nvv <- icls[1] - 1
        for (i in 1:nvv) {
          exv <- icls[i+2]
          intv <- NULL
          j1 <- 1
          j2 <- ite[exv]
          if (xx[exv] != 0.0) {
            p <- iabse[(i+1), j1:(j2+1), 2]
            nint <- length(p)
            for (k in 1:nint) 
              if (p[k] != 0) intv <- c(intv, interval[[exv]][p[k]])
            np1 <- np1 * (length(intv) - 1)
          } else if (xx[exv]==0.0) {
            intv  <- as.integer(iaa[exv, j1:j2])
            np1 <- np1 * j2
          }
          intval[[i]] <- intv
          cint <- c(cint, intval[i])
        }
        cinterval[[idx]] <- cint
        nn <- icls[1]
        ctbl.cnum[[1]] <- iby[1:np1, 2:nn, 2]
        ctbl[[1]] <- ibc[1:np1, i1:i2, 2]
        ctbl.p[[1]] <- pbc[1:np1, i1:i2, 2]
        caic[1] <- aic1[2]
      }  #  if (icl > 0) end

      idx <- idx + 1    ## The output of the Minimum AIC Model 
      cint <- NULL
      intval <- list()
      for (i in 1:nv[1]) {
        exv <- cexp[[1]][i]
        intv <- NULL
        j1 <- 1
        j2 <- ite[exv]
        if (xx[exv] != 0.0) {
          p <- iabse[(i+1), j1:(j2+1), 1]
          nint <- length(p)
          for (k in 1:nint) 
            if (p[k] != 0) intv <- c(intv, interval[[exv]][p[k]])
          np1 <- np1 * (length(intv) - 1)
        } else if (xx[exv]==0.0) {
          intv  <- as.integer(iaa[exv, j1:j2])
          np1 <- np1 * j2
        }
        intval[[i]] <- intv
      cint <- c(cint, intval[i])
      }    # for (i in 1:nv[1]) end
      cinterval[[idx]] <- cint

      iaddflg <- 0
      if (nsub > 2) { 
        aicm1 <- aic1[1]
        nbest <- 1
        for (i in 2:nsub) {
          aicm2 <- aaic[i]
          daic <- abs(aicm2 - aicm1) / max(abs(aicm2), abs(aicm1))
          if (daic > eps) break
          nbest <- nbest + 1

          if (iaddflg == 0) {

          z1 <- addaicm(i, nsamp, n, l, recode, iskip1, it, ipart, nov, item1, item2, pool, iconv, ires, iskip, isk, sk, xx, ida, da, n11, n33, pa1, ikr, jkr, ikkk, ikkkm, cexp, eps)

          if (z1$ier == 0) {
            iaddflg <- 1
            idx <- idx + 1
            np1 <- ncc[nbest]
            idflg <- is.null(dim(z1$ctable$cnum[[1]]))
            if (idflg == TRUE) ccnum <- z1$ctable$cnum[[1]][1:np1]
            if (idflg == FALSE) ccnum <- z1$ctable$cnum[[1]][1:np1, ]
            ctbl.cnum[[idx]] <- ccnum
            ctbl[[idx]] <- z1$ctable$n[[1]]
            ctbl.p[[idx]] <- z1$ctable$p[[1]]
            caic <- c(caic, aicm2)
            iabse <- z1$iabse

            cint <- NULL
            intval <- list()
            for (j in 1:nv[i]) {
              exv <- cexp[[i]][j]
              intv <- NULL
              j1 <- 1
              j2 <- ite[exv]
              if (xx[exv] != 0.0) {
                p <- iabse[(j+1), j1:(j2+1)]
                nint <- length(p)
                for (k in 1:nint) 
                  if (p[k] != 0 ) intv <- c(intv, interval[[exv]][p[k]])
              } else if (xx[exv]==0.0) {
                intv  <- as.integer(iaa[exv, j1:j2])
              }
              intval[[j]] <- intv
              cint <- c(cint, intval[j])
            } # for (j in 1:nv[i]) end
            cinterval[[idx]] <- cint
          } # if (z1$ier == 0) end
          } # if (iaddflg == 0) end 

        }  # if (i in 2:nsub) end
      }  # if (nsub > 2) end

      if (iwa != 0) ier <- c(ier, iwa)

      if (nbest > 1) {
        nr <- nbest - 1- iaddflg
        if (nr ==1 || nr > 1) ier <- c(ier, nr)
      }

      catdap2.out <- list(tway.table=cross.table, title=title, pool=pool, accuracy=xx, ires=ires, 
                          icl=icl, icls=icls,total=total, interval=interval, base.aic=baic, aic=aic,
                          aic.order=aic.order, nsub=nsub, nv=nv, ncc=ncc, aaic=aaic, cexp=cexp,
                          ctable=list(n=ctbl, p=ctbl.p, cnum=ctbl.cnum), ctable.interval=cinterval,  
                          caic=caic, ier=ier, print.level=print.level, plot=plot)

    } else { # if (ier[1] != 0)
      if (ier[1] == 2048) { # pa1
        kkj <- ier[2]
        if (kkj < ikkkm) {
          pa1n <- ceiling( sqrt(kkj / (n33 ** 3)))
          ier[2] <- pa1n
          ikkkn <- n33 ** 3 * pa1n * pa1n
          if (ikkkn > ikkkm) ier[2] <- -999    
        } else {
          ier[2] <- -999
        }  # if (kkj < ikkkm) end
      }  # if (ier[1] == 2048) end

      catdap2.out <- list(ier=ier)

    }

    class(catdap2.out) <- "catdap2"

    return(catdap2.out)
}

#================================================================

addaicm <- 
 function(iex2, nsamp, n, l, recode, iskip1, it, ipart, nov, item1, item2, pool,
 iconv, ires, iskip, isk, sk, xx, ida, da, n11, n33, pa1, ikr, jkr, ikkk, ikkkm, cexp, eps) {

#================================================================

    icl <- 1
    icl1 <- icl + 1
    nex <- length(cexp[[iex2]])
    icls <- c(nex+1, ires)
    for (i in 1:nex) icls <- c(icls, cexp[[iex2]][i])

    z <- .Call("catdap2m",
                as.integer(nsamp),
                as.integer(n),
                as.integer(l),
                as.integer(recode),
                as.integer(iskip1),
                as.integer(it),
                as.integer(ipart),
                as.integer(nov),
                as.integer(icl),
                as.integer(item1),
                as.integer(item2),
                as.integer(pool),
                as.integer(iconv),
                as.integer(ires),
                as.integer(iskip),
                as.integer(isk),
                as.double(sk),
                as.integer(icls),
                as.double(xx),
                as.integer(ida),
                as.double(da),
                as.integer(n11),
                as.integer(n33),
                as.integer(ikr),
                as.integer(jkr),
                as.integer(ikkk),
                as.double(eps))

    ier <- z[[22L]]
 
    if (ier[1] == 0 ) {
      ite <- z[[8L]]
      iby <- array(z[[16L]], dim=c(nov ,ikkk, icl1))
      ibc <- array(z[[17L]], dim=c(n11, ikkk, icl1))
      pbc <- array(z[[18L]], dim=c(n11,ikkk,icl1))
      iabse <- array(z[[20L]], dim=c(nov, n33, icl1))

      iby <- aperm(iby, c(2, 1, 3))
      ibc <- aperm(ibc, c(2, 1, 3))
      pbc <- aperm(pbc, c(2, 1, 3))

      i1 <- 1
      i2 <- ite[ires]

      np1 <- 1
      for (i in 1:(icls[1]-1)) np1 <- np1 * ite[icls[2+i]]
      nn <- icls[1]

      jby <- list(iby[1:np1, 2:nn, 2])
      jbc <- list(ibc[1:np1, i1:i2, 2])
      qbc <- list(pbc[1:np1, i1:i2, 2])
      iabse <- iabse[, , 2]

      return(list(ctable=list(n=jbc, p=qbc, cnum=jby), iabse=iabse, ier=0))

    } else {
      if (ier[1] == 2048) { # pa1
        kkj <- ier[2]
        if (kkj < ikkkm) {
          pa1n <- ceiling( sqrt(kkj / (n33 ** 3)))
          ikkkn <- n33 ** 3 * pa1n * pa1n
          ier[2] <- pa1n
          if (ikkkn > ikkkm) ier[2] <- -999    
        } else {
          ier[2] <- -999
        }  # if (kkj < ikkkm) end
        if (ier[2] == -999) cat("Error : Working area for second contingency table is too short. pa1 can no longer set a larger value.\n")
        if ( ier[2] != -999) cat(sprintf(" Error : Working area for second contingency table is too short, try pa1= %i.\n", pa1n))
      }  # if (ier[1] == 2048) end
      return(list(ier=ier[1]))
    }
}

convi <- function(data) {
    x1 <- unique(data)
    y1 <- rank(x1)
    n1 <- length(data)
    n2 <- length(y1)
    cdata <- rep(0, n1)
    for (i in 1:n1)
      for (j in 1:n2) if (data[i] == x1[j]) cdata[i] <- y1[j]
    cname <- sort(x1)
    return(list(cdata=cdata, cname=cname))
} 

#========================

print.catdap2 <- function(x,...) {

#========================

    eps <- 1.0e-10   # an error tolerance in the difference of AIC's
    ier <- x$ier[1]

    if (ier == 0) {
      n <- length(x$title)
      res <- x$ires
      lk <- x$nsub
      ktt <- max(x$nv)

      dname <- list()
      nc <- rep(0, n)
      for (i in 1:n)
        if (i != res) { 
        dname[[i]] <- dimnames(x$tway.table[[i]]$n)[[1]]
        dname[[res]] <- dimnames(x$tway.table[[i]]$n)[[2]]
      }
      for (i in 1:n)
        nc[i] <- length(dname[[i]])

#----------------------------------------------------------------------------#
      cat("\n<< List of single explanatory variables (arranged in ascending order of AIC) >>\n")
#----------------------------------------------------------------------------#
      cat(sprintf(" Response variable : %s \t(base AIC = %8.2f)\n", x$title[res], x$base.aic))
      cat("----------------------------------------------------------------------------------------\n")
      cat("\t\t\t\tNumber of\n\t\tExplanatory\tcategories\t\t\tDifference\n\t\tvariables\tof exp. var.\t  A I C     \tof AIC\t\tWeight\n")
      cat("----------------------------------------------------------------------------------------\n")
 
      nrank <- 0
      daic2 <- 111

      for (i in 1:(n-1)) {
        idx <-x$aic.order[i]
        aaaa <- x$aic[idx]
        if (i == 1) aaa1 <- x$aic[idx]
        daic <- aaaa - aaa1
        w <- exp(-1. / 2 * daic)
        if (i > 1) daic2 <- abs(aaaa - aaa2) / max(abs(aaaa), abs(aaa2))
        if (daic2 > eps) nrank <- i
        cat(sprintf("%5i\t", nrank))
        cat(sprintf("%20s\t", x$title[idx]))
        cat(sprintf("%8i    \t%8.2f  \t%8.2f  \t%8.2f\n", nc[idx], x$aic[idx], daic, w))
        aaa2 <- aaaa
      }

#-----------------------------------------------------------------------#
      cat("\n\n<< Two-way tables arranged in ascending order of AIC >>\n\n")
#-----------------------------------------------------------------------#

      i <- 1
      i1 <- 1
      i2 <- nc[res]
      nsamp <- sum(x$total[[res]][i1:i2])
      ptc <- rep(0, i2)
      for (i in i1:i2) ptc[i] <- (x$total[[res]][i] * 100.) / nsamp
      ptt <- sum(ptc[i1:i2])

      cat("<Note>\n")
      cat(sprintf(" %s\n", x$title[res]))
      if (x$accuracy[res] != 0.0) {
        cat("\tcategory    \tvalue range\n")
        for (i in i1:i2) cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", i, x$interval[[res]][i], x$interval[[res]][i+1]))
      } else {
        cat("\tcategory    \tvariable value\n")
        for (i in i1:i2) cat(sprintf("\t%8i    \t%s\n", i, dname[[res]][i]))
      }
      cat("\n\n")

      cat(sprintf("                     ( %15s )\n", x$title[res]))
      for (j in 1:(n-1)){
        j3 <- x$aic.order[j]
        j1 <- 1
        j2 <- nc[j3]
        cat(sprintf("                 ", x$title[j3]))
        for (ii in i1:i2)  cat(sprintf("%15s   ", dname[[res]][ii]))
        cat("          Total\n")
        cat(sprintf("( %15s )\n", x$title[j3]))
        for (ii in i1:(i2+2))  cat("-----------------")
        cat("\n")
        for (k in j1:j2) {
          ptr <- 0
          if (x$total[[j3]][k] != 0) ptr <- (x$total[[j3]][k] * 100.) / x$total[[j3]][k]
          cat(sprintf("%8i           ", k))
          for (ii in i1:i2) cat(sprintf("%8i ( %5.1f )", x$tway.table[[j3]]$n[k,ii], x$tway.table[[j3]]$p[k,ii]))
          cat(sprintf("%8i ( %5.1f )\n", x$total[[j3]][k], ptr))
        }  # for (k in j1:j2) end

        for (ii in i1:(i2+2))  cat("-----------------")
        cat("\n     Total         ")
        for (ii in i1:i2) cat(sprintf("%8i ( %5.1f )", x$total[[res]][ii], ptc[ii]))
        cat(sprintf("%8i ( %5.1f )\n\n", nsamp, ptt))

        cat("<Note>\n")
        cat(sprintf(" %s\n", x$title[j3]))
        if (x$accuracy[j3] != 0.0) {
          cat("\tcategory    \tvalue range\n")
          for (k in j1:j2)
            cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", k,x$interval[[j3]][k], x$interval[[j3]][k+1]))
        } else {
          cat("\tcategory    \tvariable value\n")
          for (k in j1:j2)  cat(sprintf("\t%8i    \t%s\n", k, dname[[j3]][k]))
        }
        cat("\n\n")
      }  # for (j in 1:(n-1)) end

      if (x$print.level == 0) {
#----------------------------------------------------------------------------------#
        cat("\n<< AIC's of the models with k explanatory variables (k=1,2,...) >>\n")
#----------------------------------------------------------------------------------#
        for (j in 1:ktt) {
          ijk <- 0
          nrank <- 0
          daic2 <- 111

          for (ij in 1:lk) {
            lk4 <- x$nv[ij]
            if (lk4 == j) {
              aaaa <- x$aaic[ij]
              if (ijk == 0) aaic1 <- x$aaic[ij]
              daic <- aaaa - aaic1
              w <- exp(-1. / 2 * daic)
              ijk <- ijk + 1

              if (ijk == 1) {
                cat(sprintf("\n Number of explanatory variables = %i\n", j))
                cat("----------------------------------------------------------------------------------------\n")
                cat("\t\t\t\tNumber of\n\t\t Explanatory\tcategories\t\t\tDifference\n\t\t variables\tof exp. var.\t  A I C\t\tof AIC\t\t  Weight\n")
                cat("----------------------------------------------------------------------------------------\n")
              }

              if (ijk > 1) daic2 <- abs(aaaa - aaic2) / max(abs(aaaa), abs(aaic2))
              if (daic2 > eps) nrank <- ijk
              if (lk4 == 1 || lk4 > 1)  {
                cat(sprintf("%5i ", nrank))
                lkk <- x$cexp[[ij]][1]
                cat(sprintf("%20s\t", x$title[lkk]))
                cat(sprintf("%8i\t%8.2f\t%8.2f\t%8.2f\n", x$ncc[ij], x$aaic[ij], daic, w))
              }

              if (lk4 == 2 || lk4 > 2)
                for (i in 2:lk4) {
                lkk <- x$cexp[[ij]][i]
                cat(sprintf("      %20s\n", x$title[lkk]))
              }

            }  # if (lk4 == j) end
            aaic2 <- aaaa
          }  # for (ij in 1:lk) end
        }  # for (j in 1:ktt) end
      }  # if (x$print.level == 0) end

#-------------------------------------------------------------------#
      cat("\n\n<< Summary of subsets of explanatory variables >>\n\n")
#-------------------------------------------------------------------#

      cat(sprintf(" Response variable : %s\n\n", x$title[res]))
      cat("----------------------------------------------------------------------------------------\n")
      cat("\t\t\t\tNumber of\n\t\tExplanatory\tCategories\t\t\tDifference\n\t\tvariables\tof exp. var.\t  A I C\t\tof AIC\t\t  Weight\n")
      cat("----------------------------------------------------------------------------------------\n")
      ijk <- 0
      nrank <- 0
      daic2 <- 111
      if (x$print.level == 0) lsub <- lk
      if (x$print.level == 1) lsub <- min(lk, 30)
      for (ij in 1:lsub) {
        lk4 <- x$nv[ij]
        if (lk4 > 0) {
          aaaa <- x$aaic[ij]
          if (ijk == 0 ) aaic1 <- x$aaic[ij]
          daic <- aaaa - aaic1
          w <- exp(-1. / 2 * daic)
          ijk <- ijk + 1
          if (ijk > 1) daic2 <- abs(aaaa - aaic2) / max(abs(aaaa), abs(aaic2))
          if (daic2 > eps) nrank <- ijk
          cat(sprintf("%5i ", nrank))
          lkk <- x$cexp[[ij]][1]
          cat(sprintf("%20s", x$title[lkk]))
          cat(sprintf("%8i\t\t%8.2f\t%8.2f\t%8.2f\n", x$ncc[ij], x$aaic[ij], daic, w))
          if (lk4 == 2 || lk4 > 2)
            for (i in 2:lk4) {
              lkk <- x$cexp[[ij]][i]
              cat(sprintf("      %20s\n", x$title[lkk]))
            }
          aaic2 <- aaaa
        }  # if (lk4 > 0) end
      }  # for (ij in 1:lsub) end

      iplt <- x$plot
      icl <- x$icl
      if (iplt > 0) old.par <- par(no.readonly = TRUE)
      nwindow <- 0
      ncc <- length(x$ctable$cnum)
      for (k in 1:ncc) {
        kkk <- k - icl
        if (kkk == 0) {
#---------------------------------------------------------------------------------------#
          cat("\n\n<< The output of the additional analysis >>\n\n")
#---------------------------------------------------------------------------------------#
          lk4 <- x$icls[1] - 1
        } else {
#---------------------------------------------------------------------------------------#
          cat("\n\n<< Contingency table constructed by the best subset of explanatory variables >>\n\n")
#---------------------------------------------------------------------------------------#
          lk4 <- x$nv[kkk]
        }

        cat(sprintf(" X(1) : %s\n", x$title[res]))
        j1 <- 1
        j2 <- rep(0, lk4)
        mc <- rep(0, lk4)
        etitle <- list()
        for (i in 1:lk4) {
          if (kkk == 0) {
            lkk <- x$icls[i+2]
          } else {
            lkk <- x$cexp[[kkk]][i]
          }
          cat(sprintf(" X(%i) : %s\n", (i+1), x$title[lkk]))
          mc[i] <- lkk
          etitle <- c(etitle, x$title[lkk])
          j2[i] <- length(x$ctable.interval[[k]][[i]])
          if (x$accuracy[lkk] != 0) j2[i] <- j2[i] - 1
        }
        cat("\n")

        mp <- 1
        for (i in 1:lk4) mp <- mp * j2[i]
        idm <- j2[1]
        if (lk4 > 1) for (i in 2:lk4) idm <- c(idm, j2[i])

        for (i in 1:lk4) cat(" X  ")
        cat("\t\t response variable X(1)\n")
        for (i in 1:lk4) cat(sprintf("(%i) ", (i+1)))
        for (i in i1:i2) cat(sprintf("           %3i      ", i))
        cat("     Total\n")
        for (ii in i1:(i2+2))  cat("-----------------")
        cat("\n")
        np <- dim(x$ctable$cnum[[k]])[1]
        if (is.null(np)) np <- length(x$ctable$cnum[[k]])


        idf <- is.null(dim(x$ctable$cnum[[k]])) 
        for (i in 1:np) {
            if (idf == TRUE) cat(sprintf(" %2i ", x$ctable$cnum[[k]][i]))
            if (idf == FALSE) for (j in 1:lk4) cat(sprintf(" %2i ", x$ctable$cnum[[k]][i,j]))
            for (j in i1:i2) cat(sprintf("%8i ( %5.1f )", x$ctable$n[[k]][i,j], x$ctable$p[[k]][i,j]))
            ibct <- sum(x$ctable$n[[k]][i, i1:i2])
            pbct <- sum(x$ctable$p[[k]][i, i1:i2])
            cat(sprintf("%8i ( %5.1f )\n", ibct, pbct))
        }  # end for (i in 1:np)

        for (ii in i1:(i2+2))  cat("-----------------")
        cat("\n   Total")
        if (idf == FALSE && lk4 > 2) for (j in 1:(lk4-2)) cat("    ")
        for (j in i1:i2 ) cat(sprintf("%8i ( %5.1f )", x$total[[res]][j], ptc[j]))
        cat(sprintf("%8i ( %5.1f )\n", nsamp, ptt))

        cat("\n<Note>\n")
        cat(sprintf("X(1) : %s\n", x$title[res]))
        ctable.dname <- list()

        if (x$accuracy[res] != 0.0) {
          cat("\tcategory    \tvalue range\n")
          for (i in i1:i2) cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", i, x$interval[[res]][i], x$interval[[res]][i+1]))
          ctable.dname[[1]] <- c(i1:i2)
        } else {
          cat("\tcategory    \tvariable value\n")
          for (i in i1:i2) cat(sprintf("\t%8i    \t%s\n", i, dname[[res]][i]))
          ctable.dname[[1]] <- dname[[res]]
        }

        for (i in 1:lk4) {
          if (kkk == 0) {
            lkk <- x$icls[2+i]
          } else {
            lkk <- x$cexp[[kkk]][i]
          }
          cat(sprintf("X(%i) : %s\n", (i+1), x$title[lkk]))

          if (x$accuracy[lkk] != 0.0) {
            cat("\tcategory    \tvalue range\n")
            for (j in j1:j2[i]) 
              cat(sprintf("\t%8i    \t%12.5e   -   %12.5e\n", j, x$ctable.interval[[k]][[i]][j], x$ctable.interval[[k]][[i]][j+1]))
            ctable.dname[[i+1]] <- c(j1:j2[i])
          } else {
            cat("\tcategory    \tvariable value\n")
            for (j in j1:j2[i]) cat(sprintf("\t%8i    \t%s\n", j, dname[[lkk]][j]))
            ctable.dname[[i+1]] <- dname[[lkk]]
          }
        }   #  for (i in 1:lk4) end

        cat(sprintf("\nAIC = %8.2f\n\n", x$caic[k]))

        if (iplt > 0)
          nwindow <- plot.mosaic(iplt, res, icl, k, mp, np, idm, x$title, etitle, dname, x$ctable$n, ctable.dname, nc, x$aic.order, x$aic, x$tway.table, x$total, nwindow, old.par)

      } # for (k in  1:ncc) end

      if (iplt > 0)  par(ask=old.par$ask, mfcol=old.par$mfcol, mai=old.par$mai, cex.main=old.par$cex.main)

    } else if (ier != 0) {
      print.err2(x$ier)
    } # if (ier == 0) end
 
    if (length(x$ier) > 4) {
      iwa <- x$ier[5]
      if (iwa ==1) cat("<NOTE> There is another subset with the minimum AIC.\n\n")
      if (iwa > 1) cat(sprintf("<NOTE> There are another %3i subsets with the minimum AIC.\n\n", iwa))
    }
    if (length(x$ier) > 3) {
      iwa <- x$ier[4]
      if (iwa == 3001) {
        cat(" warning : explanatory variable name in additional contingency table is wrong.\n\n")
      } else if (iwa == 3002) {
        cat(sprintf(" warning : the number of explanatory variable name in additional contingency table is less than %i.\n\n", max(n,9)))
      } else {
        if (iwa ==1) cat("<NOTE> There is another subset with the minimum AIC.\n\n")
        if (iwa > 1) cat(sprintf("<NOTE> There are another %3i subsets with the minimum AIC.\n\n", iwa))
      }
    }
}

#======================

print.err2 <- function(err) {

#======================

    ier <- err[1]
    eval <- err[2]

    if (ier == 2002) {
      cat(" Error : Working area for multi-dimensional table is not enough. Try larger pa2.\n")
    } else if (ier == 2003) {
      cat(" Error : Working area for multi-dimensional table is not enough. Try larger pa3.\n")
    } else if (ier == 2048) {
      if (eval != -999) cat(sprintf(" Error : Working area for contingency table is too short, try pa1= %i.\n", eval))
      if (eval == -999) cat("Error : Working area for contingency table is too short. pa1 can no longer set a larger value.\n")
    } else if (ier == 2037) {
      cat(" caution : the program catdap cannot deal with data sets where the number of non-zero frequency categories of the response variables is less than 2.\n")
    } else if (ier == 2035) {
      cat(" Error : the value of variable is beyond the interval specified in 'min' and 'max'.\n")
    } else if (ier == 2588) {
      cat(" Error :  lk5 > n-1 \n")
    } else if (ier == 650) {
      cat(" Error : the value of 'nvar' is too small for the additional analysis.\n")
    }
}


#====================================================

# plot function (plot.single1, plot.single2, plot.mosaic, plot.grshade)

#====================================================

plot.single1 <- function(title, dname, res, item, aic.order, aic, tway.table, total, old.par, nw) {

    n <- length(title)
    i2 <- item[res]

#---------------
    if (n < 6) {
      nc <- n
      nr <- 1
    } else if (n < 11) {
      nc <- as.integer((n+1) / 2)
      nr <- 2
    } else {
      nc <- 5
      nr <- 2
    }
    m <- nc * nr

    new.mai <- old.par$mai
    new.mai[1] <- new.mai[3] * 0.8
    new.mai[3] <- new.mai[3] * 0.5
    newcex.main <- old.par$cex.main * 0.8  
    mtitle <- "Single Explanatory Models\nin ascending order of AIC"

    nplot <- 1
    irflag <- 0

    j1 <- 1
    for (j in 1:(n-1)) {
      iex <- aic.order[j]
      ctbl <- tway.table[[iex]]$n

      if (aic[iex] > 0 && irflag == 0) {
        if (nplot%%m == 1) {
          if (j > 1 || nw > 0) dev.new()
          par(mfcol=c(nc, nr), mai = new.mai)
          nw <- nw + 1
        }
        irflag <- 1
        y <- total[[res]]
        y <- as.array(y)
        xlabel <- paste(title[res], " Total ", total[[res]][1])
        if (i2 > 1) for (ii in 2:i2) xlabel <- paste(xlabel, ":", total[[res]][ii])
        dimnames(y) <- list(dname[[res]])
        mosaicplot(y, color=TRUE, main="", ylab="AIC= 0.0\n ", xlab=xlabel,  off=0)
        nplot <- nplot + 1
      }

      if (nplot%%m == 1) {
        if (j > 1 || nw > 0) dev.new()
        par(mfcol=c(nc, nr), mai = new.mai, cex.main = newcex.main)
        nw <- nw + 1
      }

      j2 <- item[iex]
      y <- array(0, dim=c(j2, i2))
      ncc <- 0
      for (k in j1:j2)
        if (total[[iex]][k] != 0) {
          ncc <- ncc + 1
          for (ii in 1:i2)  y[ncc, ii] <- ctbl[k, ii]
        }
      dimnames(y)[1] <- list(dimnames(ctbl)[[1]])
      dimnames(y)[2] <- list(dname[[res]])
      aaic <- round(aic[iex], 2)
      mosaicplot(y, color=TRUE, main=mtitle, ylab=paste("AIC=", aaic, "\n", title[iex]), xlab=paste(title[res]),  dir=c("h", "v"), off=c(10, 0))
      mtitle <- ""
      nplot <- nplot + 1
    }  # for (j in 1:(n-1)) end

    if (irflag == 0) {
      y <- total[[res]]
      y <- as.array(y)
      xlabel <- paste(title[res], " Total ", total[[res]][1])
      if (i2 > 1) for (ii in 2:i2) xlabel <- paste(xlabel, ":", total[[res]][ii])
      dimnames(y) <- list(dname[[res]])
      mosaicplot(y, color=TRUE, main="", ylab="AIC= 0.0\n ", xlab=xlabel,  off=0)
    }
}


plot.single2 <- function(title, dname, res, item, aic.order, aic, tway.table, old.par, nw) {

    n <- length(title)
    i2 <- item[res]

    n1 <- n - 1
    if (n1 < 6) {
      nc <- n1
      nr <- 1
    } else if (n1 < 11) {
      nc <- as.integer(n / 2)
      nr <- 2
    } else {
      nc <- 5
      nr <- 2
    }
    m <- nc * nr

    new.mai <- old.par$mai
    new.mai[1] <- new.mai[3] * 0.8
    new.mai[3] <- new.mai[3] * 0.5
    newcex.main <- old.par$cex.main * 0.8  
    mtitle <- "Single Explanatory Models\nin ascending order of AIC"

    nplot <- 1
    j1 <- 1

    for (j in 1:(n-1)) {
      iex <- aic.order[j]
      ctbl <- tway.table[[iex]]$n

      if (nplot %% m == 1) {
        if (j > 1 || nw > 0) dev.new()
        par(mfcol=c(nc, nr), mai = new.mai, cex.main = newcex.main)
        nw <- nw + 1
      }

      j2 <- item[iex]
      y <- array(0, dim=c(i2, j2))
      for (ii in 1:i2) 
        for (k in j1:j2) y[ii, k] <- ctbl[k, ii]
      dimnames(y)[1] <- list(dimnames(ctbl)[[2]])
      dimnames(y)[2] <- list(dimnames(ctbl)[[1]])
      aaic <- round(aic[iex], 2)
      mosaicplot(y, color=TRUE, main=mtitle, ylab=paste("AIC=", aaic, "\n", title[iex]), xlab=paste(title[res]))
      mtitle <- ""
      nplot <- nplot + 1
    } 
}


plot.mosaic <- function(iplt, res, icl, k, mp, np, idm, title, etitle, dname, ctable, ctable.dname, nc, aic.order, aic, tway.table, total, nwindow, old.par) {

    i1 <- 1
    i2 <- nc[res]

    nex <- length(idm)
    jdm <- array(0, nex)
    for (i in 1:nex) jdm[i] <- idm[nex-i+1]
    idm <- c(jdm, i2)
    tt <- ctable[[k]]
    t <- array(tt, dim=idm)

    xtitle <- title[res]
    ytitle <- etitle[1]

    if (iplt == 1) {
      t <- aperm(t, c(nex:1, nex+1))
      for (i in 1:nex) dimnames(t)[i] <- list(ctable.dname[[i+1]])
      dimnames(t)[nex+1] <- list(ctable.dname[[1]])

      split <- "h"
      offset <- 10
      if (nex > 1)
        for (i in 2:nex) {
          ytitle <- c(ytitle, etitle[i])
          split <- c(split, "h")
          offset <- c(offset, 10)
        }
      split <- c(split, "v")
      offset <- c(offset, 0)

      if ((icl != 1 && k == 1)   || (icl == 1 && k == 2))
        plot.single1(title, dname, res, nc, aic.order, aic, tway.table, total, old.par, nwindow)

    } else if (iplt == 2) {
      t <- aperm(t, c(nex+1, nex:1))
      dimnames(t) <- ctable.dname

      if (nex > 1)
        for (i in 2:nex) {
          if (i %% 2 == 0) xtitle <- c(xtitle, etitle[i])
          if (i %% 2 == 1) ytitle <- c(ytitle, etitle[i])
        }

      split <- NULL
      offset <- NULL

      if ((icl != 1 && k == 1)   || (icl == 1 && k == 2))
        plot.single2(title, dname, res, nc, aic.order, aic, tway.table, old.par, nwindow)
    } # else if (iplt == 2) end

    if (icl == 1 && k == 1) {
      mtitle <- "Mosaicplot for additinal contingency table"
    } else {
      dev.new()
      mtitle <- "Minimum AIC Model of\n the Response Variable Distribution"
    }
    par(mfcol=c(1, 1))
    mosaicplot(t, color=TRUE, main=mtitle, xlab=paste(xtitle), ylab=paste(ytitle), dir=split,  off=offset)
    if (icl == 1 && k == 1)  nwindow <- nwindow + 1
    
    return(nwindow=nwindow)

}

plot.grshade <- function(title, aic, ires, nsamp, old.par) {

    new.mai <- old.par$mai
    new.mai[1] <- new.mai[1] * 1.5
    new.mai[2] <- new.mai[2] * 1.5
    new.mai[3] <- new.mai[3] * 2.0
    new.mai[4] <- new.mai[4] * 1.0
    par(mai=new.mai, cex.axis=0.8)

    response <- NULL
    explanatory <- NULL
    nres <- length(ires)
    for (i in nres:1) response <- c(response, title[ires[i]])
    nex <- length(title)
    for (i in 1:nex) explanatory <- c(explanatory, title[i])

    d1 <- dim(aic)[1]
    d2 <- dim(aic)[2]
    aics <- aic / nsamp
    aicsmin <- 0.0
    aicsmax <- 0.0
    for (i1 in 1:d1)
      for (i2 in 1:d2) 
        if (is.na(aics[i1,i2]) == FALSE) {
          aicsmin <- min(aicsmin, aics[i1,i2])
          aicsmax <- max(aicsmax, aics[i1,i2])
        }
    if (aicsmin  > -0.1)  aicsmin <- -0.15
    if (aicsmax < 0) aicsmax <- 0.5

    aics1 <- array(0, c(d1, d2))
    for (i in d1:1) aics1[i, ] <- aics[d1-i+1, ]
    aics1 <- aperm(aics1, c(2, 1))
    aics2 <- array(0, c(d2+1, d1))
    aics2[1:d2, ] <- aics1[1:d2, ]
    x <- (1:nrow(aics2))
    y <- (1:ncol(aics2))

    bp <- c(aicsmin, -0.1, -0.05, -0.01, 0, aicsmax) - 1.0e-10

    main.title <- paste("Gray Shading Display of All the AIC's\n\nGray level : N=Number of observations\n 5 \ : \ \ \ \ \ \ \ \ \ \ \ \ \ AIC/N < -0.10\ \n 4 \ : -0.10 =< AIC/N < -0.05\ \n 3 \ : -0.05 =< AIC/N < -0.01\ \n 2 \ : -0.01 =< AIC/N < \ 0.0\ \ \n 1 \ :\ \ 0.0\ \ \ =< AIC/N\ \ \ \ \ \ \ \ \ \ \ \ \n")

    image(x, y, aics2, col=gray(0:4/4), breaks=bp, axes = FALSE, xlab="", ylab="", main=main.title, 
sub="explanatory variable", cex.main=0.8, cex.sub=0.8)

    legend("topright", c("5","4","3","2","1"), pch=15,col=gray(0:4/4), bg="snow2", cex=0.8, pt.cex=1.1)

    axis(1, 1:d2, explanatory, las=3)
    axis(2, 1:d1, response, las=1)
    box()
}


#=====================================================================

Barplot2WayTable <- function (vname, resvar, exvar=NULL, tway.table, interval=NULL) {

#=====================================================================

    nv <- length(vname)
    if (is.null(exvar) == TRUE) exvar <- vname
    nresv <- length(resvar)
    nexv <- length(exvar)
    nc <- min(nexv, nv-1)

    if (length(tway.table) != nv) { #  catdap1
      ic <- 1
    } else if (length(tway.table[[1]]) == nv) {  # catdap1
      ic <- 1
    } else {  # catdap2
      ic <- 2
      if (is.null(interval) == TRUE) stop(" Error : Class interval is not specified.")
    }

    if (ic == 1) {  # catdap1() or catdap1c()
      nw <- 0
      for (j in 1:nresv) {
        res <- 0
        for (i in 1:nv) 
          if(vname[i] == resvar[j]) {
            res <- i
            break
          }
        if (res == 0) stop(" Error : response variable name is wrong.")

        nbar <- 0
        if (nexv > 1) par(mfcol=c(nc, 1))
        for (k in 1:nexv) {
          ex <- 0
          for (i in 1:nv)
            if(vname[i] == exvar[k]) {
              ex <- i
              break
            }
          if (ex == 0) stop(" Error : explanatory variable name is wrong.")

          if (res != ex) {
            if (nbar == 0 && nw > 0) {
              dev.new()
              if (nexv > 1) par(mfcol=c(nc, 1))
            }
            h <- tway.table[[j]][[ex]]$n
            barplot(h, names.arg=dimnames(h)[[2]], legend=dimnames(h)[[1]], args.legend=list(x="topright", title=vname[ex], cex=0.8), space=0, xlab=vname[res])
            nbar <- nbar + 1
          }
        }  # for (k in 1:nexv) end
        nw <- nw + 1
      }  # for (j in 1:nresv) end

    } else {  # catdap2()
      res <- 0
      for (i in 1:nv)
        if(vname[i] == resvar) {
          res <- i
          break
        }
      if (res == 0) stop(" Error : response variable name is wrong.")

      if (nexv > 1) par(mfcol=c(nc, 1))
      for (j in 1:nexv) {
        ex <- 0
        for (i in 1:nv)
          if(vname[i] == exvar[j]) {
            ex <- i
            break
          }
        if (ex == 0) stop(" Error : explanatory variable name is wrong.")

        if (res != ex) {
          h <- tway.table[[ex]]$n
          nr <- length(dimnames(h)[[1]])
          if (j == 1) {
            nrres <- length(dimnames(h)[[2]])
            resint <- dimnames(h)[[2]]
            resname <- tway.table[[ex]]$res
          }
          nint <- length(interval[[ex]])
          if (nr == nint) {
            barplot(h, names.arg=resint, legend=dimnames(h)[[1]], args.legend=list(x="topright", title=vname[ex], cex=0.8), space=0, xlab=vname[res])
          } else {
            itext <- NULL
            for (k in 1:nr) itext[[k]] <- paste(interval[[ex]][k], "-", interval[[ex]][k+1])
            barplot(h, names.arg=resint, legend=itext, args.legend=list(x="topright", title=vname[ex], cex=0.8), space=0, xlab=vname[res])
          }
        }
      }  # for (j in 1:nexv) end

      nint <- length(interval[[res]])
      if (nint > nrres) {
        cat("\n<Note>\n")
        cat(sprintf(" %s\n", resname))
        cat("\tcategory    \tvalue range\n")
        for (k in 1:nrres) cat(sprintf("\t%8i    \t%12.5e  -  %12.5e\n", k, interval[[res]][k], interval[[res]][k+1]))
      }
    }
    par(mfcol=c(1,1))
}

.noGenerics <- TRUE

options(warn.FPU=FALSE)
