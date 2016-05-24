# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












 A1A2A3 <- function(link = "logit",
                    inbreeding = FALSE,  # HWE assumption is the default
                    ip1 = NULL, ip2 = NULL, iF = NULL) {





  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.logical(inbreeding) || length(inbreeding) > 1)
    stop("argument 'inbreeding' must be a single logical")


  new("vglmff",
  blurb = c("G1-G2-G3 phenotype (",
            ifelse(inbreeding, "without", "with"),
            " the Hardy-Weinberg equilibrium assumption)\n\n",
            "Links:    ",
            namesof("p1", link, earg = earg, tag = FALSE), ", ", 
            namesof("p2", link, earg = earg, tag = FALSE),
            if (inbreeding) paste(",",
            namesof("f",  link, earg = earg, tag = FALSE)) else
            ""),
  deviance = Deviance.categorical.data.vgam,
  infos = eval(substitute(function(...) {
    list(Q1 = 6,
         M1 = ifelse( .inbreeding , 3, 2),
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("p1", "p2",
                if ( .inbreeding ) "f" else NULL),
         link = if ( .inbreeding )
                  c("p1" = .link , "p2" = .link , "f" = .link ) else
                  c("p1" = .link , "p2" = .link  ))
  }, list( .link = link, .inbreeding = inbreeding ))),

  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.VGAM)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("A1A1", "A1A2", "A1A3", "A2A2", "A2A3", "A3A3")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
       setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
            stop("the columns of the response matrix should have ",
                 "names (output of colnames()) ordered as ",
                 "c('A1A1','A1A2','A1A3','A2A2','A2A3','A3A3')")
    }

    predictors.names <-
     c(namesof("p1", .link , earg = .earg , tag = FALSE),
       namesof("p2", .link , earg = .earg , tag = FALSE),
       if ( .inbreeding )
       namesof("f",  .link , earg = .earg , tag = FALSE) else NULL)
    mustart <- (y + mustart) / 2


    if (is.null(etastart)) {



      mydeterminant <- weighted.mean(
                       mustart[, 2] * mustart[, 3] +
                       mustart[, 2] * mustart[, 5] +
                       mustart[, 3] * mustart[, 5], w)
      p1 <- if (is.numeric( .ip1 )) rep( .ip1 , len = n) else
            weighted.mean(mustart[, 2] * mustart[, 3], w) / mydeterminant
      p2 <- if (is.numeric( .ip2 )) rep( .ip2 , len = n) else
            weighted.mean(mustart[, 2] * mustart[, 5], w) / mydeterminant
      ff <- if (is.numeric( .iF  )) rep( .iF  , len = n) else
            weighted.mean(abs(1 - mustart[, 2] / (2 * p1 * p2)), w)
      p1 <- rep(p1, len = n)
      p2 <- rep(p2, len = n)
      ff <- rep(ff, len = n)
      p1[p1 < 0.05] <- 0.05
      p1[p1 > 0.99] <- 0.99
      p2[p2 < 0.05] <- 0.05
      p2[p2 > 0.99] <- 0.99
      ff[ff < 0.05] <- 0.05
      ff[ff > 0.99] <- 0.99

      etastart <-
        cbind(theta2eta(p1, .link , earg = .earg ),
              theta2eta(p2, .link , earg = .earg ),
              if ( .inbreeding )
              theta2eta(ff, .link , earg = .earg ) else NULL)
      mustart <- NULL  # Since etastart has been computed.

    }
  }), list( .link = link, .ip1 = ip1, .ip2 = ip2, .iF = iF,
            .inbreeding = inbreeding,
            .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    p1 <- eta2theta(eta[, 1], link = .link , earg = .earg )
    p2 <- eta2theta(eta[, 2], link = .link , earg = .earg )
    f  <- if ( .inbreeding )
          eta2theta(eta[, 3], link = .link , earg = .earg ) else 0
    p3 <- abs(1 - p1 - p2)
      cbind("A1A1" = f*p1+(1-f)*p1^2,
            "A1A2" = 2*p1*p2*(1-f),
            "A1A3" = 2*p1*p3*(1-f),
            "A2A2" = f*p2+(1-f)*p2^2,
            "A2A3" = 2*p2*p3*(1-f),
            "A3A3" = f*p3+(1-f)*p3^2)
  }, list( .link = link, .earg = earg, .inbreeding = inbreeding))),

  last = eval(substitute(expression({
    if ( .inbreeding ) {
      misc$link <-    c(p1 = .link , p2 = .link , f = .link )
      misc$earg <- list(p1 = .earg , p2 = .earg , f = .earg )
    } else {
      misc$link <-    c(p1 = .link , p2 = .link )
      misc$earg <- list(p1 = .earg , p2 = .earg )
    }

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg, .inbreeding = inbreeding ))),

  loglikelihood = function(mu, y, w, residuals = FALSE,
                           eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
      },
  vfamily = c("A1A2A3", "vgenetic"),
  deriv = eval(substitute(expression({
    p1 <- eta2theta(eta[, 1], link = .link , earg = .earg )
    p2 <- eta2theta(eta[, 2], link = .link , earg = .earg )
    p3 <- 1-p1-p2
    f  <- if ( .inbreeding )
          eta2theta(eta[, 3], link = .link , earg = .earg ) else 0
    if ( .inbreeding ) {
      dP1 <- cbind(f + 2*p1*(1-f), 2*(1-f)*p2, 2*(1-f)*(1-p2-2*p1),
                  0, -2*(1-f)*p2, -f - 2*p3*(1-f))
      dP2 <- cbind(0, 2*p1*(1-f), -2*(1-f)*p1, f+2*p2*(1-f),
                   2*(1-f)*(1-p1-2*p2), -f - 2*p3*(1-f))
      dP3 <- cbind(p1*(1-p1), -2*p1*p2, -2*p1*p3, p2*(1-p2), -2*p2*p3, 
                   p3*(1-p3))
      dl1 <- rowSums(y * dP1 / mu)
      dl2 <- rowSums(y * dP2 / mu)
      dl3 <- rowSums(y * dP3 / mu)
      dPP.deta <- dtheta.deta(cbind(p1, p2, f),
                              link = .link , earg = .earg )
      c(w) * cbind(dPP.deta[, 1] * dl1,
                   dPP.deta[, 2] * dl2, 
                   dPP.deta[, 3] * dl3)
    } else {
      dl.dp1 <- (2*y[, 1]+y[, 2]+y[, 4])/p1 -
                (2*y[,6]+y[, 4]+y[,5])/(1-p1-p2)
      dl.dp2 <- (2*y[, 3]+y[, 2]+y[,5])/p2 -
                (2*y[,6]+y[, 4]+y[,5])/(1-p1-p2)

      dp1.deta <- dtheta.deta(p1, link = .link , earg = .earg )
      dp2.deta <- dtheta.deta(p2, link = .link , earg = .earg )

      c(w) * cbind(dl.dp1 * dp1.deta,
                   dl.dp2 * dp2.deta)
    }
  }), list( .link = link, .earg = earg, .inbreeding = inbreeding ))),
  weight = eval(substitute(expression({
    if ( .inbreeding ) {
      dPP <- array(c(dP1, dP2, dP3), c(n, 6, 3))
      wz <- matrix(NA_real_, n, dimm(M))  # dimm(M)==6 because M==3
      for (i1 in 1:M)
        for (i2 in i1:M) {
          index <- iam(i1, i2, M)
          wz[, index] <- rowSums(dPP[, , i1, drop = TRUE] *
                                 dPP[, , i2, drop = TRUE] / mu) *
                                 dPP.deta[, i1] * dPP.deta[, i2]
      }
    } else {
      qq <- 1-p1-p2
      wz <- matrix(NA_real_, n, dimm(M))  # dimm(M)==3 because M==2
      ned2l.dp12  <-  2 * (1/p1 + 1/qq)
      ned2l.dp22  <-  2 * (1/p2 + 1/qq)
      ned2l.dp1dp2 <-  2 / qq
      wz[, iam(1, 1, M)] <- ned2l.dp12 * dp1.deta^2
      wz[, iam(2, 2, M)] <- ned2l.dp22 * dp2.deta^2
      wz[, iam(1, 2, M)] <- ned2l.dp1dp2 * dp1.deta * dp2.deta
    }
    c(w) * wz
  }), list( .link = link, .earg = earg, .inbreeding = inbreeding ))))
}








 MNSs <- function(link = "logit",
                  imS = NULL, ims = NULL, inS = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("MNSs Blood Group System (MS-Ms-MNS-MNs-NS-Ns phenotype)\n\n",
            "Links:    ",
            namesof("mS", link, earg = earg), ", ", 
            namesof("ms", link, earg = earg), ", ", 
            namesof("nS", link, earg = earg, tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.VGAM)

    if (length(mustart.orig))
      mustart <- mustart.orig

        ok.col.ny <- c("MS","Ms","MNS","MNs","NS","Ns")
        if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
            setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
            stop("the columns of the response matrix should have ",
                 "names (output of colnames()) ordered as ",
                 "c('MS','Ms','MNS','MNs','NS','Ns')")
    }

    predictors.names <-
       c(namesof("mS", .link , earg = .earg , tag = FALSE),
         namesof("ms", .link , earg = .earg , tag = FALSE),
         namesof("nS", .link , earg = .earg , tag = FALSE))

    if (is.null(etastart)) {
      ms <- if (is.numeric(.ims)) rep(.ims, n) else
                 c(sqrt(mustart[, 2]))
      ns <- c(sqrt(mustart[,6]))
      nS <- if (is.numeric(.inS)) rep(.inS, n) else
          c(-ns + sqrt(ns^2 + mustart[,5]))  # Solve a quadratic eqn
      mS <- if (is.numeric(.imS)) rep(.imS, n) else
              1-ns-ms-nS
      etastart <- cbind(theta2eta(mS, .link , earg = .earg ),
                       theta2eta(ms, .link , earg = .earg ),
                       theta2eta(nS, .link , earg = .earg ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .imS = imS, .ims = ims, .inS = inS, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mS <- eta2theta(eta[, 1], link = .link , earg = .earg )
    ms <- eta2theta(eta[, 2], link = .link , earg = .earg )
    nS <- eta2theta(eta[, 3], link = .link , earg = .earg )
    ns <- abs(1 - mS - ms - nS)
    cbind(MS  = mS^2 + 2*mS*ms,
          Ms  = ms^2,
          MNS = 2*(mS*nS + ms*nS + mS*ns),
          MNs = 2*ms*ns,
          NS  = nS^2 + 2*nS*ns,
          Ns  = ns^2)
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(mS = .link , ms = .link , nS = .link )

    misc$earg <- list(mS = .earg , ms = .earg , nS = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("MNSs", "vgenetic"),
  deriv = eval(substitute(expression({
    mS <- eta2theta(eta[, 1], link = .link , earg = .earg )
    ms <- eta2theta(eta[, 2], link = .link , earg = .earg )
    nS <- eta2theta(eta[, 3], link = .link , earg = .earg )
    ns <- 1-mS-ms-nS
    dP1 <- cbind(2*(mS+ms), 0, 2*(nS+ns-mS), -2*ms, -2*nS, -2*ns)
    dP2 <- cbind(2*mS, 2*ms, 2*(nS-mS), 2*(ns-ms), -2*nS, -2*ns)
    dP3 <- cbind(0, 0, 2*ms, -2*ms,  2*ns, -2*ns)  # n x 6
    dl1 <- rowSums(y * dP1 / mu)
    dl2 <- rowSums(y * dP2 / mu)
    dl3 <- rowSums(y * dP3 / mu)
    dPP.deta <- dtheta.deta(cbind(mS, ms, nS), link = .link , earg = .earg )
    c(w) * dPP.deta * cbind(dl1, dl2, dl3)
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    dPP <- array(c(dP1,dP2,dP3), c(n,6, 3))
    wz <- matrix(NA_real_, n, dimm(M))  # dimm(M)==6 because M==3
    for (i1 in 1:M)
      for (i2 in i1:M) {
        index <- iam(i1,i2, M)
        wz[,index] <- rowSums(dPP[,,i1,drop = TRUE] *
                              dPP[,,i2,drop = TRUE] / mu) *
                              dPP.deta[,i1] * dPP.deta[,i2]
    }
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}






 ABO <- function(link.pA = "logit", link.pB = "logit",
                 ipA = NULL, ipB = NULL, ipO = NULL,
                 zero = NULL) {
  link.pA <- as.list(substitute(link.pA))
  earg.pA <- link2list(link.pA)
  link.pA <- attr(earg.pA, "function.name")

  link.pB <- as.list(substitute(link.pB))
  earg.pB <- link2list(link.pB)
  link.pB <- attr(earg.pB, "function.name")


  new("vglmff",
  blurb = c("ABO Blood Group System (A-B-AB-O phenotype)\n\n",
            "Links:    ",
            namesof("pA", link.pA, earg = earg.pA, tag = FALSE), ", ", 
            namesof("pB", link.pB, earg = earg.pB, tag = FALSE)),
  deviance = Deviance.categorical.data.vgam,

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 4,
         multipleResponses = FALSE,
         parameters.names = c("pA", "pB"),
         expected = TRUE,
         zero = .zero ,
         link = c("pA" = .link.pA , "pB" = .link.pB ),
         earg = c("pA" = .earg.pB , "pB" = .earg.pB )
        )
  }, list( .link.pA = link.pA, .link.pB = link.pB,
           .earg.pA = earg.pA, .earg.pB = earg.pB,
           .zero = zero ))),

  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.VGAM)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("A","B","AB","O")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
        setequal(ok.col.ny, col.ny)) {
      if (!all(ok.col.ny == col.ny))
        stop("the columns of the response matrix should have names ",
             "(output of colnames()) ordered as c('A','B','AB','O')")
    }


    predictors.names <-
      c(namesof("pA", .link.pA , earg = .earg.pA , tag = FALSE),
        namesof("pB", .link.pB , earg = .earg.pB , tag = FALSE))
    mustart <- (y + mustart) / 2

    if (!length(etastart)) {
      pO <- if (is.Numeric( .ipO )) rep( .ipO , len = n) else
        rep(c(sqrt( weighted.mean(mustart[, 4], w)) ), len = n)
      pA <- if (is.Numeric( .ipA )) rep( .ipA , len = n) else
        rep(c(1 - sqrt(weighted.mean(mustart[, 2] + mustart[, 4], w))),
            len = n)
      pB <- if (is.Numeric( .ipB )) rep( .ipB , len = n) else
            abs(1 - pA - pO)
      etastart <- cbind(theta2eta(pA, .link.pA , earg = .earg.pA ),
                        theta2eta(pB, .link.pB , earg = .earg.pB ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link.pA = link.pA, .link.pB = link.pB,
            .ipO = ipO, .ipA = ipA, .ipB = ipB,
            .earg.pA = earg.pA, .earg.pB = earg.pB ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
      pA <- eta2theta(eta[, 1], link = .link.pA , earg = .earg.pA )
      pB <- eta2theta(eta[, 2], link = .link.pB , earg = .earg.pB )
      pO <- abs(1 - pA - pB)
      cbind(A  = pA*(pA+2*pO),
            B  = pB*(pB+2*pO),
            AB = 2*pA*pB,
            O  = pO*pO) 
  }, list( .link.pA = link.pA, .link.pB = link.pB,
           .earg.pA = earg.pA, .earg.pB = earg.pB ))),

  last = eval(substitute(expression({
    misc$link <-    c(pA = .link.pA , pB = .link.pB )
    misc$earg <- list(pA = .earg.pA , pB = .earg.pB )
    misc$expected <- TRUE
  }), list( .link.pA = link.pA, .link.pB = link.pB,
            .earg.pA = earg.pA, .earg.pB = earg.pB ))),


  loglikelihood =
  function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
        sum(dmultinomial(x = w * y, size = w, prob = mu, log = TRUE,
                         dochecking = FALSE))
    },

  vfamily = c("ABO", "vgenetic"),

  deriv = eval(substitute(expression({
    ppp <- eta2theta(eta[, 1], link = .link.pA , earg = .earg.pA )
    qqq <- eta2theta(eta[, 2], link = .link.pB , earg = .earg.pB )
    rrr <- abs(1 - ppp - qqq)


    pbar <- 2*rrr + ppp
    qbar <- 2*rrr + qqq
    naa <- y[, 1]
    nbb <- y[, 2]
    nab <- y[, 3]
    noo <- y[, 4]

    dl.dp <- (naa+nab)/ppp -   naa/pbar - 2*nbb/qbar - 2*noo/rrr
    dl.dq <- (nbb+nab)/qqq - 2*naa/pbar -   nbb/qbar - 2*noo/rrr
    dp.deta <- dtheta.deta(ppp, link = .link.pA , earg = .earg.pA )
    dq.deta <- dtheta.deta(qqq, link = .link.pB , earg = .earg.pB )

    c(w) * cbind(dl.dp * dp.deta,
                 dl.dq * dq.deta)
  }), list( .link.pA = link.pA, .link.pB = link.pB,
            .earg.pA = earg.pA, .earg.pB = earg.pB ))),

  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, dimm(M))  # dimm(M)==3 because M==2

    ned2l.dp2  <- (1 + 2/ppp + 4*qqq/qbar + ppp/pbar)
    ned2l.dq2  <- (1 + 2/qqq + 4*ppp/pbar + qqq/qbar)
    ned2l.dpdq <- 2 * (1 + qqq/qbar + ppp/pbar)

    wz[, iam(1, 1, M)] <- ned2l.dp2  * dp.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dq2  * dq.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dpdq * dp.deta * dq.deta
    c(w) * wz
  }), list( .link.pA = link.pA, .link.pB = link.pB,
            .earg.pA = earg.pA, .earg.pB = earg.pB ))))
}




 AB.Ab.aB.ab <- function(link = "logit", init.p = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("AB-Ab-aB-ab phenotype\n\n",
            "Links:    ", namesof("p", link, earg = earg, tag = TRUE)),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.VGAM)

    if (length(mustart.orig)) {
      mustart <- mustart.orig
    }

    ok.col.ny <- c("AB","Ab","aB","ab")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
       setequal(ok.col.ny, col.ny)) {
        if (!all(ok.col.ny == col.ny))
          stop("the columns of the response matrix should have ",
               "names (output of colnames()) ordered as ",
               "c('AB','Ab','aB','ab')")
    }

    predictors.names <- namesof("p", .link , earg = .earg , tag = FALSE)
    mustart <- (y + mustart) / 2

    if (is.null(etastart)) {
      p.init <- if (is.numeric( .init.p )) rep( .init.p , len = n) else
        rep(c(sqrt(4 * weighted.mean(mustart[, 4], w))), len = n)

      etastart <- cbind(theta2eta(p.init, .link , earg = .earg ))
      etastart <- jitter(etastart)
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .link = link, .init.p = init.p, .earg = earg))),
  linkinv = eval(substitute(function(eta,extra = NULL) {
    p <- eta2theta(eta, link = .link , earg = .earg )
    pp4 <- p * p / 4
    cbind(AB = 0.5 + pp4,
          Ab = 0.25 - pp4,
          aB = 0.25 - pp4,
          ab = pp4) 
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    misc$link <-    c(p = .link )

    misc$earg <- list(p = .earg )

    misc$expected <- TRUE
   }), list( .link = link, .earg = earg))),


  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
      sum(dmultinomial(x = w * y, size = w, prob = mu,
                       log = TRUE, dochecking = FALSE))
    },
  vfamily = c("AB.Ab.aB.ab", "vgenetic"),
  deriv = eval(substitute(expression({
    pp <- eta2theta(eta, link = .link , earg = .earg )

    p2 <- pp*pp
    nAB <- w * y[, 1]
    nAb <- w * y[, 2]
    naB <- w * y[, 3]
    nab <- w * y[, 4]

    dl.dp <- 8 * pp * (nAB/(2+p2) - (nAb+naB)/(1-p2) + nab/p2)

    dp.deta <- dtheta.deta(pp, link = .link , earg = .earg )

    dl.dp * dp.deta
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    ned2l.dp2 <- 4 * p2 * (1/(2+p2) + 2/(1-p2) + 1/p2)
    wz <- cbind((dp.deta^2) * ned2l.dp2)
    c(w) * wz
  }), list( .link = link, .earg = earg))))
}













 AA.Aa.aa <-
  function(linkp = "logit",
           linkf = "logit",
           inbreeding = FALSE,  # HWE assumption is the default
           ipA = NULL,
           ifp = NULL,
           zero = NULL) {
    
  linkp <- as.list(substitute(linkp))
  eargp <- link2list(linkp)
  linkp <- attr(eargp, "function.name")

  linkf <- as.list(substitute(linkf))
  eargf <- link2list(linkf)
  linkf <- attr(eargf, "function.name")

  if (!is.logical(inbreeding) || length(inbreeding) > 1)
    stop("argument 'inbreeding' must be a single logical")

  

  new("vglmff",
  blurb = c("AA-Aa-aa phenotype (",
            ifelse(inbreeding, "without", "with"),
            " the Hardy-Weinberg equilibrium assumption)\n\n",
            "Links:    ",
            namesof("pA", linkp, earg = eargp, tag = FALSE),
            if (inbreeding) paste(",",
            namesof("f",  linkf, earg = eargf, tag = FALSE)) else
            ""),
  deviance = Deviance.categorical.data.vgam,
  infos = eval(substitute(function(...) {
    list(M1 = ifelse( .inbreeding , 2, 1),
         Q1 = 3,
         multipleResponses = FALSE,
         parameters.names = c("pA",
                if ( .inbreeding ) "f" else NULL),
         expected = TRUE,
         zero = .zero ,
         link = if ( .inbreeding ) c("pA" = .linkp , "f" = .linkf ) else
                                   c("pA" = .linkp ))
  }, list( .linkp = linkp,
           .linkf = linkf, .inbreeding = inbreeding,
           .zero = zero ))),

  initialize = eval(substitute(expression({
    mustart.orig <- mustart

    delete.zero.colns <- FALSE
    eval(process.categorical.data.VGAM)

    if (length(mustart.orig))
      mustart <- mustart.orig

    ok.col.ny <- c("AA", "Aa", "aa")
    if (length(col.ny <- colnames(y)) == length(ok.col.ny) &&
      setequal(ok.col.ny, col.ny)) {
      if (!all(ok.col.ny == col.ny))
        stop("the columns of the response matrix ",
             "should have names ",
             "(output of colnames()) ordered as c('AA','Aa','aa')")
    }

    predictors.names <-
      c(namesof("pA", .linkp , earg = .eargp , tag = FALSE),
        if ( .inbreeding )
        namesof("f",  .linkf , earg = .eargf , tag = FALSE) else NULL)
    mustart <- (y + mustart) / 2
        

    if (is.null(etastart)) {
      pA <- if (is.numeric( .ipA )) rep( .ipA , len = n) else
              rep(c(sqrt( weighted.mean(mustart[, 1], w))), len = n)
      fp <- if (is.numeric( .ifp )) rep( .ifp , len = n) else
              runif(n)  # 1- mustart[, 2]/(2*pA*(1-pA))
      etastart <- cbind(theta2eta(pA, .linkp , earg = .eargp ),
                        if ( .inbreeding )
                        theta2eta(fp, .linkf , earg = .eargf ) else NULL)
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .linkp = linkp, .linkf = linkf,
            .ipA = ipA, .ifp = ifp, .inbreeding = inbreeding,
            .eargp = eargp, .eargf = eargf ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta <- as.matrix(eta)
    pA <- eta2theta(eta[, 1], link = .linkp , earg = .eargp )
    fp <- if ( .inbreeding )
          eta2theta(eta[, 2], link = .linkf , earg = .eargf ) else 0

    cbind(AA = pA^2 + pA * (1-pA) * fp,
          Aa = 2 * pA * (1-pA) * (1 - fp),
          aa = (1-pA)^2 + pA * (1-pA) * fp)
  }, list( .linkp = linkp, .linkf = linkf,
           .eargp = eargp, .eargf = eargf,
           .inbreeding = inbreeding ))),

  last = eval(substitute(expression({
    if ( .inbreeding ) {
      misc$link <-    c("pA" = .linkp, "f" = .linkf )
      misc$earg <- list("pA" = .eargp, "f" = .eargf )
    } else {
      misc$link <-    c("pA" = .linkp )
      misc$earg <- list("pA" = .eargp )
    }
    misc$expected <- TRUE
  }), list( .linkp = linkp, .linkf = linkf,
            .eargp = eargp, .eargf = eargf,
            .inbreeding = inbreeding ))),



  loglikelihood = function(mu, y, w, residuals = FALSE, eta, extra = NULL)
    if (residuals)
      stop("loglikelihood residuals not implemented yet") else {
        sum(dmultinomial(x = w * y, size = w, prob = mu,
                         log = TRUE, dochecking = FALSE))
    },
  vfamily = c("AA.Aa.aa", "vgenetic"),
  deriv = eval(substitute(expression({
    eta <- as.matrix(eta)
    pA <- eta2theta(eta[, 1], link = .linkp , earg = .eargp )
    fp <- if ( .inbreeding )
          eta2theta(eta[, 2], link = .linkf , earg = .eargf ) else 0

    if ( .inbreeding ) {
      dP1 <- cbind(2*pA*(1-fp) + fp,
                   2*(1-fp)*(1-2*pA),
                  -2*(1-pA) + fp*(1-2*pA))
      dP2 <- cbind(pA*(1-pA),
                  -2*pA*(1-pA),
                   pA*(1-pA))
      dl1 <- rowSums(y * dP1 / mu)
      dl2 <- rowSums(y * dP2 / mu)

      dPP.deta <- dtheta.deta(pA, link = .linkp , earg = .eargp )
      dfp.deta <- dtheta.deta(fp, link = .linkf , earg = .eargf )

      c(w) * cbind(dPP.deta * dl1,
                   dfp.deta * dl2)      
    } else {
      nAA <- c(w) * y[, 1]
      nAa <- c(w) * y[, 2]
      naa <- c(w) * y[, 3]
      dl.dpA <- (2*nAA+nAa)/pA - (nAa+2*naa)/(1-pA)
      dpA.deta <- dtheta.deta(pA, link = .linkp , earg = .eargp )
      dl.dpA * dpA.deta
    }  
  }), list( .linkp = linkp, .linkf = linkf,
            .eargp = eargp, .eargf = eargf,
            .inbreeding = inbreeding ))),
  weight = eval(substitute(expression({
    if ( .inbreeding ) {
      dPP <- array(c(dP1, dP2), c(n, 3, 2))
      dPP.deta <- cbind(dtheta.deta(pA, link = .linkp , earg = .eargp ),
                        dtheta.deta(fp, link = .linkf , earg = .eargf ))
      wz <- matrix(NA_real_, n, dimm(M))  # dimm(M)==3 because M==2
      for (i1 in 1:M)
        for (i2 in i1:M) {
          index <- iam(i1, i2, M)
          wz[, index] <- rowSums(dPP[,, i1, drop = TRUE] *
                                 dPP[,, i2, drop = TRUE] / mu) *
                                 dPP.deta[, i1] * dPP.deta[, i2]
        }
      c(w) * wz
    } else {
      ned2l.dp2 <- 2 / (pA * (1-pA))
      wz <- cbind(c(w) * ned2l.dp2 * dpA.deta^2)
      wz
    }
  }), list( .linkp = linkp, .linkf = linkf,
            .eargp = eargp, .eargf = eargf,
            .inbreeding = inbreeding ))))
}






