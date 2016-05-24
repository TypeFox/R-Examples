netmeta <- function(TE, seTE,
                    treat1, treat2,
                    studlab, data=NULL, subset=NULL,
                    sm,
                    level=0.95, level.comb=0.95,
                    comb.fixed=TRUE, comb.random=FALSE,
                    reference.group="",
                    all.treatments=NULL,
                    seq=NULL,
                    tau.preset=NULL,
                    title="",
                    warn=TRUE
                    ){
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch TE, treat1, treat2, seTE, studlab from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$sm <- NULL
  mf$level <- mf$level.comb <- NULL
  mf$comb.fixed <- mf$comb.random <- mf$reference.group <- NULL
  mf$all.treatments <- mf$seq <- NULL
  mf$title <- mf$warn <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  ## Catch subset (possibly) from data:
  ##
  mf2 <- match.call()
  mf2$TE <- mf2$seTE <- mf2$treat1 <- mf2$treat2 <- NULL
  mf2$studlab <- NULL
  mf2$data <- mf2$sm <- NULL
  mf2$level <- mf2$level.comb <- NULL
  mf2$comb.fixed <- mf2$comb.random <- mf2$reference.group <- NULL
  mf2$all.treatments <- mf2$seq <- NULL
  mf2$title <- mf2$warn <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$TE))) ||
        (length(mf2$subset) > length(mf$TE)))
      stop("Length of subset is larger than number of comparisons.")
    else
      mf <- mf[mf2$subset,]
  
  
  TE     <- mf$TE
  seTE   <- mf$seTE
  ##
  treat1 <- mf$treat1
  treat2 <- mf$treat2
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  ##
  if (length(mf$studlab)!=0)
    studlab <- as.character(mf$studlab)
  else{
    if (warn)
      warning("No information given for argument 'studlab'. Assuming that comparisons are from independent studies.")
    studlab <- seq(along=TE)
  }
  
  
  if (any(treat1==treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
  
  
  labels <- sort(unique(c(treat1, treat2)))
  ##
  if (!is.null(seq))
    seq <- setseq(seq, labels)
  else{
    seq <- labels
    if (is.numeric(seq))
      seq <- as.character(seq)
  }
  
  
  ##
  ## Check for correct number of comparisons
  ##
  is.wholenumber <-
    function(x, tol=.Machine$double.eps^0.5)
      abs(x - round(x)) < tol
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8*tabnarms + 1))/2)
  ##
  if (sum(sel.narms)==1)
    stop(paste("Study '", names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               "\n  Please provide data for all treatment comparisons (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep=""))
  if (sum(sel.narms)>1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep=""),
                     collapse=", "),
               "\n  Please provide data for all treatment comparisons (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep=""))
  
  
  ##
  ## Check for levels of confidence interval
  ##
  meta:::chklevel(level)
  meta:::chklevel(level.comb)
  
  
  ##
  ## Check value for reference group
  ##
  if (is.null(all.treatments))
    if (reference.group=="")
      all.treatments <- TRUE
    else
      all.treatments <- FALSE
  ##
  if (reference.group !="")
    reference.group <- setref(reference.group, labels)
  
  
  ##
  ## Check number of subgraphs
  ##
  n.subnets <- netconnection(treat1, treat2, studlab)$n.subnets
  ##
  if (n.subnets > 1)
    stop(paste("Network consists of ", n.subnets, " separate sub-networks.\n  ",
               "Use R function 'netconnection' to identify sub-networks.",
               sep = ""))
  
  
  ##
  ## Check for correct treatment order within comparison
  ##
  wo <- treat1 > treat2
  ##
  if (any(wo)){
    if (warn)
      warning("Note, treatments within a comparison have been re-sorted in increasing order.", call.=FALSE)
    TE[wo] <- -TE[wo]
    ttreat1 <- treat1
    treat1[wo] <- treat2[wo]
    treat2[wo] <- ttreat1[wo]
  }
  
  
  ##
  ## Check value for argument 'sm'
  ##
  if (missing(sm))
    if (!is.null(data) && !is.null(attr(data, "sm")))
      sm <- attr(data, "sm")
    else
      sm <- ""
  
  
  ##
  ## Generate ordered data set, with added numbers of arms per study
  ##
  p0 <- prepare(TE, seTE, treat1, treat2, studlab)
  ##
  ## Study overview
  ##
  tdata <- data.frame(studies=p0$studlab, narms=p0$narms)
  tdata <- unique(tdata[order(tdata$studies, tdata$narms),])
  studies <- tdata$studies
  narms <- tdata$narms
  ##
  ## Network meta-analysis based on prepared data set
  ##
  ## Fixed effect model
  ##
  res.f <- nma.ruecker(p0$TE, sqrt(1/p0$weights),
                       p0$treat1, p0$treat2,
                       p0$treat1.pos, p0$treat2.pos,
                       p0$narms, p0$studlab,
                       sm, level, level.comb, p0$seTE)
  ##
  ## Random effects model
  ##
  if (is.null(tau.preset))
    tau <- res.f$tau
  else
    tau <- tau.preset
  ##
  p1 <- prepare(TE, seTE, treat1, treat2, studlab, tau)
  ##
  res.r <- nma.ruecker(p1$TE, sqrt(1/p1$weights),
                       p1$treat1, p1$treat2,
                       p1$treat1.pos, p1$treat2.pos,
                       p1$narms, p1$studlab, 
                       sm, level, level.comb, p1$seTE)
  
  
  o <- order(p0$order)
  ##
  res <- list(studlab=res.f$studlab[o],
              treat1=res.f$treat1[o],
              treat2=res.f$treat2[o],
              ##
              TE=res.f$TE[o],
              seTE=res.f$seTE.orig[o],
              seTE.adj=res.f$seTE[o],
              ##
              studies=studies,
              narms=narms,
              ##
              TE.nma.fixed=res.f$TE.nma[o],
              seTE.nma.fixed=res.f$seTE.nma[o],
              lower.nma.fixed=res.f$lower.nma[o],
              upper.nma.fixed=res.f$upper.nma[o],
              ##
              leverage.fixed=res.f$leverage[o],
              w.fixed=res.f$w.pooled[o],
              ##
              TE.fixed=res.f$TE.pooled,
              seTE.fixed=res.f$seTE.pooled,
              lower.fixed=res.f$lower.pooled,
              upper.fixed=res.f$upper.pooled,
              zval.fixed=res.f$zval.pooled,
              pval.fixed=res.f$pval.pooled,
              ##
              Q.fixed=res.f$Q.pooled[o],
              ##
              TE.nma.random=res.r$TE.nma[o],
              seTE.nma.random=res.r$seTE.nma[o],
              lower.nma.random=res.r$lower.nma[o],
              upper.nma.random=res.r$upper.nma[o],
              ##
              w.random=res.r$w.pooled[o],
              ##
              TE.random=res.r$TE.pooled,
              seTE.random=res.r$seTE.pooled,
              lower.random=res.r$lower.pooled,
              upper.random=res.r$upper.pooled,
              zval.random=res.r$zval.pooled,
              pval.random=res.r$pval.pooled,
              ##
              TE.direct.fixed=res.f$TE.direct,
              seTE.direct.fixed=res.f$seTE.direct,
              lower.direct.fixed=res.f$lower.direct,
              upper.direct.fixed=res.f$upper.direct,
              zval.direct.fixed=res.f$zval.direct,
              pval.direct.fixed=res.f$pval.direct,
              ##
              TE.direct.random=res.r$TE.direct,
              seTE.direct.random=res.r$seTE.direct,
              lower.direct.random=res.r$lower.direct,
              upper.direct.random=res.r$upper.direct,
              zval.direct.random=res.r$zval.direct,
              pval.direct.random=res.r$pval.direct,
              ##
              treat1.pos=res.f$treat1.pos[o],
              treat2.pos=res.f$treat2.pos[o],
              ##
              k=res.f$k,
              m=res.f$m,
              n=res.f$n,
              Q=res.f$Q,
              df=res.f$df,
              pval.Q=res.f$pval.Q,
              I2=res.f$I2,
              tau=res.f$tau,
              tau.preset=tau.preset,                                             
              Q.heterogeneity=res.f$Q.heterogeneity,
              Q.inconsistency=res.f$Q.inconsistency,
              ##
              sm=sm,
              level=level,
              level.comb=level.comb,
              comb.fixed=comb.fixed,
              comb.random=comb.random,
              ##
              A.matrix=res.f$A.matrix,
              L.matrix=res.f$L.matrix,
              Lplus.matrix=res.f$Lplus.matrix,
              Q.matrix=res.f$Q.matrix,
              ##
              G.matrix=res.f$G.matrix[o,o],
              H.matrix=res.f$H.matrix[o,o],
              ##
              Cov.fixed = res.f$Cov,
              Cov.random = res.r$Cov,
              ##
              Q.decomp=res.f$Q.decomp,
              ##
              reference.group=reference.group,
              all.treatments=all.treatments,
              ##
              seq=seq,
              ##
              title=title,
              ##
              warn=warn,
              call=match.call(),
              version=packageDescription("netmeta")$Version
              )
  
  class(res) <- "netmeta"
  
  res
}
