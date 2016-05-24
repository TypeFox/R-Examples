`prentice.test` <-
function(
  y,                                           # data
  groups,                                      # groups
  blocks  = NULL,                              # blocks
  score   = "rank",                            # NULL: y already scored
  blkwght = "prentice",                        # block weights
  condvar = TRUE,                              # for compatibility
  alternative = "two.sided",                   # wilcox only
  mu      = 0,                                 # wilcox only
  paired  = FALSE,                             # wilcox only
  exact   = NULL,                              # wilcox only
  correct = FALSE,                             # wilcox only
  df      = -1,            # >0: as entered, =0: by design, -1: by data
  warn    = 0,
  optim     = TRUE)        # special case (kw/fm) optimization
{

  "%&%" <- function(x,y) paste(x,y,sep="") # makes shorter code
  stop.on.bdObject <- function(x, name = deparse(substitute(x)))
  {
	# simple version that works even if bdObject subclasses not defined
	if(is.element(as.character(class(x)), c("bdCharacter", "bdExclude", "bdFactor", "bdFrame", "bdInternalCache",
		"bdLogical", "bdModel", "bdNumeric", "bdObject", "bdOmit", "bdPackedObject", "bdSeries", 
		"bdSignalSeries", "bdTimeDate", "bdTimeSeries", "bdTimeSpan", "bdVector"))) {
		expr <- substitute(stop("Argument ", name, " cannot be a bdObject"), list(name = name))
		eval(expr)
	}
	invisible(NULL)
  }

  if (!is.R()) {
    if (is(y,     "bdObject")) stop.on.bdObject(y)
    if (is(groups,"bdObject")) stop.on.bdObject(groups)
    if (is(blocks,"bdObject")) stop.on.bdObject(blocks)
  }

  Norm01 <- function(x) {
    lenx <- length(x)
    if (lenx == 1) stop("can't normalize an integer")
    else return(ifelse1(all(abs(x-0.5))<0.5,x,x/(lenx+1)))
  }

  is.rank <- function(FUN) all(FUN(c(2,3,5))==c(1,2,3))
          # this condition is not sufficient, but should suffice.

  xTable <- function(blk,grp)
  ifelse1(length(blk)==0,table(0,0),
  cbind(ifelse1(min(grp)==0,NULL,0),      # guarantee grp==0 col for NA
    table(rep(blk,length.out=length(grp)),grp)))

# ----------------------------------------------------- *begin* program

  wRank <- mu.rank.nna
  Score <- NoOp

  if (is.null(score)) {                           # data already scored
    if (!condvar)
    stop("unconditional var requires score function to be known")
    wRank <- NoOp
    Score <- NoOp
  } else if (is.function(score)) {
    if (is.rank(score)) {
      Score <- NoOp
    } else {
      if (warn>=0) warning("external score function")
    }
  } else if (is.character(score)) {
    if (tolower(score) %in% c("rank","wilcoxon")) {
      Score <- NoOp
    } else if (tolower(score) %in% c("normal","vanderwaerden")) {
        Score <- MC(function(x) qnorm(Norm01(x)), list(Norm01=Norm01))
    } else {
        stop("illegal score parameter ", score)
    }
  } else stop("illegal 'score' parameter")

  if(!is.finite(mu) ) stop("argument 'mu' must be finite.")
  if(length(mu) != 1) stop("argument 'mu' must be a single value.")

  alternative <-
    char.expand(alternative,
      c("two.sided", "greater", "less"),
      stop("argument 'alternative' must match one of",
        " \"greater\", \"less\", \"two.sided\"."
          ))

  y.name      <- deparse(substitute(y))
  groups.name <- deparse(substitute(groups))

  if(substring(y.name,1,1)=="\140") {
    y.name <- substring(y.name,2,nchar(y.name)-1) %&% " of sizes "
    groups.name <- substring(groups.name,2,nchar(groups.name)-1)
  } else {
    y.name      <- "'" %&% deparse(substitute(y))      %&% "' (y)  by "
    groups.name <- "'" %&% deparse(substitute(groups)) %&% "' (groups)"
  }
  blocks.name <- "'" %&% deparse(substitute(blocks)) %&% "' (blocks)"

# ------------------------------------------------- begin handling data

  if(length(y) != length(c(groups)))
    stop(y.name, " and ", groups.name, " must have the same length.")
  groups <- unclass(groups)  # grps numbered in the order of appearance
  if (!is.null(levels(groups))) groups <- as.numeric(groups)
    if (NAgrps <- anyMissing(groups) && (warn >= 0)) warning(
      sum(is.na(groups)), " observation(s) with NA in ", groups.name)
  groups <- match(groups,unique(c(NA,groups)))-1  # NA's are group 0

  if (KW <- is.null(blocks)) {
    blocks <- 1
  } else {                   # blks numbered in the order of appearance
    if (length(y) != length(c(blocks)))
    stop(y.name, " and ", blocks.name, " must have the same length.")
    blocks <- unclass(blocks)
    if (!is.null(levels(blocks))) blocks <- as.numeric(blocks)

    if (chkMissing(blocks) && (!all(blocks.ok <- !is.na(blocks)))) {
    y      <- y     [blocks.ok]
    groups <- groups[blocks.ok]
    blocks <- blocks[blocks.ok]
    if (warn >= 0) warning( sum(!blocks.ok),
      " observation(s) with NA/NaN in ", blocks.name, " removed.")
    }
    blocks <- match(blocks,unique(blocks))
    KW <- (max(blocks)==1)  # KW by elim of all-NA blocks
  }

  m  <- xTable(blocks,groups)
  n  <- dim(m)[1]
  p  <- dim(m)[2]-1                          # planned number of groups

  if (!chkMissing(y) || (all(y.ok <- is.orderable(y)))) {
    M <- m
  } else {  
      y      <- y     [y.ok]
      groups <- groups[y.ok]
      if (!KW) blocks <- blocks[y.ok]
      if (warn >= 0) warning(
        sum(!y.ok), " non-orderable values in ", y.name, 
        " encountered.")
      M  <- xTable(blocks,groups)
      KW <- (dim(M)[1]==1) # KW by elim of NA data
  }

  if ((!KW) && (any(delBlks <- rowSums(M[,-1]>0)<2))) {
    blocks.ok <- blocks %in% which(!delBlks)
    y      <- y     [blocks.ok]
    groups <- groups[blocks.ok]
    blocks <- blocks[blocks.ok]
    if (warn >= 0) warning( sum(!blocks.ok),
      " observation(s) in blocks with <2 of ", groups.name, 
      " removed.")
    M  <- xTable(blocks,groups)
  }

  n <- dim(M)[1]
  P <- dim(M)[2]-1        # actual number of grps (w/o NA)
  NAgrps <- sum(M[,1])    # all NAgrps could have been in non-info blks

  if (P<2) {
    m <- M
  } else {
      dimnames(m)[[2]][1] <- dimnames(M)[[2]][1] <- "0"
      if (n<dim(m)[1]) {                    # fewer blocks than planned
        m <- m[dimnames(M)[[1]],,drop=FALSE]
        blocks <- match(blocks,unique(blocks))
      }
      if (P<p) {
        m <- m[,dimnames(M)[[2]],drop=FALSE]
        groups <- match(groups,sort(unique(groups)))-(NAgrps>0)
                             # now, the dimnames don't match anymore ...
  }   }
  dimnames(M) <- NULL                       # ... but it doesn't matter
  dimnames(m) <- NULL

  if (KW <- (n == 1)) blocks <- rep(1,length(y))
  WMW <- (P==2) && (n==1)
  FM  <- all(M[,-1]==1)

# ##### *end*  count planned and observed frequencies (Wittkowski 1988)

  if (P>2)
    if (!(alternative=="two.sided"))
    stop("only two-sided alternatives meaningful for ",P," groups")

  if (is.null(exact)) exact <- ((P==2) && (n==1))      # Wilcox default
  if (exact) {
    if (exact && (!(exact <- (P==2))) && (warn >= 0))
      warning("can compute exact test only for two groups")
    if (exact && (!(exact <- (n==1))) && (warn >= 0))
      warning("can compute exact test only for one block")
    if (exact && (!(exact <- is.rank(Score(wRank)))))
      stop("can compute exact p-value only with ranks as scores")
    if (exact && (!(exact <- (sum(M)<50) )) && (warn >= 0))
      warning("cannot compute exact p-value for n larger than 50")
    if (exact && (!(exact <- (len(y)==len(unique(y))))) && (warn >= 0))
      warning("cannot compute exact p-value with ties")
  }

  if (is.null(correct)) correct <- ((P==2) && (n==1))  # Wilcox default
  if (correct) {
    if (correct && (!(correct <- (n==1))) && (warn >= 0))
      warning("argument 'correction = T' ignored for ",n," blocks")
    if (correct && (!(correct <- (P==2))) && (warn >= 0))
      warning("argument 'correction = T' ignored for ",p," groups")
    if (correct && (!(correct <- (!exact))) && (warn >= 0))
      warning("argument 'correction = T' ignored for exact tests")
  }

# ------------------------------------------ ** end ** parameter checks

  qapply <- MC(function(X, INDEX, FUN=NoOp, LEVELS=1:max(INDEX)) {
    asplit <- function(data, group, l = levels(group)) {
      u <- ifelse1(is.null(l), sort(unique(group)), seq(along = l))
      n <- length(u)
      group <- match(group, u, NA, FALSE)
      return(structure(split(data, group),              # ".Internal" removed
        names = ifelse1(is.null(l), as.character(u), l) ))
    }
    as.numeric(unlist(lapply(
      if(is.R()) split(X,INDEX) else asplit(X,INDEX,LEVELS), FUN)))
  },list(NoOp=NoOp))

  FUNByIdx <- MC(function(X, INDEX, FUN, resort=TRUE, na.any=NULL ) {
    if (is.null(na.any))
      na.any <- chkMissing(X)
    if (resort || na.any) {
      r <- qapply(X, INDEX, FUN)[order(qapply(1:length(X),INDEX))]
      if (na.any) r[which.na(X)] <- NA
    } else
      r <- qapply(X, INDEX, FUN)
    return(r)
  },list(chkMissing=chkMissing,qapply=qapply))

  Centered <- MC(function(scores,blocks, Mi) {
    scores - ifelse1(length(Mi)==1,
      sum(scores)/Mi,   # Mi must exclude grp==NA
      (qapply(scores,blocks,sum)/Mi)[blocks])
  },list(qapply=qapply))

  ginv <- function(x, tol = sqrt(.Machine$double.eps)) {   # simplified
    svdX <- svd(x);
    NEZ <- (svdX$d > (tol * svdX$d[1]))
    return(structure((if (!any(NEZ)) array(0, dim(x)[2:1]) else
      svdX$v[, NEZ] %*% ((1/svdX$d[NEZ]) * t(svdX$u[, NEZ]))),
      rank = sum(NEZ)))                              # not in MASS.ginv
  }

  outer1 <- function(X) {
    nx <- length(X)                  # required as a separate line in R
    structure( X * rep(X, each=nx), dim=c(nx,nx))
  }

  # Unconditional scores have to be 'conditioned' on NAgrps
  # Ie, within each block, the remaining ranks within each tie have to
  # be 'spread' by (n+1)/(n-NA+1) - may not be the most efficient way

  spread <- function(yg) {
    yg <- matrix(yg,ncol=2)
    y1 <- yg[,1] # y1 <- y     [blocks==1]
    g1 <- yg[,2] # g1 <- groups[blocks==1]
    if (any(!g1)) {
      ny <- length(y1)
      sy <- y1[sl <- sort.list(y1)]                    # sorted by ties
      tW <- ((tm <- diff(c(0,(1:ny)[c(sy[-1] != sy[-ny],TRUE)])))+1)/
            ((tM <- qapply(g1>0,y1,sum))                         +1)
      rM <- as.numeric(
        rep(cumsum(tm)-tm,tm) +      # tie base rank + spread increments
          qapply(rep(tW,tm),rep(1:length(tm),tm),cumsum) )
      rM[rM>rep(cumsum(tm),tm)] <- 0    # re adj rk exceeding tie max rk
      y1[sl] <- rM
    }
    return(y1)
  }
  spread <- MC(function(yg) {
    yg <- matrix(yg,ncol=2)
    y1 <- yg[,1] # y1 <- y     [blocks==1]
    g1 <- yg[,2] # g1 <- groups[blocks==1]
    if (any(!g1)) {
      ny <- length(y1)
      sy <- y1[sl <- sort.list(y1)]                    # sorted by ties
      tW <- ((tm <- diff(c(0,(1:ny)[c(sy[-1] != sy[-ny],TRUE)])))+1)/
            ((tM <- qapply(g1>0,y1,sum))                         +1)
      rM <- as.numeric(
        rep(cumsum(tm)-tm,tm) +     # tie base rank + spread increments
        qapply(rep(tW,tm),rep(1:length(tm),tm),cumsum) )
      rM[rM>rep(cumsum(tm),tm)] <- 0   # re adj rk exceeding tie max rk
      y1[sl] <- rM
    }
    return(y1)
  }, list(qapply=qapply))

  ################################################### *begin* Main Body

  if ((P<2 || sum(M[,-1])<3)) {

    T1   <- 0
    V1   <- ginv(V0 <- 0) # sets rank attribute
    eff  <- 0
    W    <- 0

  } else {

    if (mu!=0) y[groups==1] <- y[groups==1] - mu           # hypothesis

    mi <- rowSums(m)
    Mi <- rowSums(M)
    Wi <- switch(tolower(blkwght),
      prentice      = (mi+1),       # 1979, see Wittkowski 1988, p.1164
      klotz         = (Mi+1),       # 1980, see Wittkowski 1998, p.1164
      skillingsmack = sqrt(Mi+1),   # 1981, see Wittkowski 1998, p.1164
      rai           = (Mi+1)/Mi,    # 1987, see Wittkowski 1998, p.1164
      stop("illegal block weight blkwght = ",blkwght)
    )                                     # '/(Mi+1)' via Normalization
    Bijk <- Wi[blocks]                                    # block wghts

    Sijk <- Score(ifelse1(optim&KW,
      wRank(y)/(Mi+1),
      FUNByIdx(y,blocks,wRank,na.any=FALSE)/(Mi[blocks]+1) ))  # Scores

    if (NAgrps <- sum(M[,1])) {                   # exclude NAgrps data
      Sijk   <- Sijk  [groups!=0]
      Bijk   <- Bijk  [groups!=0]
      blocks <- blocks[groups!=0]
      groups <- groups[groups!=0]
      Mi <- rowSums(M[,-1,drop=FALSE])
    }
    M<- M[,-1,drop=FALSE]

    Tijk <- Centered(Sijk, blocks, Mi) * Bijk         # centered scores

    Sijk <- ifelse1(condvar,         # squared scores (Wittkowski 1998)
      Tijk,                       # cond var == ''correction'' for ties
      ifelse1(
        optim&KW, ((1:Mi)/(Mi+1)-0.5) * Bijk,
        optim&FM, ((1:P )/(P +1)-0.5) * Bijk,
        Centered(
          FUNByIdx(
            ifelse1(!NAgrps,
              qapply(1:length(Tijk), blocks, wRank),
              qapply(cbind(y,groups),blocks,spread)),
              blocks,Score,na.any=FALSE),
            blocks, Mi)
          )
        )^2

      # "unconditional" is still "conditional" on the missing groups,
      # among ties including NAgrp observations, the ranks are 'spread'

    T1 <- qapply(Tijk,groups,sum)                           # avg ranks

    if (optim&KW) {

      A0i2 <- (1/(Mi-1))*sum(Sijk)
    # V0  <- A0i2 * (diag(M[1,])-(1/Mi)*outer1(M[1,])) # not needed yet
      V1   <- (1/A0i2) * diag(1/M[1,]); attr(V1,"rank") <- P-1

    } else if (optim&FM) {

      A0i2 <- (1/(P-1))*sum(matrix(Sijk,nrow=n)[1,])
                        #   V0   <- n * A0i2[1] * (diag(rep(1,P))-(1/P))
      V1   <- (1/A0i2) * diag(rep(1/n,P)); attr(V1,"rank") <- P-1

    } else {

      A0i2 <- (1/(Mi-1))*qapply(Sijk,blocks,sum)              # std dev
      V0   <- structure(dim=c(P,P), A0i2 %*% (
        t(apply(M,1,   function(x,P) diag(x))) - (1/Mi) *
        t(apply(M,1,MC(function(x) outer1(x),list(outer1=outer1))))))
      V1   <- ginv(V0)
    }

    eff <- ifelse1(p==2, (T1[2]-T1[1])/2, stdev(T1,unbiased=FALSE))

    cc <- rep(0,P)
    if ((P==2) && correct) {              # Yates continuity correction
      cc[1] <- switch(alternative,
        two.sided = 0.5 * sign(T1[1]), greater = 0.5, less = -0.5 )
      cc[2] <- -cc[1]
    }

    W  <- as.numeric((T1-cc) %*% V1 %*% (T1-cc))            # test stat

  }

  df.W <- ifelse1(     # adjust for screening/variable subset selection
            df >0,df,                              # fixed df
            df==0,p-1,                             # all groups
            attr(V1,"rank"))                       # groups with data

  if (P>2) {                                               # chi-square

    p.W <- 1 - pchisq(W, df.W)

  } else if (!exact) {

    p.W <- switch(alternative,                           # alternatives
      greater   = 1 - pnorm(sqrt(W)*sign(eff)),
      less      =     pnorm(sqrt(W)*sign(eff)),
      two.sided = 2 * pnorm(-abs(sqrt(W)))    )

  } else {                                                      # exact

    e.sum    <- M[1] * (Mi+1)/2
    r.sum  <- sum(Tijk[groups==1]) + e.sum               # Wilcoxon W_R
    if (is.R()) {# R: pWilcox == "pMannWhitney"     # Mann Whitney W_XY
      r.sum <- r.sum - M[1]*(M[1]+1)/2          # Lehmann 174:10 (1.21)
      e.sum <- e.sum - M[1]*(M[1]+1)/2
    }
    p.W <- switch(alternative,
      greater   = 1 - pwilcox(r.sum - 1, M[1], M[2]),
      less      =     pwilcox(r.sum,     M[1], M[2]),
      two.sided = min(1, 2 * if(r.sum > e.sum)
        1 -  pwilcox(r.sum - 1, M[1], M[2])
        else pwilcox(r.sum,     M[1], M[2])))

  } 
##################################################### **end** Main Body

  ret.val <- list(
    effect      = eff,
    statistic   = W,
    parameters  = df.W,
    p.value     = p.W,
    alternative = alternative,
    method      = ifelse1(exact,"Exact ","") %&%
      "Prentice (Wilcoxon/Kruskal-Wallis/Friedman) rank sum test",
    data.name   = paste(y.name,
      groups.name, if (n>1) " and " %&% blocks.name else "", sep="")
    )

  names(ret.val$statistic)  <- "statistic: chi-square" %&%
     ifelse1( !correct, "", " with correction")
  names(ret.val$parameters) <- "df"

  ret.val <- ret.val[c("effect", "statistic", "parameters", "p.value",
    "alternative", "method", "data.name")]      # R requires long names
  oldClass(ret.val) <- "htest"
  return(ret.val)

} # endfunction(prentice.test) ########################################
