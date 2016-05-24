#### Testing  groupBound() , clusterGroupBound() , etc
####          ==========     =================

## Examples in  ../man/groupBound.rd  and ../man/clusterGroupBound.rd
## should *not*
## be duplicated here.  OTOH, want strict testing (*and* be speedy !)
stopifnot(require(hdi))

## From Matrix/.../tests-tools-1.R :
showSys.time <- function(expr, ...) {
  ## prepend 'Time' for R CMD Rdiff
  st <- system.time(expr, ...)
  writeLines(paste("Time", capture.output(print(st))))
  invisible(st)
}
doExtras <- function ()
{
  interactive() || nzchar(Sys.getenv("R_hdi_check_extra")) ||
    identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
do.Extras <- doExtras()
b64 <- .Machine$sizeof.pointer == 8
cat(if(b64) "64" else "32", "bit platform; do.Extras = ", do.Extras,"\n")

##' @title Generate "Reference" Datasets where variables are in *groups*
##' @param pvec integer vector of number of variables   p[k], in k-th group, k = 1..np.
##' @param s0.v integer vector of non-zero coefficient s0[k], in k-th group, k = 1..np.
##' @param n number of observations
##' @param ... further arguments such as \code{xtype}, \code{btype}, passed to
##'  \code{\link{rXb}}
##' @param iteration see \code{\link{rXb}}
##' @param do2S ditto
##' @return a \code{\link{list}} with components \code{x} and \code{beta} similar to those
##' of \code{rXb()} and \code{p.vec} and \code{s0.vec} which are the
##' corresponding  arguments.
##' @author Martin Maechler
rXbBlock <- function(p.vec, s0.vec, n, ..., iteration = NA, do2S = FALSE) {
  stopifnot(is.numeric(p.vec), p.vec == round(p.vec),
            is.numeric(s0.vec), s0.vec == round(s0.vec),
            (np <- length(p.vec)) >= 1, length(s0.vec) == np,
            0 <= s0.vec, s0.vec <= p.vec)
  grL <-
    lapply(1:np, function(k)
      rXb(n = n, p = p.vec[k], s0 = s0.vec[k],
          ..., iteration=iteration, do2S=do2S))
  ## now concatenate the X matrices and the beta vectors:
  list(x    = do.call(cbind, lapply(grL, `[[`, "x")),
       beta = do.call(c,     lapply(grL, `[[`, "beta")),
       p.vec = p.vec, s0.vec = s0.vec)
}

##' @title From (X, beta, sigma),  generate  y = X b + eps
##' @param rXb a list, e.g. resulting from rXb (== rXb) or rXbBlock()
##' @param sigma a non-negative number
##' @return a list with components (x, y)
##' @author Martin Maechler
rXy <- function(rXb, sigma = 1) {
  stopifnot(is.list(rXb), is.matrix(X <- rXb$x), is.numeric(b <- rXb$beta),
            ncol(X) == length(b),
            length(sigma) == 1, is.finite(sigma), 0 <= sigma)
  list(x = X, y = drop(X %*% b) + rnorm(nrow(X), sd = sigma))
}

##' from a p-vector, return corresponding list of 'groups', i.e., group indices
pvec2grp <- function(pv) {
  stopifnot(is.numeric(pv), is.finite(pv), pv >= 0)
  cv <- cumsum(c(1,pv))
  mapply(seq, from=cv[-length(cv)], to = cv[-1L]-1L)
}

if(FALSE) { ## not used currently
  set.seed(47) # <-> iteration = NA to use R's seed
  d1 <- rXb(n = 64, p =  16, s0 =  7, iteration=NA, do2S = FALSE)
  d2 <- rXb(n = 64, p = 256, s0 = 12, iteration=NA, do2S = FALSE)
  stopifnot(identical(dim(d2$x), c(64L, 256L)), sum(d2$beta != 0) == 12)

  yhat.1 <- with(d1, drop(x %*% beta))
  yhat.2 <- with(d2, drop(x %*% beta))
}

set.seed(20)
d1B <- rXbBlock(p = c(7, 8, 7), s0 = c(2, 4, 1), n = 16)
p2 <- sample(c(32, 16, 16, 16, 12, 12, 5, 5, 5))
s2 <- c(rmultinom(1, size = p2, prob = ifelse(p2 > 15, 1/4, 1/8)))
rbind(p2, s2)
d2B <- rXbBlock(p = p2, s0 = s2, n = 500)
d3B <- rXbBlock(p = p2, s0 = s2, n = 128, xtype = "toeplitz", x.par = 0.95)
d4B <- rXbBlock(p = p2, s0 = s2, n = 128, xtype = "equi.corr")
d5B <- rXbBlock(p = p2, s0 = s2, n = 128, xtype = "exp.decay", x.par = c(0.25, 4))

if(interactive() && exists("db5", mode="list")) ## save when produced with "new hdi"
  save(d1B, d2B, d3B, d4B, d5B, file = "rXb_blk.rda")
if(!exists("rXb", mode="function")) ## "old hdi" -- tests should give the same (?!) -- TODO
  print(load("rXb_blk.rda"))

if(dev.interactive(TRUE)) {
  require(Matrix)
  ## if you look carefully, you recognize the three blocks of sizes (7, 8, 7):
  image(as(zapsmall(cov(d1B$x), digits=1), "sparseMatrix"))
  ## Here, have 9 blocks :
  image(as(zapsmall(cov(d2B$x), digits=1), "sparseMatrix"))
  image(as(zapsmall(cov(d3B$x), digits=1), "sparseMatrix")) # easy to see
  image(as(zapsmall(cov(d4B$x), digits=1), "sparseMatrix")) #  (ditto)
  image(as(zapsmall(cov(d5B$x), digits=1), "sparseMatrix"))
                                        # decaying quickly => basically tri-diagonal blocks
}

##- set.seed(107)
##- d1xy <- rXy(d1B, sigma = 0.5)
##- showSys.time(
##-   gb1 <- groupBound(x = d1xy$x, y = d1xy$y, group = 1:ncol(d1xy$x))
##- ) # 2.00 [lynne 2015]
##- set.seed(17)# also for random splits
##- d2xy <- rXy(d2B, sigma = 1/16)
##- showSys.time(
##-   gb2  <- groupBound(x = d2xy$x, y = d2xy$y, group = 1:ncol(d2xy$x))
##- )
##- showSys.time(
##-   gb2g <- groupBound(x = d2xy$x, y = d2xy$y, group = pvec2grp(d2B$p.vec),
##-                      silent = TRUE)
##- )
##- 
##- stopifnot(
##-   as.vector(gb1) == 0
##-  ,
##-   all.equal(as.vector(gb2),
##-             if(b64) 6.529515 else 6.598612, tol = 2e-7)
##-  ,
##-   all.equal(as.vector(gb2g),
##-             if(b64) c(0, 0, 1.924578, 0, 1.436938, 0, 0, 0, 0.3264166)
##-             else c(0, 0.01877908, 2.0447942, 0, 1.753219, 0,0,0, 0.3779607),
##-             tol = 6e-7)
##- )

if(do.Extras) { ## save time in regular  'R CMD check' ---------

  ##--------------- 3 --
  set.seed(27)# also for random splits
  d3xy <- rXy(d3B, sigma = 1/16)
  showSys.time(
    gb3  <- groupBound(x = d3xy$x, y = d3xy$y, group = 1:ncol(d3xy$x))
  )
  stopifnot(
    all.equal(as.vector(gb3), if(b64) 7.3760622 else 7.7908584, tol = 2e-7)
  )
  set.seed(19)# for random splits
  showSys.time(
    gb3g <- groupBound(x = d3xy$x, y = d3xy$y, group = pvec2grp(d3B$p.vec),
                       silent = TRUE)
  )
  stopifnot(
    all.equal(as.vector(gb3g),
              if(b64) c(0, 0.5237811, 1.354151, 0, 0, 0, 0.8847217, 0, 0)
              else    c(0, 1.1008886, 1.9249053,0, 0, 0, 0.78661001, 0, 0),
              tol = 5e-7) # 9.4e-8
  )

  ##--------------- 4 --
  set.seed(37)# also for random splits
  d4xy <- rXy(d4B, sigma = 1/16)
  showSys.time(
    gb4  <- groupBound(x = d4xy$x, y = d4xy$y, group = 1:ncol(d4xy$x),
                       silent = TRUE)
  )
  showSys.time(
    gb4g <- groupBound(x = d4xy$x, y = d4xy$y, group = pvec2grp(d4B$p.vec),
                       silent = TRUE)
  )
  stopifnot(
    all.equal(as.vector(gb4), if(b64) 5.230821 else 5.3752713,  tol = 2e-7)
    ## all.equal(as.vector(gb4), 8.475868, tol = 2e-7)
   ,
    all.equal(as.vector(gb4g),
              if(b64) c(0, 0, 3.328455, 0, 0, 0, 0, 0, 0)
              else    c(0, 0.043785673, 3.6154918, 0, 0, 0, 0, 0, 0),
              tol = 8e-7) # 1.4e-7
  )

}# only if(do.Extras)
