### This file contains all wraps to call C functions in "src/R_init_EM.c".
### Modified: Wei-Chen Chen on 2009/04/29.

# Call:
#   SEXP R_init_EM(SEXP x, SEXP n, SEXP p, SEXP nclass,
#                  SEXP short_iter, SEXP short_eps, SEXP fixed_iter,
#                  SEXP EM_iter, SEXP EM_eps,
#                  SEXP lab, SEXP labK,
#                  SEXP init_method)
# Input:
#   x: SEXP[n, p], a data matrix of n*p.
#   n: SEXP[1], a number of observations.
#   p: SEXP[1], a number of dimersions.
#   nclass: SEXP[1], a number of classes.
#   short_iter: SEXP[1], a number of short em iterations, 500 by default.	# short.iter
#   short_eps: SEXP[1], a tolerance of short em, 1e-2 by default.		# short.eps
#   fixed_iter: SEXP[1], a number of rand iterations, 1 by default.		# fixed.iter
#   n_candidate: SEXP[1], a number of candidates, 5 by default.
#   EM_iter: SEXP[1], a number of EM iterations, 1000 by default.		# EM.iter
#   EM_eps: SEXP[1], a tolerance of EM, 1e-4 by default.			# EM.eps
#   lab: SEXP[n], -1 for points with unknown clusters;
#                 0,..,(labK-1) for known.
#   labK: SEXP[1], the number of known clusters.
#   init_method: SEXP[1], initialization method. (see below for detail)
#####
#   unsupervised init_method:
#     1 = em.EM,     2 = Rnd.EM
#     11 = MBem.EM, 12 = MBRnd.EM
#     21 = acem.EM, 22 = acRnd.EM
#   semi-supervised init_method:
#     101 = em.EM,   102 = Rnd.EM
#     111 = MBem.EM, 112 = MBRnd.EM
#     121 = acem.EM, 122 = acRnd.EM
#   For iteration & time record:
#     -1 = oneRnd,      -11 = oneMBRnd
#     -101 = ss.oneRnd, -111 = ss.oneMBRnd
#####
# Output in C:
#   ret: a list contains
#     pi: SEXP[nclass], proportions of classes.
#     Mu: SEXP[nclass, p], means of MVNs.
#     LTSigma: SEXP[nclass, p * (p + 1) / 2], lower triangular sigma matrices.
#     llhdval: SEXP[1], log likelihood value.
#     nc: SEXP[nclass], numbers of observations in each class.
#     class: SEXP[n], class id's for all observations
#            starting from 0 to (nclass - 1).
# Output in R:
#     n: SEXP[1], a number of observations.
#     p: SEXP[1], a number of dimersions.
#     nclass: SEXP[1], a number of classes.
#     method: SEXP[1], a initialization method.
init.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC,
    stable.solution = TRUE, min.n = NULL, min.n.iter = 10,
    method = c("em.EM", "Rnd.EM")){
#    method = c("em.EM", "Rnd.EM", "MBem.EM", "MBRnd.EM")){
#    method = c("em.EM", "Rnd.EM", "MBem.EM", "MBRnd.EM",
#               "acem.EM", "acRnd.EM")){

#  init.code <- switch(method,
#                      "em.EM" = 1, "Rnd.EM" = 2,
#                      "MBem.EM" = 11, "MBRnd.EM" = 12,
#                      "acem.EM" = 21, "acRnd.EM" = 22,
#                      "oneRnd" = -1, "oneMBRnd" = -11)
  init.code <- switch(method[1],
                      "em.EM" = 1, "Rnd.EM" = 2,
                      "MBem.EM" = 11, "MBRnd.EM" = 12)
  n <- nrow(x)
  p <- ncol(x)

  if(! is.null(lab)){
    init.code <- init.code + 100
    labK <- max(lab)
    lab <- lab - 1
    if(length(unique(lab[lab != -1])) != labK) stop("lab is not correct.")
    if(labK > nclass) stop("lab is not correct.")
  } else{
    labK <- NULL
  }

  if(! stable.solution) min.n.iter <- 0

  total.min.n.iter <- 1
  flag <- 0
  if(is.null(min.n)) min.n <- p + 1
  repeat{
    ret <- .Call("R_init_EM",
                 as.double(t(x)),
                 as.integer(n),
                 as.integer(p),
                 as.integer(nclass),
                 as.integer(EMC$short.iter),
                 as.double(EMC$short.eps),
                 as.integer(EMC$fixed.iter),
                 as.integer(EMC$n.candidate),
                 as.integer(EMC$EM.iter),
                 as.double(EMC$EM.eps),
                 as.integer(lab),
                 as.integer(labK),
                 as.integer(init.code))
    if(all(ret$nc >= min.n)) break

    total.min.n.iter <- total.min.n.iter + 1
    if(total.min.n.iter > min.n.iter){
      flag <- 1
      break
    }
  }

  if(stable.solution && flag == 1){
    warning("A stable solution is not avaliable.")
  }
  ret$flag <- flag

  ret$pi <- ret$pi / sum(ret$pi)
  ret$Mu <- matrix(ret$Mu, nrow = nclass, byrow = TRUE)
  ret$LTSigma <- matrix(ret$LTSigma, nrow = nclass, byrow = TRUE)
  ret$class <- ret$class + 1

  ret$n <- n
  ret$p <- p
  ret$nclass <- nclass
  ret$method <- method
  class(ret) <- "emret"
  ret
}


### Original version
em.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC,
    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
  init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
          stable.solution = stable.solution, min.n = min.n,
          min.n.iter = min.n.iter, method = "em.EM")
}
rand.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC.Rnd,
    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
  init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
          stable.solution = stable.solution, min.n = min.n,
          min.n.iter = min.n.iter, method = "Rnd.EM")
}
#ac.em.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC,
#    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
#  init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
#          stable.solution = stable.solution, min.n = min.n,
#          min.n.iter = min.n.iter, method = "acem.EM")
#}
#ac.rand.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC.Rnd,
#    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
#  init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
#          stable.solution = stable.solution, min.n = min.n,
#          min.n.iter = min.n.iter, method = "acRnd.EM")
#}
mb.em.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC,
    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
  init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
          stable.solution = stable.solution, min.n = min.n,
          min.n.iter = min.n.iter, method = "MBem.EM")
}
mb.rand.EM <- function(x, nclass = 1, lab = NULL, EMC = .EMC.Rnd,
    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
  init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
          stable.solution = stable.solution, min.n = min.n,
          min.n.iter = min.n.iter, method = "MBRnd.EM")
}


### Exhausted version
exhaust.EM <- function(x, nclass = 1, lab = NULL,
    EMC = .EMControl(short.iter = 1, short.eps = Inf),
#    method = c("em.EM", "Rnd.EM", "MBem.EM", "MBRnd.EM", "abem.EM", "abRnd.EM"),
    method = c("em.EM", "Rnd.EM"),
    stable.solution = TRUE, min.n = NULL, min.n.iter = 10){
  llhdval <- -Inf

  for(i in 1:EMC$exhaust.iter){
    ret.new <- init.EM(x, nclass = nclass, lab = lab, EMC = EMC,
                       stable.solution = stable.solution, min.n = min.n,
                       min.n.iter = min.n.iter, method = method[1])
    if(ret.new$llhdval > llhdval){
      llhdval <- ret.new$llhdval
      ret <- ret.new
    }
  }

  ret
}

