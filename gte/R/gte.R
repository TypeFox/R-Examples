#' Generalized Turnbull's Estimator
#' 
#' The \code{gte} function computes the generalized Turnbull's Estimator (GTE) proposed by Dehghan and Duchesne (2011).
#' It is a nonparametric estimator of a conditional survival function given a vector of continuous covariates 
#' that can handle interval-censored lifetimes.\cr
#' The \code{print} method for objects obtained from \code{gte} only prints the output value \code{surv.summary}.\cr
#' The \code{plot} method for objects obtained from \code{gte} plots the estimate of the conditional survival
#' function, by default overlaying curves if more than one estimate is present and shading the innermost interval,
#' in which the GTE is indeterminate. 
#' 
#' For interval-censored data, the \code{\link{Surv}} function should be called with the argument 
#' \code{type="interval"} or \code{type="interval2"}. If \code{type="interval"}, the \code{event} argument
#' is mandatory. Therefore, in addition to the left and right endpoints of the censoring interval (called, respectively,
#' \code{left} and {right} for illustrative purpose), one would need a third variable (\code{status}) taking the value 0 
#' for right censored data, 1 for an event at exact time, 2 for left censored data and 3 for interval censored data. 
#' The \code{\link{Surv}} function would be called as follows:\cr
#' \code{Surv(time=left, time2=right, event=status, type="interval")}.
#' 
#' If \code{type="interval2"}, the \code{event} argument cannot be given. 
#' The value of \code{event} is derived from the \code{time} and \code{time2} argument as follows:\cr 
#' if \code{time} takes the value \code{NA}, \code{event=2} (left censored data);\cr 
#' if \code{time2} takes the value \code{NA}, \code{event=0} (right censored data);\cr
#' if \code{time=time2}, \code{event=1} (exact time);\cr 
#' otherwise, \code{event=3} (interval censored data).\cr
#' See the help page of the \code{\link{Surv}} function for more details.\cr
#'  
#' In the \code{gte} function, the data must be given through the \code{\link{Surv}} function
#' but it is internally transformed in two vectors : \code{L} and \code{R} for the left and right endpoints of 
#' the censoring interval, respectively.\cr     
#' If \code{event=0} (right censored data), then \code{L=time} and \code{R=Inf};\cr     
#' if \code{event=1} (exact time), then \code{L=time} and \code{R=time};\cr     
#' if \code{event=2} (left censored data), then \code{L=0} and \code{R=time};\cr     
#' and if \code{event=3} (interval censored data), then \code{L=time} and \code{R=time2};\cr
#' If one has vectors \code{L} and \code{R} respecting this convention, they can be given directly to \code{gte}
#' by calling \code{\link{Surv}} as follows:\cr
#' \code{Surv(L, R, type="interval2")}.     
#' 
#' @param formula A formula object with the response on the left of a ~ operator, and the covariates on the right. 
#'                The response must be a survival object as returned by the \code{Surv} function from the package 
#'                \pkg{survival} (see \bold{Details}).
#' @param data An optional data frame, list or environment  
#'             containing the variables in the model formula. If not found in data, the variables are taken from 
#'             \code{environment(formula)}, typically the environment from which \code{gte} is called.
#' @param z A matrix: each row contains the values of a covariate vector at which an estimate of the 
#'          conditional survival function is requested.
#'          If there is only one covariate, it can be a vector (possibly of length 1).  
#' @param h A vector: the values of the bandwidth parameter \eqn{h} for each covariate
#'          (default = equation 7 of Dehghan and Duchesne (2011)).
#' @param itermax maximal number of iterations for the algorithm (default=100000).
#' @param tole maximal distance between successive iterations tolerated
#'             before declaring convergence (default=0.0005). 
#' 
#' @return \item{time}{ A vector: the ordered distinct values of the left and right endpoints of the censoring 
#'                      interval (omitting the smallest value, but always including time 0).} 
#' @return \item{surv}{ A matrix: the estimates of the conditional survival function at time \code{time}. 
#'                      The \eqn{i}th column refers to the \eqn{i}th value of the covariate vector given 
#'                      in \code{z} (row \eqn{i} of \code{z}). }
#' @return \item{intmap}{A matrix : The intervals of the potential steps in the conditional survival function,
#'                       called innermost interval, over which the GTE is indeterminate. The left endpoints of
#'                       the intervals are in the first row, and the rigth endpoints in the second. The object
#'                       attribute LRin denotes whether to include each of the endpoints or not. This matrix
#'                       is computed with an internal function derived from function \code{Aintmap} of the
#'                       \pkg{interval} package.}
#' @return \item{surv.summary}{A summary of \code{surv}: the estimates of the conditional survival function 
#'                             only for the intervals of the potential steps in the function (innermost 
#'                             intervals). The row names describe the intervals, which are detailed in 
#'                             \code{intmap}. } 
#' @return \item{Call}{ The function call.}
#' 
#' @seealso \code{\link{Surv}}
#' @author Mohammad Hossein Dehghan, Thierry Duchesne and Sophie Baillargeon
#' @references Dehghan, M. H. and Duchesne, T. (2011). A generalization of Turnbull's estimator for 
#' nonparametric estimation of the conditional survival function with interval-censored data.
#' \emph{Lifetime Data Analysis}, \bold{17}, 234-255.
#' @useDynLib gte
#' @export
#' @importFrom survival is.Surv Surv
#' @examples
#' data(simul)
#' 
#' ## Calling Surv() with type="interval2"
#' Fit <- gte(Surv(L, R, type="interval2") ~ Z, data=simul, z=c(10, 20))
#' Fit
#' 
#' ## Calling Surv() with type="interval"
#' event <- ifelse(is.na(simul$R), 0,
#'                 ifelse(is.na(simul$L), 2,
#'                        ifelse(simul$R==simul$L, 1, 3)))
#' time <- ifelse(event==2, simul$R, simul$L)
#' time2 <- ifelse(event==3, simul$R, NA)
#' simul_event <- cbind(simul, time, time2, event)
#' 
#' Fit_event <- gte(Surv(time, time2, event, type="interval") ~ Z, data=simul_event, z=c(10, 20))
#' Fit_event
#' 
#' # The results are the same
#' all.equal(Fit_event$time, Fit$time)
#' all.equal(Fit_event$surv, Fit$surv)
#' 
#' ## Plotting the results
#' plot(Fit, xleg="topright")
#' 
gte <- function(formula, data, z, h=NULL, itermax=100000, tole=0.0005){
  
  call <- mfcall <- match.call()
  
  ################################################################# 
  #           To transform the given argument into the
  #           objects needed for the calculation
  ################################################################# 
  
  mfargs <- match(c("formula", "data"), names(mfcall), 0L)
  mfcall <- mfcall[c(1L, mfargs)]
  mfcall[[1L]] <- as.name("model.frame")
  mf <- eval(mfcall, parent.frame())
  
  ### Validation de la formula
  Y <- model.extract(mf, "response")
  if (!is.Surv(Y)) stop("in 'formula', the response on the left of the ~ operator must be a Surv object")
  
  ### Creation de R et L
  # A partir de mf$Surv, je vais creer L et R tel que decrit dans la section Details de la documentation
  # Ce code est une modification de la fonction SurvLR du package interval
  type <- attr(mf$Surv, "type")
  if (type == "right") {
    L <- mf$Surv[,1]
    R <- ifelse(mf$Surv[,2]==0, Inf, mf$Surv[,1] + 0.0000001)
  }
  else if (type == "counting") {
    stop("Surv object type='counting' not supported")
  }
  else if (type == "left") {
    L <- ifelse(mf$Surv[,2]==0, 0, mf$Surv[,1])
    R <- ifelse(mf$Surv[,2]==0, mf$Surv[,1], mf$Surv[,1] + 0.0000001)
  }
  else if (type == "interval") {
    L <- ifelse(mf$Surv[,3]==2, 0, mf$Surv[,1])
    R <- ifelse(mf$Surv[,3]==0, Inf, 
         ifelse(mf$Surv[,3]==3, mf$Surv[,2], 
         ifelse(mf$Surv[,3]==1, mf$Surv[,1] + 0.0000001, mf$Surv[,1])))
  }
  else {
    stop(paste("Surv obj type='", type, "' unrecognized", 
               sep = ""))
  }
# Remarque : j'ai reussi a completement remplacer les 0.0000001 par des 0 pour la censure a gauche
# et les 10000 pour des Inf pour la censure a droite. Par contre, pour les evenement a un temps exact,
# on peut fournir en entree L=R, mais dans le code on doit prendre R = L + 0.0000001, sinon les calculs
# ne sont pas bons. Mais j'enleve les lignes se referant aux R = L + 0.0000001 dans les sorties. 
  L <- pmax(L, 0)  ## Pour ramener a zero les temps negatifs
  R <- pmax(R, 0)
  N <- length(L)
              
  ### Objets relatifs aux covariables
  mm <- model.matrix(formula, data=mf) 
  # mm = matrice avec premiere colonne de 1 (intercept) et toutes les autres colonnes sont les covariables
  d <- ncol(mm) - 1 # nombre de covariables
  if (d==0) stop("at least one numeric covariate must be given in the formula")
  
  ### mise en forme et validation de z
  if(is.list(z)) z <- as.matrix(z)
  if(!is.numeric(z)) stop("'z' must be numeric")
  if (is.null(dim(z))){
    if (d==1) z <- as.matrix(z) else {
      if (length(z)!=d) stop("'z' must be of length ", d, ", the number of covariates") else z <- matrix(z, nrow=1)
    }
  } else {
    if (ncol(z) != d) stop("'z' must have ", d, " columns, i.e. the number of covariates")
  }
  
  ### traitement de h
  if (is.null(h)){
    # Valeur par defaut = equation 7 de l'article
    h <- 1.06*apply(mm[, -1, drop=FALSE], 2, sd, na.rm=TRUE)*N^(-1/5)
  } else {
    if(!is.numeric(h)) stop("'h' must be numeric")
    if (length(h)==1) h <- rep(h, d) else
      if (length(h)!=d) stop("'h' must be of length ", d, ", the number of covariates")
  }

  ### validation de itermax et tole
  if(!is.numeric(itermax) || length(itermax)!=1 || itermax<=0 || (itermax %% 1) != 0) 
    stop("'itermax' must be a single positive integer")
  if(!is.numeric(tole) || length(tole)!=1 || tole<=0) 
    stop("'tole' must be a single positive real")
  
  ################################################################# 
  #           To calculate  the indicator  matrix "Alpha1" 
  ################################################################# 
  
  LR1 <-  sort(unique(c(L,R))) 
  M <- length(LR1)-1 
  Alpha1 <-  matrix(999, ncol=M, nrow=N)  # Allocated  matrix for the indicator
  
  # Pour remplacer les valeurs Inf par un nombre plus grand que toutes les valeurs finies de L et R
  SupLR <- max(LR1[is.finite(LR1)]) + 1
  if (any(!is.finite(R))){
    R[!is.finite(R)] <- SupLR
    LR1[!is.finite(LR1)] <- SupLR
  }
  
  GetAlph.out <- .C("GetAlpha",as.double(LR1),as.double(L),as.double(R), 
      as.integer(M),  as.integer(N),alpha=as.double(Alpha1), PACKAGE="gte")
  Alpha1  <- GetAlph.out$alpha 
  
  ################################################################# 
  #           To find intervals of the potential jumps
  ################################################################# 
  
  intmap <- get.intmap(L = L, R = R)
  indic <- LR1[-1] %in% intmap[2,]
  
  ################################################################# 
  #    To calculate  the GT estimator   
  ################################################################# 
  
  S <- matrix(NA, nrow=M, ncol=nrow(z))
  colnames(S) <- if (length(z)==d) { # estimation en un seul point
      "survival|Z=z" 
    } else { 
      paste("survival|Z=z[", 1:nrow(z), ifelse(d==1, "", ", "), "]", sep="")
    }
  for (i in 1:nrow(z)){
    W <- 1    
    for (j in 1:d){
      Wj <- if(h[j]==0) { 
              rep(1, N) 
            } else {
              dnorm(mm[, -1, drop=FALSE], z[i, j], h[j])#  Using  the Gaussian kernel  function 
            }
      W <- W*Wj
    }
    W <- W/sum(W)
    
    # the intial death probability is non nul only for intervals of the potential jumps
    p <- rep(0,M)  # Initial  death probability  of B_j, j=1,...,M
    p[indic] <- 1/sum(indic)
    
    St2.out <- .C("St2",as.double(Alpha1),as.double(W),as.integer(N),as.integer(M), 
        as.integer(itermax),as.double(tole),p=as.double(p),ESP=as.double(1), PACKAGE="gte") 
    
    p <- St2.out$p 
    Ft2 <- cumsum(p)
    S[,i] <- 1 - Ft2
    S[,i] <- pmax(S[,i], 0)
  }
  

  ################################################################# 
  #    Output values adjustment   
  #################################################################
 
  # time variable for output
  time <- LR1[-1] ## this associate each survival rate to the upper bound of the corresponding interval
  time <- ifelse(time == SupLR, Inf, time)
  
  # S'il y avait des L=R (evenement a un temps exact) en entree
  id_exact <- if(type=="interval") mf$Surv[,3]==1 else mf$Surv[,2]==1
  if (sum(id_exact)>0){
    # Pour enlever les lignes des temps L et changer les temps
    # R pour leur valeur exact 
    toremove <- time %in% L[id_exact]
    time <- time[!toremove]
    S <- S[!toremove, , drop=FALSE]
    tocorrect <- time %in% R[id_exact]
    time[tocorrect] <- time[tocorrect] - 0.0000001
    
    # Correction in intmap for exact values, 
    ev <- intmap[1,] %in% L[id_exact] # indicator for exact interval
    intmap[2, ev] <- intmap[1, ev]
    attr(intmap, which = "LRin")[1, ev] <- TRUE
  }

  # Add time 0
  time <- c(0, time)
  S <- rbind(rep(1, ncol(S)), S)

  # Correction in intmap for Inf which has been replaced by max(LR1) + 1.
  if(intmap[2, ncol(intmap)] == SupLR){
    intmap[2 , ncol(intmap)] <- Inf
    attr(intmap, which = "LRin")[2 , ncol(intmap)] <- FALSE
  } 
  
  
  ################################################################# 
  #    Summary table   
  #################################################################
  surv.summary <- S[time %in% intmap[2,], ,drop=FALSE]
  rownames(surv.summary) <- paste(ifelse(attr(intmap, which = "LRin")[1,], "[", "("), 
                                  intmap[1,], ", ", intmap[2,], 
                                  ifelse(attr(intmap, which = "LRin")[2,], "]", ")"), sep="") 
  
  ################################################################# 
  #    Output   
  #################################################################
  out <- list(time=time, surv=S, intmap=intmap, surv.summary=surv.summary, call=call)
  class(out) <- "gte"
  out  
}
  
#' @rdname gte
#' @method print gte
#' @param x An object, produced by the \code{gte} function, to print or to plot.
#' @export
print.gte <- function(x, ...) {
  print.default(x$surv.summary, print.gap = 2, quote = FALSE, right=TRUE, ...)  
  invisible(x)
}

#' @rdname gte
#' @method plot gte
#' @param overlay A logical: Should the curves be overlayed when there is more than one estimate 
#'                of the conditional survival function in the \code{gte} object \code{x}? (default=\code{TRUE})
#' @param shade A logical: Should the rectangles of indeterminate NPMLE (innermost interval) 
#'              be shaded? (default=\code{TRUE})
#' @param xlab A label for the x-axis, by defaut \code{xlab = "time"}.
#' @param ylab A label for the y-axis, by defaut \code{ylab = "survival"}.
#' @param xleg x location for legend, "bottomleft" by default (see \code{\link{legend}}).
#' @param yleg y location for legend, NULL by default (see \code{\link{legend}}).
#' @param \dots Further arguments to be passed to \code{print.default} or \code{plot.default}. 
#' @export
plot.gte <- function(x, overlay = TRUE, shade = TRUE, xlab = "time", ylab = "survival", 
                     xleg="bottomleft", yleg=NULL, ...) {
  
  # Create time and S, and truncate them to the last value of intmap
  time <- x$time[x$time <= x$intmap[2, ncol(x$intmap)]]
  S <- x$surv[x$time <= x$intmap[2, ncol(x$intmap)], , drop=FALSE]
  
  # Draw the plot
  plot(x = range(time[is.finite(time)]), y = c(0,1), xlab = xlab, ylab = ylab, type="n", ...)
  op <- par()

  # Check overlay value 
  if (!overlay && ncol(S)>1) {
    warning("'overlay' was set to FALSE and 'x' contains more than one estimate of the conditional survival function: only the first estimate is plotted")
    S <- S[, 1, drop=FALSE]
  }
  
  # Add shade for rectangles of innermost intervals
  if(shade){
    interval <- x$intmap[, x$intmap[1,] != x$intmap[2,], drop=FALSE]
    if(ncol(interval)>0){
      seqcol <- if (ncol(S) == 1) 0.8 else
                if (ncol(S) == 2) c(0.7, 0.8) else
                if (ncol(S) == 3) c(0.6, 0.7, 0.8) else 
                if (ncol(S) == 4) c(0.6, 0.7, 0.8, 0.9) else
                                  seq(0.6, 0.9, length=ncol(S))
      for (j in 1:ncol(S)){
        for (i in 1:ncol(interval)){
          xInf <- interval[1, i]
          xSup <- interval[2, i] 
          yInf <- min(S[time==xSup, j])
          ySup <- max(S[time==xInf, j])
          if(!is.finite(xSup)) xSup <- op$usr[2] 
          polygon(x = c(xInf, xInf, xSup, xSup), y = c(yInf, ySup, ySup, yInf), col = gray(seqcol[j]),  border = NA)
        }
      }
    }
  }

  # If Inf is present, remove that point
  if (!is.finite(time[length(time)])){
    time <- time[-length(time)]
    S <- S[-nrow(S), , drop=FALSE]
  }
  
  # Draw lines for the survival functions
  seqlty <- 1:ncol(S)
  for (j in 1:ncol(S)) lines(x = time, y = S[,j], lty=seqlty[j])
  
  # Add a legend if more than one line was drawn
  if (ncol(S)>1){ 
    legend(x=xleg, y=yleg, legend=substr(colnames(S), 12, nchar(colnames(S))), lty = seqlty)
  }
  
}

