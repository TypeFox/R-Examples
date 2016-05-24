gvcm.cat.flex <- function(

whichCoefs,   # vector with covariates (as characters)
intercept = TRUE, # logical
data,         # dataframe with named and coded covariates
family = gaussian(), # family object

method = "REML",# for function gam
tuning = NULL,  # for function gam: argument sp

indexNrCoefs,   # vector with number of coefficients per covariate
indexPenNorm,   # vector with norm of the employed penalty (as.character)
indexPenA,      # list with the penalty matrices A_j for each covariate j, t(A_j)%*%beta_j 
indexPenWeight, # list, possible weights for the penalty terms (each entry is a vector)

#weights,
#offset, 
#start,

control = list(c=1e-05, epsilon=1e-07, gama=35, maxi=1500, nu=.5) 

#model = FALSE,
#x = FALSE,
#y = FALSE,
#plot=FALSE,

)
{
newton = list(
conv.tol = 1e-6, # relative convergence tolerance
maxNstep = 15,  #  maximum length allowed for an element of the Newton search direction (default 5)
maxSstep = 5, # is the maximum length allowed for an element of the steepest descent direction (only used if Newton fails - default 2);
maxHalf = 50 # is the maximum number of step halvings to permit before giving up (default 30).
)

#nlm <- list(
#ndigit, # number of significant digits in the GCV/UBRE score - by default this is worked out from epsilon;
#gradtol, # tolerance used to judge convergence of the gradient of the GCV/UBRE score to zero - by default set# to 10*epsilon;
#stepmax=2, # maximum allowable log smoothing parameter step
#steptol=1e-4, # the minimum allowable step length
#iterlim=200, # the maximum number of optimization steps allowed
#check.analyticals # indicates whether the built in exact derivative calculations should be checked numerically - defaults to FALSE.
#)

gamcontrol <- gam.control(newton=newton)
optimizer <- c("outer", "newton")
optimizer <- c("outer", "newton")


# definitions
  p <- length(whichCoefs)
  vonBis <- matrix(nrow=2, ncol=p)
  vonBis[1, ] <- c(0, cumsum(indexNrCoefs)[-p])+1
  vonBis[2, ] <- cumsum(indexNrCoefs)
  if(intercept) {vonBis <- vonBis+1}

# functions
  # Penalty functions
      L1 <- function(xi, control=control) sqrt((xi^2 + control$c) )^(-1)
      control$L0.log <- TRUE
      L0 <- if (control$L0.log==TRUE) {
            function(xi, control=control){p <- 1+exp(-control$gama*abs(xi))
                   2*control$gama*(sqrt(xi^2 + control$c))^(-1)*p^(-1)*(1 - 1/p)}
            } else {
            function(xi, control=control){2*(xi+control$c)^(-2)} 
            }
      L2 <- function(xi, control=control){2}
      elastic <- function(xi, control=control) {control$elastic * L1(xi, control=control) +
                       (1-control$elastic) * L1(xi, control=control)}    
      grouped <- function(beta, Aj, control=control) { # Aj matrix, dim(Aj)= q x Lj
        dfj <- ncol(Aj) 
        AjtAj <- Aj%*%t(Aj)
        rep(sqrt(dfj) / sqrt(t(beta)%*%AjtAj%*%beta + control$c), ncol(Aj) ) 
        }

  # call penalty
  A <- function(beta, PenNorm, PenA, PenWeight, control)
  {
    if (PenNorm=="grouped") {
      appro <- as.vector(PenWeight * grouped(beta, PenA, control))
    } else {
      appro <- as.vector(PenWeight * get(PenNorm)(t(PenA)%*%beta, control))
    }
    crossprod(t(PenA)*sqrt(appro))
  }

# initialize
  start <- coefs <- coefold <- rep(1, sum(indexNrCoefs) + ifelse(intercept, 1, 0))
  ff <- as.formula(paste("y ~ ", paste(whichCoefs, collapse="+"), ifelse(intercept, "+1", "-1")))
  PP <- vector("list", p)
  names(PP) <- whichCoefs
  npen <- 0
  for (j in 1:length(indexPenA)) { npen <- if (is.matrix(indexPenA[[j]])) {npen + 1} else {npen + length(indexPenA[[j]])} }
  inout <- rep(0.01, npen)
  conv <- FALSE

# iterations
  if (is.null(tuning)) {
  for (i in 1L:control$maxi) {

    for (j in 1:p) {
    
         PenAjLength <- length(indexPenA[[j]])
         if (!is.matrix(indexPenA[[j]])) {
             l <- list()
             for (k in 1:PenAjLength) {
             l[[k]] <- A(start[vonBis[1,j]:vonBis[2,j]], indexPenNorm[j],
                           indexPenA[[j]][[k]], indexPenWeight[[j]][[k]], control) 
             }
             PP[[j]] <- l
         } else {

         PP[[j]] <- list(A(start[vonBis[1,j]:vonBis[2,j]], indexPenNorm[j],
                           indexPenA[[j]], indexPenWeight[[j]], control))

         }
    }

    mj     <- gam(formula=ff, family=family, data=data,
                 # weights=NULL, subset=NULL, ?!
                 na.action=na.omit, method=method, 
                 optimizer=optimizer,
                 #optimizer=c("perf"),
                 control=gamcontrol, scale=0,
                 select=FALSE, knots=NULL, sp=tuning, 
                 min.sp=NULL, H=NULL, gamma=1,
                 fit=TRUE, paraPen=PP, in.out=list(sp=inout, scale=1)
                 , start = start
                 )

    start  <- (1-control$nu)*coefold + control$nu*mj$coefficients
    inout <- mj$sp 
    
    # converged?!
    if (any(is.na(start))) {
      coefs <- coefold
      break
    }

    #print(paste(i, ": diff = ", sum(abs(start - coefold))/sum(abs(coefold)), sep=""))
    
    if (sum(abs(start - coefold))/sum(abs(coefold)) <= control$epsilon) {
      conv  <- TRUE
      coefs <- start
      break
    } else {
      coefs <- coefold <- start
    }

  }
  } else {
  for (i in 1L:control$maxi) {

    for (j in 1:p) {
    
         PenAjLength <- length(indexPenA[[j]])
         if (!is.matrix(indexPenA[[j]])) {
             l <- list()
             for (k in 1:PenAjLength) {
             l[[k]] <- A(start[vonBis[1,j]:vonBis[2,j]], indexPenNorm[j],
                           indexPenA[[j]][[k]], indexPenWeight[[j]][[k]], control) 
             }
             PP[[j]] <- l
         } else {

         PP[[j]] <- list(A(start[vonBis[1,j]:vonBis[2,j]], indexPenNorm[j],
                           indexPenA[[j]], indexPenWeight[[j]], control))

         }
    }

    mj     <- gam(formula=ff, family=family, data=data,
                 # weights=NULL, subset=NULL, ?!
                 na.action=na.omit, method=method, # offset=offset,
                 optimizer=optimizer,
                 #optimizer=c("perf"),
                 control=gamcontrol, scale=0,
                 select=FALSE, knots=NULL, sp=tuning, # NULL,
                 min.sp=NULL, H=NULL, gamma=1,
                 fit=TRUE, paraPen=PP
                 , start = start
                 )   

    start  <- (1-control$nu)*coefold + control$nu*mj$coefficients

    # converged?!
    if (any(is.na(start))) {
      coefs <- coefold
      break
    }

    if (sum(abs(start - coefold))/sum(abs(coefold)) <= control$epsilon) {
      conv  <- TRUE
      coefs <- start
      break
    } else {
      coefs <- coefold <- start
    }

  }
  }

  if (conv == FALSE)
    warning("The algorithm did not converge.")
    
# return
  mj$converged <- conv
  mj$iter <- i
  mj

}



