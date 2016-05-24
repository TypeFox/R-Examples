# iota: Janson, H., & Olsson, U. (2001). A measure of agreement for interval or nominal multivariate observations.
#       Educational and Psychological Measurement, 61, 277-289.

iota <- function(ratings, scaledata = c("quantitative","nominal"), standardize = FALSE) {
  ratings   <- na.omit(ratings)
  scaledata <- match.arg(scaledata)
  if (is.factor(ratings[[1]][,1])) scaledata <- "nominal"
  
  ns <- nrow(ratings[[1]])
	nr <- ncol(ratings[[1]])
	nvar <- length(ratings)  # number of variables
  detail <- NULL

  if (scaledata=="quantitative") {
    if (standardize) {
      for (i in 1:nvar) {
        x <- as.numeric(ratings[[i]])
        ratings[[i]] <- (ratings[[i]]-mean(x))/sd(x)
      }
      detail <- "Variables have been z-standardized before the computation"
    }
    # Take original as new rating-structure
    ratinglist <- ratings
  }
  else if (scaledata=="nominal") {
    dummyn <- 1  # Number of dimensions in dummy list
    ratinglist <- list()

    for (i in 1:nvar) {
      # How many levels were used?
      lev <- character()
      for (j in 1:nr) {
        lev <- c(lev, levels(as.factor(ratings[[i]][,j])))
      }
      lev <- levels(factor(lev))

      # Build new rating-structure
      for (levnr in 1:length(lev)) {
        ratinglist[[dummyn]] <- matrix(0,nrow=ns,ncol=nr)
        ratinglist[[dummyn]][ratings[[i]]==lev[levnr]] <- 1/sqrt(2)
        dummyn <- dummyn + 1
      }
    }
  }

  # Compute coefficient
	nv <- length(ratinglist)  # number of variables

  doSS <- 0; deSS <- 0
  for (i in 1:nv) {
    SSt <- var(as.numeric(as.matrix(ratinglist[[i]])))*(ns*nr-1)
  	SSw <- sum(apply(ratinglist[[i]],1,var)/ns)*ns*(nr-1)
  	SSc <- var(apply(ratinglist[[i]],2,mean))*ns*(nr-1)

  	doSS <- doSS+SSw
  	deSS <- deSS+((nr-1)*SSt+SSc)
  }

  coeff  <- 1-(nr*doSS)/deSS

  method <- paste("iota for ",scaledata," data (",nvar,ifelse(nvar==1," variable)"," variables)"),sep="")
  rval <- structure(list(method = method,
                         subjects = ns, raters = nr,
                         irr.name = "iota", value = coeff,
                         detail = detail),
                    class="irrlist")
  return(rval)
}

