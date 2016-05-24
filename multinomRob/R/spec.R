#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1/
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#  $Id: spec.R,v 1.7 2005/09/27 08:04:06 wrm1 Exp $
#

# functions to interpret model formula specifications and build data to analyze

# get.xdata:
# Return model matrix corresponding to the formula in formul
# 
get.xdata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "term.labels"))==0 & attr(t1, "intercept")==0) {
    m <- NULL;  # no regressors specified for the model matrix
  }
  else {
    m <- model.matrix(formul, data=datafr);
  }
  return(m);
}

# get.ydata:
# Return response vector corresponding to the formula in formul
# 
get.ydata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "response"))==0) {
    m <- NULL;  # no response variable specified
  }
  else {
    m <- model.response(model.frame(formul, data=datafr));
  }
  return(m);
}

# get.xy:
#  Return response vectors and model matrices corresponding the formulas
#  in formlist, along with the variable names and the number of variables
#  in the model matrix for each response.  Each response is a distinct column
#  in returned matrix Y, and each model matrix is along a dimension of the
#  returned array X.  X contains 0 for columns that are not used for a
#  category of the response (corresponding to a column in Y).
#  Factor variables are expanded in X, as binary {0,1} dummy variables, reduced as
#  appropriate to match other factors and the intercept (e.g., with one factor
#  variable (facvar) and an intercept, there are length(levels(facvar))-1 dummy
#  variables, while with one factor variable and no intercept there are
#  length(levels(facvar)) dummy variables).
#  Missing data (NA) in any variable causes the entire observation to be deleted
#  (listwise deletion).
#  A matrix ypos is created to indicate which elements of the response matrix Y
#  originally had negative values:  ypos <- Y >= 0.  The ypos matrix is used to
#  indicate the number of outcome alternatives for each observation.  If there
#  are fewer than two alternatives for an observation, that observation is deleted.
#  dim(Y) == c(nobs,ncats)
#  dim(X) == c(nobs, max(c(1,xlengths)), ncats)
#  length(xlengths) == ncats
#  length(ynames) == ncats & is.character(ynames[i])
#  length(xnames) == ncats & is.list(xnames) & length(xnames[[i]]) == xlengths[i]
#     & is.character(xnames[[i]][j])
# 
get.xy <- function(formlist, datafr, print.level=0) {
  ncats <- length(formlist);
  nobs <- dim(datafr)[1];  # assume all variables are the same length
  ynames <- rep("", ncats);
  # find negative response values and missing response or regressor variable values
  # if y[j,i] is negative, ignore missing data in x[j,i]
  ypos <- matrix(FALSE, nobs, ncats);  #  in ypos, FALSE is y<0, TRUE is y>=0 or is.na(y)
  ymiss <- matrix(FALSE, nobs, ncats);  #  in ymiss, TRUE is NA, FALSE is not NA
  xmiss <- matrix(FALSE, nobs, ncats);  #  in xmiss, TRUE is NA, FALSE is not NA
  for (i in 1:ncats) {
    ti <- terms(formlist[[i]], data=datafr);
    ynames[i] <- as.character(attr(ti, "variables")[[1 + attr(ti, "response")]]) ;
    ypos[,i] <- datafr[[ ynames[i] ]] >= 0;
    ypos[,i] <- ifelse(is.na(ypos[,i]), TRUE, ypos[,i]);
    ymiss[,i] <- is.na(datafr[[ ynames[i] ]]);
    if (length(attr(ti, "term.labels")) > 0)  {
      xni <- as.character(attr(ti, "term.labels")) ;
      if (length(xni) >= 1) {
        for (ii in 1:length(xni)) {
          if (xni[ii] %in% names(datafr)) {  # skip interaction terms
            xmiss[,i] <- xmiss[,i] | is.na(datafr[[ xni[ii] ]]);
          }
        }
      }
    }
  }
  xmiss <- xmiss & ypos;
  # remove any observation for which less than two alternatives exist
  hastwo <- apply(ypos, 1, sum) >= 2;
  if (any(!hastwo)) {
    xmiss[!hastwo,] <- TRUE;
    if (print.level > 0) {
      print("some observations have fewer than two responses.");
    }
  }
  # create vector for listwise deletion
  if (print.level > 0) {
    if (any(ymiss)) {
      print("there is missing response variable data.");
    }
    if (any(xmiss)) {
      print("there is missing regressor variable data.");
    }
  }
  misslist <- apply(ymiss,1,any) | apply(xmiss,1,any);
  if (print.level > 0) {
    if (any(misslist)) {
      print("implementing listwise deletion for missing data and fewer than two responses.");
      cat("deleting observations:  ", c(1:nobs)[misslist], "\n");
    }
    # print(misslist)
  }
  # begin implementation of listwise deletion
  datafr <- datafr[!misslist,];
  ypos <- ypos[!misslist,];
#  nobs <- sum(!misslist);
  # end implementation of listwise deletion

  XYdata <- get.XYdata(formlist, datafr);
  
  dimnames(XYdata$Y) <- list(NULL, XYdata$ynames);
  dimnames(XYdata$X) <- list(NULL, NULL, XYdata$ynames);
  return(list(Y=XYdata$Y, X=XYdata$X, xnames=XYdata$xnames, ynames=XYdata$ynames,
              xlengths=XYdata$xlengths, ypos=ypos));
}

get.XYdata <- function(formlist, datafr) {
  ncats <- length(formlist);
  nobs <- dim(datafr)[1];  # assume all variables are the same length
  Y <- NULL;
  Xlist <- list();
  xlengths <- rep(0, ncats);
  ynames <- rep("", ncats);
  xnames <- list()
  for (i in 1:ncats) {
    ti <- terms(formlist[[i]], data=datafr);
    ynames[i] <- as.character(attr(ti, "variables")[[1 + attr(ti, "response")]]) ;
    Y <- cbind(Y, get.ydata(formlist[[i]], datafr));
    if (length(attr(ti, "term.labels"))==0 & attr(ti, "intercept")==0) {
      Xlist[[i]] <- 0;
      xlengths[i] <- 0;
      xnames[[i]] <- "";
    }
    else {
      Xlist[[i]] <- get.xdata(formlist[[i]], datafr);
      xlengths[i] <- dim(Xlist[[i]])[2];
      xnames[[i]] <- unlist(dimnames(Xlist[[i]])[2]);
    }
  }
  X <- array(0, dim=c(nobs, max(c(1,xlengths)), ncats));
  for (i in 1:ncats) {
    if (xlengths[i] > 0) {
      X[1:nobs, 1:xlengths[i], i] <- Xlist[[i]];
    }
  }
  return(list(Y=Y, X=X, xnames=xnames, ynames=ynames, xlengths=xlengths));
}

get.xynames <- function(formlist, datafr) {
  ncats <- length(formlist);
  ynames <- rep("", ncats);
  xnames <- list()
  for (i in 1:ncats) {
    ti <- terms(formlist[[i]], data=datafr);
    ynames[i] <- as.character(attr(ti, "variables")[[1 + attr(ti, "response")]]) ;
    if (length(attr(ti, "term.labels"))==0 & attr(ti, "intercept")==0) {
      xnames[[i]] <- "";
    }
    else {
      xnames[[i]] <- unlist(dimnames(get.xdata(formlist[[i]], datafr))[2]);
    }
  }
  return(list(xnames=xnames, ynames=ynames));
}
