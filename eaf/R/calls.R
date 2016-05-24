###############################################################################
#
#                          Copyright (c) 2011
#         Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
#             Marco Chiarandini <marco@imada.sdu.dk>
#
# This program is free software (software libre); you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, you can obtain a copy of the GNU
# General Public License at http://www.gnu.org/copyleft/gpl.html
#
# IMPORTANT NOTE: Please be aware that the fact that this program is
# released as Free Software does not excuse you from scientific
# propriety, which obligates you to give appropriate credit! If you
# write a scientific paper describing research that made substantive use
# of this program, it is your obligation as a scientist to (a) mention
# the fashion in which this software was used in the Methods section;
# (b) mention the algorithm in the References section. The appropriate
# citation is:
# 
#    Manuel Lopez-Ibanez, Luis Paquete, and Thomas Stuetzle.
#    Exploratory Analysis of Stochastic Local Search Algorithms in
#    Biobjective Optimization. In T. Bartz-Beielstein, M. Chiarandini,
#    L. Paquete, and M. Preuss, editors, Experimental Methods for the
#    Analysis of Optimization Algorithms, pages 209-222. Springer,
#    Berlin, Germany, 2010. doi: 10.1007/978-3-642-02538-9_9
# 
# Moreover, as a personal note, I would appreciate it if you would email
# manuel.lopez-ibanez@ulb.ac.be with citations of papers referencing this
# work so I can mention them to my funding agent and tenure committee.
#
################################################################################
#
# TODO:
#
#  * Follow this style for coding:
#    http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html
#
################################################################################
## for development only
##require(lattice)
#setwd("~/svn/Code/EAF")
#require(modeltools)
#source("R/common.R")
#source("R/EAF.R")
#dyn.load(file.path(".", paste("src/Reaf", .Platform$dynlib.ext, sep="")))
##paste("src/eaf",.Platform$dynlib.ext,sep=""))


### the top-level convenience function
eafplot <- function(x,...)
{
  UseMethod("eafplot");
}


eafplot.formula <-
  function(x, data = data.frame(), groups = NULL, subset = NULL,
           percentiles=NULL,include.extremes=FALSE, ...)
{
  ## formula of type time+cost~run|inst, groups=alg
  if (missing(x))
    stop("formula missing")

  formula <- x

  if ((length(formula) != 3L) || (length(formula[[2L]]) !=3L))
    stop("incorrect specification for 'formula'")

  mf <- modeltools::ModelEnvFormula(formula = formula, data = data,
                                    subset = subset, designMatrix = FALSE,
                                    responseMatrix = FALSE, ...)

  ### extract data from the ModelEnv object
  points <- mf@get("response")
  sets <- mf@get("input")
  if (length(mf@formula)==3)
    cond<-as.list(mf@get("blocks"))
  else
    cond<-NULL
  groups <- eval(substitute(groups), data, environment(formula))
  #subset <- eval(substitute(subset), data, environment(formula))

  if (!is.null(groups) && !is.factor(groups))
    stop("groups must be a factor")
  
  if (length(cond) == 0)
    {
      strip <- FALSE
      cond <- list(gl(1, length(points)))
    }

  condlevels <- lapply(cond, levels)
  cond.max.level <- unlist(lapply(cond, nlevels))

  npackets <- prod(cond.max.level)

  panel.args <- vector(mode = "list", length = npackets)

  packet.sizes <- numeric(npackets)
  if (npackets > 1)
    {
      dim(packet.sizes) <- sapply(condlevels, length)
      dimnames(packet.sizes) <- lapply(condlevels, as.character)
    }

  cond.current.level <- rep(1, length(cond))

  for (packet.number in seq_len(npackets)) {
    id <- .compute.packet(cond, cond.current.level)
    packet.sizes[packet.number] <- sum(id)

    panel.args[[packet.number]] <-
      list(points = as.matrix(points[id,]), sets = as.numeric(sets[id,]),
           groups=groups[id])
    ## MARCO: I do not think we need to care about subscripts... or do we?
    #if (subscripts)
    #   panel.args[[packet.number]]$subscripts <-
    #                   subscr[id]

    cond.current.level <-
      .cupdate(cond.current.level,
               cond.max.level)
  }
  # FIXME: I don't think this is doing the right thing.
  op <- par(mfrow=.check.layout(NULL,cond.max.level)[2:3])
  for (i in seq_len(length(panel.args))) {
    eafplot.default(panel.args[[i]]$points,
                    panel.args[[i]]$sets,
                    panel.args[[i]]$groups,
                    percentiles, ...)
  }
  par(op)
  invisible()
}



eafplot.list<-function(x,...) {
  if (!is.list(x))
    stop("'x' must be a list of data.frames with exactly three columns")

  DT<-data.frame()
  if (!is.null(names(x)))
    groups<-names(x)
  else
    groups<-c(1:length(x))
  for (i in seq_len(length(x)) ) {
    if (!is.data.frame(x[[i]]))
      stop("Each element of the list must be a data.frame with exactly three columns.")
    DT<-rbind(DT,data.frame(x[[i]],groups=groups[i]))
  }
  eafplot(as.matrix(DT[,c(1,2)]),as.numeric(DT[,3]),as.factor(DT[,4]),...)
}


eafplot.data.frame <- function(x, y = NULL, ...)
{
  eafplot.data.frame2 <- function(x, groups, main = DNAME, ...)
    eafplot(as.matrix(x[,c(1,2)]), as.numeric(x[,3]),
            groups = groups, main = main, ...)
  
  if (!is.data.frame(x))
    stop("'x' must be a data.frame with exactly three columns.\n",
         "  If you have grouping and conditioning variables, please consider using this format: 'eafplot(formula, data, ...)'")

  if (nrow(x) < 1L)
    stop("not enough (finite) 'x' observations")
  
  if (!is.null(y)) {
    if (!is.data.frame(y))
      stop("'y' must be a data.frame with exactly three columns")
    if (nrow(y) < 1L)
      stop("not enough (finite) 'y' observations")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    DT <- rbind(data.frame(x, groups = deparse(substitute(x))),
                data.frame(y, groups = deparse(substitute(y))))
    groups <- as.factor(DT[, 4])
  } else {
    DNAME <- deparse(substitute(x))
    DT <- x
    groups <- NULL
  }
  eafplot.data.frame2(DT, groups = groups, ...)
}


### Local Variables:
### ess-style: DEFAULT
### indent-tabs-mode: nil
### ess-fancy-comments: nil
### End:

