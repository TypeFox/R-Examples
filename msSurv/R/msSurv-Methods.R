####################################################################
## Accessor Functions
####################################################################

## NOTE - created an accessor function FOR EACH SLOT.

####################################################################
## State Information
####################################################################

## graphNEL object detailing MS model
setGeneric("tree",function(object) standardGeneric("tree"))
setMethod("tree",signature(object="msSurv"),
          function(object) return(object@tree))

## Number of states
setGeneric("ns",function(object) standardGeneric("ns"))
setMethod("ns",signature(object="msSurv"),
          function(object) return(object@ns))

## Event times
setGeneric("et",function(object) standardGeneric("et"))
setMethod("et",signature(object="msSurv"),
          function(object) return(object@et))

## Possible transitions
setGeneric("pos.trans",function(object) standardGeneric("pos.trans"))
setMethod("pos.trans",signature(object="msSurv"),
          function(object) return(object@pos.trans))

## Number of terminal states
setGeneric("nt.states",function(object) standardGeneric("nt.states"))
setMethod("nt.states",signature(object="msSurv"),
          function(object) return(object@nt.states))

####################################################################
## Counting Process Information
####################################################################

## Counting processes for event times
setGeneric("dNs",function(object) standardGeneric("dNs"))
setMethod("dNs",signature(object="msSurv"),
          function(object) return(object@dNs))

## At risk sets for event times
setGeneric("Ys",function(object) standardGeneric("Ys"))
setMethod("Ys",signature(object="msSurv"),
          function(object) return(object@Ys))

## Counting process for total transitions out of each state, at each time
setGeneric("sum_dNs",function(object) standardGeneric("sum_dNs"))
setMethod("sum_dNs",signature(object="msSurv"),
          function(object) return(object@sum_dNs))

## Datta-Satten counting processes for event times
setGeneric("dNs.K",function(object) standardGeneric("dNs.K"))
setMethod("dNs.K",signature(object="msSurv"),
          function(object) return(object@dNs.K))

## Datta-Satten risk sets for event times
setGeneric("Ys.K",function(object) standardGeneric("Ys.K"))
setMethod("Ys.K",signature(object="msSurv"),
          function(object) return(object@Ys.K))

## D-S counting process for total transitions out of each state, at each time
setGeneric("sum_dNs.K",function(object) standardGeneric("sum_dNs.K"))
setMethod("sum_dNs.K",signature(object="msSurv"),
          function(object) return(object@sum_dNs.K))

####################################################################
## State occupation probabilities and A-J estimators
####################################################################


## State occupation probabilities
setGeneric("ps",function(object) standardGeneric("ps"))
setMethod("ps",signature(object="msSurv"),
          function(object) return(object@ps))

## Aalen-Johansen estimates of transition probabilities at each event time
setGeneric("AJs",function(object) standardGeneric("AJs"))
setMethod("AJs",signature(object="msSurv"),
          function(object) return(object@AJs))

## Variance of state occupation probabilities
setGeneric("var.sop",function(object) standardGeneric("var.sop"))
setMethod("var.sop",signature(object="msSurv"),
          function(object) return(object@var.sop))

## Variance of A-J estimates
setGeneric("cov.AJs",function(object) standardGeneric("cov.AJs"))
setMethod("cov.AJs",signature(object="msSurv"),
          function(object) return(object@cov.AJs))

## I+dA = integrand used to calculate A-J estimate, where dA is the differential
##        for the N-A estimate of the cumulative hazard
setGeneric("I.dA",function(object) standardGeneric("I.dA"))
setMethod("I.dA",signature(object="msSurv"),
          function(object) return(object@I.dA))

## Variance of differential of N-A estimate
setGeneric("cov.dA",function(object) standardGeneric("cov.dA"))
setMethod("cov.dA",signature(object="msSurv"),
          function(object) return(object@cov.dA))

####################################################################
## State entry and exit distributions
####################################################################

## Normalized entry distribution
setGeneric("Fnorm",function(object) standardGeneric("Fnorm"))
setMethod("Fnorm",signature(object="msSurv"),
          function(object) return(object@Fnorm))

## Unconditional (subdistribution) entry distribution
setGeneric("Fsub",function(object) standardGeneric("Fsub"))
setMethod("Fsub",signature(object="msSurv"),
          function(object) return(object@Fsub))

## Normalized exit distribution
setGeneric("Gsub",function(object) standardGeneric("Gsub"))
setMethod("Gsub",signature(object="msSurv"),
          function(object) return(object@Gsub))

## Unconditional (subdistribution) exit distribution
setGeneric("Gnorm",function(object) standardGeneric("Gnorm"))
setMethod("Gnorm",signature(object="msSurv"),
          function(object) return(object@Gnorm))

## Accessor for VARIANCES

## Variance for normalized entry distribution
setGeneric("Fnorm.var",function(object) standardGeneric("Fnorm.var"))
setMethod("Fnorm.var",signature(object="msSurv"),
          function(object) return(object@Fnorm.var))

## Variance for unconditional (subdistribution) entry distribution
setGeneric("Fsub.var",function(object) standardGeneric("Fsub.var"))
setMethod("Fsub.var",signature(object="msSurv"),
          function(object) return(object@Fsub.var))

## Variance for normalized exit distribution
setGeneric("Gnorm.var",function(object) standardGeneric("Gnorm.var"))
setMethod("Gnorm.var",signature(object="msSurv"),
          function(object) return(object@Gnorm.var))

## Variance for unconditional exit distribution
setGeneric("Gsub.var",function(object) standardGeneric("Gsub.var"))
setMethod("Gsub.var",signature(object="msSurv"),
          function(object) return(object@Gsub.var))


####################################################################
## Print Method
####################################################################

setMethod("print", signature(x="msSurv"),
          function(x) {

              ## nonterminal states
              transient <- as.character(which(sapply(edges(tree(x)), function(x) length(x) > 0)))
              ## absorbing states
              absorb <- as.character(which(sapply(edges(tree(x)), function(x) length(x) == 0)))
              trans <- strsplit(colnames(dNs(x))," ")
              st.from <- sapply(trans, function(x) x[2])
              st.to <- sapply(trans, function(x) x[3])
              trans <- paste(st.from, "->", st.to, sep="", collapse=", ")

              cat("An object of class 'msSurv'.")
              cat("\n\n")

              cat("Allowable States: ")
              cat(nodes(tree(x)))
              cat("\n\n")

              cat(paste("The multistate model has", length(transient),
                        "transient state(s) and", length(absorb), "absorbing state(s).\n\n", sep = " "))

              cat("Number of timepoints with transitions: ")
              cat(length(et(x)))
              cat("\n\n")

              cat("Allowable transitions: ")
              cat(trans)
              cat("\n\n")

              cat("For a listing of available slots and methods see 'class?msSurv'\n\n")
              cat("See also: \n 'SOPt' for calculation of state occupation probabilities \n 'Pst' for calculation of transition probabilities and \n 'EntryExit' for calculation of state entry / exit distributions. \n")

          })

####################################################################
## Show Method
####################################################################

setMethod("show", signature(object="msSurv"),
          function(object) {

              ## nonterminal states
              transient <- as.character(which(sapply(edges(tree(object)), function(x) length(x) > 0)))
              ## absorbing states
              absorb <- as.character(which(sapply(edges(tree(object)), function(x) length(x) == 0)))
              trans <- strsplit(colnames(dNs(object))," ")
              st.from <- sapply(trans, function(x) x[2])
              st.to <- sapply(trans, function(x) x[3])
              trans <- paste(st.from, "->", st.to, sep="", collapse=", ")

              cat("An object of class 'msSurv'.")
              cat("\n\n")

              cat("Allowable States: ")
              cat(nodes(tree(object)))
              cat("\n\n")

              cat(paste("The multistate model has", length(transient),
                        "transient state(s) and", length(absorb), "absorbing state(s).\n\n", sep = " "))

              cat("Number of timepoints with transitions: ")
              cat(length(et(object)))
              cat("\n\n")

              cat("Allowable transitions: ")
              cat(trans)
              cat("\n\n")

              cat("For a listing of available slots and methods see 'class?msSurv'\n\n")

              cat("See also: \n 'SOPt' for calculation of state occupation probabilities \n 'Pst' for calculation of transition probabilities and \n 'EntryExit' for calculation of state entry / exit distributions. \n")

          })


####################################################################
## Summary Method
####################################################################


setMethod("summary", "msSurv",
          function(object, digits=3, all = FALSE, times = NULL,
                   ci.fun = "linear", ci.level = 0.95,
                   stateocc=TRUE, trans.pr=TRUE, dist=TRUE, DS=FALSE) {

              if (ci.level <= 0 | ci.level > 1) {
                  stop ("confidence level must be between 0 and 1")
              }

              ## print message if variances not available (need BS)
              if (is.null(cov.AJs(object)) & is.null(var.sop(object)) & is.null(Fnorm.var(object))) {
                  cat("\nNeed bootstrap to calculate variances \n\n")
              } else if (is.null(var.sop(object)) & is.null(cov.AJs(object))) {
                  cat("\nNeed bootstrap to calcualte variance of state occupation probabilities\n")
                  cat("and transition probabilities \n\n")
              } else if (is.null(var.sop(object)) & !is.null(cov.AJs(object))) {
                  cat("\nNeed bootstrap to calculate variance of state occupation probabilities \n\n")
                  CIs <- MSM.CIs(object, ci.level=0.95, trans=TRUE, sop=FALSE)
              } else {
                  CIs <- MSM.CIs(object, ci.level=0.95, trans=TRUE, sop=TRUE)
              }

              etimes <- et(object)

              if (!all){
                  if (is.null(times)) {
                      dt <- quantile(etimes, probs = c(0,0.25,0.5,0.75,1))
                      ind <- findInterval(dt,etimes) ## indicator of which information to print
                  } else {
                      if (any(times < 0)) stop("'times' must all be > 0")
                      ind <- findInterval(times, etimes)
                  }
              } else {
                  ind <- 1:length(etimes)
              }

              ## ##################################################################
              ## State Occupation Probability Section
              ## ##################################################################

              if (stateocc){

                  cat("State Occupation Information:", "\n", "\n")
                  for(i in seq(ns(object))) {
                      cat(paste("State ", nodes(tree(object))[i], "\n"))
                      if (!is.null(var.sop(object))) {
                          sop.sum <- data.frame(time=etimes[ind],
                                                estimate=ps(object)[ind,i],
                                                variance=var.sop(object)[ind,i],
                                                lower.ci=CIs$CI.p[ind,2,i],
                                                upper.ci=CIs$CI.p[ind,3,i])
                      } else {
                          sop.sum <- data.frame(time=etimes[ind],
                                                estimate=ps(object)[ind,i])
                      }
                      print(sop.sum,row.names=FALSE,digits=digits)
                      cat("\n")
                  } ## end of for statement
              } ## end of if (stateocc)


              ## ##################################################################
              ## State Entry/Exit Time Distributions
              ## ##################################################################

              if (dist){

                  ## Check whether entry/exit distributions available
                  if (is.null(Fsub(object))) {
                      cat("\nEntry distributions not available \n\n")
                  } else {
                      cat("State Entry Distribution Information:", "\n", "\n")
                      F.states <- sapply(strsplit(colnames(Fsub(object)), " "), function(x) x[2])
                      for(i in seq(ncol(Fsub(object)))) {
                          cat(paste("State ", F.states[i], "\n"))
                          dist.sum <- data.frame(time=etimes[ind],
                                                 entry.sub=Fsub(object)[ind,i],
                                                 entry.norm=Fnorm(object)[ind,i])
                          print(dist.sum,row.names=FALSE,digits=digits)
                          cat("\n")
                      } ## end of for statement
                  }
                  if (is.null(Gsub(object))) {
                      cat("\nExit distributions not available \n\n")
                  } else {
                      cat("State Exit Distribution Information:", "\n", "\n")
                      G.states <- sapply(strsplit(colnames(Gsub(object)), " "), function(x) x[2])
                      for(i in seq(ncol(Gsub(object)))) {
                          cat(paste("State ", G.states[i], "\n"))
                          dist.sum <- data.frame(time=etimes[ind],
                                                 exit.sub=Gsub(object)[ind,i],
                                                 exit.norm=Gnorm(object)[ind,i])
                          print(dist.sum,row.names=FALSE,digits=digits)
                          cat("\n")
                      } ## end of for statement
                  }

              } ## end of if (dist)

              ## ##################################################################
              ## Transition Probability Matrix Section
              ## ##################################################################

              if (trans.pr){

                  cat("Transition Probability Information:", "\n", "\n")

                  ## browser()
                  lt <- length(pos.trans(object))
                  tts <- strsplit(pos.trans(object), split = " ")
                  for (i in seq_along(pos.trans(object))) {

                      cat(paste("Transition", tts[[i]][1], "->", tts[[i]][2], "\n", sep = " "))

                      ## code to add number at risk and number transitions
                      ## print # events for transitions out of stage, else print # left
                      dns.name <- ifelse(tts[[i]][1] == tts[[i]][2],
                                         paste("dN", tts[[i]][1], ".", sep=" "),
                                         paste("dN", pos.trans(object)[i], sep=" "))

                      ifelse(dns.name %in% colnames(dNs(object)),
                             n.event <- dNs(object)[, dns.name],
                             n.event <- sum_dNs(object)[, dns.name])

                      if (DS==TRUE) {
                          ifelse(dns.name %in% colnames(dNs.K(object)),
                                 n.event.K <- dNs.K(object)[, dns.name],
                                 n.event.K <- sum_dNs.K(object)[, dns.name])
                      }

                      ys.name <- paste("y", tts[[i]][1], sep=" ")
                      n.risk <- Ys(object)[ ,ys.name]
                      if (DS==TRUE) {
                          n.risk.K <- Ys.K(object)[ ,ys.name]
                      }

                      if (dns.name %in% colnames(dNs(object))) {
                          if (!is.null(cov.AJs(object))) {
                              if (DS==TRUE) {
                                  tp.sum <- data.frame(time = etimes[ind], estimate = CIs$CI.trans[ind,1,i],
                                                       variance = CIs$CI.trans[ind,4,i],
                                                       lower.ci = CIs$CI.trans[ind,2,i],
                                                       upper.ci = CIs$CI.trans[ind,3,i],
                                                       n.risk = n.risk[ind], n.event = n.event[ind],
                                                       n.risk.K = n.risk.K[ind], n.event.K = n.event.K[ind])
                              } else {
                                  tp.sum <- data.frame(time = etimes[ind], estimate = CIs$CI.trans[ind,1,i],
                                                       variance = CIs$CI.trans[ind,4,i],
                                                       lower.ci = CIs$CI.trans[ind,2,i],
                                                       upper.ci = CIs$CI.trans[ind,3,i],
                                                       n.risk = n.risk[ind], n.event = n.event[ind])
                              }
                          } else {
                              if (DS==TRUE) {
                                  tp.sum <- data.frame(time = etimes[ind],
                                                       estimate = AJs(object)[tts[[i]][1], tts[[i]][2], ind],
                                                       n.risk  =  n.risk[ind], n.event = n.event[ind],
                                                       n.risk.K = n.risk.K[ind], n.event.K = n.event.K[ind])
                              } else {
                                  tp.sum <- data.frame(time = etimes[ind],
                                                       estimate = AJs(object)[tts[[i]][1], tts[[i]][2], ind],
                                                       n.risk  =  n.risk[ind], n.event = n.event[ind])
                              }
                          }

                      } else {
                          if (!is.null(cov.AJs(object))) {
                              if (DS==TRUE) {
                                  tp.sum <- data.frame(time = etimes[ind], estimate = CIs$CI.trans[ind,1,i],
                                                       variance = CIs$CI.trans[ind,4,i],
                                                       lower.ci = CIs$CI.trans[ind,2,i],
                                                       upper.ci = CIs$CI.trans[ind,3,i],
                                                       n.risk = n.risk[ind],
                                                       n.remain = n.risk[ind]-n.event[ind],
                                                       n.risk.K = n.risk.K[ind],
                                                       n.remain.K = n.risk.K[ind]-n.event.K[ind])
                              } else {
                                  tp.sum <- data.frame(time = etimes[ind], estimate = CIs$CI.trans[ind,1,i],
                                                       variance = CIs$CI.trans[ind,4,i],
                                                       lower.ci = CIs$CI.trans[ind,2,i],
                                                       upper.ci = CIs$CI.trans[ind,3,i],
                                                       n.risk = n.risk[ind],
                                                       n.remain = n.risk[ind]-n.event[ind])
                              }

                          }  else {
                              if (DS==TRUE) {
                                  tp.sum <- data.frame(time = etimes[ind],
                                                       estimate = AJs(object)[tts[[i]][1], tts[[i]][2], ind],
                                                       n.risk = n.risk[ind],
                                                       n.remain = n.risk[ind]-n.event[ind],
                                                       n.risk.K = n.risk.K[ind],
                                                       n.remain.K = n.risk.K[ind]-n.event.K[ind])
                              } else {
                                  tp.sum <- data.frame(time = etimes[ind],
                                                       estimate = AJs(object)[tts[[i]][1], tts[[i]][2], ind],
                                                       n.risk = n.risk[ind],
                                                       n.remain = n.risk[ind]-n.event[ind])
                              }
                          }
                      }
                      print(tp.sum, row.names = FALSE, digits = digits)
                      cat("\n")
                  } ## end of for statement
              } ## end of if trans.pr
          } ## end of summary function
          ) ## ends setMethod


####################################################################
##                  Plot Method                                   ##
####################################################################


setMethod("plot", signature(x="msSurv", y="missing"),
          function (x, states="ALL", trans="ALL", plot.type="stateocc",
                    CI=TRUE, ci.level=0.95, ci.trans="linear", ...) {


              plot.type <- match.arg(plot.type, c("stateocc", "transprob","entry.sub",
                                                  "entry.norm","exit.sub","exit.norm"))


              ####################################################################
              ## State occupation probabilities plot
              ####################################################################

              if (plot.type=="stateocc") {

                  if (states[1]=="ALL") states <- nodes(tree(x))

                  ## f.st <- factor(states)
                  ## ls <- length(states)
                  sl <- which(nodes(tree(x))%in%states) ## location of states in the matrix

                  if (CI==TRUE & !is.null(var.sop(x))) {
                      CIs <- MSM.CIs(x, ci.level=0.95) ## Calling CIs
                      rd <- CIs$CI.p
                      dimnames(rd)$state <- gsub("p", "State", dimnames(rd)$state)
                      y <- as.vector(rd[,1,sl])
                      y2 <- as.vector(rd[,2,sl]) ## lower limit
                      y3 <- as.vector(rd[,3,sl]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(rd)$state[sl], each=dim(rd)[1]))
                      st.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                        type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      st.plot <- update(st.plot, main="Plot of State Occupation Probabilites",
                                        xlab="Event Times", ylab="State Occupation Probabilities",
                                        key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                        text=list(c("Est", "Lower CI", "Upper CI")),
                                        columns=3))
                      print(st.plot)
                  } else {
                      if (CI==TRUE)
                          cat("Warning: 'var.sop'  is NULL and therefore CIs not plotted. \n")
                      ## rd <- CIs$CI.p
                      state.names <- gsub("p", "State", colnames(ps(x)))
                      y <- as.vector(ps(x)[,sl])
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(state.names[sl], each=length(et(x))))
                      st.plot <- xyplot(y ~ xvals | f.st, type="s",col=1, ...)
                      st.plot <- update(st.plot, main="Plot of State Occupation Probabilites",
                                        xlab="Event Times", ylab="State Occupation Probabilities",
                                        key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                      print(st.plot)
                  } ## end of no CIs

              } ## end of state occ plot

              ####################################################################
              ## Transition probabilities plot
              ####################################################################

              if (plot.type=="transprob") {

                  all.trans <- pos.trans(x)
                  if (trans[1] =="ALL") trans <- all.trans
                  tr <- which(all.trans%in%trans) ## location of states in the matrix

                  if (CI==TRUE & !is.null(cov.AJs(x))) {

                      CIs <- MSM.CIs(x, ci.level, ci.trans, sop=FALSE) ## Calling CIs
                      rd <- CIs$CI.trans
                      dimnames(rd)$trans <- gsub(" ", " -> ", dimnames(rd)$trans)
                      y <- as.vector(rd[,1,tr])
                      y2 <- as.vector(rd[,2,tr]) ## lower limit
                      y3 <- as.vector(rd[,3,tr]) ## upper limit
                      xvals <- rep(et(x), length(tr))
                      f.tp <- as.factor(rep(dimnames(rd)$trans[tr], each=dim(rd)[1]))
                      tr.plot <- xyplot(y + y2 + y3 ~ xvals | f.tp, allow.multiple=TRUE,
                                        type="s", lty=c(1,2,2), col=c(1,2,2),...)
                      tr.plot <- update(tr.plot, main="Plot of Transition Probabilites",
                                        xlab="Event Times", ylab="Transition Probabilites",
                                        key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                        text=list(c("Est", "Lower CI", "Upper CI")),
                                        columns=3))
                      print(tr.plot)

                  } else {
                      if (CI==TRUE)
                          cat("Warning: 'cov.AJs'  is NULL and therefore CIs not plotted. \n")
                      ## NOTE - copied from CIs function and modified accordingly
                      rd <- array(dim=c(length(et(x)), 1, length(trans)),
                                        dimnames = list(rows=et(x), cols=c("est"),
                                        trans=trans))
                      for (j in 1:length(trans)) { ## loop through possible transitions
                          idx <- unlist(strsplit(trans[j], " "))
                          rd[ , 1, j] <- AJs(x)[idx[1], idx[2] , ]
                      }
                      dimnames(rd)$trans <- gsub(" ", " -> ", dimnames(rd)$trans)
                      y <- as.vector(rd[,1,])
                      xvals <- rep(et(x), length(trans))
                      f.tp <- as.factor(rep(dimnames(rd)$trans, each=length(et(x))))
                      tr.plot <- xyplot(y ~ xvals | f.tp, type="s", lty=1, col=1, ...)
                      tr.plot <- update(tr.plot, main="Plot of Transition Probabilites",
                                        xlab="Event Times", ylab="Transition Probabilities",
                                        key = list(lines=list(col=1, lty=1), text=list("Est"), columns=1))
                      print(tr.plot)

                  }

              } ## end of 'transprob' plot


              ## ##################################################################
              ## Entry subdistribution function
              ## ##################################################################

              if (plot.type=="entry.sub") {

                  ## Check whether entry/exit distributions available at all
                  if (is.null(Fsub(x))) {
                      stop("\nEntry distributions not available \n\n")
                  }

                  ## Check whether entry/exit distributions available for requested state(s)
                  ## states for which entry dist is available
                  F.states <- sapply(strsplit(colnames(Fsub(x)), " "), function(x) x[2])
                  if (states[1]=="ALL") states <- F.states
                  if (any(!(states %in% F.states))) {
                      not.avail <- states[!(states %in% F.states)]
                      stop("\nState(s) ", not.avail, " do not have entry distributions available \n\n")
                  }


                  if (CI==TRUE & !is.null(Fsub.var(x))) {

                      ## Calling CIs for subdistribution
                      CIs <- Dist.CIs(x, ci.level, ci.trans, type = plot.type, states)
                      dimnames(CIs)$states <- paste("State", states)
                      y <- as.vector(CIs[,1,])
                      y2 <- as.vector(CIs[,2,]) ## lower limit
                      y3 <- as.vector(CIs[,3,]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(CIs)$states, each=length(et(x))))
                      ent.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                         type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      ent.plot <- update(ent.plot, main="Plot of State Entry Time Subdistributions",
                                         xlab="Event Times", ylab="State Entry Time Subdistributions",
                                         key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                             text=list(c("Est", "Lower CI", "Upper CI")),
                                             columns=3))
                      print(ent.plot)
                  }  else {
                      if (CI==TRUE)
                          cat("Warning: 'Fsub.var'  is NULL and therefore CIs not plotted. \n")
                      ## Here select appropriate columns ...
                      Dist <- Fsub(x)[,paste("F", states),drop = FALSE]
                      dimnames(Dist)[[2]]=gsub("F", "State", dimnames(Dist)[[2]])
                      y <- as.vector(Dist)
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(Dist)[[2]], each=length(et(x))))
                      ent.plot <- xyplot(y ~ xvals | f.st, type="s", col=1, ...)
                      ent.plot <- update(ent.plot, main="Plot of State Entry Time Subdistributions",
                                         xlab="Event Times", ylab="State Entry Time Subdistributions",
                                         key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                      print(ent.plot)
                  }
              } ## end of entry subdistribution plot



              ####################################################################
              ## State entry distribution (normalized)
              ####################################################################

              if (plot.type=="entry.norm") {


                  ## Check whether entry/exit distributions available at all
                  if (is.null(Fnorm(x))) {
                      stop("\nNormalized entry distributions not available \n\n")
                  }

                  ## Check whether entry/exit distributions available for requested state(s)
                  ## states for which entry dist is available
                  F.states <- sapply(strsplit(colnames(Fnorm(x)), " "), function(x) x[2])
                  if (states[1]=="ALL") states <- F.states
                  if (any(!(states %in% F.states))) {
                      not.avail <- states[!(states %in% F.states)]
                      stop("\nState(s) ", not.avail, " do not have normalized entry distributions available \n\n")
                  }


                  if (CI==TRUE & !is.null(Fnorm.var(x))) {

                      ## Calling CIs for subdistribution
                      CIs <- Dist.CIs(x, ci.level, ci.trans, type = plot.type, states)
                      dimnames(CIs)$states <- paste("State", states)
                      y <- as.vector(CIs[,1,])
                      y2 <- as.vector(CIs[,2,]) ## lower limit
                      y3 <- as.vector(CIs[,3,]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(CIs)$states, each=length(et(x))))
                      ent.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                         type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      ent.plot <- update(ent.plot, main="Plot of Normalized State Entry Time Distributions",
                                         xlab="Event Times", ylab="Normalized State Entry Time Distributions",
                                         key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                             text=list(c("Est", "Lower CI", "Upper CI")),
                                             columns=3))
                      print(ent.plot)
                  }  else {
                      if (CI==TRUE)
                          cat("Warning: 'Fnorm.var'  is NULL and therefore CIs not plotted. \n")
                      ## Here select appropriate columns ...
                      Dist <- Fnorm(x)[,paste("F", states),drop = FALSE]
                      dimnames(Dist)[[2]]=gsub("F", "State", dimnames(Dist)[[2]])
                      y <- as.vector(Dist)
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(Dist)[[2]], each=length(et(x))))
                      ent.plot <- xyplot(y ~ xvals | f.st, type="s", col=1, ...)
                      ent.plot <- update(ent.plot, main="Plot of Normalized State Entry Time Distributions",
                                         xlab="Event Times", ylab="Normalized State Entry Time Distributions",
                                         key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                      print(ent.plot)
                  }
              } ## end of normalized entry subdistribution plot


              ####################################################################
              ## State exit subdistribution
              ####################################################################

              if (plot.type=="exit.sub") {

                  ## Check whether entry/exit distributions available at all
                  if (is.null(Gsub(x))) {
                      stop("\nExit distributions not available \n\n")
                  }

                  ## Check whether entry/exit distributions available for requested state(s)
                  ## states for which exit dist is available
                  G.states <- sapply(strsplit(colnames(Gsub(x)), " "), function(x) x[2])
                  if (states[1]=="ALL") states <- G.states
                  if (any(!(states %in% G.states))) {
                      not.avail <- states[!(states %in% G.states)]
                      stop("\nState(s) ", not.avail, " do not have exit distributions available \n\n")
                  }


                  if (CI==TRUE & !is.null(Gsub.var(x))) {

                      ## Calling CIs for subdistribution
                      CIs <- Dist.CIs(x, ci.level, ci.trans, type = plot.type, states)
                      dimnames(CIs)$states <- paste("State", states)
                      y <- as.vector(CIs[,1,])
                      y2 <- as.vector(CIs[,2,]) ## lower limit
                      y3 <- as.vector(CIs[,3,]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(CIs)$states, each=length(et(x))))
                      ent.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                         type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      ent.plot <- update(ent.plot, main="Plot of State Exit Time Subdistributions",
                                         xlab="Event Times", ylab="State Exit Time Subdistributions",
                                         key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                             text=list(c("Est", "Lower CI", "Upper CI")),
                                             columns=3))
                      print(ent.plot)
                  }  else {
                      if (CI==TRUE)
                          cat("Warning: 'Gsub.var'  is NULL and therefore CIs not plotted. \n")
                      ## Here select appropriate columns ...
                      Dist <- Gsub(x)[,paste("G", states),drop = FALSE]
                      dimnames(Dist)[[2]]=gsub("G", "State", dimnames(Dist)[[2]])
                      y <- as.vector(Dist)
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(Dist)[[2]], each=length(et(x))))
                      ent.plot <- xyplot(y ~ xvals | f.st, type="s", col=1, ...)
                      ent.plot <- update(ent.plot, main="Plot of State Exit Time Subdistributions",
                                         xlab="Event Times", ylab="State Exit Time Subdistributions",
                                         key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                      print(ent.plot)
                  }

              } ## end of exit.sub


              ####################################################################
              ## State exit distribution (normalized)
              ####################################################################

              if (plot.type=="exit.norm") {

                  ## Check whether entry/exit distributions available at all
                  if (is.null(Gnorm(x))) {
                      stop("\nNormalized exit distributions not available \n\n")
                  }

                  ## Check whether entry/exit distributions available for requested state(s)
                  ## states for which exit dist is available
                  G.states <- sapply(strsplit(colnames(Gnorm(x)), " "), function(x) x[2])
                  if (states[1]=="ALL") states <- G.states
                  if (any(!(states %in% G.states))) {
                      not.avail <- states[!(states %in% G.states)]
                      stop("\nState(s) ", not.avail, " do not have normalized exit distributions available \n\n")
                  }


                  if (CI==TRUE & !is.null(Gnorm.var(x))) {

                      ## Calling CIs for subdistribution
                      CIs <- Dist.CIs(x, ci.level, ci.trans, type = plot.type, states)
                      dimnames(CIs)$states <- paste("State", states)
                      y <- as.vector(CIs[,1,])
                      y2 <- as.vector(CIs[,2,]) ## lower limit
                      y3 <- as.vector(CIs[,3,]) ## upper limit
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(CIs)$states, each=length(et(x))))
                      ent.plot <- xyplot(y + y2 + y3 ~ xvals | f.st, allow.multiple=TRUE,
                                         type="s", lty=c(1,2,2), col=c(1,2,2), ...)
                      ent.plot <- update(ent.plot, main="Plot of Normalized State Exit Time Distributions",
                                         xlab="Event Times", ylab="Normalized State Exit Time Distributions",
                                         key = list(lines=list(col=c(1, 2, 2), lty=c(1, 2, 2)),
                                             text=list(c("Est", "Lower CI", "Upper CI")),
                                             columns=3))
                      print(ent.plot)
                  }  else {
                      if (CI==TRUE)
                          cat("Warning: 'Gnorm.var'  is NULL and therefore CIs not plotted. \n")
                      ## Here select appropriate columns ...
                      Dist <- Gnorm(x)[,paste("G", states),drop = FALSE]
                      dimnames(Dist)[[2]]=gsub("G", "State", dimnames(Dist)[[2]])
                      y <- as.vector(Dist)
                      xvals <- rep(et(x), length(states))
                      f.st <- as.factor(rep(dimnames(Dist)[[2]], each=length(et(x))))
                      ent.plot <- xyplot(y ~ xvals | f.st, type="s", col=1, ...)
                      ent.plot <- update(ent.plot, main="Plot of Normalized State Exit Time Distributions",
                                         xlab="Event Times", ylab="Normalized State Exit Time Distributions",
                                         key = list(lines=list(col=c(1), lty=c(1)), text=list(c("Est"))))
                      print(ent.plot)
                  }

              } ## end of exit.norm

          } ## end of function

          )
