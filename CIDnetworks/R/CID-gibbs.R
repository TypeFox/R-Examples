
#####################################################################################
#
# Gibbs Sampler collection for CID, given the collection of input terms.

#source("COV-reference.R"); source("SBM-reference.R"); source("LSM-reference.R"); source("SR-reference.R"); library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R");

.onAttach <- function (...) {
  packageStartupMessage("CIDnetworks v 0.8.0")
}

# All subclasses available to CID.

LSM <- function(...) LSMcid$new(...)
LVM <- function(...) LVMcid$new(...)
SBM <- function(...) SBMcid$new(...)
SR <- function(...) SRcid$new(...)
BETA <- function(...) BETAcid$new(...)


EdgeCOV <- function(..., cov.type="Edge") COVARIATEcid$new(..., cov.type=cov.type)
SenderCOV <- function(...) COVARIATEcid$new(..., cov.type="Sender")
ReceiverCOV <- function(...) COVARIATEcid$new(..., cov.type="Receiver")
SendRecCOV <- function(...) COVARIATEcid$new(..., cov.type="SendRec")
IdenticalCOV <- function(...) COVARIATEcid$new(..., cov.type="Identical")


MMSBM <- function(...) MMSBMcid$new(...)
HBM <- function(...) HBMcid$new(...)


unwrap.CID.Gibbs <- function (gibbs.out) list.output.to.matrices(gibbs.out)

#{  elements <- lapply (1:length(gibbs.out[[1]]), function(el) if (grepl("out", class(gibbs.out[[1]][[el]]))) {    r1 <- lapply(1:length(gibbs.out[[1]][[el]]), function (el2) sapply(gibbs.out, function(it) it[[el]][[el2]]))    names(r1) <- names(gibbs.out[[1]][[el]])    r1  } else sapply(gibbs.out, function(it) it[[el]]))  names(elements) <- names(gibbs.out[[1]])  elements}



CIDnetwork <-
  setRefClass (
    "CIDnetwork",
    fields = list(
      n.nodes="numeric",
      edge.list="matrix",
      edge.list.rows="list",
      outcome="numeric",
      node.names="character",

      is.directed="logical",
      class.outcome="character",   #"ordinal" or "gaussian" -- if "binary", it will make it "ordinal"
      ordinal.count="numeric",
      ordinal.cutoffs="numeric",   #current draw. Every edge has the same cutoff value with respect to the outcome.

      int.outcome="numeric",  #intermediate outcome, which is Gaussian. If outcome is binary, this is what is truncated.

      intercept="numeric",    #current draw
      intercept.m="numeric",  #prior mean
      intercept.v="numeric",  #prior variance

      residual.variance="numeric",   #current draw -- always 1 if binary/ordinal
      residual.variance.ab="numeric",  #inverse gamma prior parameters -- length 2

      reciprocity.present="logical",
      reciprocal.match="integer",   #which edge has the other as its counterpart? Use for reciprocity.
      int.correlation="numeric",
      int.correlation.ab="numeric",  #symmetric beta prior.

      ## robit.augment="numeric",

      log.likelihood="numeric",   #current log likelihood of data

      components="list",   # updated with every draw -- contains COV, LSM, etc.
      comp.values="matrix" # after each draw, updates the value as expressed of that object. "value" is the component of the linear predictor for that edge ; each row corresponds to each edge; each column corresponds to each component.
      ),

    methods=list(
      initialize = function (

        edge.list,    #specs: this should be numbered from 1 to the max edge number.
#        sociomatrix,

        edge.list.rows,
        n.nodes,      #this should come from the data.
        node.names=character(),

        intercept=0,
        intercept.m=0,
        intercept.v=1000000,

        residual.variance=1,
        residual.variance.ab=c(0.001, 0.001),

        outcome=numeric(0),
        generate=FALSE,

        is.directed,
        class.outcome="ordinal",
        ordinal.count=2,
        ordinal.cutoffs=sort(rexp(ordinal.count-2, rate=50)),

        int.correlation=0,          #For now, generate nothing.
        int.correlation.ab=c(1,1),  #Symmetric beta prior.
        include.reciprocity=FALSE,  #No reciprocity unless asked.

        reinit=TRUE,
        components=list(),
        verbose=2
        )
      {

        ## Initial checks. This is largely redundant right now, but that's OK.
        ## if (is.null(class.outcome)) class.outcome <- "ordinal"

        if (!(class.outcome %in% c("binary","gaussian","ordinal"))) stop (paste("class.outcome",class.outcome,"is not supported. Must be one of", paste(c("binary","gaussian","ordinal"), collapse=",")))

        if (int.correlation != 0) include.reciprocity <- TRUE

###        if (!missing(sociomatrix)) {
###          if (!(class(sociomatrix) %in% c("array", "matrix", "data.frame"))) stop ("Sociomatrix must be a matrix or two-dimensional array.")
###          if (nrow(sociomatrix) != ncol(sociomatrix)) stop ("Sociomatrix must be square.")
###          new.nodes <- nrow(sociomatrix)
###          .self$n.nodes <<- new.nodes
###
###          if (missing(is.directed)) {
###            if (all(t(sociomatrix) == sociomatrix, na.rm=TRUE)) {
###              message ("Auto-detected a symmetric sociomatrix; assuming this is an undirected network.")
###              .self$is.directed <<- FALSE
###            } else {
###              message ("Auto-detected an asymmetric sociomatrix; assuming this is a directed network.")
###              .self$is.directed <<- TRUE
###            }
###          } else .self$is.directed <<- is.directed
###
###
###          if (!.self$is.directed) {
###            .self$edge.list <<- make.edge.list (new.nodes)
###            .self$outcome <<- sociomatrix[l.diag(new.nodes)]
###          } else {
###            .self$edge.list <<- make.arc.list (new.nodes)
###            .self$outcome <<- t(sociomatrix)[non.diag(new.nodes)]   #note: arc list changes receiver first, so it's row-based.
###          }
###
###        } else {
        if (missing(edge.list)) {
          if (missing(n.nodes)) stop ("CIDnetwork: Not detected: n.nodes, edge.list, or sociomatrix") else {
            .self$n.nodes <<- n.nodes
            if (!missing(is.directed))
              .self$is.directed <<- is.directed
            else {
              .self$is.directed <<- include.reciprocity | (int.correlation != 0)
              if(verbose > 0){
                if (.self$is.directed)
                  message ("CIDnetwork: Reciprocity was specified: assuming a DIRECTED network.")
                else
                  message ("CIDnetwork: Directedness was unspecified: assuming an undirected network.")
              }
            }

            if (.self$is.directed) {
              .self$edge.list <<- make.arc.list (n.nodes)
              .self$outcome <<- outcome
            } else {
              .self$edge.list <<- make.edge.list (n.nodes)
              .self$outcome <<- outcome
            }
          }
        } else {
          if (missing(n.nodes))
            .self$n.nodes <<- max(c(edge.list))
          else
            .self$n.nodes <<- n.nodes

          if (missing(is.directed)) {

            e1 <- unique(edge.list)

            e1[e1[,1]>e1[,2],] <- e1[e1[,1]>e1[,2], 2:1]
            e2 <- unique(e1)
            if (nrow(e1) != nrow(e2)) {
              if(verbose > 0)
                message ("CIDnetwork: Reciprocal potential edges detected; assuming this is a directed network.")
              .self$is.directed <<- TRUE
            } else {
              if(verbose > 0)
                message ("CIDnetwork: No reciprocal potential edges detected; assuming this is an undirected network.")
              .self$is.directed <<- FALSE
            }
          }else{
            .self$is.directed <<- is.directed
          }

          .self$edge.list <<- edge.list
          .self$outcome <<- outcome
        }


        if (!missing(edge.list.rows))
          .self$edge.list.rows <<- edge.list.rows
        else
          .self$edge.list.rows <<- row.list.maker(.self$edge.list)

        if (length(node.names) != .self$n.nodes){
          if(verbose > 0)
            message("Warning:  Length of node.names differs from n.nodes.  Using default node names")
          .self$node.names <<- as.character(1:.self$n.nodes)
        }else{
          .self$node.names <<- node.names
        }
        if (class.outcome == "binary") {
          .self$class.outcome <<- "ordinal"
          .self$ordinal.count <<- 2
          .self$ordinal.cutoffs <<- numeric()
        } else {
          .self$class.outcome <<- class.outcome
          .self$ordinal.count <<- ordinal.count
          .self$ordinal.cutoffs <<- ordinal.cutoffs
        }

        .self$intercept <<- intercept
        .self$intercept.m <<- intercept.m
        .self$intercept.v <<- intercept.v

        if (class.outcome == "gaussian") {
          .self$residual.variance <<- residual.variance
        } else {
          .self$residual.variance <<- 1
        }
        .self$residual.variance.ab <<- residual.variance.ab


        ##reinitialize components, just in case.
        if (length(components)>0) {
          if (class(components) != "list") components.t <- list(components) else components.t <- components
          if(reinit){
            for (kk in 1:length(components.t)) {
              components.t[[kk]]$reinitialize(n.nodes=.self$n.nodes,
                                              edge.list=.self$edge.list,
                                              node.names=.self$node.names)

              if (class(components.t[[kk]]) %in% c("SBMcid", "MMSBMcid", "HBMcid")) {
                .self$intercept.m <<- 0
                .self$intercept.v <<- 0.000000000001
              }
            }
          }
          .self$components <<- components.t
#          if(length(which(sapply(components ,class) == "COVARIATEcid")) >1){
#            .self$components <- combine.covariates(components)
#          }
        } else .self$components <- list()

        ##message("Component initialization complete.")



        ##Let's do reciprocity now! ACT, 5-21-14
        .self$int.correlation <<- int.correlation
        .self$int.correlation.ab <<- int.correlation.ab
        if(verbose > 1) message("Reciprocity: ",include.reciprocity)

        reciprocity.present <<- include.reciprocity
        if (reciprocity.present) {
          ## find matches.
          reciprocal.match <<- apply(.self$edge.list, 1, function(ee) min(which (.self$edge.list[,1] == ee[2] & .self$edge.list[,2] == ee[1])))
        } else {
          reciprocal.match <<- rep(as.integer(NA), length(outcome))
        }
        ##int.correlation=0,          #For now, generate nothing.
        ##int.correlation.ab=c(1,1),  #Symmetric beta prior.
        ##include.reciprocity=FALSE,  #No reciprocity unless asked.



        if (generate) .self$generate()

        ## Control statements currently force .self$outcome to have already
        ## been set anyways.
        ##        else if (length(outcome) == nrow(.self$edge.list)) {
        ##          .self$outcome <<- outcome
        ##        }

        if (length(.self$outcome) != nrow(.self$edge.list))
          stop (paste0("In CIDnetwork$initialize: Outcome variable length: ",length(.self$outcome),"; Edges: ", nrow(.self$edge.list),". Did you forget to add a required variable?"))

        ## Check for empty categories, currently disallowed.
        if (class.outcome == "ordinal") {
          counts <- tabulate (.self$outcome + 1, .self$ordinal.count)
          if (any(counts == 0))
            stop ("The following ordinal categories have no outcomes: ", paste(which(counts == 0)-1, collapse=" "), "; required non-zero counts for ", paste (0:(.self$ordinal.count-1), collapse=" "))
        }

        ##message("Generation complete.")
        update.intermediate.outcome()
        ##message("CID Initialization complete.")

      },

      reinitialize = function (n.nodes=NULL, edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          edge.list.rows <<- row.list.maker(edge.list)
        }

        if (!is.null(node.names)) {
            if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)


        if (length(components)>0) for (kk in 1:length(components)) {
          components[[kk]]$reinitialize (.self$n.nodes, .self$edge.list, .self$node.names)
          if (class(components[[kk]]) %in% c("SBMcid", "MMSBMcid", "HBMcid")) {
            intercept.m <<- 0
            intercept.v <<- 0.000000000001
          }
        }


        if (reciprocity.present) {
          ## find matches.
          reciprocal.match <<- apply(edge.list, 1, function(ee) min(which (edge.list[,1] == ee[2] & edge.list[,2] == ee[1])))
        } else {
          reciprocal.match <<- rep(as.integer(NA), length(outcome))
        }

      },

      pieces = function (include.name=FALSE) {
        if (length(components)>0) {
          c(list(intercept=intercept,
                 residual.variance=residual.variance,
                 log.likelihood=log.likelihood,
                 ordinal.cutoffs=ordinal.cutoffs,
                 int.correlation=int.correlation),
            lapply(components, function(cc) cc$pieces(include.name)))
        } else {
          list(intercept=intercept,
               residual.variance=residual.variance,
               log.likelihood=log.likelihood,
               ordinal.cutoffs=ordinal.cutoffs,
               int.correlation=int.correlation)
        }
      },

      show = function () {
        message("CIDnetwork object properties:")

        message(paste("class.outcome:", class.outcome))
        if (class.outcome == "ordinal") message(paste("Ordinal groups:", ordinal.count))
        message(paste("Nodes:", n.nodes))
        if (is.directed) message(paste("Potential Directed Arcs:", nrow(edge.list))) else message(paste("Potential Undirected Edges:", nrow(edge.list)))


        message("Intercept: ", intercept)
        message("Variance: ", residual.variance)
        if (length(components)>0) for (kk in 1:length(components)) {
          message(class(components[[kk]])); components[[kk]]$show()
        }

      },
#      combine.covariates <- function(components.arg){
#        cov.ind <- which(sapply(components,class) == "COVARIATEcid")
#
#    },
      plot = function (coefs=coef.cov, names=1:length(coefs), sd=NULL, interval=NULL, ...) {
        if (length(components) > 0) for (cc in 1:length(components)) components[[cc]]$plot()
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },

      value = function (redo=FALSE) {
        if (redo) update.comp.values()
        rowSums(comp.values) + intercept
      },
      value.ext = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        if (length(components)>0){
            comp.values.here <- sapply (components, function(cc) cc$value())
        }else{
            comp.values.here <- matrix(0, nrow=nrow(edge.list))
        }
        (rowSums(comp.values.here) + parameters$intercept)[edges]
      },

      summary = function () show(),


      generate = function () {

        # Built 5-21-14. Tested: ----
        if (reciprocity.present) {
          int.deviation <- matrix(rnorm(2*nrow(edge.list), 0, sqrt(residual.variance)), ncol=2)
          int.deviation[,2] <- int.correlation*int.deviation[,1] + sqrt(1-int.correlation^2)*int.deviation[,2]

          swaps <- which(edge.list[,1]>edge.list[,2] & !is.na(reciprocal.match))

          int.deviation[swaps,1] <- int.deviation[reciprocal.match[swaps],2]
          int.outcome <<- value(redo=TRUE) + int.deviation[,1]
        } else {
          int.outcome <<- rnorm(nrow(edge.list),
                                value(redo=TRUE),
                                sqrt(residual.variance))
        }

        ## Generate user-observed outcome from intermediate outcome.
        if (class.outcome == "ordinal") {
          temp.out <- 0*int.outcome
          ordinal.steps <- c(0, ordinal.cutoffs)
          for (ii in 1:length(ordinal.steps)) temp.out <- temp.out + 1*(int.outcome > ordinal.steps[ii])
          outcome <<- temp.out
        }

        if (class.outcome == "gaussian") outcome <<- int.outcome

      },


      rem.values = function(kk) {if (kk>0) value() - comp.values[,kk] else value() - intercept},

      update.comp.values = function () {
        if (length(components)>0) comp.values <<- sapply (components, function(cc) cc$value()) else comp.values <<- matrix(0, nrow=nrow(edge.list))
      },



      ## Update Z_ij.
      update.intermediate.outcome = function () {
          value.hold <- value(redo=TRUE)
          if(class.outcome == "gaussian"){
              int.outcome <<- outcome
          }
          if(any(is.na(outcome))){
              io.temp <- int.outcome
              io.temp[is.na(outcome)] <- rnorm(sum(is.na(outcome)),
                                               value.hold,1)
              int.outcome <<- io.temp

          }
          if (class.outcome == "ordinal") {
              ##first, update the values themselves.
              io.temp <- int.outcome

              ##matched.value <- value.hold[reciprocal.match]
              lohi <- edge.list[,1] < edge.list[,2]

              breaker.lower <- c(-Inf, 0, ordinal.cutoffs)
              breaker.upper <- c(0, ordinal.cutoffs, Inf)

              breakers.lower <- breaker.lower[outcome + 1]
              breakers.upper <- breaker.upper[outcome + 1]


              ##are any going to be trouble?
              p.gen <- rep(0.5, length(outcome))
              for (ii in 1:ordinal.count) {
                  p.gen[!is.na(outcome) & outcome == ii-1] <-
                      pnorm(breaker.upper[ii], value.hold[!is.na(outcome) & outcome == ii-1], 1) - pnorm(breaker.lower[ii], value.hold[!is.na(outcome) & outcome == ii-1], 1)
              }
              p.gen[is.na(outcome)] <- 0

              ##for (ii in 1:ordinal.count) {
              ## unmatched edges.
              condition <- p.gen > 1e-10 & is.na(reciprocal.match)
              io.temp[condition] <-
                  rtnorm(sum(condition),
                         value.hold[condition], 1,
                         lower=breakers.lower[condition],
                         upper=breakers.upper[condition])

                                        #matched edges, (lower, higher).
              condition <- p.gen > 1e-10 & !is.na(reciprocal.match) & lohi
              io.temp[condition] <-
                  rtnorm(sum(condition),
                         value.hold[condition] + int.correlation*(io.temp[reciprocal.match[condition]] - value.hold[reciprocal.match[condition]]),
                         1 - int.correlation^2,
                         lower=breakers.lower[condition],
                         upper=breakers.upper[condition])

                                        #matched edges, (higher, lower).
              condition <- p.gen > 1e-10 & !is.na(reciprocal.match) & !lohi
              io.temp[condition] <-
                  rtnorm(sum(condition),
                         value.hold[condition] + int.correlation*(io.temp[reciprocal.match[condition]] - value.hold[reciprocal.match[condition]]),
                         1 - int.correlation^2,
                         lower=breakers.lower[condition],
                         upper=breakers.upper[condition])

              p.gen[is.na(outcome)] <- 0.5
              io.temp[p.gen <= 1e-10] <- (breakers.lower[p.gen <= 1e-10]+breakers.upper[p.gen <= 1e-10])/2
              io.temp[p.gen <= 1e-10 & outcome == 0 ] <- (breakers.upper[outcome == 0 & p.gen <= 1e-10])
              io.temp[(p.gen <= 1e-10) & (outcome == ordinal.count - 1)] <- (breakers.lower[outcome == ordinal.count - 1 & p.gen <= 1e-10])

              int.outcome <<- io.temp

              ##now, change the cutoff values, which lie between the Z values for each one. Assume a flat prior for now.
              if (length(ordinal.cutoffs)>0) for (kk in 1:length(ordinal.cutoffs)) {
                  effective.range <- c(max(int.outcome[outcome <= kk]), min(int.outcome[outcome >= kk+1]))
                  if (is.na(effective.range[1])) effective.range[1] <- 0
                  if (is.na(effective.range[2])) effective.range[1] <- 10000
                  ordinal.cutoffs[kk] <<- runif(1, effective.range[1], effective.range[2])
              }

          }
      },



      ## Simple slice sampler for autocorrelation term \rho. Pretty optimized! Act, 6-3-14
      int.correlation.prior = function (pp = int.correlation) {
        dbeta((1+pp)/2,  int.correlation.ab[1], int.correlation.ab[2], log=TRUE) - log(2)
      },

      #Note: this is being drawn from the main likelihood, not the intermediate likelihood, because it's more efficient for the whole chain.
      draw.int.correlation = function () {

        these.pieces <- pieces()
        current.value <- int.correlation

        #First: draw uniform at current value and establish the vertical cut.
        cut.1 <- log(runif(1)) +
          log.likelihood.by.value(value.ext(these.pieces), these.pieces) + #, use.intermediate=TRUE) +  BD 12/31
              int.correlation.prior(current.value)
        limits <- c(-0.9999,0.9999)

        emergency.count <- 0
        repeat {
          prop.value <- runif(1, limits[1], limits[2])   #Draw proposal value between the bounds.
          these.pieces$int.correlation <- prop.value
          if (log.likelihood.by.value(value.ext(these.pieces), these.pieces) +  #, use.intermediate=TRUE) + BD
              int.correlation.prior(prop.value) > cut.1) {   #is the current density above the threshold?
            int.correlation <<- prop.value; break    #keep it and save.
          } else limits[1 + 1*(prop.value > current.value)] <- prop.value    #Trim down.

          emergency.count <- emergency.count+1; if (emergency.count > 1000) stop ("In draw.int.correlation, way too many narrows-down were made. This is either at a peak or it's broken.")
        }

      },

      #What does the grid say?
      test.int.correlation = function (rho=seq(-19,19)/20) {
        these.pieces <- pieces()
        sapply (rho, function (rr) {
          these.pieces$int.correlation <- rr
          log.likelihood.by.value(value.ext(these.pieces), these.pieces) + #, use.intermediate=TRUE) + BD 12/31
              int.correlation.prior(rr)
        })
      },




      draw.intercept = function (verbose=2) {

        outcomeresid <- int.outcome - rem.values(0);
        varpiece <- solve(nrow(edge.list)/residual.variance + 1/intercept.v)
        meanpiece <- varpiece*(sum(outcomeresid)/residual.variance + intercept.m/intercept.v)
        if (verbose > 2) message ("Intercept ",meanpiece," ",sqrt(varpiece))
        intercept <<- rnorm(1, meanpiece, sqrt(varpiece))

      },

      ## Note: getting this to work for 2D pmvnorm. ACT, 5-23-14
      log.likelihood.by.value = function (value.this=value(),
        pieces.this=pieces(),
        sumup=TRUE,
        use.intermediate=FALSE)
      {
        cm <- function(pp) matrix(c(1,pp,pp,1), nrow=2)
        output <- 0*int.outcome
        if (class.outcome == "gaussian" | use.intermediate) {
          outcomeresid <- int.outcome - value.this
          if (reciprocity.present) {
            picks <- which(is.na(reciprocal.match))
            output[picks] <- dnorm(outcomeresid[picks], 0, sqrt(pieces.this$residual.variance), log=TRUE)

            picks2 <- which (!is.na(reciprocal.match) & edge.list[,1] < edge.list[,2])
            output[picks2] <- dmvnorm(cbind(outcomeresid[picks2], outcomeresid[reciprocal.match[picks2]]),
                                      rep(0, 2),
                                      pieces.this$residual.variance*cm(pieces.this$int.correlation),
                                      log=TRUE)
          } else output <- dnorm(outcomeresid, 0, sqrt(pieces.this$residual.variance), log=TRUE)

          if (sumup) output <- sum(output[!is.na(outcome)])
        } else if (class.outcome == "ordinal") {

          breaker.lower <- c(-Inf, 0, pieces.this$ordinal.cutoffs)
          breaker.upper <- c(0, pieces.this$ordinal.cutoffs, Inf)

          ##  Getting ranges for non-NA outcomes
          breakers.lower <- breaker.lower[outcome[!is.na(outcome)] + 1]
          breakers.upper <- breaker.upper[outcome[!is.na(outcome)] + 1]
          value.this <- value.this[!is.na(outcome)]

          ##  BD Added 12/31
          reciprocal.match.na <- reciprocal.match[!is.na(outcome)]

          # first: unmatched edges.
          if (reciprocity.present) {
            picks <- which(is.na(reciprocal.match.na))
            output[picks] <- log(
              pnorm(breakers.upper[picks],
                    value.this[picks],
                    sqrt(pieces.this$residual.variance)) -
              pnorm(breakers.lower[picks],
                    value.this[picks],
                    sqrt(pieces.this$residual.variance)))

            ## now, matched edges.  BD added !is.na(outcome)
            picks2 <- which (!is.na(reciprocal.match.na) & edge.list[!is.na(outcome),1] < edge.list[!is.na(outcome),2])

            ## ACT was working on this. See basefunctions line 51
            output[picks2] <- my.pmvnorm(lower=matrix(breakers.lower[c(picks2, reciprocal.match.na[picks2])], ncol=2),
                                         upper=matrix(breakers.upper[c(picks2, reciprocal.match.na[picks2])], ncol=2),
                                         meanval=matrix(value.this[c(picks2, reciprocal.match.na[picks2])], ncol=2),
                                         sigma=pieces.this$residual.variance,
                                         rho=pieces.this$int.correlation)

            #this.output <- sapply(picks2, function (pp) {
            #  pmvnorm (breakers.lower[c(pp, reciprocal.match.na[pp])],
            #          breakers.upper[c(pp, reciprocal.match.na[pp])],
            #           mean=value.this[c(pp, reciprocal.match.na[pp])],
            #           sigma=pieces.this$residual.variance*matrix(c(1, pieces.this$int.correlation, pieces.this$int.correlation, 1), nrow=2))})
            #output[picks2] <- log(this.output)

          } else {

            output <- log(
              pnorm(breakers.upper,
                    value.this,
                    sqrt(pieces.this$residual.variance)) -
              pnorm(breakers.lower,
                    value.this,
                    sqrt(pieces.this$residual.variance)))
          }
          if(sumup) output <- sum(output)
        }

        return(output)
      },

      update.log.likelihood = function () {
        log.likelihood <<- log.likelihood.by.value ()
      },

      draw.variance = function (verbose=2) {
        outcomeresid <- int.outcome - value();

        if (verbose > 2) message ("Variance: ",nrow(edge.list), " ", sum(outcomeresid^2))
        residual.variance <<-
          1/rgamma(1,
                   residual.variance.ab[1] + nrow(edge.list)/2,
                   residual.variance.ab[2] + sum(outcomeresid^2)/2)

      },

      draw = function (verbose=FALSE) {

        if (reciprocity.present) draw.int.correlation()
        if (class.outcome != "gaussian") update.intermediate.outcome()

        if (length(components)>0) for (kk in 1:length(components)) {
          update.comp.values()
          components[[kk]]$outcome <<- .self$int.outcome - rem.values(kk)
          components[[kk]]$residual.variance <<- residual.variance
          components[[kk]]$draw()

          if (exists("shift", components[[kk]])) {   ### - where does this currently arise?
            intercept <<- intercept + components[[kk]]$shift
            components[[kk]]$shift <<- 0
          }
        }

        #variance and intercept.
        update.comp.values()
        draw.intercept(verbose)

        if (class.outcome == "gaussian") {
          update.comp.values()
          draw.variance(verbose)
          update.comp.values()
        }

        ## update correlation!!! ACT, 5-23-14

        update.log.likelihood()

      },

      random.start = function () {
        intercept <<- rnorm (1, 0, 1)
#        draw.intercept()
        if (length(components)>0) for (kk in 1:length(components)) components[[kk]]$random.start()
        if (class.outcome == "gaussian") draw.variance()
        if (reciprocity.present) draw.int.correlation()

        update.log.likelihood()

      },

      gibbs.full = function (
        report.interval=100, draws=100,
        burnin=-1, thin=10,   #auto-cut.
        auto.burn.cut=TRUE,
        auto.burn.count=200,
        make.random.start=TRUE,
        auto.converge=FALSE,
        extend.max = 10,
        extend.count = 100,
        verbose = 2) {

        #if (auto.burn.cut) {thin <- 10; burnin <- 0}
        if (thin < 1) stop ("thin must be an integer greater than zero.")
        if (burnin >= 0) {
          if(verbose > 1)
            message ("Overriding auto burn-in and cut. Initially discarding ",burnin," iterations.")
          auto.burn.cut <- FALSE
        } else if (burnin < 0 & !auto.burn.cut) {
          if(verbose > 0)
            message ("Negative burn-in inputted -- defaulting to auto-burn-in.")
          auto.burn.cut <- TRUE
        }


        out <- list()
        if (make.random.start) random.start()

        if (auto.burn.cut) {
          if(verbose > 1)
            message ("CID Auto-Burning In")
          repeat {
            loglikes <- rep(NA, auto.burn.count)
            theseruns <- 1:auto.burn.count
            time.one <- proc.time()[3]
            for (kk in theseruns) {
              draw()
              loglikes[kk] <- log.likelihood
            }
            time.two <- proc.time()[3] - time.one

          #  obj1 <- lm(loglikes ~ theseruns)
          #  slopes <- summary(obj1)$coef[2, c(1,4)]
          #  if (slopes[1]<0 | slopes[2] > 0.05) break else message ("CID Still Auto-Burning In")
          #  if (any(is.na(loglikes))) {print (loglikes)} else {
            slopes <- sum(loglikes*(theseruns - mean(theseruns)))
            if (is.na(slopes) | slopes >= 0){
              if(verbose > 1) message ("CID Still Auto-Burning In")
            }else ##if (slopes<0)
              break
            ##else message ("CID Still Auto-Burning In")
            ##}
          }

          if(verbose > 0)
            message ("CID Auto-Burned In. Estimated Seconds Remaining: ",
                     round(time.two/auto.burn.count*thin*draws))
          burnin <- 0
        }


        if(burnin > 0) time.one <- proc.time()[3]
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {

            out[[index]] <- pieces(include.name=TRUE)
            if (report.interval > 0 & index %% report.interval == 0) {
              if(verbose > 1) message("CID ",index)
            }
          } else if (round(index)==index) {
            if (report.interval > 0) if (index %% report.interval == 0){
              if(index != 0){
                if(verbose > 1) message("CID burnin ",index)
              }else{
                time.two <- proc.time()[3] - time.one
                if(verbose > 1) message ("CID Burned In. Estimated Seconds Remaining: ", round((time.two/burnin)*thin*draws))
              }
            }
          }
        }

        if(auto.converge){
          extend.iter <- 0
          extend.count <- min(extend.count,draws)

          while(!convergence.test(out) & extend.iter < extend.max) {
            if(verbose > 0)
              message("CID Convergence not detected.  Extending chain")
            extend.iter <- extend.iter + 1

            ##  Shifting Chain
            if(draws > extend.count){
              for(ii in 1:(draws - extend.count)){
                out[[ii]] <- out[[ii + extend.count]]
              }
            }

                ##  Getting new draws
            for(kk in 1:(extend.count*thin)){
              draw();
              index <- draws - extend.count + kk/thin
              if(round(index) == index){
                out[[index]] <- pieces(include.name=TRUE)
              }
            }
          }

          if(!convergence.test(out)){
            if(verbose > 0)
              message("CID Convergence not detected after maximum extensions")
          }else{
            if(verbose > 0)
              message("CID Converged")
          }
        }

        return(out)
      },

        convergence.test = function(gibbs.out){
            log.lik <- unlist(gibbs.switcheroo(gibbs.out)$log.lik)
            draws <- length(log.lik)
            chain.1 <- log.lik[1:(draws/10)]
            chain.2 <- log.lik[(draws/2 + 1):draws]
            sd.est <- sqrt(var(chain.1)/length(chain.1) +
                           var(chain.2)/length(chain.2))

            return((mean(chain.2) - mean(chain.1)) <= qnorm(0.999) * sd.est)
        },

      #Removes magic number later. This should be the length of (intercept, residual, loglik, cutoffs, int.correlation).
      non.comp.count = function () 5,


      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),


      gibbs.mean = function(gibbs.out){

        switched <- gibbs.switcheroo(gibbs.out)
        intercept.mean <- mean(unlist(switched$intercept))
        components.mean <- NULL

        if(length(components) > 0){
          for(cc in 1:length(components)){
            components.mean <- c(components.mean,components[[cc]]$gibbs.mean(switched[[cc+non.comp.count()]]))
          }
        }
        if(ordinal.count > 2){
          ordinal.cuts.mean <- mean(unlist(switched$ordinal.cuts))
        }else{
          ordinal.cuts.mean <- ordinal.cutoffs
        }

        mean.object <- CID(edge.list,outcome,
                           intercept=intercept.mean,
                           components=components.mean,
                           intercept.m=intercept.m,intercept.v=intercept.v,
                           residual.variance=residual.variance,
                           residual.variance.ab=residual.variance.ab,
                           #ordinal.count=ordinal.count,
                           ordinal.cutoffs=ordinal.cuts.mean,
                           int.correlation=int.correlation,
                           int.correlation.ab=int.correlation.ab,
                           include.reciprocity=reciprocity.present,
                           reinit=FALSE,
                           verbose=0)
        mean.object$update.log.likelihood()
        mean.object$update.comp.values()
        return(mean.object)


      },


      gibbs.switcheroo = function (gibbs.out) {   #returns the chain for each component/parameter.
        out <- lapply(1:length(gibbs.out[[1]]), function(el)
                      lapply(1:length(gibbs.out), function(el2) gibbs.out[[el2]][[el]]))
        out.names <- c("intercept","res.variance","log.lik",
                        "ordinal.cuts","reciprocity")
        s.total <- length(out)
        if(s.total > non.comp.count()){
            sapply(out,FUN=function(x)return(class(x[[1]])))
            comp.classes <- sapply(out,FUN=function(x) return(class(x[[1]])))
            out.names <- c(out.names,
                           comp.classes[(non.comp.count()+1):s.total])
        }
        names(out) <- out.names
        return(out)
      },




#      gibbs.summary = function (gibbs.out) {
#        switched <- gibbs.switcheroo (gibbs.out)
#        s.sum <- function (int1) c(min=min(int1), max=max(int1), mean=mean(int1), sd=sd(int1), quantile(int1, c(0.025, 0.975)))
#        out <- list()

#        out$intercept <- s.sum(unlist(switched[[1]]))
#        message ("Intercept:"); print(out$intercept)

#        out$residual.variance <- s.sum(unlist(switched[[2]]))
#        if (out$residual.variance[4] != 0) {  #SD
#          message ("Residual variance:")
#          print (out$residual.variance)
#        }

#        out$log.likelihood <- s.sum(unlist(switched[[3]]))
#        message ("Log likelihood:"); print(out$log.likelihood)

#        if (class.outcome == "ordinal") if (ordinal.count > 2) {
#          o1 <- matrix(unlist(switched[[4]]), nrow=ordinal.count-2)
#          out$ordinal.cutoffs <- t(apply(o1, 1, s.sum))
#          rownames(out$ordinal.cutoffs) <- paste("Cutoff", 1:(ordinal.count-2), 2:(ordinal.count-1), sep="-")
#          message ("Ordinal cutoffs:")
#          print (out$ordinal.cutoffs)
#        }

#        if (length(components) > 0) for (cc in 1:length(components)) {
#          message ("Component ", class(components[[cc]]),":")
#          out[[cc+non.comp.count()]] <- components[[cc]]$print.gibbs.summary(switched[[cc+non.comp.count()]])
#          names(out)[cc+non.comp.count()] <- class(components[[cc]])
#        }

#        out <- structure(out,class="summary.CID.Gibbs")
#        return(invisible(out))
#      },

    gibbs.summary = function (gibbs.out) {
        switched <- gibbs.switcheroo (gibbs.out$results)
        s.sum <- function (int1) {c(min=round(min(int1),3), max=round(max(int1),3), estimated.mean= round(mean(int1),3), estimated.sd=round(sd(int1),3), round(quantile(int1, c(0.025, 0.975)),3))
}
        out <- list()

        out$intercept <- s.sum(unlist(switched[[1]]))

        out$residual.variance <- s.sum(unlist(switched[[2]]))

        out$log.likelihood <- s.sum(unlist(switched[[3]]))

        if (class.outcome == "ordinal") if (ordinal.count > 2) {
          o1 <- matrix(unlist(switched[[4]]), nrow=ordinal.count-2)
          out$ordinal.cutoffs <- t(apply(o1, 1, s.sum))
          rownames(out$ordinal.cutoffs) <- paste("Cutoff", 1:(ordinal.count-2), 2:(ordinal.count-1), sep="-")
        }else{
            out$ordinal.cutoffs <- NULL
        }

        if (length(components) > 0) for (cc in 1:length(components)) {
            ncc <- non.comp.count()
            out[[cc+ncc]] <- components[[cc]]$gibbs.summary(switched[[cc+ncc]])
            names(out)[cc+non.comp.count()] <- class(components[[cc]])
        }

        out$CID.object <- gibbs.out$CID.object
        out <- structure(out,class="summary.CID.Gibbs")
        return(out)
      },

#      print.gibbs.summary = function (gibbs.out) {
#        switched <- gibbs.switcheroo (gibbs.out)
#        s.sum <- function (int1) {c(min=round(min(int1),3), max=round(max(int1),3), estimated.mean = round(mean(int1)), estimated.sd=round(sd(int1),3), round(quantile(int1, c(0.025, 0.975)),3))}
#        out <- list()
#
#        out$intercept <- s.sum(unlist(switched[[1]]))
#        message ("Intercept:"); print(out$intercept)
#
#        out$residual.variance <- s.sum(unlist(switched[[2]]))
#        if (out$residual.variance[4] != 0) {  #SD
#          message ("Residual variance:")
#          print (out$residual.variance)
#        }
#
#        out$log.likelihood <- s.sum(unlist(switched[[3]]))
#        message ("Log likelihood:"); print(out$log.likelihood)
#
#        if (class.outcome == "ordinal") if (ordinal.count > 2) {
#          o1 <- matrix(unlist(switched[[4]]), nrow=ordinal.count-2)
#          out$ordinal.cutoffs <- t(apply(o1, 1, s.sum))
#          rownames(out$ordinal.cutoffs) <- paste("Cutoff", 1:(ordinal.count-2#), 2:(ordinal.count-1), sep="-")
#          message ("Ordinal cutoffs:")
#          print (out$ordinal.cutoffs)
#        }
#
#        if (length(components) > 0) for (cc in 1:length(components)) {
#          message ("Component ", class(components[[cc]]),":")
#          out[[cc+non.comp.count()]] <- components[[cc]]$print.gibbs.summary(#switched[[cc+non.comp.count()]])
#          names(out)[cc+non.comp.count()] <- class(components[[cc]])
#        }
#
#        out <- structure(out,class="summary.CID.Gibbs")
#        return(invisible(out))
#      },


        print.gibbs.summary = function(gibbs.sum){

            message ("Intercept:"); print(gibbs.sum$intercept[-(1:2)])
            if (gibbs.sum$residual.variance[4] != 0) {  #SD
		message ("Residual variance:")
		print (gibbs.sum$residual.variance)
            }

            message ("Log likelihood:"); print(gibbs.sum$log.likelihood[-(1:2)])

            if(!is.null(gibbs.sum$residual.variance[-(1:2)])){
              if(gibbs.sum$residual.variance["estimated.sd"] != 0){
                message ("Residual variance:")
                print (gibbs.sum$residual.variance)
              }
            }

            if (class.outcome == "ordinal") if (ordinal.count > 2) {
        	message ("Ordinal cutoffs:")
        	print (gibbs.sum$ordinal.cutoffs)
            }

            if(!is.null(gibbs.sum$ordinal.cutoffs)){
                message ("Ordinal cutoffs:")
                print (gibbs.sum$ordinal.cutoffs)
            }

            if (length(components) > 0){
              cov.count <- 0
              for (cc in 1:length(components)) {
                ##		message ("Component ", class(components[[cc]]),":")
                ix <- which(names(gibbs.sum) == class(components[[cc]]))
                if(length(ix) > 1){
                  cov.count <- cov.count + 1
                  ix <- ix[cov.count]
                }
                for(ii in 1:length(ix)){
                  components[[cc]]$print.gibbs.summary(gibbs.sum[[ix]])
                }

#		print(round(gibbs.sum$class(components[[cc]])),3)
              }
            }
            return()#invisible(gibbs.sum))
        },

      gibbs.plot = function (gibbs.out, DIC=NULL, which.plots=1:(length(components)+5), auto.layout=TRUE) {

        switched <- gibbs.switcheroo (gibbs.out)

        if (auto.layout) {
          if (intercept.v < 0.0001) which.plots <- which.plots[which.plots != 1]
          if (class.outcome != "gaussian") which.plots <- which.plots[which.plots != 2]
          if (ordinal.count == 2) which.plots <- which.plots[which.plots != 4]
          if (!reciprocity.present) which.plots <- which.plots[which.plots != 5]
          cols <- ceiling(length(which.plots)/2)
          par(mfrow=c(2, cols))
        }

        ## 1: Grand Intercept
        if (1 %in% which.plots & intercept.v >= 0.0001) plot.default (unlist(switched[[1]]),
                                   main="Grand Intercept",
                                   xlab="Iteration",
                                   ylab="Intercept",
                                   type="l")

        ## 2: Residual variance. Skipped if not gaussian.
        if (2 %in% which.plots & class.outcome=="gaussian")
          plot.default (unlist(switched[[2]]), main="Residual Variance",
                                              xlab="Iteration",
                        ylab="Residual Variance")
        main.label <- "Log-likelihood"; if (!is.null(DIC)) {
          main.label <- paste0(main.label, ": DIC = ",signif(DIC[1], 5))
          if (length(DIC)>=2) main.label <- paste0(main.label, "\nDeviance of Average = ",signif(DIC[2], 5))
          if (length(DIC)>=3) main.label <- paste0(main.label, "\nEffective Parameter Count = ",signif(DIC[3], 5))
        }

        ## 3: Log-likelihood
        if (3 %in% which.plots) {
          plot.default (unlist(switched[[3]]), main=main.label,
                        xlab="Iteration",ylab="Log Likelihood",
                        type="l")
        }

        ## 4: Ordinal cutoffs.
        if (4 %in% which.plots & class.outcome=="ordinal" & ordinal.count>2) {
          draws <- length(unlist(switched[[1]]))
          xx <- sort(rep(1:draws, ordinal.count-2))
          plot.default (c(1,xx), c(0,unlist(switched[[4]])), col=c(0, rep(1:(ordinal.count-2), draws)),
                        main="Ordinal Cutoff Values",
                        xlab="Iteration",
                        ylab="Cutoff Value")
          abline(h=0, col=8)
        }

        ## 5: Reciprocity
        if (5 %in% which.plots & reciprocity.present) {
          main.label <- "Reciprocal Correlation"
          plot.default (unlist(switched[[5]]), main=main.label,
                        xlab="Iteration",
                        ylab="Reciprocal Correlation")


        }

        #6 to (5 + length(components)): components!
        if (length(components) > 0) for (cc in 1:length(components)) if ((non.comp.count()+cc) %in% which.plots) components[[cc]]$gibbs.plot(switched[[cc+non.comp.count()]])

      },


      DIC = function (gibbs.out, add.parts=FALSE) {
        #all.values <- gibbs.value(gibbs.out)
        #deviance.of.average <- -2*log.likelihood.by.value (apply(all.values, 1, mean))
        deviance.of.average <- -2 * gibbs.mean(gibbs.out)$log.likelihood
        average.deviance <- -2 * mean(unlist(gibbs.switcheroo(gibbs.out)$log.lik))
                                 #2*apply(all.values, 2, log.likelihood.by.value))
        output <- c(DIC=2*average.deviance - deviance.of.average)
        if (add.parts) output <- c(output,
                                   deviance.of.average=deviance.of.average,
                                   effective.parameters=average.deviance - deviance.of.average,
                                   average.deviance=average.deviance)
        return(output)
      },

      marginal.loglikelihood = function (gibbs.out) {
        all.values <- gibbs.value(gibbs.out)
        model.log.likelihoods <- apply(all.values, 2, log.likelihood.by.value)
        1/mean(1/model.log.likelihoods)
      },

      pseudo.CV.loglikelihood = function (gibbs.out) {
        all.values <- gibbs.value(gibbs.out)
        model.log.likelihoods <- apply(all.values, 2, log.likelihood.by.value, sumup=FALSE)
        each.like <- log(1/apply(1/exp(model.log.likelihoods), 1, mean))
        sum(each.like)
      }

    )
    )

### print.CIDnetwork <- function (x, ...) {}
### plot.CIDnetwork <- function (x, ...) {}
### summary.CIDnetwork <- function (object, ...) {}




###############################################################
####################  USER FUNCTIONS  #########################
###############################################################


#CID <- function (...) CIDnetwork$new(...)
CID.generate <- function (...) CID (..., generate=TRUE)


CID <- function (input,    # Must be edge.list, sociomatrix, CID.object or CID.Gibbs.object
                 outcome,  # Only used when input is an edge.list
                 n.nodes,
                 node.names,
                 intercept = 0,# Only used if input is missing.
                 components,
                 class.outcome="ordinal",
                 fill.in.missing.edges=missing(outcome),
                 generate=FALSE,
                 verbose=2,
                 ...) {
#  browser()
  if (missing(components)) components <- list()

  if (missing(input) & missing(n.nodes)) stop ("CID: no input was provided. Requires at least one of: sociomatrix, edge list, n.nodes.")
  if (!missing(input)){
    if(!(is(input,"matrix") | is(input,"data.frame") | is(input,"array"))) stop ("CID: input must be an edge.list or sociomatrix.")
  }
  if (missing(input)) {
    return(CIDnetwork$new(n.nodes=n.nodes, intercept=intercept,components=components, class.outcome=class.outcome, generate=TRUE, verbose=verbose, ...))
  } else { #get started!

    if ((is(input,"matrix") | is(input,"data.frame") | is(input,"array")) && ncol(input) == 2) { ##if we have an edgelist
      edge.list <- input

#this is an impossible state      if (!(class(edge.list) %in% c("matrix", "array", "data.frame"))) stop ("The edge.list provided must be an array, matrix or data frame with two columns; this is not an array, matrix or data frame.")
#also impossible      if (ncol(edge.list) != 2) stop ("The edge.list provided must have two columns.")

      edge.list <- cbind(as.character(edge.list[,1]), as.character(edge.list[,2]))
      if (any(edge.list[,1] == edge.list[,2])) {
        if(verbose > 0){
          message ("Removing self edges.")
        }
        nonselfies <- which(edge.list[,1] != edge.list[,2])
        edge.list <- edge.list[nonselfies,]
        if (!missing(outcome)) outcome <- outcome[nonselfies]
      }

      node.names.tmp <- unique(c(as.character(edge.list[,1]), as.character(edge.list[,2])))
      is.num.list <- FALSE
      if(!missing(node.names)) {
        if(length(node.names) >= length(node.names.tmp)){
          if(any(is.na(match(node.names.tmp,node.names)))){ ##node names don't match
            suppressWarnings(node.nums <- as.numeric(node.names.tmp))
            is.num.list <- !any(is.na(node.nums))
            if(!is.num.list){  ## If edge.list has non-numbers, use edge.list names
              if(verbose > 0)
                message("Warning:  edge.list contains nodes not named in node.names.")
              node.names <- node.names.tmp
            }else{
              node.names.tmp <- as.character(sort(node.nums))  ##  sorting the node name numbers
              if(max(node.nums) > length(node.names)){
                if(verbose > 0)
                  message("Warning:  edge.list contains node numbers larger than node.names vector.")
                node.names <- node.names.tmp
              }
            }
          }
        }else{
          if(verbose > 0)
            message("Warning:  Length of node.names less than n.nodes.  Using default node names")
          node.names <- node.names.tmp
        }
      }else{
        node.names <- node.names.tmp
      }

      n.nodes <- length(node.names)
      if(is.num.list){
        numbered.edge.list <- cbind(match(edge.list[,1], node.names.tmp),
                                    match(edge.list[,2], node.names.tmp))
      }else{
        numbered.edge.list <- cbind(match(edge.list[,1], node.names),
                                    match(edge.list[,2], node.names))
      }

      if (missing(outcome) | fill.in.missing.edges) {
        if(verbose > 0){
          if (missing(outcome))
            message("Assuming that this is a complete network with specified edges as binary ties.")
          else
            message("Filling in unspecified edges as zeroes.")
        }

        ##  OLD CODE THAT ASSUMED UNDIRECTED NETWORK
###        new.edge.list <- make.edge.list(n.nodes)
###        rowmatch <- sapply (1:nrow(edge.list),
###                            function(rr) min(which((numbered.edge.list[rr,1] == new.edge.list[,1] &
###                                                    numbered.edge.list[rr,2] == new.edge.list[,2]) |
###                                                   (numbered.edge.list[rr,2] == new.edge.list[,1] &
###                                                    numbered.edge.list[rr,1] == new.edge.list[,2]))))

        new.edge.list <- make.edge.list(n.nodes)
        new.edge.list <- rbind(new.edge.list,new.edge.list[,2:1])
        new.edge.list <- new.edge.list[order(new.edge.list[,1],new.edge.list[,2]),]


        rowmatch <- sapply (1:nrow(edge.list),
                            function(rr) which((numbered.edge.list[rr,1] == new.edge.list[,1] &
                                                numbered.edge.list[rr,2] == new.edge.list[,2])))

        rowmatch <- rowmatch[is.finite(rowmatch)]

        temp.outcome <- rep(0, nrow(new.edge.list))
        if (missing(outcome)) {
          temp.outcome[rowmatch] <- 1
        } else {
          temp.outcome[rowmatch] <- outcome[is.finite(rowmatch)]
        }

        outcome <- temp.outcome
        edge.list <- new.edge.list

      } else {
        edge.list <- numbered.edge.list
      }
    } else { ## we have a sociomatrix
      sociomatrix <- input
      if (nrow(sociomatrix) != ncol(sociomatrix)) stop ("The provided sociomatrix is not square.")

      n.nodes <- nrow(sociomatrix)
      if (any(sociomatrix != t(sociomatrix))) {
        if(verbose > 0)
          message ("CID: Auto-detected an asymmetric sociomatrix; assuming this is a directed network.")
        edge.list <- make.arc.list (n.nodes)
        outcome <- t(sociomatrix)[non.diag(n.nodes)]
      } else {
        if(verbose > 0)
          message ("CID: Auto-detected a symmetric sociomatrix; assuming this is an undirected network.")
        edge.list <- make.edge.list (n.nodes)
        outcome <- sociomatrix[u.diag(n.nodes)]  #just to be clear.
      }

      ## Inferring Node Names
      if(missing(node.names) || length(node.names) != ncol(sociomatrix)) {
        if(!missing(node.names)){
          if(verbose > 0) message("Warning: Length of node.names differs from nodes in sociomatrix.  Using default node.names")
        }
        if(!is.null(colnames(sociomatrix))){
          node.names <- colnames(sociomatrix)
        }else{
          node.names <- 1:n.nodes
        }
      }
    }

    ordinal.count <- max(outcome,na.rm=TRUE)+1

##    if (is.null(class.outcome)) {

    if (class.outcome=="ordinal") {
        outcome.na <- outcome[!is.na(outcome)]
      if (!all(round(outcome.na)==outcome.na)) {
        if(verbose > 0)
          message ("Detected non-integer values for outcomes; fitting this model with Gaussian outcomes.")
        class.outcome <- "gaussian"
      } else {
        if(verbose > 1)
          message ("Fitting: ordinal outcome with ",ordinal.count," states.")
      }
    } else {
      if(verbose > 1)
        message ("Fitting: Gaussian outcome.")
    }

    CID.object <- CIDnetwork$new(edge.list, n.nodes=n.nodes,
                                 intercept=intercept,
                                 components=components, outcome=outcome,
                                 class.outcome=class.outcome,
                                 ordinal.count=ordinal.count,
                                 node.names=as.character(node.names),
                                 generate=generate,
                                 verbose=verbose, ...)

    return(CID.object)
  }

}


CID.Gibbs <- function (input, # Must be edge.list, sociomatrix, CID.object or
                              # CID.Gibbs.object
                       outcome, # Only used when input is an edge.list
                       node.names,
                       ##edge.list,
                       ##outcome,
                       ##sociomatrix,

                       ##CID.object,
                       ##CID.Gibbs.object,

                       ##n.nodes=max(edge.list),

                       components, #=list(),
                       class.outcome=NULL,
                       fill.in.missing.edges=missing(outcome),
                       new.chain=FALSE,

                       draws=100,
                       burnin=-1,  ## burnin = -1 performs auto-burnin.
                       thin=10,
                       report=100,
                       auto.converge=FALSE,
                       extend.max = 10,
                       extend.count = 100,
                       verbose = 2,

                       ...) {
  #edge.list=n0$edge.list; outcome=n0$outcome; components=list(); n.nodes=max(edge.list)
  #edge.list=dolphins; components=list(LSM(2)); class.outcome=NULL
  #data(prison); sociomatrix=prison

  if(!is(input,"CID.Gibbs")) {  #missing(CID.Gibbs.object)){
    if(!is(input,"CIDnetwork")) {  #then it's a sociomatrix/edge.list/outcome combo.
      CID.object <- CID (input,
                         outcome,
                         components=components,
                         class.outcome,
                         node.names=node.names,
                         fill.in.missing.edges=fill.in.missing.edges,
                         verbose=verbose, ...)
    } else {
      CID.object <- input
      if (!missing(components)) CID.object$components <- components
      CID.object$reinitialize()
    }
    res <- CID.object$gibbs.full(draws=draws, burnin=burnin, thin=thin,
                                 report=report,
                                 auto.converge=auto.converge,
                                 extend.max = extend.max,
                                 extend.count = extend.count,
                                 verbose=verbose, ...)
  } else {
    if (!missing(components)) warning ("You cannot change the components of a CID.Gibbs object. Call this using the appropriate CID object if you wish to change the components.")
    CID.Gibbs.object <- input
    CID.object <- CID.Gibbs.object$CID.object
    res.2 <- CID.object$gibbs.full(draws=draws, burnin=burnin, thin=thin,
                                   report=report,
                                   auto.converge=auto.converge,
                                   extend.max = extend.max,
                                   extend.count = extend.count, ...)
    ##  Concatenating Chains
    if(new.chain){
      res <- res.2
    } else {
        res <- c(CID.Gibbs.object$results,res.2)
    }
  }

  # At this point we should have a CID object. Hurrah!


  DIC <- CID.object$DIC(res, add.parts=TRUE)

  #marg.ll <- CID.object$marginal.loglikelihood(res)
  #pseudo.CV.ll <- CID.object$pseudo.CV.loglikelihood(res)

  output <- list(results=res,  #unwrap.CID.Gibbs(
                 CID.object=CID.object,
                 CID.mean=CID.object$gibbs.mean(res),
                 DIC=DIC)
  class(output) <- "CID.Gibbs"

  return(output) #,
   #           marg.ll=marg.ll,
   #           pseudo.CV.ll=pseudo.CV.ll))

}


print.CID.Gibbs <- function (x, ...) {
    #x$CID.object$gibbs.summary(x$results, ...)
   x$CID.object$show(...)
}

summary.CID.Gibbs <- function (object, ...) {
  object$CID.object$gibbs.summary(object, ...)

}

print.summary.CID.Gibbs <- function(x, ...){
    x$CID.object$print.gibbs.summary(x)
}


plot.CID.Gibbs <- function (x, ...) {
  x$CID.object$gibbs.plot (x$results, x$DIC, ...)
}


