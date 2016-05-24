
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R");
# Latent Vector Model: Reference Class
# Y_ij =  q_i'q_j + e_ij

#This version: no multiplicative factor, just yet. Too much bother at this point to get it right.

LVMcid <-
  setRefClass(
    "LVMcid",
    fields = list(
      dimension="numeric",

      latent.vector.pos="matrix",
      #mult.factor="numeric",
      #mult.factor.m="numeric",
      #mult.factor.v="numeric",

      latent.vector.pos.m="numeric",
      latent.vector.pos.V="matrix",
      latent.vector.pos.P="matrix",
      #latent.vector.tune="numeric",

      ##inherited from main. Must fix later, but OK for now.

      node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      edge.list.rows="list"    #,

      ),

    methods=list(

      initialize = function (

        dimension=1,

        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        edge.list.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),

        latent.vector.pos=matrix(rnorm(dimension*n.nodes), nrow=n.nodes),
        #mult.factor=-1,
        #mult.factor.m=0,
        #mult.factor.v=10000,

        latent.vector.pos.m=rep(0, dimension),
        latent.vector.pos.V=diag(10000, dimension),
        #latent.vector.tune=0.1,

        generate=FALSE

        ) {

        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$edge.list.rows <<- edge.list.rows
        .self$residual.variance <<- residual.variance
        .self$node.names <<- as.character(1:.self$n.nodes)

        .self$dimension <<- dimension
        .self$latent.vector.pos <<- latent.vector.pos
        .self$latent.vector.pos.m <<- latent.vector.pos.m
        .self$latent.vector.pos.V <<- latent.vector.pos.V
        .self$latent.vector.pos.P <<- solve(latent.vector.pos.V)

        #.self$mult.factor <<- mult.factor
        #.self$mult.factor.m <<- mult.factor.m
        #.self$mult.factor.v <<- mult.factor.v
        #.self$latent.vector.tune <<- latent.vector.tune

        #adjust.lsp()

        if (generate) .self$generate() else .self$outcome <<- outcome

      },

      #adjust.lsp = function (mult.up=TRUE) {
      #  mft <- mean(edge.list.distance(latent.vector.pos, edge.list))
      #  mult.factor <<- mult.factor*mft
      #  latent.vector.pos <<- latent.vector.pos/mft
      #},

      reinitialize = function (n.nodes=NULL,
        edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes  #.self$
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          edge.list.rows <<- row.list.maker(edge.list)
        }
        if (nrow(latent.vector.pos) != .self$n.nodes) {
          message ("Reinitializing LVM Vectors")
          latent.vector.pos <<- matrix(rnorm(dimension*n.nodes), nrow=n.nodes)
          #adjust.lsp()
        }
        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)
      },

      pieces = function (include.name=FALSE) {
        out <- list (latent.vector.pos=latent.vector.pos) #, mult.factor=mult.factor)
        class(out) <- "LVMout"
        out
      },

      show = function () {
        message("t(latent.vector.pos):"); print(t(latent.vector.pos))
#       message("mult.factor:"); print(mult.factor)
      },
      plot = function (pos=latent.vector.pos, ...) {
        latent.space.plot (pos, arrowlines=TRUE, labels=node.names, ...)
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },


      value = function () {cosine.closeness(latent.vector.pos, edge.list)},
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
        cosine.closeness(parameters[[1]], rbind(edge.list[edges,])) },


      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },



      random.start = function () {
        latent.vector.pos <<- matrix(rnorm(dimension*n.nodes), nrow=n.nodes)
        #mult.factor <<- rnorm(1, mult.factor.m, sqrt(mult.factor.v))
      },

      draw = function (verbose=0) {  # tune=latent.vector.tune
      #d1 <- LSMcid$new(); latent.vector.pos <- d1$latent.vector.pos; mult.factor <- d1$mult.factor; edge.list <- d1$edge.list; edge.list.rows <- d1$edge.list.rows; n.nodes <- d1$n.nodes; mult.factor.m=0; mult.factor.v=10000; latent.vector.tune=0.1

        lsdim <- dim(latent.vector.pos)[2]
        latent.vector.pos.hold <- latent.vector.pos

        for (dd in 1:n.nodes) {

          #Gibbs draw for one node given others. Direct!
          row1 <- which(edge.list[,1] == dd)
          row2 <- which(edge.list[,2] == dd)

          #get the counterpart.
          xx.mat <- rbind(matrix(latent.vector.pos.hold[edge.list[row1,2],], ncol=lsdim),
                          matrix(latent.vector.pos.hold[edge.list[row2,1],], ncol=lsdim))

          cls.VV <- solve(t(xx.mat)%*%xx.mat/residual.variance + latent.vector.pos.P)
          cls.mean <- cls.VV%*%(t(xx.mat)%*%outcome[c(row1,row2)]/residual.variance + latent.vector.pos.P%*%latent.vector.pos.m)

          latent.vector.pos.hold[dd,] <- c(rmvnorm(1, cls.mean, cls.VV))

        }

  #Rotate back.
        latent.vector.pos.hold <-
          postprocess.latent.positions(latent.vector.pos.hold, recenter=FALSE)
        rownames(latent.vector.pos.hold) <- node.names

        latent.vector.pos <<- latent.vector.pos.hold

      },

      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1, make.random.start=FALSE) {
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("LVM ",index)
          } else if (round(index)==index) {
            if (report.interval > 0) if (index %% report.interval == 0) message("LVM burnin ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        lsp.all <- sapply(gibbs.out, function(gg) gg$latent.vector.pos)
        output <- matrix(apply(lsp.all, 1, mean), nrow=n.nodes)
        rownames(output) <- node.names
        colnames(output) <- paste0("pos",1:ncol(output))
        return(output)
      },
      print.gibbs.summary = function (gibbs.out) {
        get.sum <- gibbs.summary(gibbs.out)
        message ("Mean Latent Vector Positions:")
        print(get.sum)
        return(invisible(get.sum))
      },

      gibbs.mean = function(gibbs.out){
        get.sum <- gibbs.summary(gibbs.out)

        return(LVM(dimension=dimension,n.nodes=n.nodes,
                   edge.list=edge.list,
                   edge.list.rows=edge.list.rows,
                   residual.variance=residual.variance,
                   outcome=outcome,
                   latent.vector.pos=get.sum,
                   latent.vector.pos.m=latent.vector.pos.m,
                   latent.vector.pos.V=latent.vector.pos.V))
      },

      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum, main = "Mean Latent Vector Positions from Gibbs Sampler", ...)
      },

      gibbs.node.colors = function (gibbs.out) {
        rep("#DDDDFF", n.nodes)
      }

      )
    )




