

# Nodal Intercept Network Model: Reference Class.
# Note: This is now directional. The undirected case is now BETA().

BETAcid <-
  setRefClass(
    "BETAcid",
    fields = list(

      intercept.sr="matrix",
      intercept.sr.var="numeric",
      intercept.sr.var.ab="numeric",

      node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      edge.list.rows="list"    #,
      ),

    methods = list(
      initialize = function (

        intercept.sr.var=1,

        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        edge.list.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),

        intercept.sr.var.ab=c(0.001, 0.001),

        intercept.sr=cbind(rnorm(n.nodes, 0, sqrt(intercept.sr.var))),

        generate=FALSE

        ) {

        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$edge.list.rows <<- edge.list.rows
        .self$residual.variance <<- residual.variance
        .self$node.names <<- as.character(1:.self$n.nodes)

        .self$intercept.sr.var <<- intercept.sr.var
        .self$intercept.sr.var.ab <<- intercept.sr.var.ab
        .self$intercept.sr <<- intercept.sr

        if (generate) .self$generate() else .self$outcome <<- outcome

      },

      reinitialize = function (n.nodes=NULL,
        edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          edge.list.rows <<- row.list.maker(edge.list)
        }
        if (length(intercept.sr) != n.nodes) {
          message ("Reinitializing BETA Intercepts")
          intercept.sr <<- cbind(rnorm(n.nodes, 0, sqrt(intercept.sr.var)))
        }
        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)
      },

      pieces = function (include.name=FALSE) {
        out <- list(intercept.sr=intercept.sr, intercept.sr.var=intercept.sr.var)
        class(out) <- "BETAout"
        #if (include.name) out <- c("SR", out)
        out
      },

      show = function () {
        message("intercept.sr.var:"); print(intercept.sr.var)
        message("intercept.sr:"); print(intercept.sr)
      },
      plot = function (coefs=intercept.sr, names=node.names, sd=NULL, interval=NULL, ...) {
        dotchart.coef (coefs, names, sd, interval, ...)
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },


      value = function () {intercept.sr[edge.list[,1]] + intercept.sr[edge.list[,2]]},
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
#      value.ext = function (int.sr=intercept.sr, edges=1:nrow(edge.list)) {
        parameters[[1]][edge.list[edges,1]] + parameters[[1]][edge.list[edges,2]]
      },


      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },



      random.start = function () {
        intercept.sr.var <<- rgamma(1, 1) #1/rgamma(1, intercept.sr.var.ab[1], intercept.sr.var.ab[1])
        intercept.sr <<- cbind(rnorm (n.nodes, 0, sqrt(intercept.sr.var)))
        names(intercept.sr) <<- node.names
      },

      draw = function(verbose=0) {

        intercept.block <- xtx(edge.list, n.nodes)
        intercept.outcome <- xtyc(edge.list, outcome, n.nodes)

        varblock <- solve(intercept.block/residual.variance + diag(1/intercept.sr.var, n.nodes))
        meanblock <- varblock%*%(intercept.outcome/residual.variance + 0)   #mean zero.

        intercept.sr <<- cbind(c(rmvnorm(1, meanblock, varblock)))
        names(intercept.sr) <<- node.names

        intercept.sr.var <<-
          1/rgamma(1,
                   intercept.sr.var.ab[1] + n.nodes/2,
                   intercept.sr.var.ab[2] + sum(intercept.sr^2)/2)

      },

      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1, make.random.start=FALSE) {
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("SR ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        coef.cov.mat <- sapply(gibbs.out, function(gg) gg$intercept.sr)
        ob1 <- cbind(mean=apply(coef.cov.mat, 1, mean),
                     sd=apply(coef.cov.mat, 1, sd),
                     q2.5=apply(coef.cov.mat, 1, quantile, 0.025),
                     q97.5=apply(coef.cov.mat, 1, quantile, 0.975))
        rownames(ob1) <- node.names
        ob1
      },
      print.gibbs.summary = function (gibbs.out) {
        get.sum <- gibbs.summary(gibbs.out)
        message ("Sender-Receiver Coefficients:")
        print (get.sum)
        return(invisible(get.sum))
      },

      gibbs.mean = function(gibbs.out){
        get.sum <- gibbs.summary(gibbs.out)

        return(BETA(intercept.sr.var=1,n.nodes=10,
                    edge.list=edge.list,
                    edge.list.rows=edge.list.rows,
                    residual.variance=residual.variance,
                    outcome=outcome,
                    intercept.sr.var.ab=intercept.sr.var.ab,
                    intercept.sr=cbind(get.sum[,"mean"])))
      },


      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum[,1], interval = get.sum[,3:4], main = "BETA-Intercept Summary from Gibbs Sampler", ...)
      },

      gibbs.node.colors = function (gibbs.out) {
        get.sum <- gibbs.summary(gibbs.out)
        value <- (get.sum[,1]-min(get.sum[,1]))/(max(get.sum[,1])-min(get.sum[,1]))
        rgb (rep(1, n.nodes), 1-value, 1-value)
      }

      )
    )



