

# Nodal Intercept Network Model: Reference Class.
# Note: This is now directional. The undirected case is now BETA().

SRcid <-
  setRefClass(
    "SRcid",
    fields = list(

      intercept.sr="matrix",
      intercept.sr.Var="matrix",     #2x2 covariance matrix. Capitalized because it's a matrix.

      intercept.sr.Var.a="numeric",  # Wishart prior.
      intercept.sr.Var.b="matrix",   # Wishart prior.

      node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      edge.list.rows="list"    #,
      ),

    methods = list(
      initialize = function (

        intercept.sr.Var=diag(10,2),

        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        edge.list.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),

        intercept.sr.Var.a=0.001,
        intercept.sr.Var.b=diag(0.001, 2),

        intercept.sr = rmvnorm(n.nodes, rep(0,2), intercept.sr.Var),

        generate=FALSE

        ) {

        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$edge.list.rows <<- edge.list.rows
        .self$residual.variance <<- residual.variance
        .self$node.names <<- as.character(1:.self$n.nodes)

        .self$intercept.sr.Var <<- intercept.sr.Var
        .self$intercept.sr.Var.a <<- intercept.sr.Var.a
        .self$intercept.sr.Var.b <<- intercept.sr.Var.b
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
          message ("Reinitializing SR Intercepts")
          intercept.sr <<- rmvnorm(n.nodes, rep(0,2), intercept.sr.Var)
        }
        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)
      },

      pieces = function (include.name=FALSE) {
        out <- list(intercept.sr=intercept.sr, intercept.sr.Var=intercept.sr.Var)
        class(out) <- "SRout"
        #if (include.name) out <- c("SR", out)
        out
      },

      show = function () {
        message("intercept.sr.Var:"); print(intercept.sr.Var)
        message("intercept.sr:"); print(intercept.sr)
      },

      ## change to scatterplot.
      plot = function (coefs=intercept.sr, names=node.names,
        sd.a=NULL, sd.b=NULL,
        interval.a=NULL, interval.b=NULL,
        ...) {
        dotchart.two (coefs[,1], coefs[,2],
                      names, sd.a, interval.a, sd.b, interval.b,
                      ...)
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },


      value = function () {intercept.sr[edge.list[,1],1] + intercept.sr[edge.list[,2],2]}
      ,
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
#      value.ext = function (int.sr=intercept.sr, edges=1:nrow(edge.list)) {
        parameters[[1]][edge.list[edges,1],1] + parameters[[1]][edge.list[edges,2],2]
      },


      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },



      random.start = function () {

        intercept.sr.Var <<- solve(matrix(c(rWishart(1, 3, diag(0.01, 2))), nrow=2))
        ##intercept.sr.Var <<- rgamma(1, 1) #1/rgamma(1, intercept.sr.Var.ab[1], intercept.sr.Var.ab[1])
        intercept.sr <<- matrix(c(rmvnorm (n.nodes, rep(0,2), intercept.sr.Var)), ncol=2)
        rownames(intercept.sr) <<- node.names
      },


      draw = function(verbose=0) {

        edge.list.plus <- cbind(edge.list[,1], edge.list[,2]+n.nodes)
        intercept.block <- xtx(edge.list.plus, 2*n.nodes)
        intercept.outcome <- xtyc(edge.list.plus, outcome, 2*n.nodes)

        inv.sr.Var.big <- kronecker (solve(intercept.sr.Var), diag(1, n.nodes))

        varblock <- solve(intercept.block/residual.variance + inv.sr.Var.big)
        meanblock <- varblock%*%(intercept.outcome/residual.variance + 0)   #mean zero.

        intercept.sr <<- matrix(c(rmvnorm(1, meanblock, varblock)), ncol=2)
        rownames(intercept.sr) <<- node.names

        intercept.sr.Var <<-
          solve(matrix(c(rWishart(1,
                                  intercept.sr.Var.a + n.nodes/2,
                                  solve(intercept.sr.Var.b + t(intercept.sr)%*%intercept.sr/2))),
                       nrow=2))

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

      gibbs.mean = function(gibbs.out){
        get.sum <- gibbs.summary(gibbs.out)

        return(SR(intercept.sr.Var=intercept.sr.Var,
                 n.nodes=n.nodes,
                  edge.list=edge.list,
                  edge.list.rows=edge.list.rows,
                  residual.variance=residual.variance,
                  outcome=outcome,
                  intercept.sr.Var.a=intercept.sr.Var.a,
                  intercept.sr.Var.b=intercept.sr.Var.b,
                  intercept.sr=matrix(get.sum[,"mean"],ncol=2)))

      },

      gibbs.summary = function (gibbs.out) {
        coef.cov.mat <- sapply(gibbs.out, function(gg) c(gg$intercept.sr))
        ob1 <- cbind(mean=apply(coef.cov.mat, 1, mean),
                     sd=apply(coef.cov.mat, 1, sd),
                     q2.5=apply(coef.cov.mat, 1, quantile, 0.025),
                     q97.5=apply(coef.cov.mat, 1, quantile, 0.975))
        rownames(ob1) <- c(paste0(node.names,"S"), paste0(node.names,"R"))
        ob1
      },
      print.gibbs.summary = function (gibbs.sum) {
        message ("Sender-Receiver Coefficients:")
        print (gibbs.sum)
        return()
      },


      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (cbind(get.sum[1:n.nodes, 1], get.sum[n.nodes + 1:n.nodes,1]),

              interval.a = get.sum[1:n.nodes,3:4],
              interval.b = get.sum[n.nodes + 1:n.nodes,3:4],

              main = "SR-Intercept Summary from Gibbs Sampler", ...)
      },

      gibbs.node.colors = function (gibbs.out) {
        get.sum <- gibbs.summary(gibbs.out)
        value <- (get.sum[,1]-min(get.sum[,1]))/(max(get.sum[,1])-min(get.sum[,1]))
        rgb (rep(1, n.nodes), 1-value, 1-value)
      }

      )
    )



