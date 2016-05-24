
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R");

# Covariates component: Reference Class

COVARIATEcid <-
  setRefClass(
    "COVARIATEcid",
    fields = list(
      covariates="array",
      coef.cov="numeric",
      coef.cov.m="numeric",
      coef.cov.V="matrix",
      coef.cov.P="matrix",

      cov.type="character",

      cov.names="character",

      cov.block="matrix",

      node.names="character",
      cov.node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      edge.list.rows="list"    ,
      center="logical"

      ),

    methods=list(
      initialize = function (

        covariates=array(1, c(nrow(edge.list),1)),
        coef.cov=c(), #rep(0, ncol(covariates)),

        cov.names=character(), #as.character(colnames(covariates)),

        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        edge.list.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),

        coef.cov.m=0, #ncol(covariates)),
        coef.cov.V=diag(1000000, 1), #rep(ncol(covariates), 2)),

        cov.type=c("Edge","Sender","Receiver","SendRec","Identical"),
        center=TRUE,


        generate=FALSE
        ) {

        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$edge.list.rows <<- edge.list.rows
        .self$residual.variance <<- residual.variance
        .self$node.names <<- as.character(1:.self$n.nodes)
        .self$center <<- center
        .self$cov.type <<- match.arg(cov.type)
        .self$cov.node.names <<- as.character(NA)

        #handle covariates.
        cov.arg <- covariates
        ##        if(cov.type == "Edge"){
        if (!(class(cov.arg) %in% c("array", "matrix","data.frame","numeric"))) stop ("COV: Covariates must be a matrix, array or data.frame object.")

        if(class(cov.arg) %in% c("array","matrix","numeric")){
          if(!is(cov.arg,"array")){ #if given a vector, convert to column
            cov.arg <- as.matrix(cov.arg)
          }
          if(cov.type != "Edge"){
            if(!is.null(rownames(cov.arg)))
              .self$cov.node.names <<- rownames(cov.arg)
          }
        }else{
          if(cov.type != "Edge"){
            if(!is.numeric(cov.arg[,1])){
              if(ncol(cov.arg) > 1){
                .self$cov.node.names <<- as.character(cov.arg[,1])
                cov.arg <- cov.arg[,-1]
              }else{
                stop("COV:  Covariates must have at least one numeric column")
              }
            }else{
              if(!is.null(rownames(cov.arg))){
                message("COV:  Assuming rownames correspond to node names")
                .self$cov.node.names <<- rownames(cov.arg)
              }
            }
          }
          cov.arg <- as.matrix(cov.arg)
          if(!is.numeric(cov.arg))
            stop("COV:  Covariates must be numeric")
        }


        .self$covariates <<- cov.arg
        rm(covariates)  #  removing the passed version of covariates
        if (nrow(covariates) == ncol(covariates)) { #If matrix form covariates
          if (length(dim(covariates)) != 3) covariates <<- array (covariates, c(dim(covariates), 1))
          temp.dim <- dim(covariates)[3]
          temp.names <- dimnames(covariates)[[3]]
          if (length(coef.cov) == temp.dim) .self$coef.cov <<- coef.cov else .self$coef.cov <<- rep(0, temp.dim)
        } else {
          temp.dim <- ncol(covariates)   #works for nodal or edge list covariates.
          temp.names <- colnames(covariates)
          if (length(coef.cov) == temp.dim) .self$coef.cov <<- coef.cov else .self$coef.cov <<- rep(0, temp.dim)
        }

        if(length(cov.names) == 0){
          if(is.null(temp.names)){
            .self$cov.names <- paste("COV",1:length(.self$coef.cov),sep=".")
          }else{
            .self$cov.names <- temp.names
          }
        }

#        if (length(cov.names) != length(.self$coef.cov)) {
#          warning ("Covariate matrix columns: ", ncol(.self$covariates), ", covariate names ", length(cov.names), ";replacing with defaults.")
#          .self$cov.names <- paste0 ("cov", 1:length(.self$coef.cov))
#        } else {
#          .self$cov.names <<- cov.names
#        }

        if (length(coef.cov.m) != temp.dim) .self$coef.cov.m <<- rep(coef.cov.m, temp.dim)[1:temp.dim] else .self$coef.cov.m <<- coef.cov.m
        if (nrow(coef.cov.V) != temp.dim) .self$coef.cov.V <<- diag(diag(coef.cov.V), temp.dim) else .self$coef.cov.V <<- as.matrix(coef.cov.V)

        #.self$cov.block <<- t(covariates)%*%covariates

        .self$coef.cov.P <<- solve(.self$coef.cov.V)

        if (generate) .self$generate() else .self$outcome <<- outcome

      },

      reinitialize = function (n.nodes=NULL, edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          edge.list.rows <<- row.list.maker(edge.list)
        }
        if (!is.null(node.names) & length(node.names) == .self$n.nodes) {
          node.names <<- node.names
        }else
          node.names <<- as.character(1:.self$n.nodes)


        if (cov.type == "Edge") {

          if (nrow(covariates) == n.nodes & ncol(covariates) == n.nodes) {
            message("Detected a sociomatrix-style covariate array. Adjusting to match the edge list.")
            covariates <<- as.matrix(mat.cov.to.edge.list.cov (covariates, arc.list=edge.list))
          }

        }else{

          if(!any(is.na(cov.node.names))) {
            names.match <- match(node.names,cov.node.names)
            if(!any(is.na(names.match))){
              covariates <<- covariates[names.match,,drop=FALSE]
            }else{
              message("COV:  node.names for covariates don't match those for edge.list.  Inaccurate coefficients may result")
            }
          }

          if (cov.type == "Sender")
            covariates <<- covariates[edge.list[,1],,drop=FALSE]

          if (cov.type == "Receiver")
            covariates <<- covariates[edge.list[,2],,drop=FALSE]

          if (cov.type == "SendRec")
            covariates <<- covariates[edge.list[,1],,drop=FALSE] + covariates[edge.list[,2],,drop=FALSE]

          if (cov.type == "Identical")
            covariates <<- 1*(covariates[edge.list[,1],,drop=FALSE] == covariates[edge.list[,2],,drop=FALSE])
        }

        #Now, we can correct! Check dimensions of covariate matrix against edge list.

        if (nrow(covariates) < nrow(edge.list)) {
          stop (paste0("The number of rows in the COV matrix, ",nrow(covariates)," is less than the number of edges, ",nrow(edge.list),"."))
        }
        if (nrow(covariates) > nrow(edge.list)) {
          stop (paste0("The number of rows in the COV matrix, ",nrow(covariates)," is greater than the number of edges, ",nrow(edge.list),"."))
        }

        ## .self$cov.block <<- t(covariates)%*%covariates
        if(center){
            covariates.new  <- t( t(covariates) - colMeans(covariates))
            if(!identical(covariates.new, covariates)){
                message("Centering covariates to speed sampler.  To disable centering pass center = FALSE")
                covariates <<- covariates.new
            }

        }
        cov.block <<- t(covariates)%*%covariates  ##



      },

      pieces = function (include.name=FALSE) {
        out <- list (coef.cov=coef.cov)
        class(out) <- "COVout"
        out
      },

      show = function (show.cov=FALSE) {
        message("coef.cov:"); print(coef.cov)
        if (show.cov) {message("t(covariates):"); print(t(covariates))}
      },
      plot = function (coefs=coef.cov, names=cov.names, sd=NULL, interval=NULL, main = "Covariate Summary from Class Values", ...) {
        dotchart.coef (coefs, names, sd, interval, main=main, ...)
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },


      value = function () {covariates%*%coef.cov},
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {cbind(covariates[edges,])%*%parameters[[1]]},



      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },



      random.start = function () {coef.cov <<- c(rmvnorm (1, coef.cov.m, coef.cov.V))},   #ncol(covariates)

      draw = function (verbose=0) {
        varblock <- solve(cov.block/residual.variance + coef.cov.P)
        meanblock <- varblock%*%(t(covariates)%*%outcome/residual.variance +
                                 coef.cov.P%*%coef.cov.m)
        coef.cov <<- c(rmvnorm(1, meanblock, varblock))
      },

      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1, make.random.start=FALSE) {
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("COV ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        coef.cov.mat <- rbind(sapply(gibbs.out, function(gg) gg$coef.cov))
        ob1 <- data.frame(estimated.mean=round(apply(coef.cov.mat, 1, mean),3),
                          estimated.sd=round(apply(coef.cov.mat, 1, sd),3),
                          q2.5=round(apply(coef.cov.mat, 1, quantile, 0.025),3),
                          q97.5=round(apply(coef.cov.mat, 1, quantile, 0.975),3))
        rownames(ob1) <- cov.names
        return(ob1)
      },
      print.gibbs.summary = function (gibbs.sum) {
        message (cov.type," Coefficients:")
        print (gibbs.sum)
        return()
      },

      gibbs.mean = function(gibbs.out){
        get.sum <- gibbs.summary(gibbs.out)

        return(EdgeCOV(covariates=covariates,
                       coef.cov=get.sum[,"estimated.mean"],
                       cov.names=cov.names,
                       n.nodes=n.nodes,
                       edge.list=edge.list,
                       edge.list.rows=edge.list.rows,
                       residual.variance=residual.variance,
                       outcome=outcome,
                       coef.cov.m=coef.cov.m,
                       coef.cov.V=coef.cov.V,
                       cov.type=cov.type,
                       center=center))
      },


      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum[,1], interval = get.sum[,3:4], main = "Covariate Summary from Gibbs Sampler", ...)
      },

      gibbs.node.colors = function (gibbs.out) {
        rep("#DDDDFF", n.nodes)
      }

      )
    )



#SenderCOVcid <-
#  setRefClass(
#    "SenderCOVcid",
#    contains=c("EdgeCOVcid"),

#    methods=list(

#      reinitialize = function (n.nodes=NULL, edge.list=NULL)
#      {
#        if (!is.null(n.nodes)) n.nodes <<- n.nodes
#        if (!is.null(edge.list)) {
#          edge.list <<- edge.list
#          edge.list.rows <<- row.list.maker(edge.list)
#        }

        #reinitialize covariates to match.
#        if (!(class(covariates) %in% c("array", "matrix"))) stop ("COV: Covariates must be a matrix or array object.")
#        covariates <<- as.array(cbind(covariates))
#        if (nrow(covariates) != n.nodes) stop ("SenderCOV: number of nodes is not the number of rows of covariates.")

#        covariates <<- covariates[edge.list[,1],]
#      }


#      )
#    )


