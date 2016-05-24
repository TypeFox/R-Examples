setMethod(
          'coef',
          'browntree',
          function (object, ...) {
            list(
                 sigma=object@sigma,
                 theta=object@theta,
                 sigma.sq.matrix=sym.par(object@sigma)
                 )
          }
          )

setMethod(
          'coef',
          'hansentree',
          function (object, ...) {
            list(
                 sqrt.alpha=object@sqrt.alpha,
                 sigma=object@sigma,
                 theta=object@theta,
                 alpha.matrix=sym.par(object@sqrt.alpha),
                 sigma.sq.matrix=sym.par(object@sigma)
                 )
          }
          )

# coerce methods allow for coercion to a data.frame
setAs(
      from='ouchtree',
      to='data.frame',
      def = function (from) {
        df <- data.frame(
                         nodes=from@nodes,
                         ancestors=from@ancestors,
                         times=from@times,
                         labels=from@nodelabels,
                         row.names=from@nodes
                         )
        rownames(df) <- from@nodes
        df
      }
      )

setAs(
      from='browntree',
      to='data.frame',
      def = function (from) {
        cbind(
              as(as(from,'ouchtree'),'data.frame'),
              as.data.frame(from@data)
              )
      }
      )

setAs(
      from='hansentree',
      to='data.frame',
      def = function (from) {
        cbind(
              as(as(from,'ouchtree'),'data.frame'),
              as.data.frame(from@regimes),
              as.data.frame(from@data)
              )
      }
      )


setMethod('logLik','browntree',function(object)object@loglik)

setMethod(
          "summary",
          "browntree",
          function (object, ...) {
            cf <- coef(object)
            dof <- object@nchar+length(object@sigma)
            deviance=-2*logLik(object)
            aic <- deviance+2*dof
            aic.c <- aic+2*dof*(dof+1)/(object@nterm*object@nchar-dof-1)
            sic <- deviance+log(object@nterm*object@nchar)*dof
            list(
                 call=object@call,
                 sigma.squared=cf$sigma.sq.matrix,
                 theta=cf$theta,
                 loglik=logLik(object),
                 deviance=deviance,
                 aic=aic,
                 aic.c=aic.c,
                 sic=sic,
                 dof=dof
                 )
          }
          )

setMethod('logLik','hansentree',function(object)object@loglik)

setMethod(
          "summary",
          "hansentree",
          function (object, ...) {
            cf <- coef(object)
##            if (length(object@hessian)>0)
##              se <- sqrt(diag(solve(0.5*object@hessian)))
            dof <- length(object@sqrt.alpha)+length(object@sigma)+sum(sapply(object@theta,length))
            deviance=-2*logLik(object)
            aic <- deviance+2*dof
            aic.c <- aic+2*dof*(dof+1)/(object@nterm*object@nchar-dof-1)
            sic <- deviance+log(object@nterm*object@nchar)*dof
            list(
                 call=object@call,
                 conv.code=object@optim.diagn$convergence,
                 optimizer.message=object@optim.diagn$message,
                 alpha=cf$alpha.matrix,
                 sigma.squared=cf$sigma.sq.matrix,
                 optima=cf$theta,
                 loglik=logLik(object),
                 deviance=deviance,
                 aic=aic,
                 aic.c=aic.c,
                 sic=sic,
                 dof=dof
                 )
          }
          )

setMethod(
          'print',
          'ouchtree',
          function (x, ...) {
            print(as(x,'data.frame'),...)
            invisible(x)
          }
          )

setMethod(
          'print',
          'browntree',
          function (x, ...) {
            cat("\ncall:\n")
            print(x@call)
            print(as(x,'data.frame'),...)
            sm <- summary(x)
            cat('\nsigma squared:\n')
            print(sm$sigma.squared)
            cat('\ntheta:\n')
            print(sm$optima)
            print(unlist(sm[c("loglik","deviance","aic","aic.c","sic","dof")]))
            invisible(x)
          }
          )

setMethod(
          'print',
          'hansentree',
          function (x, ...) {
            cat("\ncall:\n")
            print(x@call)
            print(as(x,'data.frame'),...)
            if (length(x@optim.diagn)>0) {
              if (x@optim.diagn$convergence!=0)
                cat("\n",sQuote("optim")," convergence code: ",x@optim.diagn$convergence)
              if (!is.null(x@optim.diagn$message))
                cat("\n",sQuote("optim")," diagnostic message: ",x@optim.diagn$message)
            }
            sm <- summary(x)
            cat('\nalpha:\n')
            print(sm$alpha)
            cat('\nsigma squared:\n')
            print(sm$sigma.squared)
            cat('\ntheta:\n')
            print(sm$optima)
            print(unlist(sm[c("loglik","deviance","aic","aic.c","sic","dof")]))
            invisible(x)
          }
          )

setMethod(
          'show',
          'ouchtree',
          function (object) {
            print(as(object,'ouchtree'))
            invisible(NULL)
          }
          )

setMethod(
          'show',
          'browntree',
          function (object) {
            print(as(object,'browntree'))
            invisible(NULL)
          }
          )

setMethod(
          'show',
          'hansentree',
          function (object) {
            print(as(object,'hansentree'))
            invisible(NULL)
          }
          )
