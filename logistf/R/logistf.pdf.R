logistf.pdf <-
function (x, y, pos,
 firth = TRUE, weight, control, plcontrol, offset=NULL, b, beta=NULL, loglik=NULL, fit=NULL, old=FALSE)
{

    # Georg Heinze, 3 Jan 2012
    # georg.heinze@meduniwien.ac.at
    # This function can be called similarly as logistf.fit, but one has to provide a beta value to evaluate and a variable

    if (missing(control))
        control <- logistf.control()
    if (missing(plcontrol))
        plcontrol <- logistpl.control()
    if (is.null(offset))
        offset <- rep(0,nrow(x))
    res<-matrix(0,1,3)
    k <- ncol(x)
    n <- nrow(x)

    if (is.null(fit) & is.null(beta) & is.null(loglik)) {
      if(old) fit <- logistf.fit.old(x, y, weight = weight, offset = offset, firth = firth, control = control)
      else fit <- logistf.fit(x, y, weight = weight, offset = offset, firth = firth, control = control)
      beta<-fit$coefficients
      loglik<-fit$loglik[2]
    }
#    else  if(is.null(beta) | is.null(loglik)) stop("logistf.pdf: Please specify either beta and loglik or fit object.\n")
    
#      std.pos <- diag(fit$var)[pos]^0.5
#      coefs <- fit$beta
#      covs <- fit$var
      init<-beta
      init[pos]<-b
#      loglik.fit<-fit$loglik[2]
#      alltimes <<- alltimes + system.time(
#       vorzeit<-Sys.time()
        if(old) xx <- logistf.fit.old(x, y, weight = weight, offset = offset,
                firth = firth, col.fit <- (1:k)[-pos], init = init,
                control = control)
        else  xx <- logistf.fit(x, y, weight = weight, offset = offset,
                firth = firth, col.fit <- (1:k)[-pos], init = init,
                control = control)
#                  )
#       alltimes<<- alltimes+ Sys.time()-vorzeit
#        ncalls <<- ncalls + 1
#        init<-xx$beta
        res[1,1]<-b
        res[1,2]<-2*(loglik-xx$loglik)
        res[1,3]<-1-(1-pchisq(res[,2],1))/2
        res[1,3][res[,1]<beta[pos]]<-1-res[,3][res[,1]<beta[pos]]
        results<-list(beta=res[1,1], chisq=res[1,2], pdf=res[1,3])

     results
   }

