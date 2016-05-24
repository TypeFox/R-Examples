glmc.control <- function(epsilon.glm = 1e-8, maxit.glm= 100, trace.glm= FALSE,trace.optim= 0,fnscale.optim=-1, parscale.optim = rep.int(1,1), ndeps.optim = rep.int(0.001,1),maxit.optim = 100,abstol.optim = -Inf, reltol.optim= sqrt(.Machine$double.eps), alpha.optim = 1, beta.optim = 0.5, gamma.optim = 2, REPORT.optim= 10, type.optim = 1, lmm.optim = 5,factr.optim= 1e+07, pgtol.optim = 0, tmax.optim = 10, temp.optim =10,maxit.weights = 25, gradtol.weights = 1e-07, svdtol.weights = 1e-09,itertrace.weights = FALSE)
{

#glm

    if(!is.numeric(epsilon.glm) || epsilon.glm <= 0)
	stop("value of epsilon.glm must be > 0")
    if(!is.numeric(maxit.glm) || maxit.glm <= 0)
	stop("maximum number of iterations must be > 0")

#optim

    if(!is.numeric(reltol.optim) || reltol.optim <= 0)
	stop("value of reltol.optim must be > 0")
    if(!is.numeric(maxit.optim) || maxit.optim <= 0)
	stop("maximum number of iterations in optim must be > 0")
    if(trace.optim < 0){
        stop("trace.weights should be a non-negative integer. If positive, tracing information on\n the progress of the optimization is produced. Higher values\n may produce more tracing information.")
    }

#weights

   if(!is.numeric(gradtol.weights) || gradtol.weights <= 0)
	stop("value of gradtol.weights must be > 0")
    if(!is.numeric(maxit.weights) || maxit.weights <= 0)
	stop("maximum number of iterations for weight step must be > 0")
    if(itertrace.weights < 0){
        stop("itertrace.weights should be a non-negative integer. If positive, tracing information on\n the progress of the optimization is produced. Higher values\n may produce more tracing information.")
    }

glm.control<-list(epsilon=epsilon.glm, maxit=maxit.glm, trace=trace.glm)

weights.control<-list(maxit=maxit.weights, gradtol=gradtol.weights, svdtol=svdtol.weights, itertrace=itertrace.weights)

optim.control<-list(trace=trace.optim,fnscale=fnscale.optim,parscale=parscale.optim , ndeps = ndeps.optim, maxit=maxit.optim,abstol=abstol.optim, reltol=reltol.optim, alpha=alpha.optim, beta=beta.optim, gamma=gamma.optim , REPORT=REPORT.optim, type=type.optim , lmm=lmm.optim,factr=factr.optim, pgtol=pgtol.optim, tmax=tmax.optim, temp=temp.optim)
   
#    list(epsilon = epsilon, maxit = maxit, trace = trace,reltol.weights = reltol.weights,maxit.weights = maxit.weights, trace.weights = trace.weights)

list(glm.control=glm.control,optim.control=optim.control,weights.control=weights.control)
}
