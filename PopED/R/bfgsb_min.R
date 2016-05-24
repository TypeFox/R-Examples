#' Nonlinear minimization using BFGS with box constraints
#' 
#' This is the implementation of a Broyden Fletcher Goldfarb Shanno (BFGS) 
#' method for nonlinear minimization with box constraints.
#' 
#' @param f_name A function name (as a text string) that returns an objective function and the gradient of that objective function, in that order. 
#' See \code{\link{calc_ofv_and_grad}} as used in \code{\link{Doptim}}.
#' @param f_options Options for the f_name argument.
#' @param x0 the intial values to optimize
#' @param l the lower bounds
#' @param u the upper bounds 
#' @param options a list of additional settings arguments
#'  
#' @return A list containing:
#' \item{x_k}{The objective function.}
#' \item{f_k}{The gradient.}
#' \item{B_k}{The hessian.}
#' @family Optimize
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_optimize.R
#' @example tests/testthat/examples_fcn_doc/examples_bfgsb_min.R
#' 
#' @export
#' @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

bfgsb_min <- function(f_name,f_options, x0, l, u, options=list()){
    ## This is the implementation of a Broyden Fletcher Goldfarb Shanno (BFGS)
    ## method for nonlinear minimization with box constraints.
    if(any(x0>u | x0<l)){
        stop(sprintf('Starting value violates given bounds!'))
    }
    tp=0
    if((size(x0,2)>1)){
        tp=1
        x0=t(x0)
        l=t(l)
        u=t(u)
    }
    x_k=x0
    ##initial hessian is the identity
    B_k=diag(rep(1,length(x0)))
    ##calc f at x_k
    f_options[[1]] <- x_k
    returnArgs <- do.call(f_name,f_options) 
    f_k <- returnArgs[[1]]
    gf_k <- returnArgs[[2]]
    
    ##options
    if(is.null(options$pgtol)) options$pgtol=1e-4
    if(is.null(options$factr)) options$factr=1e-8
    if(is.null(options$update_fhandle)) options$update_fhandle="print_information"
    
    do.call(options$update_fhandle,list(0, f_k, x_k))
    iter=1
    ##as long as gradient is bigger than certain value
    while(projected_gradient_norm(l,u,x_k,gf_k)>options$pgtol){
        ##calc cauchy point
        x_cp=cauchy_point(x_k, l, u, gf_k, B_k)
        ##determine search direction ignoring the bounds
        xk1=subspace_min(x_k, l, u, x_cp, gf_k, B_k)
        dk=xk1-x_k
        ##do line seach along dk
        returnArgs <- line_search(f_name,f_options,l, u, x_k, f_k, gf_k, dk,options) 
        f_k1 <- returnArgs[[1]]
        gf_k1 <- returnArgs[[2]]
        x_k1 <- returnArgs[[3]]
        s_k=x_k1-x_k
        ##stop if change in f is too small
        if((f_k-f_k1)<=options$factr*max(matrix(c(abs(f_k1),abs(f_k),1),nrow=1,byrow=T))){
            break
        }
        do.call(options$update_fhandle,list(iter, f_k1, x_k1))
        y_k=gf_k1-gf_k
        x_k=x_k1
        gf_k=gf_k1
        f_k=f_k1
        ##update hessian with BFGS formula if Wolfe conditions are met
        if((t(y_k)%*%s_k)/(-t(gf_k)%*%s_k) > .Machine$double.eps){
            den1 <- t(y_k)%*%s_k
            den2 <- t(s_k)%*%B_k%*%s_k
            B_k=B_k+y_k%*%t(y_k)/den1[,]-B_k%*%s_k%*%t(B_k%*%s_k)/den2[,]
        }
        iter=iter+1
    }
    if(tp){
        x_k=t(x_k)
    }
    return(list(  x_k=  x_k, f_k= f_k, B_k  = B_k  )) 
}

print_information <- function(iter, f, x){
    fprintf('%d\t%g\t', iter,f)
    fprintf('%6.2f\t', x)
    fprintf('\n')
    return() 
}

