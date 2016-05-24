## R interface function that calls C-function 'minQuad' and passes arguments by 'reference' (DUP = FALSE):
## Any variable changed in C must be initialized in C, therefore arguments to C-function are partitioned as:
## (i)   READ-ONLY-C  or INPUT variables : can  be initialized in R and unchanged in C
## (ii)  WRITE-ONLY-C or OUTPUT variables: must be initialized in C (R's may be 'overwritten') and are changed in C
## NOTE: R-OUTPUT variables: make an extra copy of the C-result to be returned to R
##       Arguments that are changed in C should be locally declared in this function
##       Arguments that are changed in C must be initialized in C(!), as once the 'last' call to this function
##       set the value of the local variable, our next call from R may be able to reset it, hence
##       no such variable can be an input to the C function
##       Hence the two sets of arguments 'INPUT' and 'OUTPUT'.



minQuad <- function(H,b,C = 1.0,n1=0,n2=0, mem.efficient = FALSE,alpha = NULL,
        lower = NULL,upper = NULL,mat.constr = NULL,lhs.constr = NULL,rhs.constr = NULL,
        control = list(DUP = TRUE,maxit = 1e4, tol = 1e-04, 
            verbose = FALSE, ret.ws = FALSE,ret.data = FALSE,
            rank = 0,
            method = c("default","tron","loqo","exhaustive","x"),
            optim.control = list(),
            q = 2, 
            ws =  c("v", "v2","greedy", "rv2wg", "rvwg", "rv", "rv2")
        )
)   
{   
    if (is.null(b)) stop("b has to be specified")
    n1n2 <- length(b)
    
    if (!length(H)) stop("matrix 'H' has to be specified")
    if (!mem.efficient) {
        if ((nrow(H)!=ncol(H)) || (nrow(H)!=n1n2)) stop("dim of H is not right: "%+%nrow(H)%+%" by "%+%ncol(H)%+%"should be "%+%n1n2%+%" by "%+%n1n2) 
    } else{
        if (is.null(n1) || is.null(n2)) stop("n1 and n2 have to be specified when memory.efficient = TRUE")
        if ((nrow(H)!=ncol(H)) || (nrow(H)!=n1+n2)) stop("dim of H is not right: "%+%nrow(H)%+%" by "%+%ncol(H)%+%"should be "%+%n1+n2%+%" by "%+%n1+n2)
    } 

    # constraint in case of weights
    n.constr = as.integer(0)
    if(length(mat.constr)){
        if(!is.matrix(mat.constr)) mat.constr <- matrix(mat.constr,nrow = 1)
        n.constr <- nrow(mat.constr)
        if(length(lhs.constr)!= n.constr || length(rhs.constr)!= n.constr)stop('Length of left or right constraint vector does not match nrow(constraint matrix)')
        mat.constr <- as.double(mat.constr)
        lhs.constr <- as.double(lhs.constr)
        rhs.constr <- as.double(rhs.constr)
        control$method <- "loqo"
        if(!length(lower)) lower <- double(n1n2)
        if(!length(upper)) upper <- as.double(rep(C,length.out = n1n2))     
    }else{
        mat.constr <- double()
        lhs.constr <- double()
        rhs.constr <- double()
    }    
    control <- do.call('control.minQuad',control)                 
    
  
    
    
    ## Pass matrix H by reference to C, without making any extra copy in memory
    ## H has to be a 'matrix' or a simple vector and its R-numeric-type ('double' or 'integer') has to conform to C's argument 'double *' or
    #   'int
    #   *'.
    ## H is passed as a 'double*' to a vector that is a row_major_order repr. of R's nxn matrix
    ## to conform to C's argument 'double *' or 'int *':
    ## H <- as.double(H): is NOT a 'free' operation even if H is a matrix of R-type 'double'
    ## H <- as.matrix(H): is a 'free' operation as long as is.matrix(H) = TRUE
    ## is.matrix(H) = FALSE: if H is a data.frame an extra copy will be made by as.matrix().
    
    if (!is.double(H))H <- as.double(H) #this check is free but we must make
    
    ## Input arguments: these are unchanged in C and therefore may be initialized in R
    if(length(C) < n1n2) C <- rep(C,length.out = n1n2)
    if (!is.double(C)) 
        C <- as.double(C)
    if (is.null(alpha)) 
        alpha <- 0.5*C   
    if (!is.double(alpha)) 
        alpha <- as.double(alpha) 
    if (!is.double(b)) 
        b <- as.double(b)
    if (!is.integer(n1)) 
        n1 <- as.integer(n1)
    if (!is.integer(n2)) 
        n2 <- as.integer(n2)

            
    ## Output arguments - must be initialized in C because they are changed in C and subsequent calls to this
    # The last call to C-function FIXES(!) the initial values of these for the next call
    # Thus these R-initializations may be overwritten in C for all except that of the first call to minQuad() in C.
    alpha1 <- double(n1n2)
    value <- 0.0 
    epsilon <- 0.0
    convergence <- 0
    iterations <- 0
    ws <- integer(1);if(control$ret.ws) ws <- integer(control$q * control$maxit) 
    
    x <- NULL    
    if((control$method == 'd')){
        if(control$q == 1){
            if(mem.efficient) warning("Memory efficient input for 'H' is not supported for method solveQuad()")
            control$ret.ws <- FALSE
#       Solvequad solves f = a'Qa - b'a = 2{0.5a'Qa-0.5b'a} = 2{0.5a'Qa+c'a} hence pass b_solveQuad = -2*c = -2.0*b_minQuad
#       Also we need to divide the returned function value
            x <- .C("solve_Quad",
            C = C,n1n2 = n1n2, b = -2.0*b, 
            big = as.integer(!mem.efficient),Q = H, 
            n1=n1, n2=n2,K=H, 
            alpha0 = alpha, 
            maxit = as.integer(control$maxit), tol = as.double(control$tol), 
            trace = as.integer(control$verbose),
            alpha = alpha1,value = as.double(value), iterations = as.integer(iterations), convergence = as.integer(convergence), epsilon = as.double(epsilon), 
            DUP = TRUE, NAOK = control$NAOK,PACKAGE = "aucm")
            x$b <- -.5*x$b # corresponds to 'b' in 0.5a'Q' + b'a
            x$value = x$value / 2.0 # corresponds to 0.5a'Q' + b'a
        }else{
            on.exit('Free_minQuad_2')
            x <- .C("minQuad_2", 
            n1n2 = n1n2, b = -2.0*b, C = C,H = H,
            mem_efficient = as.integer(mem.efficient),
            n1=n1, n2=n2,
            maxit = as.integer(control$maxit), 
            tol = as.double(control$tol),machine.double.eps = as.double(.Machine$double.eps), 
            alpha0 = alpha, trace = as.integer(control$verbose),
            as.integer(control$ret.ws),
            # output variables
            alpha = alpha1,value = as.double(value), iterations = as.integer(iterations), convergence = as.integer(convergence), 
            epsilon = as.double(epsilon), 
            ws = ws,
            # .C variables
            DUP = TRUE, NAOK = control$NAOK,PACKAGE = "aucm");
            x$b <- -.5*x$b # corresponds to 'b' in 0.5a'Q' + b'a
            x$value = x$value / 2.0 # <=> 0.5a'Q' + b'a
        }
    }else{
        on.exit('Free_minQuad')
        
        # the order of c("t","l","e","x") must match enumeration in C code {0..3}
        # it is to match the 'case' statements of the 'switch' stmt. in minquad.c  
        method <- as.integer(match(control$method,c("t","l","e","x"),nomatch = 1) - 1)    
        # the order of  c("rvwg","rv2wg","rv","v2","v","greedy") must match enumeration in C code {0..5}
        ws.method <- as.integer(match(control$ws[1],c("rvwg","rv2wg","rv","rv2","v","v2","greedy"),nomatch = 0) - 1)    
        x <- .C("minQuad", 
            n1=n1, n2=n2,
            n1n2 = n1n2, 
            mem_efficient = as.integer(mem.efficient),
            H = H,
            b = b,         
            C = C,
            alpha0 = alpha,
            n.constr = n.constr,
            mat.constr = mat.constr,
            lhs.constr = lhs.constr,
            rhs.constr = rhs.constr,
            rank = as.integer(control$rank),
            s = as.integer(0),#as.integer(control$s),
            q = as.integer(control$q),
            ws.method = ws.method,
            ret.ws = as.integer(control$ret.ws),        
            method = method,
            optim.control = as.double(unlist(control$optim.control)),
            maxit = as.integer(control$maxit), 
            tol = as.double(control$tol),machine.double.eps = as.double(.Machine$double.eps), 
            trace = as.integer(control$verbose),
            # output variables
            alpha = alpha1,value = as.double(value), iterations = as.integer(iterations), convergence = as.integer(convergence), epsilon = as.double(epsilon), 
            ws = ws,
            # .C variables
            DUP = TRUE, NAOK = control$NAOK,PACKAGE = "aucm"
        )
    }   
    t1 <- Sys.time()
    
    ix <- match(c("convergence", "alpha", "value", "iterations", "epsilon"), names(x), nomatch = FALSE)
    ix <- ix[ix != 0]
    
    ## MUST COPY results back for returning to R because the result of 'return(x)'
    ## could be overwritten by the next call to this function
    res <- vector(mode = "list", length = length(ix))
    names(res) <- names(x)[ix]
    for (i in 1:length(ix)) res[[i]] <- x[[ix[i]]]
    
    res$n    <- n1n2
    res$nSV  <- sum(0.0 < res$alpha)
    res$nBSV <- sum(res$alpha == C)
    res$nFSV <- sum((0.0 < res$alpha) & (res$alpha < C))
    
    if(control$ret.data)res<-c(res,(list(H=H,b=b,n1=n1,n2=n2))) #0.5*a'Ha+b'a

        
##  dim = (2 x maxit), no extra memory is used (indexing in C naturally results in this)    
    if(control$ret.ws){
        res$ws <- x$ws;
        dim(res$ws) <- c(control$q,control$maxit)
        res$ws <- res$ws[,1:res$iterations,drop = FALSE]
    }
    res$control <- control
    class(res) = c("minQuad",class(res))
    
    if (res$convergence!=0) warning(paste("convergence failed with code ",res$convergence))
    return(res)
}

print.minQuad=function(x, ...){
    fit=x
    fit$alpha=paste("a numeric vector of size",length(x$alpha))
    class(fit)="list"
    print(fit)
}

# arguments to minQuad() that are passed through rauc()
control.minQuad <- function(
    maxit = 1e4, tol = 1e-04, 
    q = 0,    # size of working set
    ws =  c("v","v2","greedy","rv2wg","rvwg","rv","rv2"),
    method = c("default","tron","loqo","exhaustive","x"),
    optim.control = list(),
    rank = 0,
    DUP = FALSE,
    NAOK = FALSE,
    verbose = FALSE, 
    ret.ws = FALSE,
    ret.data = FALSE
)
{

    default <- list(
        maxit = 1e4, tol = 1e-04, 
        q = 0,
        ws = "v",
        method = "d",
        optim.control = list(maxfev = 1000,fatol = tol,frtol = tol,cgtol=tol,gtol=tol,fmin = -.Machine$double.xmax),
        rank = 0,
        DUP = FALSE,
        NAOK = FALSE,
        verbose = FALSE, 
        ret.ws = FALSE,
        ret.data = FALSE)

    default.optim.control <- list(
        t = list(maxfev = 1000,fatol = tol,frtol = tol,cgtol=tol,gtol=tol,fmin = -.Machine$double.xmax),
        l = list(bound = 10,margin=0.05,maxiter=40,sigfig = 7,inf = 1e6)
    )

    default$optim.control <- default.optim.control[default$method] 
    nmd <- names(default)
    
    control <- list(maxit = maxit, tol = tol, q = q,rank = rank,
                ws = ws,method = method,optim.control = optim.control,
                DUP = DUP,NAOK = NAOK, verbose = verbose, ret.ws = ret.ws,ret.data = ret.data)
                
    if(!length(control$method)) stop("Invalid method parameter ",control$method)
    control$method <- substr(control$method[1],1,1)
    if(!any(control$method == c("e","t","d","x","l")))
        stop("Invalid method parameter ",control$method)

    if(any(control$method == names(default.optim.control))){
        if(!length(control$optim.control)){
            control$optim.control <- default.optim.control[[control$method]]
        }else{
            default.control.optim.control <- default.optim.control[[control$method]]
            nmdcoc <- names(default.control.optim.control)
            control$optim.control <- control$optim[match(nmdcoc,names(control$optim.control),nomatch = 0)]
            unset.names <- setdiff(nmdcoc,names(control$optim.control))
            control$optim.control[unset.names] <- default.control.optim.control[unset.names]    
        }
        control$optim.control <- control$optim.control[order(names(control$optim.control))]
    }else{
        control$optim.control <- NULL
    }
 

#   if(!control$s) control$s <- control$maxit + 1 #so that (iter %% s == 0) is never TRUE
    if(control$rank < 0) control$rank <- 0
      
# if method is 'tron' then q = 0 in C <-> dynamic
# if method is any other then q = 2 by default  
    control$q <- as.integer(control$q)      
    if(control$q < 0) 
        stop("Size of working set 'q' should be a positive integer or possibly zero for method 'tron'.")
    if((q == 0) && (control$method != 't')) q <- 2
    
# rename "working.set" to "ws"
    # wx <- which(names(control)=="working.set")
    # if(length(wx)) names(control)[wx] <- "ws"
    
    if(!length(control$ws))
         control$ws <- default$ws
    control$ws <- control$ws[1]
    if(!any(control$ws ==  c("rvwg","rv2wg","rv","rv2","v","v2","greedy")))
        stop("Invalid working set parameter ",control$ws)

    control
}

# free memory allocated in C
Free_minQuad_2 <- function(){.C("free_minQuad_2",package = "aucm")}
Free_minQuad <- function(){.C("free_minQuad",package = "aucm")}


# objfun <- function(H,b,a,mem.efficient = FALSE,n1=0,n2=0)
# {
    # .C('R_objective',as.integer(n1),as.integer(n2),as.integer(ifelse(mem.efficient,n1*n2,nrow(H)))
    # ,as.integer(mem.efficient),H,a,b,result = double(1),package = "aucm")$result
# }
#objfun(fit2$K,fit2$b,fit2$alpha,mem.efficient = TRUE,n1 = fit2$n1,n2=fit2$n2)

    
    
