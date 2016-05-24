setClass("hyperdirichlet",
         representation = representation("numeric" , NC="numeric", pnames="character", validated="logical"),
         prototype      = prototype(NC=NA_real_ , pnames=character() , validated=FALSE)
         )

setGeneric("NC",function(x){standardGeneric("NC")})
setMethod("NC","hyperdirichlet",
          function(x){x@NC}
          )

setGeneric("params",function(x){standardGeneric("params")})
setMethod("params","hyperdirichlet",
          function(x){x@.Data}
          )

setGeneric("pnames",function(x){standardGeneric("pnames")})
setMethod("pnames","hyperdirichlet",
          function(x){x@pnames}
          )

setGeneric("validated",function(x){standardGeneric("validated")})
setMethod("validated","hyperdirichlet",
          function(x){x@validated}
          )

                                        # There are no occurrences of
                                        # "@" below this line.

setGeneric("pnames<-",function(x,value){standardGeneric("pnames<-")})
setMethod("pnames<-","hyperdirichlet",
          function(x,value){
            hyperdirichlet(params(x),NC=NC(x),pnames=value,validated=TRUE)
          }
          )



# dim() is a primitive function so calling setGeneric("dim") is
# currently unnecessary.

setMethod("dim","hyperdirichlet",
          function(x){round(log(length(params(x)))/log(2))}
          )

".mean_hyperdirichlet" <- function(x, normalize=TRUE, ...){
  f <- function(i){
    jj <- rep(0,dim(x))
    jj[i] <- 1
    return(mgf(x,jj, ...))
  }
  out <- sapply(seq_len(dim(x)),f)
  if(normalize){
    out <- out/sum(out)
  }
  return(out)
}
  
setGeneric("mean")
setMethod("mean" , "hyperdirichlet", .mean_hyperdirichlet)

".hd_add" <- function(e1,e2){
  stopifnot(is.hyperdirichlet(e1) & is.hyperdirichlet(e2))
  stopifnot(dim(e1) == dim(e2))
  if( .Generic == "+"){
    return(hd_add(e1, e2, assume_validated=FALSE))
  } else {
    return(.hd_binary_op_error(e1,e2))
  }
}

"hd_add" <- function(e1, e2, assume_validated = FALSE){
  stopifnot(is.hyperdirichlet(e1) & is.hyperdirichlet(e2))
  n <- dim(e1)
  stopifnot(n == dim(e2))
  e1z <- length(pnames(e1))==0
  e2z <- length(pnames(e2))==0
  if(e1z & !e2z){
    jj.p <- pnames(e2)
  } else if (!e1z & e2z){
    jj.p <- pnames(e1)
  } else if (!e1z & !e2z){
    jj.p <- pnames(e1)
  } else {
    jj.p <- character(0)
  }
  if(assume_validated){
    return(hyperdirichlet(powers(e1) + params(e2) , pnames = jj.p, validated = TRUE))
  } else {
    return(hyperdirichlet(powers(e1) + params(e2) , pnames = jj.p, validated = validated(e1) & validated(e2)))
  }
}

".hd_binary_op_error" <- function(e1,e2){
  stop("The only binary operations allowed on hyperdirichlet objects are:  to add two of the same dimensions, as in 'a+b'; and to multiply a hyperdirichlet object by a (positive) numeric scalar")
}

".hd_no_unary" <- function(e1,e2){stop("No unary operators are defined on hyperdirichlet objects")}

".hd_prod_scalar" <- function(e1,e2){
  if(! (.Generic == "*")){
    return(.hd_binary_op_error(e1,e2))
  } 
  if(is.hyperdirichlet(e1) & is.numeric(e2)){
    return(.hd_prod_scalar_worker(e1,e2))
  } else if(is.numeric(e1) & is.hyperdirichlet(e2)){
    return(.hd_prod_scalar_worker(e2,e1))  # arguments e1 and e2 transposed
  } else {
    stop("this cannot happen")
  }
}

".hd_prod_scalar_worker" <- function(e1,e2){ # this function does the work
  if( !(is.hyperdirichlet(e1) & is.numeric(e2))){
    stop("e1 must be hyperdirichlet and e2 must be numeric")
  }
  if(length(e2)>1){
    stop("numeric operand must be length one")
  }
  if(e2<0){
    stop("numeric operand must not be negative")
  }
  hyperdirichlet(x=params(e1)+powers(e1)*(e2-1) , NC=NA, pnames=pnames(e1), validated=validated(e1))
}

setMethod("Arith",signature(e1 = "hyperdirichlet", e2="missing"), .hd_no_unary )
setMethod("Arith",signature(e1 = "hyperdirichlet" , e2="hyperdirichlet") , .hd_add)
      
setMethod("Arith",signature(e1 = "hyperdirichlet" , e2="numeric")   , .hd_prod_scalar)
setMethod("Arith",signature(e1 = "numeric" , e2="hyperdirichlet")   , .hd_prod_scalar)

setMethod("Arith",signature(e1 = "hyperdirichlet" , e2="ANY")   , .hd_binary_op_error)
setMethod("Arith",signature(e1 = "ANY" , e2="hyperdirichlet")   , .hd_binary_op_error)

"is.hyperdirichlet" <- function(x){is(x,"hyperdirichlet")}

".hd_valid" <- function(object){

  x <- params(object)
  normc <- NC(object)
  pn <- pnames(object)

  if(validated(object)){return(TRUE)}

  if(!is.numeric(x)){
    return("x not numeric")
  } else if (!  (is.numeric(normc))){
    return("NC must be numeric or NA")
  } else if (length(normc) != 1){
    return("NC must be length 1")
  } else if (abs(2^dim(object)-length(x)) > 0.5){
    return("params must be a vector of length 2^n for integer n")
  } else if ( (length(pn) > 0)  &  (length(pn) != dim(object))  ){
    return(paste("pnames not correct length---should be length ",dim(object)))
  } else if (!is.proper(object)){
    return("object unnormalizable")
  } else {
    return(TRUE)
  }
}

setValidity("hyperdirichlet", .hd_valid)

setGeneric("as.hyperdirichlet",
           function(x, calculate_NC=FALSE, ...){standardGeneric("as.hyperdirichlet")}
           )

".as_hyperdirichlet" <- function(x, calculate_NC, validated=FALSE,  ...){
  if(is.hyperdirichlet(x)){
    if(!is.na(NC(x))){
      return(x)
    }
    pn <- pnames(x)
  } else {
    pn <- NULL
  }
  if(calculate_NC){
    return(hyperdirichlet(x, NC = calculate_B(x, ...), pnames=pn, validated=validated))
  } else {
    return(hyperdirichlet(x, NC = NA_real_, pnames=pn, validated=validated))
  }
}

"matrix_to_HD" <- function(x, calculate_NC = FALSE, bernoulli = NULL, ...){
    if(is.null(bernoulli)){
      jj <- x[!is.na(x)]
      if(all(  (jj==0)|(jj==1) )){
        bernoulli <- TRUE
      } else {
        bernoulli <- FALSE
      }
    }
    if(bernoulli){
      return(bernoulli_matrix_to_HD(x, calculate_NC = calculate_NC, ...))
    } else {
      return(multinomial_matrix_to_HD(x, calculate_NC = calculate_NC, ...))
    }
  }

"bernoulli_matrix_to_HD" <- function(x, calculate_NC = FALSE, ...){
  stopifnot(is.matrix(x))
  jj <- x[!is.na(x)]
  stopifnot(all((jj==0) | (jj==1)))
  d <- ncol(x)
  out <- uniform(d)
  for(i in seq_len(nrow(x))){
    set <- x[i,,drop=TRUE]

    jjj <- set[!is.na(set)]
    if(!any(jj == 0)){
      warning(paste("row ", i,"does not have  losers"))
    }
    if(!any(jj == 1)){
      warning(paste("row ", i,"does not have winners"))
    }
    
    out <- hd_add(out , single_bernoulli_obs(d, win=which(set==1), lose=which(set==0)), assume_validated=TRUE)
  }
  out <- as.hyperdirichlet(x=out, calculate_NC=calculate_NC, validated=TRUE, ...)
  pnames(out) <- colnames(x)
  return(out)
}

"multinomial_matrix_to_HD" <- function(x, calculate_NC = FALSE, ...){
  stopifnot(is.matrix(x))
  n <- ncol(x)
  out <- uniform(n)
  for(i in seq_len(nrow(x))){
    jj <- x[i,,drop=TRUE]
    wanted <- !is.na(jj)
    stopifnot(all(jj[wanted] >= 0))
    out <- hd_add(out , mult_restricted_obs(n, which(wanted), jj[wanted]), assume_validated = TRUE)
  }
  out <- as.hyperdirichlet(out, calculate_NC = calculate_NC, ...)
  pnames(out) <- colnames(x)
  return(out)
}

setMethod("as.hyperdirichlet","hyperdirichlet", .as_hyperdirichlet )
setMethod("as.hyperdirichlet","numeric", .as_hyperdirichlet)
setMethod("as.hyperdirichlet","matrix", matrix_to_HD)

setAs("numeric", "hyperdirichlet",
      function(from){ hyperdirichlet(from)}
      )
setAs("matrix" , "hyperdirichlet",
      function(from){ matrix_to_HD(from)}
      )

"hyperdirichlet" <- function(x, NC, pnames=character(), validated=FALSE){
  if(length(x)==1){ return(uniform(x))}
  if(missing(NC)){NC <- NA_real_}
  if(is.na(NC)){ NC <- NA_real_  }
  if(is.null(pnames)){ pnames <- character()}
  new("hyperdirichlet" , x , NC=NC , pnames=pnames, validated=validated)  # This is the only time new() is called
}

".head_hyperdirichlet" <- function(x, n=6, ...){
  .print_hyperdirichlet(x, n=n, do_head = TRUE, ...)
}

".tail_hyperdirichlet" <- function(x, n=6, ...){
  .print_hyperdirichlet(x, n=n, do_head = FALSE, ...)
}

setGeneric("head")
setGeneric("tail")

setMethod("head" , signature="hyperdirichlet" , .head_hyperdirichlet)
setMethod("tail" , signature="hyperdirichlet" , .tail_hyperdirichlet)

".print_hyperdirichlet_worker" <- 
function(x, n=0, do_head = NA, ...){ # This function does the work
  out <- cbind(binmat(dim(x),pnames=pnames(x)),params=params(x),powers=powers(x))
  rownames(out) <- paste("[", seq_len(nrow(out)),"]",sep="")
  if(n==0){
    return(out)
  } else {
    if(do_head){
      return(head(out, n=n, ...))
    } else {
      return(tail(out, n=n, ...))
    }
  }
}
    
".print_hyperdirichlet" <- function(x, n=0, do_head = TRUE, ...){
  jj <- .print_hyperdirichlet_worker(x, n=n, do_head=do_head, ...)
  print(jj)
  cat("\n")
  if(is.na(NC(x))){
    cat("Normalizing constant not known")
    cat("\n")
  } else {
    cat(paste("Normalizing constant: ", NC(x),"\n"))
  }
  return(invisible(jj))
}

setMethod("show", "hyperdirichlet",
          function(object){.print_hyperdirichlet(object)}
          )

".get_logical" <- function(x){   
  sum((2^((length(x)-1):0))[x])
}

setMethod("[", "hyperdirichlet",
          function(x, i, j, drop){
            x <- params(x)
            if(!missing(j)){
              warning("second argument to extractor function ignored")
            }
            
            if(is.logical(i)){
              if(is.matrix(i)){
                return(x[1 + apply(i,1, .get_logical)])
              } else {
                return(x[1 + .get_logical(i)])
              }
            } else {
              return(x[i])
            }
          }
          )


setReplaceMethod("[", signature(x="hyperdirichlet"),
                 function(x,i,j, ...,  recalculate=FALSE, validated = FALSE, value){
                   jjn <- NC(x)
                   jjx <- params(x)
                   if(!missing(j)){
                     warning("second argument to extractor function ignored")
                   }
                   if(is.logical(i)){
                     if(is.matrix(i)){
                       jjx[1 + apply(i,1, .get_logical)] <- value
                     } else {
                       jjx[1 + .get_logical(i)] <- value
                     }
                   } else {
                     jjx[i] <- value
                   }
                   if(recalculate){
                     return(hyperdirichlet(jjx, NC=calculate_B(jjx, ...), pnames=pnames(x), validated = TRUE))
                   } else {
                     return(hyperdirichlet(jjx, NC=NA, pnames=pnames(x), validated = validated))
                   }
                 }
                 )
                 
"calculate_B" <- function(x, disallowed=NULL, give=FALSE, ...){
  if(!is.hyperdirichlet(x)){
    x <- as.hyperdirichlet(x, calculate_NC=FALSE)
  }
  n <- dim(x)-1  # "-1" because this is the dimension of the integrand, not the number of p's.

  ## First, special dispensation for the Dirichlet:
  if(is.null(disallowed) & is.dirichlet(x)){
    return(diri_norm(params(x)[.pow2(dim(x))]))
  }
  
  if(is.null(disallowed)){
    f <- function(e){
      dhyperdirichlet_e(e,HD=x, include.Jacobian=TRUE)
    }
  } else {
    f <- function(e){
      jj <- e_to_p(c(1,e))
      return(
             as.numeric(!disallowed(jj)) *
             dhyperdirichlet_e(e,HD=x, include.Jacobian=TRUE)
             )
    }
  }
  out <- adaptIntegrate(f,lowerLimit=rep(0,n),upperLimit=rep(1,n), ...)
  if(give){
    return(out)
  } else {
    return(out$integral)
  }
}

"B" <- function(x, ...){
  stopifnot(is.hyperdirichlet(x))
  out <- NC(x)
  if(is.na(out)){
    return(calculate_B(x, ...))
  } else {
    return(out)
  }
}

"Jacobian" <-
function(e){
  n <- length(e)
  e1 <- e[1]
  e <- e[-c(1,length(e))]
  prod(cumprod(1-e))*e1^n
}

"binmat" <-
function(n, alternatives=NULL, pnames=NULL){
  if(is.null(alternatives)){
    alternatives <- 0:1
  }
  out <-  as.matrix(expand.grid(rep(list(alternatives),n)))
  out <- out[,rev(seq_len(n)),drop=FALSE]
  if(length(pnames) == 0){
    colnames(out) <- paste("p",seq_len(n),sep="")
  } else {
    colnames(out) <- pnames
  }
  return(out)
}

"dhyperdirichlet_e" <-
function(e, HD, include.Jacobian=TRUE){

  stopifnot(is.hyperdirichlet(HD))
  parameters <- params(HD)
  e <- c(1,e)
  p <- e_to_p(e)
  out <- dhyperdirichlet(p,parameters)
  if(include.Jacobian){
    return(Jacobian(e)*out)
  } else {
    return(out)
  }
}

"dhyperdirichlet" <-
function(p, HD, include.NC = FALSE, TINY = 1e-10, log=FALSE){
  if(is.matrix(p)){
    return(apply(p,1, match.fun(sys.call()[[1]]),
                 HD=HD, include.NC=include.NC, TINY=TINY, log=log))
  }
  
  HD <- as.hyperdirichlet(HD, validated=TRUE)

  parameters <- params(HD)
  n <- length(p)
  stopifnot(n == dim(HD))
  p <- pmax(p, TINY)
  p <- p/sum(p)
  jj <- binmat(n,c(FALSE,TRUE))
  pp <- apply(jj,1,function(x){sum(p[x])})

  out <- sum(log(pp[-1])*parameters[-1]) - sum(log(p))


  jj.nc <- NC(HD)
  if(log){ # Return log of prob
    if(include.NC){
      return(out - log(jj.nc))
    } else {
      return(out)
    }
  } else {  # Return prob (not log of prob)
    out <- exp(out)
    if(include.NC){
      return(out/jj.nc)
    } else {
      return(out)
    }
  }
}

"dirichlet" <-
function(params, powers, pnames){
  if ( !xor(missing(params),missing(powers))){
    stop("supply exactly one of params or powers")
  }

  if(missing(params)){
    x <- powers+1
  } else {
    x <- params
  }
  
  if(missing(pnames)){
    pnamesx <- names(x)
  } else {
    pnamesx <- pnames
  } 
    
  if(any(x<0)){
    stop("parameters must be positive; powers must be >1")
  }
  
  lx <- length(x)
  out <- rep(0,2^lx)
  out[.pow2(lx)] <- rev(x)
  return(hyperdirichlet(out, NC=diri_norm(x),validated=TRUE, pnames=pnamesx))
}

".pow2" <- function(n){
  1+2^(0:(n-1))
}

"diri_norm" <- function(x){
  ##  prod(gamma(x))/gamma(sum(x))
  exp(sum(lgamma(x))-lgamma(sum(x)))
}

"gd_norm" <- function(a,b){
  ## prod(beta(a,b))
  exp(sum(lbeta(a,b)))
}

"is.dirichlet" <- function(x){
  p <- params(x)
  if(all(p[-.pow2(dim(x))] == 0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

"dirichlet_params" <- function(x){
  out <- params(x)[rev(.pow2(dim(x)))]
  jj <- pnames(x)
  if(length(jj) > 0){
    names(out) <- pnames(x)
  }
  return(out)
}

"dirichlet_params<-" <- function(x,value){
  jj <- params(x)
  jj[rev(.pow2(dim(x)))] <- value
  return(hyperdirichlet(jj , pnames=pnames(x)))
}

"uniform" <-
function(n, pnames=NULL){
  jj <- rep(1,n)
  names(jj) <- pnames
  return(dirichlet(jj))
}

"single_obs" <-
function(d,n){
  jj <- rep(0,d)
  jj[n] <- 1
  return(obs(jj))
}

"obs" <- function(x){
  dirichlet(x+1)
}

"single_multi_restricted_obs" <- function(d,n,x){  # smro(4, 2, 1:3) means 4 players, winner number 2, restricted to 1,2,3
  out <- hyperdirichlet(d)
  out[seq_len(d) %in% x] <- -1
  return(out + single_obs(d,n))
}

"mult_restricted_obs" <- function(d, a, nobs){  #mro(5, a=c(2,4), nobs=c(10,20)) means 5 players, 30 games between p2 and p4 of which p2 won 10 and p4 won 20.
  stopifnot(length(a) == length(nobs))
  out <- params(uniform(d))
  out[1+2^(d-a)] <- nobs+1
  jj <- seq_len(d) %in% a
  out[1+sum(2^(which(rev(jj))-1))]  <- -sum(nobs)
  return(hyperdirichlet(out, validated = TRUE))
}
  
"single_bernoulli_obs" <- function(d,win,lose){   #sbo(5, c(1,2),3) means 5 players, 1&2 vs 3, and 1&2 won.
  out <- params(uniform(d))
  if(any(intersect(win,lose))){
    stop("'win' and 'lose' must have empty intersection")
  }

  winners <- seq_len(d) %in% win
  players <- seq_len(d) %in% union(win,lose)

  jj <- 1+sum(2^(which(rev(winners))-1))
  out[jj] <- out[jj] + 1

  jj2 <- 1+sum(2^(which(rev(players))-1))
  out[jj2] <- -1

  return(hyperdirichlet(out, validated=TRUE))
}

"mult_bernoulli_obs" <- function(d,team1,team2,wins1,wins2){ #mbo(5, 1:2,3:4, 10,11) means 5 players, 1&2 vs 3&4, 21 (=10+11) games, with 1&2 winning 10 and 3&4 winning 11.
  out <- uniform(d)
  if(any(intersect(team1,team2))){
    stop("team1 and team2 must have empty intersection")
  }
  out[seq_len(d) %in% team1] <- out[seq_len(d) %in% team1] + wins1
  out[seq_len(d) %in% team2] <- out[seq_len(d) %in% team2] + wins2
  out[seq_len(d) %in% union(team1,team2)] <-
  out[seq_len(d) %in% union(team1,team2)] - (wins1 + wins2)
  return(out)
}

"bernoulli_obs" <- function(d, winners, losers){
  stopifnot(length(winners) == length(losers))
  out <- uniform(d)
  for(i in seq_along(winners)){
    out <- out + single_bernoulli_obs(d, winners[[i]], losers[[i]])
  }
  return(out)
}

"p_to_e" <-
function(p){
  c(
    sum(p),
    p[-length(p)]/rev(cumsum(rev(p))[-1])
    )
}

"e_to_p" <-
function(e){
  e[1]*
    cumprod(c(1,1-e[-1]))*
      c(e[-1],1)
}

"gd" <-
function(a, b, b0=0, pnames=NULL){
  stopifnot(length(a) == length(b))
  k <- length(a)+1
  
  out <- numeric(2^k)
  for(i in 1:(k-1)){
    out[1+2^i] <- a[k-i]   # ie {p_i}^{a_i-1}  (sic: 'out' holds the parameters, not the powers)
  }
  
  out[2] <- b[k-1]         # ie {p_k}^{b_{k-1}-1}

  for(i in seq(from=2, to=k-1)){
    out[2^i] <- b[k-i] - (a[k-i+1] + b[k-i+1])    # ie (sum_j=i^k p_j)^{b_{i-1}=(a_i+b_i)}

  }
  out[2^k] <- b0 - (a[1]+b[1])
  return(hyperdirichlet(out, gd_norm(a,b), validated=TRUE,pnames=pnames))
}

"is.proper" <-
function(x, irregardless=FALSE){
  if(!irregardless){
    if(validated(x)){return(TRUE)}
  }
  if(is.hyperdirichlet(x)){x <- params(x)}
  n <- round(log(length(x))/log(2))
  jj <- binmat(n)
  for(i in seq_len(2^n-1)){
    ref <- as.logical(jj[i,,drop=TRUE])
    wanted <- apply(jj,1,function(a){identical(ref,ref|as.logical(a))})
    if(sum(x[wanted]) < 0){
      print(i)
      return(FALSE)
    }
  }
  return(TRUE)
}

"justpairs" <-
function(x){
  stopifnot(is.matrix(x))
  stopifnot(nrow(x)==ncol(x))
  
  n <- nrow(x)
  f <- function(x){sum(2^(n-x))}
  
  out <- numeric(2^n)
  
  two1s <- which(upper.tri(x),arr.ind=TRUE)
  two0s <- which(lower.tri(x),arr.ind=TRUE)
  
  out[1 + apply(two1s,1,f)] <- x[two1s]
  
  if(n==4){
    warning("lower triangular entries discarded")
  } else {
    out[2^n - apply(two0s,1,f)] <- x[two0s]
  }
  return(dirichlet(rep(2,nrow(x))) + hyperdirichlet(out))
}

"paircomp" <- function(x){
  diag(x) <- 0
  jj <- x+t(x)
  jj[lower.tri(jj)] <- 0
  
  out <-
    hyperdirichlet(params(dirichlet(powers=colSums(x))) + params(uniform(nrow(x))) - params(justpairs(jj)))

  pnames(out) <- rownames(x)
  return(out)
}

"mgf" <-
function(x, powers, ...){
    B(x+dirichlet(powers=powers), ...) / B(x, ...)
  }

"maximum_likelihood"<- function(HD, start_p = NULL, give=FALSE, disallowed = NULL, zero=NULL, ...){
  if(is.null(zero)){
    return(mle(HD=HD, start_p=start_p , give=give , disallowed=disallowed , ...))
  } else {
    return(mle_restricted(HD=HD , start_p=start_p , give=give, disallowed=disallowed , zero=zero, ...))
  }
}

"mle" <- function(HD, start_p = NULL, give=FALSE, disallowed = NULL, ...){
  stopifnot(is.hyperdirichlet(HD))
  parameters <- params(HD)
  
  n <- dim(HD)
  if(is.null(start_p)){
    start_p <- rep(1/n,n)
  }
  
  if(length(start_p) == n){
    start_p <- start_p/sum(start_p)
    start_p <- start_p[-(n-1)]
  } 
  
  f <- function(pm1){
    # for "pm1" read "p[-1]" or "p of minus 1"
    p <- c(pm1,1-sum(pm1))
    if(any(p<0) | sum(p)>1){ return(Inf) }
    if(!is.null(disallowed)){
      if(disallowed(p)){ return(Inf) }
    }
    
    out <- -dhyperdirichlet(p,HD=parameters, include.NC = FALSE , log = TRUE)
    return(out)
  }


  ## Now chop off the last element of start_p so that the optimization
  ## can work with independent variables:

  start_p <- start_p[-n]
  
  out <- optim(par=start_p, fn=f, ...)
  if(give){
    return(out)
  } else {
    jj <- out$par
    jj <- c(jj,1-sum(jj))

    if(length(pnames(HD))==0){
      names(jj) <- paste("p",seq_len(dim(HD)),sep="")
    } else {
      names(jj) <- pnames(HD)
    }
    
    support <- dhyperdirichlet(p=jj , HD = HD ,log=TRUE)
    return(list(
                MLE        = jj,
                likelihood = exp(support),
                support    = support
                ))
  }
}

"mle_restricted" <- function(HD, start_p = NULL, give=FALSE, disallowed = NULL,  zero=NULL, ...){
  parameters <- params(as.hyperdirichlet(HD, FALSE))
  
  n <- dim(HD)
  if(is.null(start_p)){
    start_p <- rep(1/n,n)
  }
  
  if(length(start_p) == n){
    start_p <- start_p/sum(start_p)
    start_p <- start_p[-(n-1)]
  } 

  if(!is.logical(zero)){
    zero <- seq_len(n) %in% zero
  }
  nonzero <- !zero
  m <- max(which(nonzero))
  l <- sum(nonzero)
  
  f <- function(pshort){
    p <- rep(0,n)
    p[nonzero][-l] <- pshort
    p[nonzero][l] <- 1-sum(p)
    if(any(p[nonzero]<0) | sum(p[nonzero])>1){ return(Inf) }
    if(!is.null(disallowed)){
      if(disallowed(p)){ return(Inf) }
    }
    
    out <- -dhyperdirichlet(p,HD=parameters, include.NC = FALSE , log = TRUE)
    return(out)
  }
    
  ## Now chop off the last element of start_p so that the optimization
  ## can work with independent variables:
  start_p <- start_p[nonzero][-l]
  out <- optim(par=start_p, fn=f, ...)
  
  
  if(give){
    return(out)
  } else {
    jj <- rep(0,n)
    jj[nonzero][-l] <- out$par
    jj[nonzero][l] <- 1-sum(jj)

    if(length(pnames(HD))==0){
      names(jj) <- paste("p",seq_len(dim(HD)),sep="")
    } else {
      names(jj) <- pnames(HD)
    }

    support <- dhyperdirichlet(p=jj , HD = HD ,log=TRUE)
    return(list(
                MLE        = jj,
                likelihood = exp(support),
                support    = support
                ))
  }
}

"powers" <- function(x){
  stopifnot(is.hyperdirichlet(x))
  out <- params(x)
  jj <- 1+2^(seq_len(dim(x))-1)
  out[jj] <- out[jj]-1
  return(out)
}

"triplot" <- function(HD, l=100,  do_image = TRUE, do_contour = TRUE, discard=0.05, labels=NULL, ...){
  maxy <- sin(pi/3)
  x <- seq(0,1,len=l)
  y <- seq(0, maxy, len=l)

  p1 <- outer(x,y, function(x,y) x - y/maxy/2)
  p2 <- outer(x,y, function(x,y) y/maxy)
  p3 <- 1-p1-p2

  p <- abind(p1,p2,p3,along=3)

  if(is.hyperdirichlet(HD)){
    stopifnot(dim(HD) == 3)
    f <- function(p){ ifelse(0 < p[1] & 0 < p[3], dhyperdirichlet(p=p , HD = HD, log=TRUE), NA) }
    if(is.null(labels)){labels <- pnames(HD)}
  } else {
    f <- function(p){ ifelse(0 < p[1] & 0 < p[3], HD(p), NA) }
  }
  
  z <- apply(p,1:2,f)

  jj <- z[!is.na(z)]
  jjmin <- quantile(jj,discard)
  z[ !is.na(z) & z<jjmin] <- jjmin

  if(do_image){
    image(z, asp=sqrt(3)/2, axes=FALSE, ...)
  }
  if(do_contour){
    contour(z,add=TRUE, ...)
  }

  segments(0   , 0 , 0.5 , 1 , lwd=4)
  segments(0.5 , 1 , 1   , 0 , lwd=4)
  segments(1   , 0 , 0   , 0 , lwd=4)

    if(length(labels) == 3){
    
      jj <- 0.05
      d1 <- jj*cos(pi/6)
      d2 <- jj*sin(pi/6)
      
      text(0-d1  ,0-d2,labels[1])
      text(0.5,   1+jj,labels[2])
      text(1+d1  ,0-d2,labels[3])
  }
    
    return(invisible(z)) 
}

"rhyperdirichlet" <- 
function (n, HD, start=NULL, sigma=NULL) {
    d <- dim(HD)
    if(is.null(sigma)){sigma <- 1/d}
    out <- matrix(NA, n, d)
    if(is.null(start)){
      start <- rep(1,d)
    }
    start <- start/sum(start)
    
    out[1, ] <- start
    for (i in 2:n) {
      ## following lines suggested by Simon Byrne via gmail, 30 Oct 2012
      z <- rmvnorm(n = 1, sigma = sigma * diag(nrow = d))
      z <- z - mean(z)
      proposed <- out[i - 1, ] + z
      if(all(proposed>0)){
        num <- dhyperdirichlet(proposed  , HD=HD , log=FALSE)
        den <- dhyperdirichlet(out[i-1,] , HD=HD , log=FALSE)
        if ((num == 0) & (den == 0)) {
          print("this cannot happen")
          alpha <- 0
        } else {
          alpha <- min(1, num/den)
        }
        if (runif(1) < alpha) {
          out[i, ] <- proposed
        } else {
          out[i, ] <- out[i - 1, ]
        }
      } else {
        out[i, ] <- out[i - 1, ]
      }
    }
    return(out)
  }
"probability" <- function(x , disallowed, ...){
  return(calculate_B(x , disallowed=disallowed , give=FALSE, ...) / B(x, ...))
}


"fitz" <- function(dat, include.missing=TRUE, validated=NULL){

  stopifnot(length(dat) == 26)

  if(is.null(validated)){
    if(any(dat<0)){
      warning("some data negative.  Setting validated to FALSE")
      validated <- FALSE
    } else {
      validated <- TRUE
    }
  }
  
  n <- 3
  b <- uniform(2^n)
  pnames(b) <- 
    apply(as.matrix(expand.grid(rep(list(c("T","F")),n)))[,rev(seq_len(n))],1,paste,collapse="")
  
  f <- function(...){
    out <- rep(FALSE,2^n)
    out[c(list(...),recursive=TRUE)] <- TRUE
    return(out)
  }
  
  b <- "[<-"(b,f(1), value = dat[1] + 1, validated=validated)
  b <- "[<-"(b,f(2), value = dat[2] + 1, validated=validated)
  b <- "[<-"(b,f(3), value = dat[3] + 1, validated=validated)
  b <- "[<-"(b,f(4), value = dat[4] + 1, validated=validated)
  b <- "[<-"(b,f(5), value = dat[5] + 1, validated=validated)
  b <- "[<-"(b,f(6), value = dat[6] + 1, validated=validated)
  b <- "[<-"(b,f(7), value = dat[7] + 1, validated=validated)
  b <- "[<-"(b,f(8), value = dat[8] + 1, validated=validated)
  
  if(include.missing){
    b <- "[<-"(b,f(1,5), value = dat[09], validated=validated)
    b <- "[<-"(b,f(2,6), value = dat[10], validated=validated)
    b <- "[<-"(b,f(3,7), value = dat[11], validated=validated)
    b <- "[<-"(b,f(4,8), value = dat[12], validated=validated)
    b <- "[<-"(b,f(1,3), value = dat[13], validated=validated)
    b <- "[<-"(b,f(2,4), value = dat[14], validated=validated)
    b <- "[<-"(b,f(5,7), value = dat[15], validated=validated)
    b <- "[<-"(b,f(6,8), value = dat[16], validated=validated)
    b <- "[<-"(b,f(1,2), value = dat[17], validated=validated)
    b <- "[<-"(b,f(3,4), value = dat[18], validated=validated)
    b <- "[<-"(b,f(5,6), value = dat[19], validated=validated)
    b <- "[<-"(b,f(7,8), value = dat[20], validated=validated)
    
    b <- "[<-"(b,f(1,3,5,7), value = dat[21], validated=validated)
    b <- "[<-"(b,f(2,4,6,8), value = dat[22], validated=validated)
    b <- "[<-"(b,f(1,2,5,6), value = dat[23], validated=validated)
    b <- "[<-"(b,f(3,4,7,8), value = dat[24], validated=validated)
    b <- "[<-"(b,f(1,2,3,4), value = dat[25], validated=validated)
    b <- "[<-"(b,f(5,6,7,8), value = dat[26], validated=validated)
  }
  
  return(b)
}

"maxmult" <- function(M, start_a=NULL, give=FALSE, method="nlm", ...){

  jj <- apply(M,1,sum)
  stopifnot(max(jj)==min(jj))

  
  if(is.null(start_a)){
    k <- ncol(M)
    n <- mean(jj)
    phat <- (colSums(M)/sum(M))[-k] # sum(phat) != 1
    Sigma_m <- -n*outer(phat,phat)
    diag(Sigma_m) <- n*phat*(1-phat) # equation 17
    Sigma_m <- Sigma_m
    Sigma_cm <- cov(M[,-k])
    C_hat <- (det(Sigma_cm)/det(Sigma_m))^(1/(k-1))  # equation 25

    start_a <- colMeans(M) * (n-C_hat)/(n*(C_hat-1))
    names(start_a) <- colnames(M)
  }
  
  func <- function(x,l){ # notation "l" from Mosimann 1962
    ifelse(any(l<0),Inf,
           lfactorial(sum(x)) -sum(lfactorial(x))
              +lgamma(sum(l)) +sum(lgamma  (x+l))
              -sum(lgamma(l)) -lgamma(sum  (x+l))
           )
  }
  
  f <- function(l){ -sum(apply(M, 1, FUN=func,l=l))}
  
  switch(method,
         "nlm" = {
           jj <- nlm(f, start_a, ...)
           if(!give){ # awkard test used for consistency with 'give' in mle()
             out <- jj$estimate
             names(out) <- colnames(M)
             return(out)
           }
         },
         "optim" = {
           jj <- optim(start_a, f, ...)
           if(!give){
             out <- jj$par
             return(out)
           }
         },
         { return(start_a)
         }
         )  
  return(jj)  # at this point 'give' *must* be TRUE
}

  
