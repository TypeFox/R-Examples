.DBS <-
function( y , m , groups=NULL )
           {
             if( is.null(groups) ) groups = rep(1,ncol(y))
             groups <- as.factor(groups)

             if( any(y<0) ) stop( "Count data has negative values" )

             z <- y/m
             z[is.na(z)]<-0

             yy = .DBSplitIntoGroups( y , groups )[ as.character(unique(groups)) ]
             mm = .DBSplitIntoGroups( m , groups )[ as.character(unique(groups)) ]

             M <- rowSums(m)
             U <- (1/M) * rowSums(m*.rho(z))
             Mk = sapply( mm , rowSums )
             Tk = mapply( FUN = function(y,m) rowSums(y)/rowSums(m) , yy , mm )
             R <- (1/M) * rowSums(Mk*.rho(Tk))

             S <- M*(R - U)
             S = pmax(S,sqrt(.Machine$double.eps))
             return( S )
           }
.DBSplitIntoGroups <-
function (y, groups)
  {
    if( is.null(groups) ) groups <- rep(1,ncol(y))
    groups <- as.factor( groups )
    nrows <- nrow(y)
    yy <- lapply(split(t(y), groups), FUN = function(u) matrix(u, nrow = nrows, byrow = TRUE))[ unique(as.character(groups)) ]
    for (i in 1:length(unique(groups))) rownames(yy[[i]]) <- rownames(y)
    return( yy )
  }
.eta <-
function( vec ) log(vec)-log(1-vec)
.gr.nlgbp <-
function( pars , a , S )
           {
             b = pars[1]
             q = pars[2]
             N = length(S)
             S2 = log( S + q )
             S3 = 1/(S+q)
             dg = digamma( b ) - digamma( a + b )

             derivs = c( -sum(S2) + N*log(q) - sum(dg) , sum(-(a+b)*S3) + N*b/q )

             gr = c( derivs[1] , derivs[2] )
             return( -gr )
           }
.Hessian <-
function( pars , a , S )
        {
          b = pars[1]
          q = pars[2]
          N = length(S)
          S2 = log( S + q )
          S3 = 1/(S+q)
          S4 = 1/(S+q)^2
          tg = trigamma(b) - trigamma(a+b)

          mtx = matrix( c( -sum(tg) , -sum(S3) + N/q , -sum(S3) + N/q , sum((a+b)*S4) - N*b/q^2 ) , 2 , 2 )
          return( -mtx )
        }
.nll.gbp <-
function( pars , a , S )
           {
                 b = pars[1]
                 q = pars[2]
                 N = length(S)
                 S1 = log( S )
                 S2 = log( S + q )
                 lB = lbeta( a , b )
                 llval = sum((a-1)*S1) - sum((a+b)*S2) + N*b*log(q) - sum(lB)
                 return( -llval )
           }
.nll.gbp2 <-
function( pars , a , S )
           {
             nllval <- .nll.gbp( pars , a , S )
             attr(nllval, "gradient") <- .gr.nlgbp( pars = pars , a = a, S = S )
             attr(nllval, "hessian") <- .Hessian( pars = pars , a = a, S = S )
             return( nllval )
           }
.nll.gbp.delta2 <-
function( gamma , a , S )
           {
             delta = gamma / (1 - gamma)
             b = delta*a
             m.S = tapply( S , names(S) , mean )
             m.S = m.S[match( names(S) , names(m.S) )]
             q = delta*m.S
             N = length(S)
             S1 = log( S )
             S2 = log( S + q )
             lB = lbeta( a , b )
             llval = sum((a-1)*S1 - (a+b)*S2 + b*log(q) - lB)
             return( -llval )
           }
.psi <-
function( vec ) -log(1-vec)
.rho <-
function( vec )
        {
          out = .psi(vec) - vec * .eta(vec)
          out[is.na(out)] = 0
          return(out)
        }


.nll.gbp.gamma.bq <- function( pars , a , S )
           {
                 gamma.b <- pars[1]
                 gamma.q <- pars[2]
                 b = gamma.b/(1-gamma.b)
                 q = gamma.q/(1-gamma.q)
                 N = length(S)
                 S1 = log( S )
                 S2 = log( S + q )
                 lB = lbeta( a , b )
                 llval = sum((a-1)*S1) - sum((a+b)*S2) + N*b*log(q) - sum(lB)
                 return( -llval )
           }

.gr.nlgbp.gamma.bq <- function( pars , a , S )
           {
             gamma.b <- pars[1]
             gamma.q <- pars[2]
             b = gamma.b/(1-gamma.b)
             q = gamma.q/(1-gamma.q)
             N = length(S)
             S2 = log( S + q )
             S3 = 1/(S+q)
             dg = digamma( b ) - digamma( a + b )

             db1.gamma.b1 = 1/(1-gamma.b)^2
             dq1.gamma.q1 = 1/(1-gamma.q)^2

             derivs = c( -sum(S2) + N*log(q) - sum(dg) , sum(-(a+b)*S3) + N*b/q ) * c( db1.gamma.b1 , dq1.gamma.q1 )

             gr = c( derivs[1] , derivs[2] )
             return( -gr )
           }

.Hessian.gamma.bq <- function( pars , a , S )
        {
          gamma.b <- pars[1]
          gamma.q <- pars[2]
          b = gamma.b/(1-gamma.b)
          q = gamma.q/(1-gamma.q)
          N = length(S)
          S2 = log( S + q )
          S3 = 1/(S+q)
          S4 = 1/(S+q)^2
          dg = digamma(b) - digamma(a+b)
          tg = trigamma(b) - trigamma(a+b)

          db1.gamma.b1 = 1/(1-gamma.b)^2
          dq1.gamma.q1 = 1/(1-gamma.q)^2
          db2.gamma.b2 = 2/(1-gamma.b)^3
          dq2.gamma.q2 = 2/(1-gamma.q)^3

          d2.b2 = -sum(tg)
          d2.q2 = sum((a+b)*S4) - N*b/q^2
          d2.b1q1 = -sum(S3) + N/q

          derivs1.bq = c( -sum(S2) + N*log(q) - sum(dg) , sum(-(a+b)*S3) + N*b/q )
          d1.b1 = derivs1.bq[1]
          d1.q1 = derivs1.bq[2]

          d2.gamma.b2 = d2.b2 * (db1.gamma.b1)^2 + d1.b1 * db2.gamma.b2
          d2.gamma.q2 = d2.q2 * (dq1.gamma.q1)^2 + d1.q1 * dq2.gamma.q2
          d2.gamma.b1.gamma.q1 = d2.b1q1 * db1.gamma.b1 * dq1.gamma.q1

          mtx = matrix( c( d2.gamma.b2 , d2.gamma.b1.gamma.q1 , d2.gamma.b1.gamma.q1 , d2.gamma.q2 ) , 2 , 2 )
          return( -mtx )
        }
