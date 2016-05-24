#' @title Binomial option pricing
#'
#' @description \code{binomopt} using the binomial pricing algorithm
#'     to compute prices of European and American calls and puts.
#'
#' @name binom
#'
#' @aliases binomopt binomplot binomial
#'
#' @return By default, \code{binomopt} returns the option price. If
#'     \code{returnparams=TRUE}, it returns a list where \code{$price}
#'     is the binomial option price and \code{$params} is a vector
#'     containing the inputs and binomial parameters used to compute
#'     the option price. Optionally, by specifying
#'     \code{returntrees=TRUE}, the list can include the complete
#'     asset price and option price trees, along with trees
#'     representing the replicating portfolio over time. The  current
#'     delta, gamma, and theta are also returned. If
#'     \code{returntrees=FALSE} and \code{returngreeks=TRUE}, only the
#'     current price, delta, gamma, and theta are returned. The function
#'     \code{binomplot} produces a visual representation of the
#'     binomial tree.
#'
#' @usage
#'
#' binomopt(s, k, v, r, tt, d, nstep = 10, american = TRUE,
#'     putopt=FALSE, specifyupdn=FALSE, crr=FALSE, jarrowrudd=FALSE,
#'     up=1.5, dn=1.5, returntrees=FALSE, returnparams=FALSE,
#'     returngreeks=FALSE)
#' 
#' binomplot(s, k, v, r, tt, d, nstep, putopt=FALSE, american=TRUE,
#'     plotvalues=FALSE, plotarrows=FALSE, drawstrike=TRUE,
#'     pointsize=4, ylimval=c(0,0),
#'     saveplot = FALSE, saveplotfn='binomialplot.pdf',
#'     crr=FALSE, jarrowrudd=FALSE, titles=TRUE, specifyupdn=FALSE,
#'     up=1.5, dn=1.5)
#'
#'
#' @param s Stock price
#' @param k Strike price of the option
#' @param v Volatility of the stock, defined as the annualized
#'     standard deviation of the continuously-compounded return
#' @param r Annual continuously-compounded risk-free interest rate
#' @param tt Time to maturity in years
#' @param d Dividend yield, annualized, continuously-compounded
#' @param nstep Number of binomial steps. Default is \code{nstep = 10}
#' @param american Boolean indicating if option is American
#' @param putopt Boolean \code{TRUE} is the option is a put
#' @param specifyupdn Boolean, if \code{TRUE}, manual entry of the
#'     binomial parameters up and down. This overrides the \code{crr}
#'     and \code{jarrowrudd} flags
#' @param up,dn If \code{specifyupdn=TRUE}, up and down moves on the
#'     binomial tree
#' @param crr \code{TRUE} to use the Cox-Ross-Rubinstein tree
#' @param jarrowrudd \code{TRUE} to use the Jarrow-Rudd tree
#' @param returntrees If \code{returntrees=TRUE}, the list returned by
#'     the function includes four trees: for the price of the
#'     underlying asset (stree), the option price (oppricetree), where
#'     the option is exercised (exertree), and the probability of
#'     being at each node. This parameter has no effect if
#'     \code{returnparams=FALSE}, which is the default.
#' @param returnparams Return the vector of inputs and computed
#'     pricing parameters as well as the price
#' @param returngreeks Return time 0 delta, gamma, and theta in the
#'     vector \code{greeks}
#' @param plotvalues display asset prices at nodes
#' @param plotarrows draw arrows connecting pricing nodes
#' @param drawstrike draw horizontal line at the strike price
#' @param pointsize CEX parameter for nodes
#' @param ylimval \code{c(low, high)} for ylimit of the plot
#' @param saveplot boolean; save the plot to a pdf file named
#'     \code{saveplotfn}
#' @param saveplotfn file name for saved plot
#' @param titles automatically supply appropriate main title and x-
#'     and y-axis labels
#' 
#' @importFrom graphics lines plot par points abline arrows mtext text 
#' @importFrom grDevices dev.off pdf
#' 
#' @details Returns an option price, a vector of the parameters used
#'     to compute the price.  Optionally returns
#'     the following matrices, all but but two of which have
#'     dimensionality \eqn{(\textrm{nstep}+1)\times (\textrm{nstep}+
#'     1)}{(nstep+1)*(nstep+1)}:
#' 
#' \describe{
#' 
#' \item{stree}{the binomial tree for the price of the underlying
#'     asset.}
#' 
#' \item{oppricetree}{the binomial tree for the option price at each
#'     node}
#' 
#' \item{exertree}{the tree of boolean indicators for whether or not
#'     the option is exercisd at each node}
#' 
#' \item{probtree}{the probability of reaching each node}
#'
#' \item{delta}{at each node prior to expiration, the number of units
#'     of the underlying asset in the replicating portfolio. The
#'     dimensionality is \eqn{(\textrm{nstep})\times
#'     (\textrm{nstep})}{nstep*nstep}}
#'
#' \item{bond}{at each node prior to expiration, the bond position in
#'     the replicating portfolio. The dimensionality is
#'     \eqn{(\textrm{nstep})\times (\textrm{nstep})}{nstep*nstep}}
#' 
#' }
#' 
#' \code{binomplot} plots the stock price lattice and shows
#' graphically the probability of being at each node (represented as
#' the area of the circle at that price) and whether or not the option
#' is optimally exercised there (green if yes, red if no), and
#' optionally, ht, depending on the inputs
#'
#' @note By default, \code{binomopt} computes the binomial tree using
#'     up and down moves of \deqn{u=\exp((r-d)*h + \sigma\sqrt{h})}{u
#'     = exp((r-d)*h + v*h^(0.5))} and \deqn{d=\exp((r-d)*h -
#'     \sigma\sqrt{h})}{d = exp((r-d)*h - v*h^(0.5))} You can use
#'     different trees: There is a boolean variable \code{CRR} to use
#'     the Cox-Ross-Rubinstein pricing tree, and you can also supply
#'     your own up and down moves with \code{specifyupdn=TRUE}. It's
#'     important to realize that if you do specify the up and down
#'     moves, you are overriding the volatility parameter.
#'
#' @examples
#' s=40; k=40; v=0.30; r=0.08; tt=0.25; d=0; nstep=15
#' 
#' binomopt(s, k, v, r, tt, d, nstep, american=TRUE, putopt=TRUE)
#' 
#' binomopt(s, k, v, r, tt, d, nstep, american=TRUE, putopt=TRUE,
#'     returnparams=TRUE)
#'
#' ## matches Fig 10.8 in 3rd edition of Derivatives Markets
#' x <- binomopt(110, 100, .3, .05, 1, 0.035, 3, american=TRUE,
#'     returntrees=TRUE, returnparams=TRUE)
#' print(x$oppricretree)
#' print(x$delta)
#' print(x$bond)
#' 
#' binomplot(s, k, v, r, tt, d, nstep, american=TRUE, putopt=TRUE)
#' 
#' binomplot(s, k, v, r, tt, d, nstep, american=FALSE, putopt=TRUE)
#' 
#' 


## this matches fig 10.8 in the 3rd edition:
## binomopt(110, 100, .3, .05, 1, 0.035, 3, american=TRUE,
##          returntrees=TRUE, returnparams=TRUE)

#' @export
binomopt <- function(s, k, v, r, tt, d,
                     nstep=10, american = TRUE, putopt=FALSE,
                     specifyupdn=FALSE, crr=FALSE, jarrowrudd=FALSE,
                     up=1.5, dn=1.5, returntrees=FALSE,
                     returnparams=FALSE, returngreeks=FALSE) {
    ## set up the binomial tree parameters
    
    h <- tt/nstep
    if (!specifyupdn) {
        if (crr) {
            up <- exp(sqrt(h)*v)
            dn <- exp(-sqrt(h)*v)
        } else if (jarrowrudd) {
            up <- exp((r-d-0.5*v^2)*h + sqrt(h)*v)
            dn <- exp((r-d-0.5*v^2)*h - sqrt(h)*v)
        } else {
            up <- exp((r-d)*h + sqrt(h)*v)
            dn <- exp((r-d)*h - sqrt(h)*v)
        }
    }
    p <- (exp((r-d)*h) - dn)/(up - dn)
    nn <- 0:nstep
    payoffmult <- ifelse(putopt,-1,1)
    ##
    ## Construct stock price matrix
    ##
    stree <- matrix(0, nstep+1, nstep+1)
    stree[1, ] <- s*up^nn
    for (i in 2:(nstep+1))
        stree[i, i:(nstep+1)] <- stree[i-1, (i-1):nstep]*dn
    Vc <- matrix(0, nstep+1, nstep+1)  ## Initialize opt price matrix
    Vc[, nstep+1] <- pmax(0, (stree[, nstep+1] - k)*payoffmult)
    ##
    ## recurse
    ##
    for (i in (nstep):1) {
        Vnc <- exp(-r*h)*(p*Vc[1:i,i+1] + (1-p)*Vc[2:(i+1),i+1])
        if (american) {
            Vc[1:i,i] <- pmax(Vnc, (stree[1:i,i] - k)*payoffmult)
        } else {
            Vc[1:i,i] <- Vnc
        }
    }
    price <- c(price=Vc[1, 1])
    if (returntrees | returngreeks) {
        deltatree <- matrix(0, nrow=nstep, ncol=nstep)
        bondtree <- matrix(0, nrow=nstep, ncol=nstep)
        for (i in 1:(nstep)) {
            deltatree[1:i, i] <- exp(-d*h)*(Vc[1:i, i+1] - Vc[2:(i+1), i+1])/
                (up-dn)/stree[1:i, i]
            bondtree[1:i, i] <- exp(-r*h)*(up* Vc[2:(i+1), i+1] -
                                      dn*Vc[1:i, i+1])/(up-dn)
        }
       
        delta <- deltatree[1, 1]
        if (nstep >= 2) {
            gamma <- (deltatree[1, 2] - deltatree[2, 2])/
                (stree[1, 2] - stree[2, 2])
            epsilon <- (up*dn - 1)*s
            theta <- (Vc[2, 3] - epsilon*delta - 0.5*epsilon^2*gamma
                - Vc[1, 1])/(2*h)/365
        } else {
            theta <- NA
            gamma <- NA
        }
        greeks <- c(delta=delta, gamma=gamma, theta=theta)
    }
    if (!returnparams & !returntrees & !returngreeks) return(price=price)
    params=c(s=s, k=k, v=v, r=r, tt=tt, d=d,
             nstep=nstep, p=p, up=up, dn=dn, h=h)
    if (returntrees) {
        probtree <- matrix(0, nstep+1, nstep+1)
        A <- matrix(p^nn, nstep+1, nstep+1, byrow=TRUE)
        B <- matrix((1-p)^nn/p^nn, nstep+1, nstep+1)
        for (i in 0:nstep) { probtree[,i+1] <- choose(i,nn) }
        exertree <- (payoffmult*(stree - k) == Vc)
        probtree <- A*B*probtree
        return(list(price=price, greeks=greeks, params=params,
                    oppricetree=Vc, stree=stree, probtree=probtree,
                    exertree=exertree, deltatree=deltatree,
                    bondtree=bondtree)
               )
    } else if (returngreeks) {
        return(c(price, greeks, params))
    } else {
        return(c(price, params))
    }
}


#' @export
binomplot <- function(s, k, v, r, tt, d, nstep, putopt=FALSE,
                      american=TRUE, plotvalues=FALSE,
                      plotarrows=FALSE, drawstrike=TRUE, 
                      pointsize=4, ylimval=c(0,0),
                      saveplot = FALSE, saveplotfn='binomialplot.pdf',
                      crr=FALSE, jarrowrudd=FALSE, titles = TRUE,
                      specifyupdn=FALSE, up=1.5, dn=1.5) {
    ## see binomopt for more details on tree
    ## construction. "plotvalues" shows stock price values;
    ## "drawstrike" if true draws a line at the strike price; "probs"
    ## makes pointsizes proportional to probability of that point,
    ## times "pointsize"
    ##
    ## If no value given for ylimval and setylim=TRUE, there will be
    ## an error
    ## 
    
    setylim <- ifelse((sum(ylimval^2)==0), FALSE, TRUE)
    y <- binomopt(s, k, v, r, tt, d, nstep, american, putopt,
                  specifyupdn, crr, jarrowrudd, up, dn,
                  returnparams=TRUE, returntrees=TRUE)
    h <- tt/nstep
    for (i in c('up', 'dn', 'p')) assign(i, y$params[i])
    for (i in c('stree', 'exertree', 'oppricetree', 'probtree'))
        assign(i, y[['i']])
    nn <- 0:nstep
    payoffmult <- ifelse((putopt),-1,1)
    stree <- y$stree
    exertree <- y$exertree
    oppricetree <- y$oppricetree
    probtree <- y$probtree
    
    ## The rep command replicates each entry in nn a different number of
    ## times (1st entry once, second, entry twice, etc. Need to add 1
    ## because the first entry in nn is zero, which implies zero reps. The
    ## point of the stree restriction is not to plot zeros.
    plotcolor <- ifelse(exertree,"green3","red")
    if (saveplot) pdf(saveplotfn)
    plot(rep(nn, nn+1)*h, stree[stree>0]
        ,ylim=ifelse(c(setylim, setylim),ylimval,
                     c(min(stree[stree>1]-2),max(stree)*1.03))
        ,col=plotcolor[stree>0]
        ,pch=21
         ## ifelse returns an object with the size of the first
         ## argument. So in order for it to pass an array the first
         ## argument is an array filled with the boolean "probs". 
        ,cex=ifelse(stree[stree>0], sqrt(probtree[stree>0])*pointsize, 1)
        ,bg=plotcolor[stree>0] ## only matters for pch 21-25
        ,xlab=ifelse(titles, "Binomial Period", "")
        ,ylab=ifelse(titles,  "Stock Price", "")
        ,main=if (titles) paste(ifelse(american,"American","European"),
                                ifelse(putopt,"Put","Call"))
         )
    if (titles)
        mtext(paste0("Stock = ",format(s, digits=3),
                     ", Strike = ",format(k, digits=3),
                     ", Time = ",format(tt,digits=4),
                     ifelse(tt==1," year,"," years,")
                    ," Price = ",format(oppricetree[1,1],digits=5)))
    if (drawstrike) abline(h=k)
#    yoffset <- ifelse(setylim, 0.075*ylimval[1], 0.03*max(stree))
    yoffset <- ifelse(setylim, 0.03*(ylimval[2]-ylimval[1]),
                      0.03*max(stree))
    if (plotarrows) {
        for (i in 1:nstep) {
            for (j in 1:i) {
                arrows((i-1)*h, stree[j,i],c(i,i)*h,
                       c(stree[j,i+1],stree[j+1,i+1]), length=0.06)
            }
        }
    }
    if (plotvalues) {
        for (i in 1:(nstep+1)) {
            text((i-1)*h,stree[1:i,i]+yoffset,format(stree[1:i,i], digits=3),
                 cex=0.7)
        }
    }
    if (saveplot) dev.off()
}
