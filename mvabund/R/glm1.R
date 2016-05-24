glm1 = function(y, X, lambda, family="negative.binomial", weights = rep(1, length(y)), b.init = NA, phi.init=NA,
phi.method="ML", tol=c(1.e-8,.Machine$double.eps), n.iter=100, phi.iter=1)
{
    counter = 1
    zero.counter=0

    family.old=family

    # turn family into a family argument
    if (is.character(family)) 
    {
      if (family == "negbinomial" || family == "negative.binomial")
          family = negative.binomial(1/tol[2])
      else    
        family = get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) 
        family = family()

    # if it is negative binomial, set is.nb=T
    if(pmatch("Negative Binomial",family$family,nomatch=0)==0)
      is.nb=F
    else
    {
      is.nb=T
      if(is.na(phi.init))
      {
        phi.init = family$var(1)-1 #extract overdispersion parameter by evaluating variance at mu=1 and solving for phi.
        if(phi.init>2*tol[2])
          phi.iter=n.iter*2 #if phi was specified in family argument, make sure it stays fixed in estimation.
      }
    }    
    
    if (length(lambda) == 1)
        lambda = c( 0, rep(lambda, dim(X)[2]-1) )

    likes = c()
    phis = c()
    mult = 0

    #GETTING INITIAL ESTIMATE B.lasso as the mean of the response vector as follows: (mean,0,0,0,.....,0)
    if (any(is.na(b.init)) ) #Testing whether the lasso/glm vector has been initialised to NAs
        b.lasso = c( family$linkfun(mean(y)), rep(0, dim(X)[2] - 1))
    if (any(is.na(b.init)) == FALSE)
        b.lasso = b.init

    #initialise is.in (the active set,  parameters in the model (and this occurs only if the derivatives are greater than lambda)) and is.different vectors
    is.in        = abs(b.lasso) > tol[2] | lambda==0
    sign.change  = rep(-1, dim(X)[2])
    is.different = is.in  # Should initially be true if the is.in variable is true

    signs = sign(b.lasso)

    res = mu.update(b.lasso,y,X,weights,family,is.nb,phi.init,lambda,init=T,tol=tol)


#BEGINNING THE LASSO ESTIMATION ALGORITHM.
    while(any(is.different == TRUE))
    {

      likes = c(likes, res$like)
      phis = c(phis, res$phi)
      
### UPDATE BETA ###

        #coeff from previous iteration
        b.old = b.lasso
        phi.old = res$phi

        xw   = t(as.vector(res$weii) * t(t(X[,is.in])))
        xwx   = solve(xw %*% X[,is.in] + diag( sum(is.in) )*min(tol[1]/100,sqrt(tol[2])) )
        betaa = xwx %*% ( xw %*% res$z - as.matrix(lambda[is.in] * signs[is.in]))
#        betaa = qr.solve(xw %*% X[,is.in] + diag( sum(is.in) )*sqrt(tol), xw %*% z - as.matrix(lambda[is.in] * signs[is.in]))

# end update beta


#DW modification: move deletion step


### UPDATE MU, WEIGHTS, etc ###
        b.lasso[is.in] = betaa #changing the lasso estimate to the most recent estimate
        res      = mu.update(b.lasso,y,X,weights,res$family,is.nb,phi.old,lambda,init=F,phi.step=F,phi.method,tol=tol)
        #at this stage, not updating phi, hence phi.step=F
# end updating mu etc

### HALF-STEP IF REQ'D TO ENSURE INCREASING LL ###

        dev.change = (res$like - likes[length(likes)])/abs(likes[length(likes)]+0.1)
        half.step = 0
        while(dev.change< (-tol[1]) & half.step<4)
        {
            #do a half-step
            b.lasso[is.in] = 0.5 * ( b.old[is.in] + b.lasso[is.in] )

            half.step = half.step + 1

            res = mu.update(b.lasso,y,X,weights,res$family,is.nb,phi.old,lambda,init=F,phi.step=F,phi.method,tol=tol)
            #still not updating phi
            dev.change = (res$like - likes[length(likes)])/abs(likes[length(likes)]+0.1)
        }
# end half-step


### DELETION STEP ###

        #Check signs
        sign.change                   = rep(1,dim(X)[2])
        sign.change[is.in]            = signs[is.in]*sign(b.lasso[is.in])

        sign.change[lambda==0]        = 1
        sign.change[sign.change == 0] = 1

        if (any(sign.change != 1) == TRUE)
        {
            delta                       = b.lasso[is.in]-b.old[is.in]
            prop                        = min(abs(b.old[is.in]/delta)[sign.change[is.in] != 1])
            #print(paste("exclusion step, variable", "...","excluded"))
            b.lasso[is.in]              = b.old[is.in] + prop*delta
            mult                        = 2
            res = mu.update(b.lasso,y,X,weights,res$family,is.nb,phi.old,lambda,init=F,phi.step=F,phi.method,tol=tol)
            dev.change = (res$like - likes[length(likes)])/abs(likes[length(likes)]+0.1)
        }
        is.in[abs(b.lasso) < tol[2]] = FALSE
        b.lasso[is.in == FALSE]     = 0
        signs          = sign(b.lasso)
# end deletion step


### INCLUSION STEP ###

#need to update mu and IRLS weights first to get correct score eqns.
#but note that nuisance params should not be updated yet - otherwise score eqns won't be satisfied
        xw  = t(as.vector(res$weii) * t(t(X)))
        score    =  t( t(xw) / as.vector(res$deriv) ) %*% (y - res$mu)
        score.lamb = lambda
        score.lamb[lambda == 0] = abs(score[lambda == 0]) + max(lambda) + 100000
        viol = as.vector(abs(score))/as.vector(score.lamb)
        bigviol = max(viol)
        viol.out = (viol==bigviol & is.in == FALSE)

        if (any(sign.change != 1) == FALSE)
        {
            #Check to see if any variables need to be added
            if ( any( viol.out ) & bigviol > (1+tol[1]))
            {
                {
                    to.add               = which(viol.out)[1] #only add one at a time
                    is.in[to.add]        = TRUE
                    signs[to.add]        = sign(score[to.add])
                    mult = 1
                }
            }
        }
# end inclusion step

#UPDATE PHI (if required)
        if(is.nb)
        {
            phi.step = counter%%phi.iter==0
            res = mu.update(b.lasso,y,X,weights,res$family,is.nb,phi.old,lambda,init=F,phi.step=phi.step,phi.method,tol=tol)

            phis[length(phis)] = res$phi

            dev.change = (res$like - likes[length(likes)])/abs(likes[length(likes)]+0.1)
        }
        likes[length(likes)] = res$like

        diff         = mult+abs(dev.change)
        is.different = rep(diff > tol[1], dim(X)[2]) #To stay in the loop if not converged
        mult         = 0
        if ( counter>n.iter)
        {
          is.different = FALSE
            warning('non-convergence!!!')
            counter = Inf
        }

        counter = counter + 1
        df = sum(abs(b.lasso)>tol[2])

#        summ=c(counter,max(lambda),df,bigviol,res$phi,diff,res$like)
#        names(summ)=c("counter","lam","df","bigviol","phi","diff","like")
#        print(summ)
        if(df==0) #add some random terms to model if none there
        {
            zero.counter = zero.counter + 1
            if(zero.counter>3)
                break
            is.in[sample(1:dim(X)[2],5)]=T
        }
    }
    likes = likes + sum(as.vector(lambda)*abs(b.lasso))
    
    check = any( abs(score[is.in==F])>lambda[is.in==F]+tol[1] ) || any( abs( abs(score[is.in])-lambda[is.in] ) > tol[1] )    
    #end IRLS/LASSO algorithm

    names(b.lasso) = dimnames(X)[[2]]   
    dimnames(score)[[1]]=dimnames(X)[[2]]   
    dimnames(res$mu)=dimnames(y)

    return(list(coefficients=b.lasso, fitted.values=res$mu, logLs = likes, phis=phis, phi=res$phi, score=score, counter = counter, check=check, family=res$family, weights=weights))
}


mu.update = function(b.lasso, y, X, weights, family, is.nb, phi, lambda, init=F, phi.step=T, phi.method="ML",tol=c(1.e-8,.Machine$double.eps))
{
    eta   = X %*% b.lasso
    mu = family$linkinv(eta)

    if(is.nb==F)
        phi = NA #this is not needed anywhere, included to avoid error when storing phis at end

    if (is.nb == T)
    {
        if(init == T) #initialising phi
        {
            disp = sum( (y - mu)^2/(mu )    ) / (length(y) - dim(X)[2])
            if (is.na(phi)) #if initial phi is NA, choose a starting value
            {
            #Initialising phi like in chi-sq dampening (using the dispersion, not the likelihood):
                if (disp > 1)
                    phi = 1
                else
                    phi <- tol[1]
            }
        }
        if(init == F)
        {
            if( phi.step==T )
            { #do phi-step only if required
                if (phi.method=="ML") #do mL estimation of phi
                {
                    #now update phi (via k):
                    k <- 1/phi

                    # the first derivative of the likelihood wrt phi;
                    dl.dk     = sum( - (y-mu)/(mu+k) + digamma(y+k) - log(mu+k) ) + length(y) * ( log(k) - digamma(k) )
                    dl.dphi   = - k^2 * dl.dk
                    like.num  = dl.dphi

                    #the second derivative; will involve the trigamma
                    d2l.dk2    = sum( (y-mu)/(mu+k)^2 + trigamma(y+k) - 1/(mu+k) ) + length(y) * ( 1/k - trigamma(k) )
                    d2l.dphi2  = 2*k^3*dl.dk + k^4 * d2l.dk2
                    like.den   = d2l.dphi2

                    phi = phi + like.num/abs(like.den) #plus ./abs(like.den) to get out of concave down areas
                }
                else
                {
                    like.num = sum( (y - mu)^2/(mu + phi*mu^2)  ) / (length(y) - dim(X)[2])  -1
                    like.den = - sum(mu^2*(y - mu)^2/(mu + phi*mu^2)^2) / (length(y) - dim(X)[2])
                    phi <- phi - like.num/like.den
                }
            }
            # keep phi finite
            phi <- min(phi,1/tol[1])     
        }
        # if phi is not positive, use poisson for logL calculations
        if(phi<=tol[1])
        {
          phi = tol[1]
          family = poisson()
        }
        else
        {
          k = 1/phi
          family = negative.binomial(k) #update overdispersion in family object
        }
    }

    if( family$family=="binomial" || family$family=="poisson" || pmatch("Negative Binomial",family$family,nomatch=0)==1 )
      dev = NA #to save a little computation time, don't get dev when it is not needed by the aic function.
    else
      dev = family$dev(y,mu,weights) #only actually needed by gaussian, gamma, inverse.gaussian
    like  = -sum(family$aic(y,1,mu,weights,dev))/2 - sum(as.vector(lambda)*abs(b.lasso))
    vari = family$variance(mu)
    deriv = family$mu.eta(eta)
    z     = eta + (y - mu)/deriv
    weii  = weights*deriv^2/vari
    return(list(eta=eta,mu=mu,like=like,vari=vari,deriv=deriv,z=z,weii=weii,phi=phi,family=family))
}
