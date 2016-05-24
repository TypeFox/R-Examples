# simple loglikelihood
loglik.sim.dw <-  function(lik.par , x  ,    lq.extra ,  lbeta.extra , logit)
{
  l.l = length(lik.par)
  if( l.l != 2) stop('error in the number of parameters! ',l.l,' <> T = ',2)
  q    = lik.par[1] * lq.extra
  beta = lik.par[2] * lbeta.extra

  if (all(beta >= 0) && all(q<= 1) && all(q >= 0) ){
    loglik = sum(log(q^(x^beta)-q^((x+1)^beta)))
  }else{
    #if(any(q==0) || any(q==1) || any(beta==0)){
    #loglik = -Inf
    #}else{
    loglik = NA
    #cat('\r NA')
    #}
  }
  return(loglik)
}
# loglikelihood - regression on q
loglik.q.dw  <-  function(lik.par , x , y , lq.extra ,  lbeta.extra  , logit )
{
  l.like.par = length(lik.par)
  if( l.like.par != ncol(x)+1) stop('error in the number of parameters!',l.like.par,'<> T=',ncol(x)+1)
  beta = lik.par[ l.like.par] * lbeta.extra
  theta= lik.par[-l.like.par] * lq.extra
  if(logit){
    q    = exp (x %*% theta-log(1+exp(x %*% theta)) )
  }else{
    q    = exp (- exp ( x %*% theta ))
  }

  if (all(beta >= 0) && all(q <= 1) && all(q >= 0) ){
    loglik = sum(log(q^(y^beta)-q^((y+1)^beta)))
  }else{
    #if(any(q==0) || any(q==1) || any(beta==0)){
    # loglik = -Inf
    #}else{
    loglik = NA
    #cat('\r NA')
    #}
  }
  return(loglik)
}
#  loglikelihood - regression on beta
loglik.beta.dw <- function(lik.par , x , y , lq.extra ,  lbeta.extra , logit)
{
  l.l = length(lik.par)
  if(l.l  != ncol(x)+1) stop('error in the number of parameters!You:',l.l,' <> T = ',ncol(x)+1)
  q      = lik.par[ 1] * lq.extra    #par[ length(par)]
  theta  = lik.par[-1] * lbeta.extra #par[-length(par)]
  beta   = exp(x %*% theta)
  #if(q %in% theta) print('Q in theta!')

  if (all(beta >= 0) && all(q <= 1) && all(q >= 0) ){
    loglik = sum(log(q^(y^beta)-q^((y+1)^beta)))
  }else{
    #     if(any(q == 0) || any(q == 1) || any(beta ==0)){
    #       loglik = -Inf
    #     }else{
    loglik = NA
    #cat('\r NA')
    # }
  }
  return(loglik)
}
#  loglikelihood - regression on both beta and q
loglik.qbeta.dw<- function(lik.par , x , y , lq.extra ,  lbeta.extra , logit)
{
  ncx        = ncol(x)
  l.l        = length(lik.par)
  if( l.l != 2*ncx) stop('error in the number of parameters!',l.l,' <> T = ',2*ncol(x) )
  theta.q    = head(lik.par , ncx ) * lq.extra
  theta.beta = tail(lik.par , ncx ) * lbeta.extra

  x1 = x #+ xxx1
  x2 = x #+ xxx2

  if(all(x[,1]==1)){
    x1[,1]=1
    x2[,1]=1
  }

  if(logit){
    q          = exp(x2 %*% theta.q-log(1+exp(x %*% theta.q)))
  }else{
    q          = exp(-exp(x2 %*% theta.q))
  }
  beta         = exp(x1 %*% theta.beta)

  if (all(beta >= 0) && all(q <= 1) && all(q >= 0) ){
    loglik = sum(log(q^(y^beta)-q^((y+1)^beta)))
  }else{
    #if((any(q==0) || any(q==1)) || any(beta==0)){
    # loglik = -Inf
    #}else{
    loglik = NA
    # cat('\r NA')
    # }
  }
  return(loglik)
}
###########   end of likelihoods ##################
#  prior distribution
prior.tot.dw <-  function(pr.par   , q.par , b.par      ,
                          dist.q   , dist.b             ,
                          lq , lb  , l.par , dist.l     ,
                          penalized , fixed.l)
{
  #------ Lambda prior parameters if needed
  par.length   = length(pr.par)

  if((penalized == TRUE) && (fixed.l <= 0) ){
    if(par.length != (lb+lq+1)) stop('error in the number of parameters!',par.length,' <> T = ',(lb+lq+1))
    #lp           = pr.par[(lb+lq+1):par.length]
    lp           = pr.par[par.length]
    l.par1       = l.par[1]
    l.par2       = l.par[2]
    lambda.prior =  sum (dist.l(lp ,  l.par1  ,   l.par2  ,   log = TRUE ))
  }else if ( (penalized == TRUE) && (fixed.l > 0) ){
    if(par.length != (lb+lq)) stop('error in the number of parameters!',par.length,' <> T = ' ,(lb+lq))
    lambda.prior = 0
    lp           = fixed.l
  }else{
    if(par.length != (lb+lq)) stop('error in the number of parameters!',par.length,' <> T = ',(lb+lq))
    #lp           = rep(1,lp+lq)
    lp           = 1
    lambda.prior = 0
  }
  #------ Q prior parameters
  q.par1 =               q.par[1]
  q.par2 =    sqrt(lp) * q.par[2]
  #q.par2 =  head(lp,lq)* q.par[2]
  #------ Beta prior parameters
  b.par1 =                 b.par[1]
  b.par2 =     sqrt(lp) *  b.par[2]
  #b.par2 =  tail(lp,lb) *  b.par[2]
  #------ modedl parameters
  qp     = head(pr.par , lq            )
  bp     =     (pr.par[ (lq+1):(lq+lb)])
  #------ priors
  theta.prior  = sum (dist.q(qp ,  q.par1  ,   q.par2 ,    log = TRUE))
  gamma.prior  = sum (dist.b(bp ,  b.par1  ,   b.par2  ,   log = TRUE))
  prior =  ( theta.prior + gamma.prior + lambda.prior)
  return(prior)
}
#  posteriot distribution
posterior.tot.dw <- function(par     ,
                             x       ,  y      ,
                             lq.extra          ,
                             lbeta.extra       ,
                             para.q  , para.b  ,
                             prior.function    ,  ####### END OF PRIOR ####
                             dist.q  , q.par   ,
                             dist.b  , b.par   ,
                             dist.l  , l.par   ,
                             penalized  , fixed.l ,
                             logit ,
                             ...
)
{
  if ((para.b == FALSE) && (para.q == TRUE)) {
    loglik   = loglik.q.dw   (lik.par = head(par,ncol(x) + 1),
                              x = x , y = y ,
                              lq.extra = lq.extra ,
                              lbeta.extra = lbeta.extra,
                              logit = logit)
    prior    = prior.function(pr.par = par,
                              q.par  = q.par   , b.par  = b.par   ,
                              dist.q = dist.q  , dist.b = dist.b  ,
                              #lb = 1           , lq = ncol(x)-1   ,
                              lb = 1           , lq = ncol(x)     ,
                              l.par = l.par    , dist.l = dist.l  ,
                              penalized =  penalized , fixed.l = fixed.l)
    post     = loglik + prior
  }else if ((para.b == TRUE ) && (para.q == FALSE)) {
    loglik   = loglik.beta.dw(lik.par = head(par,ncol(x) + 1), x = x, y = y ,
                              lq.extra = lq.extra ,lbeta.extra = lbeta.extra,
                              logit = logit)
    prior    = prior.function(pr.par = par,
                              q.par  = q.par   , b.par  = b.par  ,
                              dist.q = dist.q  , dist.b = dist.b ,
                              dist.l = dist.l  , l.par  = l.par  ,
                              #lb = ncol(x)-1   , lq = 1  ,
                              lb = ncol(x)     , lq = 1          ,
                              penalized = penalized , fixed.l = fixed.l)
    post     = loglik + prior
  }else  if ((para.b == TRUE) && (para.q == TRUE)) {
    loglik   =  loglik.qbeta.dw(lik.par  = head(par,2*ncol(x)), x = x, y = y ,
                                lq.extra = lq.extra ,lbeta.extra = lbeta.extra,
                                logit = logit)
    prior    =  prior.function (pr.par    = par    ,
                                q.par  = q.par  , b.par  = b.par  ,
                                dist.q = dist.q , dist.b = dist.b ,
                                dist.l = dist.l , l.par  = l.par  ,
                                lb = ncol(x)    , lq = ncol(x)    ,
                                penalized = penalized ,  fixed.l = fixed.l)
    post   = loglik + prior
  }else{
    loglik   =  loglik.sim.dw (lik.par  = head(par,2), x = y  ,
                               lq.extra = lq.extra ,lbeta.extra = lbeta.extra,
                               logit = logit)
    prior    =  prior.function(pr.par    = par                ,
                               q.par  = q.par  , b.par  = b.par  ,
                               dist.q = dist.q , dist.b = dist.b ,
                               dist.l = dist.l , l.par = l.par   ,
                               lb = 1          , lq = 1 ,
                               penalized = penalized ,  fixed.l = fixed.l )
    post   = (loglik + prior)
  }
  return(post)

}

#  Proposal distribution
proposalfunction.tot.dw <- function(par   , v.scale  ,
                                    para.q, para.b   ,
                                    s , i , cov.m    ,
                                    lb , lq ){
  length.par = length(par)
  if ((cov.m == 1) &&
      (i > 100) &&
      (i %% 5 == 0)
  ) {
    s   =    dw.cov.matrix(x = s , i = i )
    ncs =    ncol(s)
    #b1  =    mvnfast::rmvn(n = 1,mu = par*0,sigma = (2.38^2)/ncs *s )
    b1  =    MASS::mvrnorm(n = 1   , mu   = par*0 ,  Sigma = (2.38^2)/ncs *s )
    b2  =    rnorm  (n = 1   , mean = 0     ,  sd    = rep(1/ncs , ncs))
    prop.d = par + b1 + v.scale*b2
  }else if (cov.m == 2) {
    prop.d     =    par + runif(length.par ,  -abs(v.scale)/2 , abs(v.scale)/2 )
  }else if (cov.m == 3) {
    prop.d = par + rlaplace(n = length.par , mu =  0 ,sigma = v.scale )
  }else{
    prop.d     =    par + rnorm(length.par,  0               , v.scale )
  }
  if (length.par > lb + lq) {
    lp = (lb + lq + 1):(length.par)
    prop.d[lp] = abs(prop.d[lp])
  }
  return(prop.d)
}




# Reversible jumps & Metropolis-Hasting sampling
MH.tot.dw <- function(formula, data, startvalue      ,
                      para.q  , para.b               ,
                      dist.q  , dist.b               ,
                      q.par   , b.par                ,
                      dist.l  , l.par                ,
                      fixed.l , logit                ,
                      cov.m   , penalized            ,
                      v.scale , iterations           ,
                      lq.extra, lbeta.extra          ,
                      prior.function                 ,
                      RJ = FALSE
)
{
  if ( para.q  || para.b ) {
    call <- match.call()
    mf <- match.call(expand.dots = TRUE)
    m  <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y  <- model.response(mf, "numeric")
    x  <- model.matrix(mt, mf, contrasts)
    #up  =  .0000000005
    #down= 0
    #xxx1 <<- rnorm(length(x),down,up)
    #xxx2 <<- rnorm(length(x),down,up)
    if (para.b == TRUE) {lb = ncol(x)}else{lb = 1}
    if (para.q == TRUE) {lq = ncol(x)}else{lq = 1}
  }else{
    y  = as.vector(data)
    x  = NA
    lb = lq = 1
    RJ = FALSE
  }

  #   message(paste('MH - Penalized = ',penalized,' and qReg = ',para.q, ', betaReg = ',para.b,
  #                 '\n  and ', (lb+lq),
  #                 ' parameters and ', length(startvalue) , 'initial values applied.' ))

  len.par = length(startvalue)
  #### This is for RJ
  if(RJ){
    #message('Reversible-Jump in action ... ')
    penalized = FALSE ; fixed.l = -1
  }
  total.a.par = lb + lq
  prob   = rep(.5 , total.a.par) # initial  probability of existing variables
  rjmu   = rep( 0 , total.a.par)
  rjsig  = rep(.5 , total.a.par)
  #pparam = startvalue[1:total.a.par]
  #pparam[pparam!=0] = 1
  model.chain       = c()
  chain             = array(dim = c(iterations + 1 , len.par));
  chain[1,]         = startvalue
  ####### Start iterations  #######
  #pb    = txtProgressBar(min = 0 , max = iterations,initial = 0 , style = 3)
  #pb    = winProgressBar(title="MCMC progress - close me by pressing Alt+F4 if necessary.", label="0% done", min=0, max=iterations, initial=0)
  pr2   = posterior.tot.dw( par = chain[1 , ]  ,
                            para.q = para.q , para.b = para.b ,
                            q.par = q.par   , b.par = b.par   ,
                            dist.q = dist.q , dist.b = dist.b ,
                            dist.l =  dist.l, l.par  = l.par  ,
                            x = x   , y = y , logit = logit   ,
                            lq.extra    = lq.extra            ,
                            lbeta.extra = lbeta.extra         ,
                            penalized   = penalized           ,
                            prior.function = prior.function   ,
                            fixed.l     = fixed.l
  )
  state         = startvalue
  functionat    = pr2
  # ---- CORE
  i = 1;  j = 1;  k = 0 ; ReAc = 0; error = 0 ;
  it.interval = round(iterations/20) ; interval = 0
  while (i <= iterations) {
    if (i %% it.interval == 0) {
      acceptance = apply(chain,2, function(x){1-mean(duplicated(x))} )
      acceptance = round(max (acceptance) * 100, 2)
      info <- paste(round(i/iterations*100),'% done, Acceptance = ',round(acceptance,2),'%')
      #setWinProgressBar(pb, i,label = info )
      #setTxtProgressBar(pb, i,label = acceptance)
      cat('\r ',  info)
    }
    if (i < 15) { #just to make sure that the process is not stuck
      vs = 10 + v.scale
    }  else  {
      vs =      v.scale
    }


    proposal = proposalfunction.tot.dw(par = chain[i , ],
                                       para.q = para.q ,para.b = para.b,
                                       v.scale = vs   ,
                                       s = chain , i = i   ,
                                       cov.m  = cov.m      ,
                                       lb = lb , lq = lq)

    # To make sure that proposal is not acting on zero parmeters.
    if(RJ){
      zeros = which (chain[i , ] == 0)
      if(length(zeros) > 0) {
        proposal[zeros] = 0
      }
    }
    pr1 = posterior.tot.dw( par = proposal ,
                            para.q = para.q , para.b = para.b ,
                            q.par  = q.par  , b.par = b.par   ,
                            dist.q = dist.q , dist.b = dist.b ,
                            dist.l = dist.l , l.par  = l.par  ,
                            x = x   , y = y,  logit = logit   ,
                            lq.extra = lq.extra               ,
                            lbeta.extra = lbeta.extra         ,
                            penalized = penalized             ,
                            prior.function = prior.function   ,
                            fixed.l = fixed.l
    )

    probab = min(1,exp( pr1 - pr2 )) #* correction
    if (!is.na(probab)) {
      if ((runif(1) <= probab ) ) {
        chain[i + 1 , ]  = proposal
        ReAc [i + 1   ]  = 1
        #pr4              = pr2
        pr2              = pr1
      }else{
        chain[i + 1 , ] = chain[i , ]
        ReAc [i + 1   ] = 0
        #pr4             = pr1
      }

      ###### ----- This is a maximizing procedure ---- ###
      if (tail(functionat,1)  <  pr1 )  {
        state            = rbind(state , proposal)
        functionat       =  c(functionat,pr1)
        #cat('\r OK')
      }
      ###### END of Maximizing

      #----- START Riversible jump!
      sz = 1  # one sample each step!
      if(lb==1 && lq>1){
        interval = 1:lq
        r  = sample(x = interval , size = sz)
      }else if(lb>1 && lq==1){
        interval = 2:total.a.par
        r  = sample(x = interval , size = sz)
      }else if(lb>1 && lq>1){
        interval = 1:total.a.par
        r  = sample(x =  interval, size = sz)
      }else{
        RJ = FALSE
      }
      if(RJ){
        pr4 = pr2
        tr = chain[i + 1 , r]
        newchain    = as.vector(chain[i + 1 , ])
        if ( all(tr == 0) ){
          rv          = rnorm(length(r) , rjmu[r] , rjsig[r])
          newchain[r] = rv
          pr3 = posterior.tot.dw( par = newchain ,
                                  para.q = para.q , para.b = para.b ,
                                  q.par  = q.par  , b.par = b.par   ,
                                  dist.q = dist.q , dist.b = dist.b ,
                                  dist.l = dist.l , l.par  = l.par  ,
                                  x = x , y = y   , logit = logit   ,
                                  lq.extra = lq.extra               ,
                                  lbeta.extra = lbeta.extra         ,
                                  penalized = penalized             ,
                                  prior.function = prior.function   ,
                                  fixed.l = fixed.l
          )
          #        num   =  pr3 + dnorm(rv , chain[i + 1 , r]    , chain[i + 1 , len.par], log = TRUE) + log(prob[r])
          num   =  pr3 + sum(dist.b(rv  , b.par[1], b.par[2], log = TRUE)) + sum(log(  prob[r]))
          den   =  pr4 + sum(dnorm (rv  , rjmu[r] , rjsig[r], log = TRUE)) + sum(log(1-prob[r]))
        }else{
          rv          = 0
          newchain[r] = rv
          pr3 = posterior.tot.dw( par = newchain ,
                                  para.q = para.q , para.b = para.b ,
                                  q.par  = q.par  , b.par = b.par   ,
                                  dist.q = dist.q , dist.b = dist.b ,
                                  dist.l = dist.l , l.par  = l.par  ,
                                  x = x  , y = y  , logit = logit   ,
                                  lq.extra = lq.extra               ,
                                  lbeta.extra = lbeta.extra         ,
                                  penalized = penalized             ,
                                  prior.function = prior.function   ,
                                  fixed.l = fixed.l
          )
          #den   =  pr4 + dnorm(tr , chain[i + 1 , r]  , chain[i + 1 , len.par] , log = TRUE) + log(prob[r])
          num   =  pr3 + sum(dnorm (tr , rjmu[r], rjsig[r], log = TRUE)) + sum(log(1-prob[r]))
          den   =  pr4 + sum(dist.b(tr ,b.par[1], b.par[2], log = TRUE)) + sum(log(  prob[r]))
        }
        if(runif(1) <= min(1,exp(num-den)) ){
          chain[i + 1 , r] = rv
          pr2              = pr3
        }else{
          chain[i + 1 , r] = tr
        }
        if(chain[i + 1 , r]==0) adre = -1 else adre = 1
        #pparam[r]        = pparam[r] + adre
        model.chain[i]   = adre*r
      }
      #----- END Riverseible jump!
      i = i + 1; j = 1
    }else{
      j = j + 1; k = k + 1; #cat('\r E:',k)
      if ( j > 1500 ) stop('error in model! Maybe you need a different set of initial values or model parameters.',call. = 0)
    }

  }
  # ---
  #close(pb)
  cat('\r\n')

  if (k > 0) {cat('\n There are ',k,' ignored values in the process!')} ; error = k
  return(list( chain = chain     , RejAcc = ReAc    ,
               error = error     , cov.m = cov.m    ,
               minf = functionat , minState = state ,
               x = x             , y = y            ,
               lb = lb           , lq = lq          ,
               logit = logit     , RJ = RJ          ,
               model.chain = checkin(model.chain,interval ) )
  )
}



# semi final function (after is the final algorithm)
par.bayesian.tot.dw <- function(formula = NA,
                                data,
                                initials    , burn.in         ,
                                q.par       , b.par           ,
                                dist.q      , dist.b          ,
                                para.q      , para.b          ,
                                dist.l      , l.par           ,
                                fixed.l     , logit           ,
                                v.scale     , iterations      ,
                                cov.m       ,
                                sampling = c('indp','syst','bin'),
                                ###
                                penalized                        ,
                                lq.extra = 1 , lbeta.extra = 1   ,
                                prior.function = prior.tot.dw,
                                RJ = FALSE,
                                ...
)
{

  #library('MASS')
  if(RJ & (para.q || para.b)){
    #message('Reversible-Jump in action ... ')
    penalized = FALSE ; fixed.l = -1
  }else{
    RJ = FALSE
  }
  cat('\n\n============================== Sampler configuration ============================== \n')
  cat('Iterations:',iterations,'\t|', 'Data:',is.data.frame(data),	 '\t \t|','Length of Initials:',length(initials),'\n')
  cat('RegQ:',para.q==1,	 '\t \t|','RegB:',para.b==1,				       '\t \t|','Formula:',is(formula,"formula"),'\n')
  # cat('Q.pars:',q.par,	 '\t \t|','Q.Dist:',is.function(dist.q),	 '\t \t|','Lambda.Dist:',deparse(quote(dist.l)),'\n')
  # cat('B.pars:',b.par,	 '\t \t|','B.Dist:',deparse(quote(dist.b)),'\t \t|','Lambda.par:' ,l.par,'\n')
  cat('Logit:',logit==1,	 '\t \t|','Scale:',v.scale,	               '\t\t| Rev.Jumps:',RJ==1,'\n')
  cat('Penalized:',penalized==1,	 '\t|','Fixed.penalty:',fixed.l>0,' |\t','\n')
  # cat('\n Q.prior:',substr(deparse(dist.b),0,30),'...\n')
  # cat('\n B.prior:',substr(deparse(dist.l),0,30),'...\n')
  # cat('\n Lambda hyper prior:',substr(deparse(dist.q),0,30),'...\n')
  cat('----------------------------------------------------------------------------------   \n' )
  cat('Proposal (1=Including covariates,2=Uniform,3=Laplace,>3=Gaussian):',cov.m,'\n' )
  cat('----------------------------------------------------------------------------------   \n' )
  cat('Chain summary (bin=Burn-in, syst=Systematic, indp=Independent):',sampling,'\n' )
  cat('----------------------------------------------------------------------------------   \n' )
  cat('* If Penalized=TRUE then you need to set all distributions.                \n' )
  cat('----------------------------------------------------------------------------------   \n' )
  cat('* If RJ=TRUE then Penalized is automatically set to FALSE and fixed.l diactivates. \n' )
  cat('__________________________________________________________________________________  \n' )


  # Start the clock!
  ptm <- proc.time()

  m.chain=MH.tot.dw(formula = formula , data=data,
                    startvalue = initials,
                    para.b = para.b , para.q = para.q  ,
                    dist.q = dist.q , dist.b = dist.b  ,
                    q.par  = q.par   , b.par = b.par   ,
                    dist.l =  dist.l, l.par  = l.par   ,
                    fixed.l = fixed.l                  ,
                    cov.m = cov.m                      ,
                    v.scale = v.scale                  ,
                    iterations = iterations            ,
                    penalized =  penalized             ,
                    lq.extra = lq.extra                ,
                    lbeta.extra = lbeta.extra          ,
                    prior.function = prior.function    ,
                    logit = logit , RJ = RJ
  )


  unused.sample = round(iterations * burn.in)
  #----- SAMPLING
  if(sampling =='syst'){
    #cat('\n Systematic sampling applied.\n')
    samp.seq   = seq(2,iterations,length.out = iterations-unused.sample)
    samp.seq   = round(samp.seq)
    chain      = m.chain$chain      [samp.seq,]
    model.chain= m.chain$model.chain[samp.seq ]
  }else if(sampling == 'indp'){
    #cat('\n Independent sampling applied. \n')
    samp.seq   = sample(x = 2:iterations,size = iterations-unused.sample)
    chain      = m.chain$chain      [samp.seq,]
    model.chain= m.chain$model.chain[samp.seq ]
  }else{
    #cat('\n Burn-in applied. \n')
    chain      = m.chain$chain      [-(1:unused.sample),]
    model.chain= m.chain$model.chain[-(1:unused.sample) ]
  }

  acceptance = apply(chain,2, function(x){1-mean(duplicated(x[x!=0]))}) * 100

  # Stop the clock
  cat('\n Procedure finished in ',( proc.time() - ptm  )[1],' seconds. \n')

  out = list(chain           =  chain              ,
             formula         = formula                    ,
             acceptance.rate = acceptance                 ,
             data = data        , iterations = iterations ,
             q.par = q.par      , b.par = b.par           ,
             initials = initials, burn.in = burn.in       ,
             dist.q =  dist.q   , dist.b =  dist.b        ,
             para.b = para.b    , para.q = para.q         ,
             dist.l =  dist.l   , l.par  = l.par          ,
             fixed.l = fixed.l                            ,
             v.scale = v.scale                            ,
             sampling = sampling                          ,
             RejAccChain = m.chain$RejAcc                 ,
             cov.m = cov.m                                ,
             error    = m.chain$error                     ,
             x  =m.chain$x , y = m.chain$y                ,
             minf = m.chain$minf , minState = m.chain$minState,
             lq.extra = lq.extra               ,
             lbeta.extra = lbeta.extra         ,
             penalized = penalized             ,
             prior.function = prior.function   ,
             lb = m.chain$lb , lq = m.chain$lq ,
             model.chain = model.chain         ,
             logit = logit, RJ = m.chain$RJ    ,
             duration =  proc.time() - ptm
  )
  return(out)
}

plot.bdw = function(x,est = Mode, prob = 0.95,
                    adj = 2, r.outliers = TRUE,
                    density = FALSE, exc.tun = FALSE ,
                    ...){
  dw.object = x
  bdw.plot( dw.object$chain,
            estimate.statistic = est , adj = adj ,
            remove.outliers = r.outliers, sampling = TRUE ,
            ...)
  cat('\n =======',round(prob*100),'% Confidence interval ======= \n')
  dw.HPDinterval   ( dw.object,sub='Parameters',prob = prob,
                     bw = adj, remove.outliers = r.outliers,
                     density = density,
                     exclude.tuning = exc.tun,...)
}


# Draw a plot for bayesian object
bdw.plot <- function(dw.object,
                     estimate.statistic = Mode,
                     adj = 1 , truepar = NA, remove.outliers = TRUE,
                     sampling = TRUE  ,...)
{
  if(requireNamespace('coda') == FALSE){
    stop ('Package Coda is needed in order to execute this function!' ,call. = FALSE)
  }
  #dw.object = dw.object$chain
  truepar= as.matrix(truepar)
  chain  = dw.object$chain
  if(dw.object$fixed.l > 0) {chain = chain[,-ncol(chain)]}
  #### initials ####
  acc    = dw.object$acceptance.rate; acc=round(acc,2)
  ncc    = ncol(chain)
  nrc    = nrow(chain)
  nctp   = ncol(truepar)
  col.names = colnames(truepar)
  names     = chain.name(ncc,dw.object$lq,dw.object$lb)
  names2    = extract.vars(dw.object,reg=TRUE)
  symb.p    = names$symb.p
  symb.p2   = names2$symb.p
  scl    = dw.object$v.scale;  scl=round(scl,3)
  lw     = 3
  col    = 3
  lty    = 3
  lsize  = 2
  #int.shift = 1.1

  if(all(col.names!="")) {cname.chi =col.names}else{cname.chi='E/T'}
  if(length(scl) == 1) scl = rep(scl,ncc)
  if(exists('n.repeat',where = dw.object) == TRUE) {nrep=dw.object$n.repeat}else{nrep=1}

  par(mfrow = c(2,3),cex=.75)
  for(i in 1:ncc){
    chi = chain[,i]
    chi = na.omit(chi)
    if((sampling == TRUE) && (nrep>1) ){
      schi = sample(chi,nrc/nrep)
      pchi = chi[seq(1,nrc,length.out = nrc/nrep )]
      rnd.rep = round(runif(1,0,nrep-1)) /nrep * nrc
      achi =  (rnd.rep+1):(rnd.rep + nrc/nrep)
    }else{
      schi = pchi = achi = chi
    }
    if(dw.object$RJ){
      schi = schi[schi!=0]
      pchi = pchi[pchi!=0]
      achi = achi[achi!=0]
      chi  =chi  [chi !=0]
    }
    m = estimate.statistic(schi);#m=round(m,2)
    if(remove.outliers == TRUE) {
      schi = remove_outliers(schi)
      schi = na.omit(schi)
    }
    #schi = c(min(schi)*4/5,schi,max(schi)*6/5)
    #pchi = c(min(pchi)*4/5,pchi,max(pchi)*6/5)
    ac=acc[i]
    hist(schi, add=0, probability = TRUE,
         breaks = 30,
         #main = paste('Posterior of  parameter ',i,'/',ncc, ' acc = ',ac, '%',sep = '') ,
         main = symb.p[i],
         xlab=paste('Bayesian estimation = ',round(m,2), ', Scale = ',scl[i], ', AcR = ', round(ac,5))
         #xlim=c(min(schi)[1]*int.shift,max(schi)[1]*int.shift)
    )
    lines(density(schi,adjust  = adj),col='black')
    if(nrow(truepar) >= i && all(is.na(truepar) == FALSE) ){
      abline(v=truepar[i,],
             col=col:(col+nctp-1),
             lty=lty:(lty+nctp-1),
             lw=lw
      )
      legend('topleft',
             legend = c(
               paste(1:nctp, '.',cname.chi,':',round(truepar[i,],2)),
               paste('Bayesian:',round(m,2))
             ),
             seg.len = lsize,
             #fill = c(col:(col+nctp-1),2),
             lty  = c(lty:(lty+nctp-1),2),
             col  = c(col:(col+nctp-1),2),
             bg   = 'transparent' ,
             lwd=lw-1)
    }
    abline(v = m , col = 2 , lw=lw , lty = 2)
    abline(...)
    if (nrc > 50000){ ptype = 'b'} else {ptype='l'}
    plot(pchi,xlab='Index',
         #ylim = c(min(pchi)[1]*int.shift,max(pchi)[1]*int.shift),
         ylab='Estimation', type = ptype ,
         #main=paste('Convergence of parameter ',i ,'/', ncc,sep = ''),
         main = symb.p[i],
         col='gray' ,
         pch ='.' )
    if(nrow(truepar) >= i && all(is.na(truepar[i,]) == FALSE) ){
      abline(h=truepar[i,],
             col=col:(col+nctp-1),
             lty=lty:(lty+nctp-1),
             lw=lw)
      legend('topleft',
             legend = c(
               paste(1:nctp, '.',cname.chi,':',round(truepar[i,],2)),
               paste('Bayesian:',round(m,2))
             ),
             seg.len = lsize,
             #fill = c(col:(col+nctp-1),2),
             lty  = c(lty:(lty+nctp-1),2),
             col  = c(col:(col+nctp-1),2),
             bg   = 'transparent' , lwd=lw-1)
    }
    abline(h = m , col=2 , lw=lw,lty=2)
    #ACF must be computed from the ORIGINAL chain!
    acf=acf(achi,plot = FALSE , lag.max = ceiling(sqrt(nrc/nrep)))
    if(any(is.nan(acf$acf)== TRUE) ) acf=rep(0,ceiling(sqrt(nrc/nrep)))
    #plot(acf,type='h',xlab='Lag',ylab='ACF',main=paste('ACF of parameter ',i ,'/', ncc,sep = ''),col='gray')
    plot(acf,type='h',xlab='Lag',ylab='ACF',main = symb.p[i],col='gray')
    cat('\r',i,' of ',ncc,' plot completed.')
  }
  par(mfrow=c(1,1))
  if(dw.object$RJ){
    model.chain = dw.object$model.chain
    h=table(model.chain)
    h2=length(h)/2
    p=(h[(h2+1):(2*h2)]/(h[h2:1]+h[(h2+1):(2*h2)]))
    bplot=  barplot(-h[h2:1],col=2,ylim=c(-1.5*max(abs(h)),2.5*(max(abs(h)))),xpd=1, xaxt='n',density = 80)
    barplot(h[(h2+1):(2*h2)],add=1,col=3,xpd=1, xaxt='n')
    #axis(1, at=bplot ,labels=paste('V.',1:(h2),' | ',round(p,2),sep=''),las=3,cex.axis=.8)
    abline(v=bplot,col='gray',lty=3)
    axis(1, at=bplot ,labels=symb.p2,las=3,cex.axis=.8)
    legend('topright', horiz = 1, legend = c('Accepted','Rejected'),col = c(3,2),fill = c(3,2),density = c(100,80))
    text (x=bplot,y=1.5*max(abs(h)),labels = paste('+',round(p*100,1),'%'))
    title(main = 'Reversible-Jumps MH model selection')
  }

  #crosscorr.plot(as.mcmc(chain),main='Chain correlation diagram')
}


#############################
bay.laplace.calc <- function (logpost, mode, ...)
{
  options(warn = -1)
  fit = optim(mode, logpost, gr = NULL, ..., hessian = TRUE,
              control = list(fnscale = -1) )
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  p = length(mode)
  int = p/2 * log(2 * pi) + 0.5 * log(det(h)) + logpost(mode,
                                                        ...)
  stuff = list(mode = mode, var = h, int = int, converge = fit$convergence ==
                 0)
  return(stuff)
}
########BF
dw.BF <- function (dw.object,  est.stat = Mode , ...)
{
  options(warn = -1)
  if (dw.object$chain$RJ){
    mode.v = apply (dw.object$chain$chain , 2 , function(x) { x=x[x!=0] ; est.stat(x)} )
  }else{
    mode.v = apply (dw.object$chain$chain , 2 ,  est.stat )
  }
  fit = optim( par = mode.v , fn = posterior.tot.dw ,
               gr = NULL, ..., hessian = TRUE      ,
               control = list(fnscale = -1)        ,
               para.b  = dw.object$chain$para.b    ,
               para.q  = dw.object$chain$para.q    ,
               dist.q  = dw.object$chain$dist.q    ,
               dist.b  = dw.object$chain$dist.b    ,
               q.par   = dw.object$chain$q.par     ,
               b.par   = dw.object$chain$b.par     ,
               dist.l  = dw.object$chain$dist.l    ,
               l.par   = dw.object$chain$l.par     ,
               x       = dw.object$chain$x         ,
               y       = dw.object$chain$y         ,
               logit   = dw.object$chain$logit     ,
               lbeta.extra = dw.object$chain$lbeta.extra,
               lq.extra = dw.object$chain$lq.extra   ,
               penalized = dw.object$chain$penalized ,
               prior.function = dw.object$chain$prior.function,
               fixed.l = dw.object$chain$fixed.l
  )
  options(warn = 0)
  mode = fit$par
  h = -solve(fit$hessian)
  p = length(mode)
  int = p/2 * log(2 * pi) + 0.5 * log(det(h)) +
    posterior.tot.dw(par = mode,
                     para.b  = dw.object$chain$para.b    ,
                     para.q  = dw.object$chain$para.q    ,
                     dist.q  = dw.object$chain$dist.q    ,
                     dist.b  = dw.object$chain$dist.b    ,
                     q.par   = dw.object$chain$q.par     ,
                     b.par   = dw.object$chain$b.par     ,
                     dist.l  = dw.object$chain$dist.l    ,
                     l.par   = dw.object$chain$l.par     ,
                     x       = dw.object$chain$x         ,
                     y       = dw.object$chain$y         ,
                     logit   = dw.object$chain$logit     ,
                     lbeta.extra = dw.object$chain$lbeta.extra,
                     lq.extra = dw.object$chain$lq.extra ,
                     penalized = dw.object$chain$penalized,
                     fixed.l = dw.object$chain$fixed.l,
                     prior.function = dw.object$chain$prior.function)

  stuff = list(mode = mode,
               var = h,
               bf = int,
               converge = fit$convergence == 0
  )
  return(stuff)
}

##############################################################################
# Estimate Deviance Information Criterion (DIC)
#
# References:
#   Bayesian Data Analysis.
#   Gelman, A., Carlin, J., Stern, H., and Rubin D.
#   Second Edition, 2003
#
#   Bayesian predictive information criterion for the evaluation of
#     hierarchical Bayesian and empirical Bayes models.
#   Ando, T.
#   Biometrika, 2007
#
# Input:
#   x       : matrix of posterior samples
#   lik     : vector of the likelihood of the posterior samples
#   lik.fun : function that calculates the likelihood
#   ...     : other parameters that are passed to 'lik.fun'
#
# Output:
#   list()
#     DIC   : Deviance Information Criterion
#     IC    : Bayesian Predictive Information Criterion
#     pD    : Effective number of parameters (pD = Dbar - Dhat)
#     pV    : Effective number of parameters (pV = var(D)/2)
#     Dbar  : Expected value of the deviance over the posterior
#     Dhat  : Deviance at the mean posterior estimate
##############################################################################
# calc.dic <- function(x,lik,lik.fun,...) {
#   D.bar <- -2*mean(lik)
#   theta.bar <- summary(x)$statistics[,"Mean"]
#   D.hat <- -2*lik.fun(theta.bar,...)
#   pD <- D.bar - D.hat
#   pV <- var(-2*lik)/2
#   list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
# }
dw.info.Criteria = function (dw.object,  est.stat = Mode ,
                             prob = .95, sample.p = TRUE , ...){
  lik   =       l.tot.dw(dw.object , est.stat = NA   , chain.values = TRUE  , sample.p = sample.p)
  D.bar = -2 * mean(lik)
  D.hat = -2 * (l.tot.dw(dw.object , est.stat = mean , chain.values = FALSE , sample.p = sample.p))
  pD    = D.bar - D.hat
  pV    = var(-2*lik)/2
  #############################
  L  = ( l.tot.dw(dw.object , est.stat = est.stat , sample.p = sample.p) )
  k  = dw.object$chain$lb + dw.object$chain$lq
  if (dw.object$chain$penalized || dw.object$chain$RJ){
    zr = dw.HPDinterval(dw.object = dw.object ,prob =  prob ,draw = FALSE , density = FALSE )
    df = k - sum(zr[1:k,4])
  }else{
    df =k
  }

  n  = length(dw.object$chain$y)
  BIC  = -2*L + k *  log(n)
  CAIC = -2*L + k * (log(n)+1)
  AIC  = -2*L + k *2
  AICc = AIC + 2*k*(k+1)/(n-k+1)
  QIC  = -2 * (L - k*log(k))/n
  stuff = list(    AIC  = AIC ,
                   AICc = AICc,
                   BIC  = BIC ,
                   QIC  = QIC ,
                   CAIC = CAIC ,
                   df   = df   ,
                   DIC=pD+D.bar,
                   BPIC=2*pD+D.bar,
                   pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat,
                   loglik = L
  )
  return(stuff)

}
summary.bdw = function(object,est = Mode ,  prob =.95 , samp = TRUE , ...){
  dw.object = object
  cat('\r Please wait ...')
  BF = dw.BF (dw.object, est.stat = est , ... )
  cr = dw.info.Criteria ( dw.object ,est.stat = est, prob = prob ,sample.p = samp , ...)
  info = dw.object$chain
  cat('\r\n',
      '============================== Sampler ================================ \n' ,
      '    Iterations : ' , info$iterations    , '\t Logit : ' , info$logit==1 , '\t Scale : ' ,round(info$v.scale,5) , '\n' ,
      '    Rev.Jump   : ' , info$RJ==1         , '\t RegQ  : ' , info$para.q==1 ,'\t RegB  : ' ,info$para.b==1 , '\n' ,
      '    Penalized  : ' , info$penalized==1  , '\t Fixed.penalty : ' , info$fixed.l>0, '\n' ,
      '============================ Model Summary ============================ \n'  ,
      '    AIC : ' , cr$AIC , '\t AICc : ' , cr$AICc , '\t BIC     : ' ,  cr$BIC , '\n' ,
      '    QIC : ' , cr$QIC , '\t CAIC : ' , cr$CAIC , '\t LogPPD  : ' ,   BF$bf , '\n' ,
      '    DIC : ' , cr$DIC , '\t PBIC : ' , cr$BPIC , '\t df      : ' ,   cr$df , '\n' ,
      '=======================================================================  \n'
  )
  #NextMethod('t')


}
######### likelihood
l.tot.dw <- function(dw.object , est.stat = Mode , chain.values = FALSE , sample.p )
{
  para.q   = dw.object$chain$para.q
  para.b   = dw.object$chain$para.b
  l.par    = dw.object$chain$lb + dw.object$chain$lq
  x        = dw.object$chain$x
  y        = dw.object$chain$y
  logit    = dw.object$chain$logit
  lq.extra = dw.object$chain$lq.extra
  lbeta.extra   = dw.object$chain$lbeta.extra


  if(chain.values == FALSE){
    par      = apply(dw.object$chain$chain[,1:l.par] , 2 , est.stat)
    par      = matrix(par,nrow=1)
  }else{
    par = dw.object$chain$chain[,1:l.par]
    ny  = nrow(par)
    smp = sample(1:ny,round(ny*sample.p))
    par = par[smp,]
  }

  if(      (para.b == FALSE) && (para.q == TRUE)){
    loglik = apply(X =  par , MARGIN = 1 , FUN = loglik.q.dw     , x=x , y=y , lq.extra=lq.extra , lbeta.extra=lbeta.extra , logit = logit)
  }else if((para.b == TRUE ) && (para.q == FALSE)){
    loglik = apply(X =  par , MARGIN = 1 , FUN = loglik.beta.dw  , x=x , y=y , lq.extra=lq.extra , lbeta.extra=lbeta.extra , logit = logit)
  }else  if((para.b == TRUE) && (para.q == TRUE)) {
    loglik = apply(X =  par , MARGIN = 1 , FUN = loglik.qbeta.dw , x=x , y=y , lq.extra=lq.extra , lbeta.extra=lbeta.extra , logit = logit)
  }else{
    loglik = apply(X =  par , MARGIN = 1 , FUN = loglik.sim.dw   , x=y       , lq.extra=lq.extra , lbeta.extra=lbeta.extra , logit = logit)
  }
  return(  loglik   )
}
########### HPDinterval
dw.HPDinterval <- function (dw.object , prob = 0.95 , bw = 2.5 ,
                            remove.outliers = TRUE ,
                            draw = TRUE , density = FALSE , exclude.tuning = FALSE,
                            var.lab = NA  ,fixedcol = 1 , ...){
  if(requireNamespace('coda') == FALSE){
    stop ('Package Coda is needed in order to execute this function!' ,call. = FALSE)
  }
  lq = dw.object$chain$lq
  lb = dw.object$chain$lb
  chain = dw.object$chain$chain
  nc.chain = ncol(chain)
  if((nc.chain> lq+lb) && (exclude.tuning==TRUE)){
    chain    = chain[,1:(lb+lq)]
    nc.chain = nc.chain -1
  }

  if (remove.outliers == TRUE){
    chi = apply(chain,2,remove_outliers)
  }else{
    chi = chain
  }
  if(dw.object$chain$RJ){
    est   = apply(chi,2,function(x)   {x=x[x!=0] ; x=na.omit(x); Mode(x = x)})
    res   = apply(chi,2,function(x)   {x=x[x!=0] ; x=na.omit(x); if(length(x)>5) {coda::HPDinterval(coda::as.mcmc(x),prob = prob)}else{c(-10^-8,10^-8)} })
    if(density==TRUE) {d     = apply(chi,2,function(x,dw){x=x[x!=0] ; x=na.omit(x); density(na.omit(x),adjust  = bw)})}
  }else{
    est   = apply(chi,2,function(x)   {x=na.omit(x); Mode(x = x)})
    res   = apply(chi,2,function(x)   {x=na.omit(x); if(length(x)>5) {coda::HPDinterval(coda::as.mcmc(x),prob = prob)}else{c(-10^-8,10^-8)} })
    if(density==TRUE) {d     = apply(chi,2,function(x,dw){x=na.omit(x); density(na.omit(x),adjust =  bw)})}
  }
  #acc.r = round(dw.object$chain$acceptance.rate,2)
  chname = chain.name(nc.chain,lq,lb,var.lab)
  symb   = chname$symb
  symb.p = chname$symb.p
  las    = chname$las
  marb   = chname$marb

  colnames(res)           = symb
  rownames(res)  = c('lower','upper')
  res=t(res)
  lower = res[,1]
  upper = res[,2]
  Zero.included = apply(res,1,function(x){prod(sign(x))<=0})


  if (draw){
    if(lb+lq > 20){cexv= 1/log(nrow(res)+1)+.5}else{cexv=1}
    par(mar =c(marb,4,3,2),cex=cexv)
    nx = nrow(res)
    matplot( cbind ( res , est)                     ,
             col=c(2,3,1) , lty = c ( 3 , 3 , 1 )   ,
             type = c('n' , 'n' , 'p')              ,
             xaxt = 'n'                             ,
             lwd = 1                               ,
             pch = c( 0 , 0 , 1 )                   ,
             ylim = c(min(res,est,na.rm = 1)-.5,max(res,est,na.rm = 1)+.5),
             xlim = c(1, ncol(chain)+1)             ,
             #xlab =  ' \r\n Variables'                   ,
             ylab = paste(round(prob*100,3),'% HPD interval ')
             #xlab = paste('Variables, Acceptance rate = ', round(dw.object$chain$acceptance.rate[1],2),'%' )
    )
    abline ( h = 0 , v =1:nx, col = 'cornsilk2' , lty = 6)
    title(...)
    col.count = 1
    for(k in 1:nx){
      if(abs(res[k,1]-res[k,2]) > 10^-6){
        arrows(k                     ,
               res[k,1]              ,
               k                     ,
               res[k,2]              ,
               length = 0.1          ,
               angle  = 90           ,
               code   = 3            ,
               #col   = c ( col.count%%2+2 , col.count%%2+2 )    ,
               col    =  c(fixedcol*Zero.included[k]+2,fixedcol*Zero.included[k]+2) ,
               lwd    = 2,
               lty    = fixedcol*Zero.included[k]*2+1
        )
        col.count = col.count + 1
      }
    }
    if(density == TRUE){
      for(i in 1:ncol(chain)){
        #with(d,expr = lines(d[[i]]$y/2.65+i,(d[[i]]$x) , lwd =1 ,col=(i)%%2+2 ,lty = i%%2+1))
        with(d,expr = lines(d[[i]]$y/sd(d[[i]]$y)/4+i,(d[[i]]$x) , lwd =1 ,col= fixedcol*Zero.included[i]+2 ,lty = fixedcol*Zero.included[i]*2+1 ) )
        abline(v = i  , lty = 3 , col = 'gray')
      }
    }
    axis(1, 1:nx , labels = symb.p , las = las )
    if(dw.object$chain$fixed.l <= 0  && dw.object$chain$penalized==TRUE){abline (v = ncol(chain)  , lty = 5 , col = 'gray') }
  }

  res=cbind(lower,est,upper,Zero.included)
  return(res)
}


#### Plot density for a fixed x
dw.bay.plot.density <- function (x, dw.object , xmax = 15 ,
                                 smooth = .5 , draw = TRUE , var.lab = NA , ...){

  ncx = ncol(dw.object$chain$x)
  logit = dw.object$chain$logit
  if(     dw.object$chain$para.q == TRUE && dw.object$chain$para.b == FALSE){
    a = 1:ncx
    b = ncx+1
    q.t    = cbind(1 , x) %*% as.matrix(dw.object$res[a])
    b.t    = dw.object$res[b]
    if(logit){
      q      = exp(q.t - log ( 1 + exp(q.t) ) )
    }else{
      q      = exp(-exp(q.t ))
    }
    beta   = rep(b.t , 1)
  }else if(dw.object$chain$para.q == FALSE && dw.object$chain$para.b == TRUE){
    a = 1
    b = 2:(ncx+1)
    q.t    = dw.object$res[a]
    b.t    = cbind( 1 , x ) %*% as.matrix(dw.object$res[b])
    q      = rep( q.t , 1)
    beta   = exp(b.t)
  }else if(dw.object$chain$para.q == TRUE && dw.object$chain$para.b == TRUE){
    a = 1:ncx
    b = (ncx+1):(2*ncx)
    q.t    = cbind( 1 , x ) %*% as.matrix(dw.object$res[a])
    b.t    = cbind( 1 , x ) %*% as.matrix(dw.object$res[b])
    if(logit){
      q      = exp(q.t - log ( 1 + exp(q.t) ) )
    }else{
      q      = exp(-exp(q.t ))
    }
    beta   = exp(b.t)
  }else{
    a = 1
    b = 2
    q     = rep (dw.object$res[a] , 1)
    beta  = rep (dw.object$res[b] , 1)
  }

  res = c()
  d = seq(from = 0 , to = xmax , by = 1 )
  for(i in d ){
    res[i+1] = ddw(i , q = q , beta = beta )
  }
  res = smooth.spline (res , spar = smooth )
  if(draw == TRUE){
    plot( res$x , res$y , type = 'l' ,main = 'PMF' , xlab ='x' , ylab = 'pmf' , ...)
  }
  return(list ( y = res$y , range = res$x , q = q , beta = beta ) )
}




# multicore
bdw.mc <- function(dw.object     ,
                          n.repeat = 10 ,
                          cores = 0   )
{
  if(requireNamespace('foreach')    == FALSE |
     requireNamespace('doParallel') == FALSE |
     requireNamespace('BDWreg')        == FALSE ) stop ('Package needed!')
  #require('foreach')
  # Register all cores
  cat('\r Processing ...')
  if (cores<1) cores = parallel::detectCores()
  doParallel::registerDoParallel(cores = cores)
  ini   = dw.object$chain$initials
  dw.temp.object = dw.object
  ###################
  stime = proc.time()
  i=1
  ch4   = foreach(i = 1:n.repeat , .inorder = FALSE  ) %dopar% {
    dw.temp.object$chain$initials   =  ini/i
    do.call(par.bayesian.tot.dw  ,  dw.temp.object$chain)
  }
  message('Algorithm executed in ',round((proc.time()-stime)[3]) ,' seconds.')


  chain           = ch4[[1]]$chain
  acceptance.rate = ch4[[1]]$acceptance.rate
  minf            = ch4[[1]]$minf
  minState        = ch4[[1]]$minState
  RejAccChain     = ch4[[1]]$RejAccChain
  error           = ch4[[1]]$error
  initials        = ch4[[1]]$initials
  model.chain     = ch4[[1]]$model.chain
  ###########
  foreach::foreach(i = 2:n.repeat,.inorder = TRUE) %do% {
    chain           = rbind(chain        , ch4[[i]]$chain)
    acceptance.rate = (acceptance.rate   + ch4[[i]]$acceptance.rate)
    minf            = c(minf             , ch4[[i]]$minf)
    minState        = rbind(minState     , ch4[[i]]$minState)
    RejAccChain     = c(RejAccChain      , ch4[[i]]$RejAccChain)
    error           = c(error            + ch4[[i]]$error)
    initials        = rbind(initials     , ch4[[i]]$initials)
    model.chain     = c(model.chain      , ch4[[i]]$model.chain)
    cat('\r merging ...',i,' of ' ,n.repeat)
  }
  ### close connection here
  message('Closing connections ...')
  base::closeAllConnections()
  message('End.')

  dw.object$chain$chain=chain
  dw.object$chain$acceptance.rate=acceptance.rate/n.repeat
  dw.object$chain$minf=minf
  dw.object$chain$minState=minState
  dw.object$chain$RejAccChain=RejAccChain
  dw.object$chain$error=error
  dw.object$chain$initials=initials
  dw.object$chain$n.repeat=n.repeat
  dw.object$chain$model.chain=model.chain
  dw.object$chain$iterations=dw.object$chain$iterations*n.repeat
  dw.object$chain$duration = (proc.time()-stime)
  dw.object$chain$cores=cores
  dw.object$all = ch4
  class(dw.object) = 'bdw'
  return(dw.object  )

}


