SEIRRV = function(t,y,p,more=NULL)
{
    beta = p['mean.beta'] + p['ampl.beta']*(1 - 2*(t%%1>0.6 & t%%1<0.75))
    lambda = (beta*y[,'I'] + p['imports'])/p['popsize']

    r = y
    names(r) = names(y)

    r[,'S'] = p['birthrate']*p['popsize']*(1-p['p']) + p['alpha.3']*y[,'R2'] +
                p['polar.prob']*p['gamma']*y[,'I'] - (lambda+p['deathrate'])*y[,'S']

    r[,'E'] = lambda*y[,'S'] + (1-p['boost.prob'])*p['foi.mod']*y[,'R2']- (p['sigma']+p['deathrate'])*y[,'E']

    r[,'I'] = p['sigma']*y[,'E'] - (p['gamma']+p['deathrate'])*y[,'I']

    r[,'R1'] = (1-p['polar.prob'])*p['gamma']*y[,'I'] + p['boost.prob']*p['foi.mod']*lambda*y[,'R2'] -
              (p['deathrate']+p['alpha.1'])*y[,'R1']

    r[,'R2'] = p['alpha.1']*y[,'R1'] + p['alpha.2']*y[,'V'] -
              (p['deathrate']+p['alpha.3']+p['foi.mod']*lambda)*y[,'R2']

    r[,'V'] = p['deathrate']*p['popsize']*p['p'] - (p['deathrate']+p['alpha.2'])*y[,'V']

    return(r)
}


SEIRR = function(t,y,p,more=NULL)
{
    beta = p['mean.beta'] + p['ampl.beta']*(1 - 2*(t%%1>0.6 & t%%1<0.75))
    lambda = (beta*y[,'I'] + p['imports'])/p['popsize']

    r = y
    names(r) = names(y)

    r[,'S'] = p['birthrate']*p['popsize'] + p['alpha.3']*y[,'R2'] +
                p['polar.prob']*p['gamma']*y[,'I'] - (lambda+p['deathrate'])*y[,'S']

    r[,'E'] = lambda*y[,'S'] + (1-p['boost.prob'])*p['foi.mod']*y[,'R2']- (p['sigma']+p['deathrate'])*y[,'E']

    r[,'I'] = p['sigma']*y[,'E'] - (p['gamma']+p['deathrate'])*y[,'I']

    r[,'R1'] = (1-p['polar.prob'])*p['gamma']*y[,'I'] + p['boost.prob']*p['foi.mod']*lambda*y[,'R2'] -
              (p['deathrate']+p['alpha.1'])*y[,'R1']

    r[,'R2'] = p['alpha.1']*y[,'R1'] -
              (p['deathrate']+p['alpha.3']+p['foi.mod']*lambda)*y[,'R2']

    return(r)
}

SEIR = function(t,y,p,more=NULL)
{
    beta = p['mean.beta'] + p['ampl.beta']*(1 - 2*(t%%1>0.6 & t%%1<0.75))
    lambda = (beta*y[,'I'] + p['imports'])/p['popsize']

    r = y
    names(r) = names(y)

    r[,'S'] = p['birthrate']*p['popsize'] + p['alpha.3']*y[,'R2'] +
                p['polar.prob']*p['gamma']*y[,'I'] - (lambda+p['deathrate'])*y[,'S']

    r[,'E'] = lambda*y[,'S'] - (p['sigma']+p['deathrate'])*y[,'E']

    r[,'I'] = p['sigma']*y[,'E'] - (p['gamma']+p['deathrate'])*y[,'I']

    r[,'R1'] = (1-p['polar.prob'])*p['gamma']*y[,'I']  -
              (p['deathrate']+p['alpha.1'])*y[,'R1']

    return(r)
}



#SEIRRV_old = function(t,y,p,more)
#{
#    params = more$p.fun(t,more$p.more)
#    beta = more$beta.fun(t,p,more$betadef)
#    lambda = (beta*I + p['i'])/params[,'N']
#
#    r = y
#    names(r) = names(y)
#    
#    r[,'S'] = params[,'nu']*params[,'N']*(1-p['p']) + p['alpha3']*y[,'R2'] + p['phi']*params[,'gamma']*y[,'I'] -
#                lambda*y[,'S'] - params[,'mu']*y[,'S']
#                
#    r[,'E'] = lambda*y[,'S'] + (1-p['xi'])*p['epsilon']*y[,'R2']- (params[,'sigma']+params[,'mu']*y[,'E']
#
#    r[,'I'] = params[,'sigma']*y[,'E'] - (params[,'gamma']+params[,'mu']*y[,'I']
#
#    r[,'R1'] = (1-p['phi'])*params[,'gamma']*y[,'I'] + p['xi']*p['epsilon']*lambda*y[,'R2'] -
#              (params[,'mu']+p['alpha1']*y[,'R1']
#
#    r[,'R2'] = p['alpha1']*y[,'R1'] + p['alpha2']*y[,'V'] -
#              (params[,'mu']+p['alpha3']+p['epsilon']*lambda)*y[,'R2']
#    
#    r[,'V'] = params[,'mu']*params[,'N']*p['p'] - (params[,'mu']+p['alpha2'])*y[,'V']
#    
#    return(r)
#}
#
#
#SEIRV = function(t,y,p,more)
#{
#    params = more$p.fun(t,more$p.more)
#    beta = more$beta.fun(t,p,more$betadef)
#
#    r = y
#    names(r) = names(y)
#    
#    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) + p['alpha']*y[,'R'] + p['phi']*params[,'gamma']*y[,'I'] + 
#                p['alpha2']*y[,'V'] - y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S'] 
#    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
#    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
#    r[,'R'] = (1-p['phi'])*params[,'gamma']*y[,'I']  - params[,'mu']*y[,'R'] - p['alpha']*y[,'R'] 
#    r[,'V'] = params[,'mu']*(1-p['p'])*params[,'N'] - p['alpha2']*y[,'V'] - params[,'mu']*y[,'V']
#
#    return(r)
#}
#
#
#
#SEIRS = function(t,y,p,more)
#{
#    params = more$p.fun(t,more$p.more)
#    beta = more$beta.fun(t,p,more$betadef)
#
#    r = y
#    names(r) = names(y)
#    
#    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) + p['alpha']*y[,'R'] + 
#                 - y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S'] 
#    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
#    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
#    r[,'R'] = params[,'gamma']*y[,'I'] - params[,'mu']*y[,'R'] - p['alpha']*y[,'R'] 
#
#    return(r)
#}
#
#
#SEIR = function(t,y,p,more)
#{
#    params = more$p.fun(t,more$p.more)
#    beta = more$beta.fun(t,p,more$betadef)
#
#    r = y
#    names(r) = names(y)
#
#    r[,'S'] = params[,'mu']*params[,'N']*(1-p['p']) +
#                 - y[,'S']*(p['i']+beta*y[,'I']) - params[,'mu']*y[,'S']
#    r[,'E'] = y[,'S']*(p['i']+beta*y[,'I']) - params[,'sigma']*y[,'E'] - params[,'mu']*y[,'E']
#    r[,'I'] = params[,'sigma']*y[,'E'] - params[,'gamma']*y[,'I'] - params[,'mu']*y[,'I']
#    r[,'R'] = params[,'gamma']*y[,'I'] - params[,'mu']*y[,'R']
#
#    return(r)
#}