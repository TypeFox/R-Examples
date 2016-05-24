chemo.fun = function(times,y,p,more=NULL)
{
  r = y;

  p = exp(p)
  
#  if(any( pars < 0)){ return(0*r) }
#  else{
  Q = p['p1']*y[,'C1'] + p['p2']*y[,'C2'];
#  Qs = Q*(Q>=p['Qstar']) + p['Qstar']*(Q<p['Qstar'])
  Qs = Q*exp(100*(Q-p['Qstar']))/(1+exp(100*(Q-p['Qstar']))) + p['Qstar']/(1+exp(100*(Q-p['Qstar'])))  
  
  
  r[,'N']  = p['delta']*(p['NI']-y[,'N']) - p['rho']*y[,'C1']*y[,'N']/(p['KC1']+y[,'N']) - p['rho']*y[,'C2']*y[,'N']/(p['KC2']+y[,'N']);
  r[,'C1'] = y[,'C1']*(p['XC']*p['rho']*y[,'N']/(p['KC1']+y[,'N']) - p['p1']*p['G']*(y[,'B']+y[,'S'])/(p['KB']+Qs) - p['delta'] );
  r[,'C2'] = y[,'C2']*(p['XC']*p['rho']*y[,'N']/(p['KC2']+y[,'N']) - p['p2']*p['G']*(y[,'B']+y[,'S'])/(p['KB']+Qs) - p['delta'] );
  r[,'B']  = y[,'B']*( p['XB']*p['G']*Q/(p['KB']+Qs) - (p['delta']+p['m']+p['lambda']));
  r[,'S']  = p['lambda']*y[,'B'] - (p['delta']+p['m'])*y[,'S'];
#}
  return(r)
}


chemo.ode = function(times,y,p)
{

  p = exp(p)
  y = exp(y)
  
  r = y;
  C = y['C1']+y['C2'];
  Q = p['p1']*y['C1'] + p['p2']*y['C2'];
#  Qs = Q*(Q>=p['Qstar']) + p['Qstar']*(Q<p['Qstar'])
  Qs = Q*exp(10*(Q-p['Qstar']))/(1+exp(10*(Q-p['Qstar']))) + p['Qstar']/(1+exp(10*(Q-p['Qstar'])))  
      
  r['N']  = p['delta']*(p['NI']-y['N']) - p['rho']*y['C1']*y['N']/(p['KC1']+y['N']) - p['rho']*y['C2']*y['N']/(p['KC2']+y['N']);
  r['C1'] = y['C1']*(p['XC']*p['rho']*y['N']/(p['KC1']+y['N']) - p['p1']*p['G']*(y['B']+y['S'])/(p['KB']+Qs) - p['delta'] );
  r['C2'] = y['C2']*(p['XC']*p['rho']*y['N']/(p['KC2']+y['N']) - p['p2']*p['G']*(y['B']+y['S'])/(p['KB']+Qs) - p['delta'] );
  r['B']  = y['B']*( p['XB']*p['G']*Q/(p['KB']+Qs) - (p['delta']+p['m']+p['lambda']));
  r['S']  = p['lambda']*y['B'] - (p['delta']+p['m'])*y['S'];

  return(list(r/y))
}

# Variable names:  N, C1, C2, B, S

# Paramter names: p1, p2, NI, delta, m, lambda, XC, KC1, KC2, rho, G, XB, KB, Qstar

