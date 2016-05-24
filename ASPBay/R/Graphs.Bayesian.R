Graphs.Bayesian <-
function(M, burn=0, xbins=200, ORlim=c(1,5), conf.int=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95), print=TRUE)
{
  PSI <- M$OR
  Fa <- M$f_ab+M$f_aB
  Fb <- M$f_ab+M$f_Ab
  D <- M$f_ab*M$f_AB-M$f_Ab*M$f_aB
  R2 <- D^2/(Fa*(1-Fa)*Fb*(1-Fb))
  
  psi <- PSI*(PSI>=1) + 1/PSI*(PSI<1)
  fa <- Fa*(PSI>=1) + (1-Fa)*(PSI<1)

  if (burn!=0)
  {
    psi2 <- psi[-(1:burn)]
    r22 <- R2[-(1:burn)]
	fa2 <- fa[-(1:burn)]
	fb2 <- Fb[-(1:burn)]
  } else {
    psi2 <- psi
    r22 <- R2
    fa2 <- fa
	fb2 <- Fb
  }

  # r² and OR graph  
  int <- cred_int(r22,psi2,xbins,conf.int=conf.int)
  hex_r2_OR <- hexbinplot(y~x,data.frame(x=r22,y=psi2),aspect=1, xbins=xbins, style = "colorscale",colorcut=unique(int$colorcut),colramp= function(n){colorhex(n,beg=15,end=225)}, ylim=ORlim, ylab='OR', xlab=expression(r^2))
  hex_r2_OR <- graph_cred_int(hex_r2_OR,conf.int=conf.int,count=int$count)
  if (print==TRUE)
  {
    dev.new()
    print(hex_r2_OR)
  }
  
  # frequencies graph  
  int <- cred_int(fa2,fb2,xbins,conf.int=conf.int)
  hex_fa_fb <- hexbinplot(y~x,data.frame(x=fa2,y=fb2),aspect=1, xbins=xbins, style = "colorscale",colorcut=unique(int$colorcut),colramp= function(n){colorhex(n,beg=15,end=225)}, xlab='fa', ylab='fb')
  hex_fa_fb <- graph_cred_int(hex_fa_fb,conf.int=conf.int,count=int$count)
  if (print==TRUE)
  {
    dev.new()
    print(hex_fa_fb)
  }
  
  return(list(hex_r2_OR=hex_r2_OR, hex_fa_fb=hex_fa_fb))
}