`brune.search` <-
function(infreq, inspec, f1, f2, omega0, fcorn, tstar0, gamma  )
{

  ##  print(paste(sep=' ', "Bsearch", f1, f2, omega0, fcorn, tstar0, gamma))
flim = infreq>=f1&infreq<=f2
x = infreq[flim]
y = inspec[flim]
n = length(x)
fcorn = fcorn
omega0  = omega0

dgam = c(gamma-0.1*abs(gamma),gamma+2*abs(gamma), 0) 
ngam = 20
tst1 = tstar0-0.3*abs(tstar0)
tst2 = tstar0+0.3*abs(tstar0)
nstar = 20

dstar = c(tst1, tst2,0)

  sear<-.C("CALL_DGAMMA", PACKAGE = "RSEIS",
    as.double(x),
    as.double(y),
    as.integer(n) ,
    as.double(fcorn),
    as.double(omega0),
    as.double(dgam), 
    as.integer(ngam),		  
    as.double(dstar), 
    as.integer(nstar)
    )
gam3 =  sear[[6]]
tstar3 =  sear[[8]]

bruney = sear[[2]]

if(all(bruney>0) & all(y>0))
  {
    chisqrd = sum ( (log10(y) - log10(bruney))*(log10(y) - log10(bruney)) )
    
  }
else
  {
    chisqrd =NA
    
  }

return(list(omega0=omega0,tstar0=tstar3[3]  , fc=fcorn,  alpha=0, gamma=gam3[3], chisqrd=chisqrd ) )

}

