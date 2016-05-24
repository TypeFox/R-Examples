`get.corner` <-
function( INfreq, INspec, dt, f1, f2, PLOT=FALSE, VERBOSE = FALSE)
  {
####  calculate the best fit Omega0, Corner freq and Tstar
    ##  for a Brune Model
    ## computations are done in the LOG-LOG domain
    ## so we fit a flat line for Omega and a linear fall off for tstar
    
    ## example: mc = get.corner(Spec$f, lspec, dt, 0.01, 15.0, PLOT=TRUE)

    
    if(missing(PLOT)) { PLOT=FALSE }
    if(missing(VERBOSE)) {  VERBOSE = FALSE}

    
    n = length(INfreq);
    corn = 0
    ave = 0
    slope = 0
    interc = 0
    K = 0

    flag = INfreq>=f1 & INfreq<=f2

    freqlim  =  INfreq[flag]
    speclim  =  INspec[flag]

    Lfreqlim = log10(freqlim)
    Lspeclim  = log10(speclim)

    n = length(freqlim);
    
  barf<-.C("CALL_DCORN", PACKAGE = "RSEIS",
    as.double(Lfreqlim),
    as.double(Lspeclim),
    as.integer(n) ,
    as.integer(K) ,
    as.double(ave) ,
    as.double(slope) ,
    as.double(interc)
    )

    K = unlist(barf[[4]])+1
    ave=unlist(barf[[5]]);
 
    slope = unlist(barf[[6]]);
    interc  = unlist(barf[[7]]);

    

    ###  now get the corner frequency

    newx = (ave-interc)/slope

    ## fc = freqlim[K]

    fc = 10^newx
    
    PI = 3.14159265358979;
    gamma = 2.0;
   
   
    
    f = fc;

   
    
    ftem = f*(PI * log10(2.718281828));
    fa = slope*log10(fc);
    fb = 0.5*log10( 1+ (f/fc)^(2*gamma));
    tem = -(fa+fb)/ftem;
    tstar0 = tem/10.;
    omega0 = 10^ave;

    if(VERBOSE)
      {
     print(paste(sep=' ', "Gcorn:",
                 formatC(ave) ,
                 formatC(omega0),
                 formatC(fc),
                 formatC(slope),
                 formatC(interc),
                 formatC(tstar0) ))

   }
    if(PLOT==TRUE)
      {

        opar = par("usr")
        par(mfrow=c(2,1))
        
        plot(Lfreqlim, Lspeclim, type='l', main="Restricted Freq", xlab="Log Freq", ylab="Log Displacement")
        abline(h=ave, col=4)
        abline(interc, slope, col=2)

        plot(log10(INfreq), log10(INspec), type='l', main="Full Freq", xlab="Log Freq", ylab="Log Displacement")
        abline(h=ave, col=4)
        abline(interc, slope, col=2)

        
       invisible( par(opar))

      }
    ret = list(corn=fc,ave=ave,slope=slope,interc=interc,tstar0=tstar0,omega0=omega0 )

    return(ret)

    ## example: mc = get.corner(Spec$f, lspec, dt, 0.01, 15.0, PLOT=TRUE)


}

