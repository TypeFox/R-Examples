`deconinst` <-
function(data, sintr=0.008, KAL = PreSet.Instr(), key=1, Calibnew=c(1,1.0, 0.0 ) , waterlevel=1.e-8)
  {

    ##  x = EJ10$dat[,2]

    #### data = deconinst(x, 0.008, KAL,1, waterlevel=1.e-8)

    
##% deconvolve instrument response
##% Usage: data=deconinst(data,Calibold,Calibnew,sintr,key,waterlevel);
##%
##% Deconvolve old instrument response from data and convolve new 
##% instrument response.  Stabalize using a waterlevel (eg = 1e-6).
##%
##%INPUT PARAMETERS:
##% data        = real matrix of time series stored by columns
##% Calibold   = complex matrix of instrument responses stored by columns
##%               There are 62 rows in a format described below.
##% Calibnew   = complex column vector describing the new instrument response
##%               If there are 62 rows it contains poles and zeros in the 
##%               format described below.
##%               If Calibnew has fewer than 62 rows, the first
##%               element is a key to the description of the new response
##%               If Calibnew(1)==3 then use a cos (hanning) taper
##%                 In this case Calibnew contains 7 numbers:
##%                 [3,waterlevel,gain,f1,f2,f3,f4]
##% sintr       = real row vector of sample intervals (s)
##%               must all be the same
##% key         = row vector containing the data columns to deconvolve
##%               for example deconvolve all m columns if key=[1:m] 
##% waterlevel  = optional parameter used to stabalize the deconvolution
##%               default value is 1.e-8
##%OUTPUT parameters:
##% data        = deconvolved data
##%

    ########  sample interval
    if(missing(sintr))
      {
        
        sintr=0.008
      }

    #########  list of decon intrument responses (poles &  zeros)
    if(missing(KAL))
      {
        KAL = PreSet.Instr()
#########> names(KAL)
#########[1] "40T"     "3T"      "L28"     "LE3D20s" "GEOSP1"

        
      }

    #########   which instrument in KAL to choose
    if(missing(key))
      {
        key=1
      }                                      

 #########
    if(missing(waterlevel))
      {
        waterlevel=1.e-8
      }                                           ##% default waterlevel

 #########  instrument to convolve to (usually set to nothing)    
    if(missing("Calibnew"))
      { 
        Calibnew = c(1,1.0, 0.0 )      
      }

    
    n=length(data)                          ##% number of data
    nn=next2(n)                             ##% next power of 2 for FFT
    ###   could use nextn(n,2) here
    f=makefreq(nn,sintr);                     ##% define frequency vector


    instnew =1
    meandsnew = 1
    
    ##  need to zero pad the data
    ## remove the mean also 
   ##  why = c(data-mean(data),rep(0,nn-n))
    why = c(data,rep(0,nn-n))


    

    DATA=fft(why)
    ###  what are these?  appropriate for Myake-jima decon?
   ## Calibnew = c(3,1.0, 0.4882812, 0, 0.0, 0.4882812)
   ## Calibnew = c(3,1.0, 0.4882812, 0, -0.4882812, 0.0)
   ## Calibnew = c(3,1.0, 0.4882812 )
   ## Calibnew = c(3,1.0, 0.4882812 )

    ##   better to try this:
   ## Calibnew  Calibnew = c(1,1.0, 0.0 )

    
    
    calibkey = Calibnew[1]
    if(calibkey==3)
      {
        meandsnew=Re(Calibnew[2])

############  use this to filter the output in the freq domain.
        
        
    ##     g=abs(f);
     ##    fcut  = Calibnew[3]*max(g)
    ##     i3=(g<=fcut);
    ##     fk = f[i3]
    ##     f1 = fk[1]
    ##     f2 = max(fk)
    ##     f3 = min(fk)
    ##     f4 = fk[length(fk)]
     ##    i1=(f>=f1&f<=f2);
    ##     i2=(f>=f3&f<=f4);
    ##     instnew=f*0;
   ##      instnew[!i3]=instnew[!i3]+1;
    ##     instnew[i1]=0.5*(1-cos(pi*(f[i1]-f1)/(f2-f1)));
   ##      instnew[i2]=0.5*(1+cos(pi*(f[i2]-f3)/(f4-f3)));


        f1=Re(Calibnew[3]);
        f2=Re(Calibnew[4]);
        f3=Re(Calibnew[5]);
        f4=Re(Calibnew[6]);
        instnew=f*0;
        g=abs(f);
        
        i1=which(g>=f1&g<=f2);
        
        i2=which(g>=f3&g<=f4);
        
        i3= which(g>=f2&g<=f3);
        
        
        instnew[i3]=instnew[i3]+1;
        
        instnew[i1]=0.5*(1-cos(pi*(g[i1]-f1)/(f2-f1)));
        instnew[i2]=0.5*(1+cos(pi*(g[i2]-f3)/(f4-f3)));
        
        
      }
    calibkey = Calibnew[1]
    if(calibkey==1)
      {
        meandsnew=Re(Calibnew[2])
        instnew=f*0;
        instnew=instnew+1
      }

    
    
    ##  disp(['Deconvolving instrument response from trace ',int2str(k)]);

##   here need to convert to real units
    
  ##     Sense = sensitivity for converting volts to m/s


    ######   if the number of poles and zeros are 0 and the
    ######   sesnitivity is non-zero, then divide by sense
    ######    and return  
    if(KAL[[key]]$np==0 & KAL[[key]]$nz==0 & KAL[[key]]$Sense != 0 )
      {
        d = data/KAL[[key]]$Sense
        return(d)
      }


    
    instold=INSTresponse(KAL,key, f, plotkey=NULL);
    
    ###plot(instold$resp, xlim=locator(2)$x)
    ###plot(f,   Mod(instold$transfer))
    ### l = locator(2)
    ###plot(f,   Mod(instold$transfer),  xlim=l$x)

    temp1= Re(instold$transfer*Conj(instold$transfer)) ;
    
    gamma=max(temp1)*waterlevel;
    
    ###  temp= Re(instnew$transfer*Conj(instold$transfer)) / (temp1+gamma);
   ###  temp=  instnew / (temp1+gamma);
###  temp=  instnew / (instold$transfer+gamma);

    ###  tempdata=Re(fft(DATA*temp,inverse = TRUE ));


    temp=instnew*Conj(instold$transfer)/(temp1+gamma)
    
    tempdata=Re(fft(DATA*temp,inverse = TRUE )/nn);
    
    da=tempdata[1:n];

    ###########   convert to m/s using the sensitivity
    
    ###########      Sense is the sensitivity of the instrument
    meandsold=KAL[[key]]$Sense

    ##########   this is the part that converts from Volts to m/s (velocity)
    ##########    to get displacement one must also integrate
    d=da*meandsnew/meandsold;
    ### l = locator(2)
        ###  plot(d, xlim=l$x)
    
    return(d)
    
}

