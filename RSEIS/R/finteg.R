`finteg` <-
function(data, dt)
  {

    waterlevel=1.e-8
    
    n=length(data)                          ##% number of data
    nn=next2(n)                             ##% next power of 2 for FFT
    f=makefreq(nn, dt);                     ##% define frequency vector
    why = c(data-mean(data),rep(0,nn-n))
    DATA=fft(why)


    temp1 =  dt*f*complex(real=0, imaginary=1)*2*pi
    ## gamma=max(Re(temp1))*waterlevel;
    temp1[Im(temp1)==0.0 ] = dt*complex(real=0, imaginary=waterlevel)*2*pi

    
    
    temp=complex(real=1, imaginary=0.0)/(temp1)
    
    tempdata=Re(fft(DATA*temp,inverse = TRUE )/nn);

    da=tempdata[1:n];
    return(da)

    
  }

