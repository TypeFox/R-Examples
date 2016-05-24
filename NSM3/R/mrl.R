mrl<-function (data, alpha, main=NULL, ylim=NULL,xlab=NULL,...) 
{
  # arguments:    
  
  #data is a vector of survival times in any order
  #(1-alpha) is the approximate coverage probability for the confidence band
  
  n = length(data)
  data = sort(data)
  S = M = Fem = quant = a = MU = ML = PM = PMU = PML = 0
  
  # calculation of S(x), M(x), and the empirical dsn at the survival times
  S[1]= n
  M[1] = mean(data)
  Fem[1] = 1
  for(i in 2:n){
    S[i]= n-(i-1)
    M[i] = (sum(data[i:n])/(n-(i-1)))-data[i-1]
    Fem[i] = (n-(i-1))/n
  }
  
  
  # Table of critical values for Hall-Wellner confidence band.
  quant = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75)
  a = c(2.807, 2.241, 1.96, 1.534, 1.149, 0.871)
  
  for(i in 1:6){
    if(quant[i]==alpha){
      aalpha = a[i]
    }
  }
  
  Dn = (aalpha*var(data)^.5)/(n^.5)
  
  # calculation of the bands ML and MU
  for(i in 1:n){
    ML[i] = M[i] - (Dn/Fem[i])
    MU[i] = M[i] + (Dn/Fem[i])
  }
  
  # calculation of PM(x), PMU(x), PML(x) for plotting.
  PM[1] = M[1]
  PM[2] = M[1] - data[1]
  PMU[1] = MU[1]
  PMU[2] = MU[1] - data[1]
  PML[1] = ML[1]
  PML[2] = ML[1] - data[1]
  
  for(i in 2:n){
    PM[2*i-1] = M[i]
    PM[2*i] = M[i] + (data[i-1]- data[i])
    PMU[2*i-1] = MU[i]
    PMU[2*i] = MU[i] + (data[i-1]- data[i])
    PML[2*i-1] = ML[i]
    PML[2*i] = ML[i] + (data[i-1]-data[i])
  }
  
  
  if(is.null(ylim))
    ylim=range(c(PM, PMU, PML))
  
  if(is.null(main))
    main= "Plot of Mean Residual Life and bounds"
  
  if(is.null(xlab))
    xlab="Time"
  
  
  x.data=seq(from=min(data), to=max(data), length.out=length(PM))
  plot(x.data,PM, type="l", ylim=ylim, xlab=xlab,main=main, ...)
  lines(x.data,PMU, lty=2)
  lines(x.data,PML, lwd=2, lty=2)
  
  
  
  
  # Create a list of objects containing the vector mean residual life PM,
  # and the vectors PMU and PML of upper and lower bounds for the mean
  # residual life, respectively.
  
  
  # To extract PM, PMU, PML, use results [[?]].
  results<-list(PM, PMU, PML)
  results
}
