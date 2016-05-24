`autoreg` <-
function(a, numf=1024 , pord = 100, PLOT=FALSE,  f1=.01, f2=50)
  {
    if(missing(numf)) {numf=1024}
    if(missing(f1)) {f1=.01}
    if(missing(f2)) {f2=50}
    if(missing(pord)) {pord = 100}
    if(missing(PLOT)) {  PLOT=FALSE }


    dts  = leests(a)
       
    fdt = 0.5*(seq(from=0, to=(numf-1), by=1))/numf / dts$dt;
    
    kout=rep(0, length=numf)

    sig = dts$y - mean(dts$y)
    n = length(sig)

    ary = .C("CALL_ARspec",  PACKAGE = "RSEIS",
      as.double(sig ), as.double(kout), as.integer(n),  as.integer(numf),  as.integer(pord) )

    amp = ary[[2]]
    
    if(PLOT==TRUE)
      {
        flag = fdt>=f1&fdt<=f2
        plot(fdt[flag], amp[flag], type='l', xlab="Freq, Hz", ylab="Amp")
        
      }
    invisible(list(amp=amp, freq=fdt))
  }

