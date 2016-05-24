`pwlet2freqs` <-
function(noctave, nvoice,  dt, flip=TRUE, tab.FREQ, plot=FALSE,  perc = 0.85)
{
  if(missing(perc)) {  perc = 0.85 }
  
  i1 = sort(rep(c(1:noctave), times=nvoice))
  jj = rep(c(0:(nvoice-1)), times=noctave)

  sa = 2^(i1+jj/nvoice)

  efs = scal2freqs(sa, dt)

  if(flip==TRUE)
    {
      efs = rev(efs)
    }

  I1 = matrix(rep(efs, times=length(tab.FREQ)), ncol=length(tab.FREQ))
  I2 = matrix(rep(tab.FREQ, times=length(efs)), ncol=length(tab.FREQ), byrow=TRUE)

  IA = apply(abs(I1-I2), 2, which.min)


  Iat = nvoice*(log2(sa[IA])-1)/(nvoice*noctave)

   why   = RPMG::RESCALE(Iat , 0 , perc , 0, 1 )

  
  if(plot==TRUE)
    {
      abline(h=why, lty=2, col=rgb(0.5, 0.5, 0.5) )
      alabs = as.character(tab.FREQ)
      ###  alabs[length(alabs)] = paste(sep=" ", alabs[length(alabs)], "Hz")
      alabs[length(why)]  = paste(alabs[length(why)], "Hz")
      axis(side = 4, at=why, labels=alabs)
      ## mtext(side=4, line=1, at = why[length(why)] , text="Hz", adj=(-1.1))
    }

  FOUT = why
  invisible(list(why=why, Iat=Iat, efs=efs) )
}

