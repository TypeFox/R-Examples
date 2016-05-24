"sigplot" <-
function(mat=.HPROFIL, column=5, sigcol=7, tit="Title of Graph", xtit="Step", ytit="Uncertainty Coefficient", override=FALSE, lowy=0, highy=1.0) {

  # STORE THE NUMBER OF STEPS
  maxsteps <- dim(mat)[1]

  if(!override) {
    minval <- 0
    maxval <- ceiling(max(mat[,column]))
  }
  else {
    minval <- lowy
    maxval <- highy
  }
  
  dat <- mat[,column]
  
  # SIGNIFICANCES (0 = SIG, 1 = NOT SIG)
  sig <- mat[,sigcol]
  sig[sig < 0.05] <- 0
  sig[sig >= 0.05] <- 1
  
  # PREPARE PLOT AREA
  plot(1:maxsteps, mat[,column], type="n", xlim=c(1,maxsteps), ylim=c(minval,maxval), xlab=xtit, ylab=ytit, main=tit)
  par(new=T)
  
  # DRAW BLACK LINE
  plot(1:maxsteps, mat[,column], type="l", col="black", xlim=c(1,maxsteps), ylim=c(minval,maxval), xlab=xtit, ylab=ytit)
  par(new=T)
  
  # PLOT EACH POINT INDIVIDUALLY - EACH TIME CHECKING FOR SIGNIFICANCE
  for(step in 1:maxsteps) {
    if(sig[step] == 0) {
      plot(step:step, dat[step:step], type="p", col="gray", lwd=8, xlim=c(1,maxsteps), ylim=c(minval,maxval), xlab=xtit, ylab=ytit)
    
    }
    else {
      plot(step:step, dat[step:step], type="p", col="black", lwd=8, xlim=c(1,maxsteps), ylim=c(minval,maxval), xlab=xtit, ylab=ytit)
    
    } 
    par(new=T)
  }

  
  
}

