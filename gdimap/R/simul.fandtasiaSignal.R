## Based on Matlab code by
## Angelos Barmpoutis, PhD
## fanDTasia ToolBox

simul.fandtasiaSignal <-
function(g, gridsz=32, b=4000, sigma=NULL, savedir=tempdir())
{
  order <- 4 # In this demo we compute 4th-order tensors
  sz <- dim(g)[1]
  GradientOrientations <- rbind(c(1, 0, 0), g[1:sz,])
  b_value <- c(10, rep(b, sz))
  ngrad <- dim(GradientOrientations)[1]
  S <- array(1, dim=c(gridsz,gridsz,1,ngrad))
  gridsz2 <- gridsz/2
	gc()
  cat("Simulating DW-MRI field using Barmpoutis' algorithm\n")
  cat("Please be patient ...\n")
  tline <- floor(c(0.2,0.4,0.6,0.8)*ngrad)
  tperc <- c(20, 40, 60, 80)
  for (i in 2:ngrad) {
    tt <- which(tline == i)
    if(length(tt) != 0) {
      cat(paste(tperc[tt],"% ", sep="")); cflush() }
    for (x in 1:gridsz) {
      for (y in 1:gridsz) {
        f1_flag <- 0
        f2_flag <- 0
        if (((x*x+y*y) > (gridsz2*gridsz2)) & ((x*x+y*y) < (gridsz*gridsz))) {
          v <- c(y/x, -1, 0)
          v <- v / sqrt(t(v)%*%v)
          fiber_orientation1 <- v
          f1_flag <- 1
        }
        if ((x < y+10) & (x > y-10)) {
          fiber_orientation2 <- c(sqrt(2)/2, sqrt(2)/2, 0)
          f2_flag <- 1
        }
        if (f1_flag==0 & f2_flag==1) {
          fiber_orientation1 <- fiber_orientation2
        }
        else {
          if (f1_flag==1 & f2_flag==0) {
             fiber_orientation2 <- fiber_orientation1
          }
          else {
             if (f1_flag==0 & f2_flag==0) {
              fiber_orientation1 <- c(0, 0, 1)
              fiber_orientation2 <- c(0, 0, 1)
            }
          }
        }
        S[x,y,1,i] <- S[x,y,1,1] *
          (SimulateDWMRI(fiber_orientation1,GradientOrientations[i,]) +
          SimulateDWMRI(fiber_orientation2,GradientOrientations[i,]))/2
        # Riccian type noise model (typical std ~ 0.01)
        if(!is.null(sigma)) {
           S[x,y,1,i] <- sqrt((S[x,y,1,i] +
            sigma*rnorm(1))^2+(sigma*rnorm(1))^2)
        }
      }
    }
  }
  cat("100% completed\n")
  f <- paste(savedir,"/simfield",sep="")
  options(niftiAuditTrail = FALSE)
  writeNIfTI(S, filename=f)
  cat("wrote",f,"\n")
}



SimulateDWMRI <-
function(fiber,gradient)
{
  fiber <- fiber/sqrt(fiber%*%fiber)
  gradient <- gradient/sqrt(gradient%*%gradient)
  cosine <- fiber%*%gradient
  cp <- c(0.0672, 0.1521, 0.3091, 0.4859, 0.6146)
  numofcp <- length(cp)
  out <- 0
  for (k in 1:numofcp) {
      out <- out+cp[k]*basisN(abs(cosine),k,numofcp)
  }
  out <- out + cp[numofcp] * basisN(abs(cosine), numofcp+1, numofcp)
}


basisN <-
function(argument, segment_id, total_segments)
{
  #Evaluates the 2nd-order spline basis
  #argument must be a real in [0...1]
  #segment_id must be an integer in [1...total_segments]
  argument <- 1-argument;
  out <- numeric(length(argument))
  for (i in 1:length(argument)) {
     if ((argument[i] >= (segment_id-2)/total_segments) && (argument[i] < (segment_id-1)/total_segments) ) {
          t <- (argument[i]-(segment_id-2)/total_segments)/((segment_id-1)/total_segments-(segment_id-2)/total_segments)
          out[i] <- 0.5*t*t
    }
    else {
      if( (argument[i] >= (segment_id-1)/total_segments) && (argument[i] < segment_id/total_segments) ) {
          t <- (argument[i]-(segment_id-1)/total_segments)/((segment_id)/total_segments-(segment_id-1)/total_segments)
          out[i] <- -t*t+t+0.5
      }
       else {
        if ( (argument[i] >= (segment_id)/total_segments) && (argument[i] < (segment_id+1)/total_segments) ) {
          t <- (argument[i]-(segment_id)/total_segments)/((segment_id+1)/total_segments-(segment_id)/total_segments)
          out[i] <- 0.5*(1-t)*(1-t)
        }
        else {
          out[i] - 0
        }
      }
    }
  }
  return(out)
}


