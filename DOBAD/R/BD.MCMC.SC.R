
### "..." is for good form
### simMethod: 0 means Accept-Reject; 1 is our sim method.
BD.MCMC.SC <- function(Lguess, Mguess,
                       beta.immig, #not a guess, of course
                       alpha.L, beta.L, alpha.M, beta.M,
                       data,
                       burnIn=100, N=1000, n.fft=1024,
                       verbose=1, verbFile=NULL,
                       simMethod=-1,
                       ...){
  UseMethod("BD.MCMC.SC", data)
}

