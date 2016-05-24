CTMC.simulate.piecewise <-
function(rate.fns, jumpLim.fns, T.times, init.state){
  M <- length(rate.fns);
  init.state.i <- init.state;
  simPieces <- vector("list", length=M);
  for (i in 1:M){
    simPieces[i] <- CTMC.simulate(rate.fn=rate.fns[[i]], jumpLim.fn=jumpLim.fns[[i]],
                              T.time=T.times[i+1]-T.times[i], init.state=init.state.i);
    init.state.i <- getStates(simPieces[[i]])[length(getStates(simPieces[[i]]))]
  }
  combineCTMC(simPieces);
}

