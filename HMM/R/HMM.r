initHMM = function(States, Symbols, startProbs=NULL, transProbs=NULL,
    emissionProbs=NULL)
{
  nStates    = length(States)
  nSymbols   = length(Symbols)
  S          = rep(1/nStates,nStates)
  T          = 0.5*diag(nStates) + array(0.5/(nStates),c(nStates,nStates))
  E          = array(1/(nSymbols),c(nStates,nSymbols))
  names(S)   = States
  dimnames(T)= list(from=States,to=States)
  dimnames(E)= list(states=States,symbols=Symbols)
  if(!is.null(startProbs)){S[]  = startProbs[]}
  if(!is.null(transProbs)){T[,] = transProbs[,]}
  if(!is.null(emissionProbs)){E[,] = emissionProbs[,]}
  return(list(States=States,Symbols=Symbols,startProbs=S,transProbs=T,
      emissionProbs=E))
}

simHMM = function(hmm, length)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  states   = c()
  emission = c()
  states   = c(states, sample(hmm$States,1,prob=hmm$startProbs))
  for(i in 2:length)
  {
    state  = sample(hmm$States, 1, prob=hmm$transProbs[states[i-1],])
  	states = c(states, state)
  }
  for(i in 1:length)
  {
    emi      = sample(hmm$Symbols, 1, prob=hmm$emissionProbs[states[i],])
  	emission = c(emission, emi)
  }
  return(list(states=states,observation=emission))
}

viterbi = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  v          = array(NA,c(nStates,nObservations))
  dimnames(v)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    v[state,1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,observation[1]])
  }
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      maxi = NULL
      for(previousState in hmm$States)
      {
        temp = v[previousState,k-1] + log(hmm$transProbs[previousState,state]) 
        maxi = max(maxi, temp)
      }
      v[state,k] = log(hmm$emissionProbs[state,observation[k]]) + maxi
    }
  }
  # Traceback
  viterbiPath = rep(NA,nObservations)
  for(state in hmm$States)
  {
    if(max(v[,nObservations])==v[state,nObservations])
    {
      viterbiPath[nObservations] = state
      break
    }
  }
  for(k in (nObservations-1):1)
  {
    for(state in hmm$States)
    {
      if(max(v[,k]+log(hmm$transProbs[,viterbiPath[k+1]]))
          ==v[state,k]+log(hmm$transProbs[state,viterbiPath[k+1]]))
      {
        viterbiPath[k] = state
        break
      }
    }
  }
  return(viterbiPath)
}

forward = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  f          = array(NA,c(nStates,nObservations))
  dimnames(f)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    f[state,1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,observation[1]])
  }
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      logsum = -Inf
      for(previousState in hmm$States)
      {
        temp   = f[previousState,k-1] + log(hmm$transProbs[previousState,state])
		if(temp > - Inf)
		{
			logsum = temp + log(1 + exp(logsum - temp ))
		}
      }
      f[state,k] = log(hmm$emissionProbs[state,observation[k]]) + logsum
    }
  }
  return(f)
}

backward = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  b          = array(NA,c(nStates,nObservations))
  dimnames(b)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    b[state,nObservations] = log(1)
  }
  # Iteration
  for(k in (nObservations-1):1)
  {
    for(state in hmm$States)
    {
      logsum = -Inf
      for(nextState in hmm$States)
      {
        temp   = b[nextState,k+1] + log(hmm$transProbs[state,nextState]*hmm$emissionProbs[nextState,observation[k+1]])
		if(temp > - Inf)
		{
        	logsum = temp + log(1 + exp(logsum-temp))
		}
      }
      b[state,k] = logsum
    }
  }
  return(b)
}

dishonestCasino = function()
{
  # setup HMM
  nSim          = 2000
  States        = c("Fair","Unfair")
  Symbols       = 1:6 #c("1er","2er","3er","4er","5er","6er")
  transProbs    = matrix(c(.99,.01,.02,.98), c(length(States),length(States)), byrow=TRUE)
  emissionProbs = matrix(c(rep(1/6,6),c(rep(.1,5),.5)), c(length(States),length(Symbols)), byrow=TRUE)
  hmm = initHMM(States, Symbols, transProbs=transProbs, emissionProbs=emissionProbs)
  sim = simHMM(hmm,nSim)
  vit = viterbi(hmm, sim$observation)
  f   = forward(hmm, sim$observation)
  b   = backward(hmm, sim$observation)
  # todo: probObservations is not generic!
  f[1,nSim]->i
  f[2,nSim]->j
  probObservations = (i + log(1+exp(j-i)))
  posterior = exp((f+b)-probObservations)
  x = list(hmm=hmm,sim=sim,vit=vit,posterior=posterior)
  
  readline("Plot simulated throws:\n")
  mn = "Fair and unfair die"
  xlb = "Throw nr."
  ylb = ""
  plot(x$sim$observation,ylim=c(-7.5,6),pch=3,main=mn,xlab=xlb,ylab=ylb,bty="n",yaxt="n")
  axis(2,at=1:6)

  readline("Simulated, which die was used:\n")
  text(0,-1.2,adj=0,cex=.8,col="black","True: green = fair die")
  for(i in 1:nSim)
  {
    if(x$sim$states[i] == "Fair")
      rect(i,-1,i+1,0, col = "green", border = NA)
    else
      rect(i,-1,i+1,0, col = "red", border = NA)   
  }

  readline("Most probable path (viterbi):\n")
  text(0,-3.2,adj=0,cex=.8,col="black","Most probable path")
  for(i in 1:nSim)
  { 
    if(x$vit[i] == "Fair")
      rect(i,-3,i+1,-2, col = "green", border = NA)
    else
      rect(i,-3,i+1,-2, col = "red", border = NA)  
  }

  readline("Differences:\n")
  text(0,-5.2,adj=0,cex=.8,col="black","Difference")
  differing = !(x$sim$states == x$vit)
  for(i in 1:nSim)
  {
    if(differing[i])
      rect(i,-5,i+1,-4, col = rgb(.3, .3, .3), border = NA)
    else
      rect(i,-5,i+1,-4, col = rgb(.9, .9, .9), border = NA)  
  }

  readline("Posterior-probability:\n")
  points(x$posterior[2,]-3, type="l")
  #points(x$posterior[2,]-5, type="l")

readline("Difference with classification by posterior-probability:\n")
  text(0,-7.2,adj=0,cex=.8,col="black","Difference by posterior-probability")
  differing = !(x$sim$states == x$vit)
  for(i in 1:nSim)
  {
   if(posterior[1,i]>0.5)
   {
    if(x$sim$states[i] == "Fair")
      rect(i,-7,i+1,-6, col = rgb(.9, .9, .9), border = NA)
    else
      rect(i,-7,i+1,-6, col = rgb(.3, .3, .3), border = NA)
   }else
   {
    if(x$sim$states[i] == "Unfair")
      rect(i,-7,i+1,-6, col = rgb(.9, .9, .9), border = NA)
    else
      rect(i,-7,i+1,-6, col = rgb(.3, .3, .3), border = NA)
   }
  }

  readline("Difference with classification by posterior-probability > .95:\n")
  text(0,-7.2,adj=0,cex=.8,col="black","Difference by posterior-probability > .95")
  differing = !(x$sim$states == x$vit)
  for(i in 1:nSim)
  {
   if(posterior[2,i]>0.95 || posterior[2,i]<0.05)
   {
    if(differing[i])
      rect(i,-7,i+1,-6, col = rgb(.3, .3, .3), border = NA)
    else
      rect(i,-7,i+1,-6, col = rgb(.9, .9, .9), border = NA)
   }
   else
   {
     rect(i,-7,i+1,-6, col = rgb(.9, .9, .9), border = NA)
   }
  }
  invisible(x)
}

posterior = function(hmm, observation)
{
	hmm$transProbs[is.na(hmm$transProbs)]       = 0
	hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
	f = forward(hmm, observation)
	b = backward(hmm, observation)
	probObservations = f[1,length(observation)]
	for(i in 2:length(hmm$States))
	{
		j = f[i,length(observation)]
		if(j > - Inf)
		{
			probObservations = j + log(1+exp(probObservations-j))
		}
	}
	posteriorProb = exp((f+b)-probObservations)
	return(posteriorProb)
}

baumWelch = function(hmm, observation, maxIterations=100, delta=1E-9, pseudoCount=0)
{
	tempHmm = hmm
	tempHmm$transProbs[is.na(hmm$transProbs)]       = 0
	tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
	diff = c()
	for(i in 1:maxIterations)
	{
		# Expectation Step (Calculate expected Transitions and Emissions)
		bw = baumWelchRecursion(tempHmm, observation)
		T  = bw$TransitionMatrix
		E  = bw$EmissionMatrix
		# Pseudocounts
		T[!is.na(hmm$transProbs)]    = T[!is.na(hmm$transProbs)]    + pseudoCount
		E[!is.na(hmm$emissionProbs)] = E[!is.na(hmm$emissionProbs)] + pseudoCount
		# Maximization Step (Maximise Log-Likelihood for Transitions and Emissions-Probabilities)
		T = (T/apply(T,1,sum))
		E = (E/apply(E,1,sum))
		d = sqrt(sum((tempHmm$transProbs-T)^2)) + sqrt(sum((tempHmm$emissionProbs-E)^2))
		diff = c(diff, d)
		tempHmm$transProbs    = T
		tempHmm$emissionProbs = E
		if(d < delta)
		{
			break
		}
	}
	tempHmm$transProbs[is.na(hmm$transProbs)]       = NA
	tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = NA
	return(list(hmm=tempHmm,difference=diff))
}

baumWelchRecursion = function(hmm, observation)
{
	TransitionMatrix    = hmm$transProbs
	TransitionMatrix[,] = 0
	EmissionMatrix      = hmm$emissionProbs
	EmissionMatrix[,]   = 0
	f = forward(hmm,  observation)
	b = backward(hmm, observation)
	probObservations = f[1,length(observation)]
	for(i in 2:length(hmm$States))
	{
		j = f[i,length(observation)]
		if(j > - Inf)
		{
			probObservations = j + log(1+exp(probObservations-j))
		}
	}
	for(x in hmm$States)
	{
		for(y in hmm$States)
		{
			temp = f[x,1] + log(hmm$transProbs[x,y]) +
					log(hmm$emissionProbs[y,observation[1+1]]) + b[y,1+1]
			for(i in 2:(length(observation)-1))
			{
				j = f[x,i] + log(hmm$transProbs[x,y]) +
						log(hmm$emissionProbs[y,observation[i+1]]) + b[y,i+1]
				if(j > - Inf)
				{
					temp = j + log(1+exp(temp-j))
				}
			}
			temp = exp(temp - probObservations)
			TransitionMatrix[x,y] = temp
		}
	}
	for(x in hmm$States)
	{
		for(s in hmm$Symbols)
		{
			temp = -Inf
			for(i in 1:length(observation))
			{
				if(s == observation[i])
				{
					j = f[x,i] + b[x,i]
					if(j > - Inf)
					{
						temp = j + log(1+exp(temp-j))
					}
				}
			}
			temp = exp(temp - probObservations)
			EmissionMatrix[x,s] = temp
		}
	}
	return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=EmissionMatrix))
}

viterbiTraining = function(hmm, observation, maxIterations=100, delta=1E-9, pseudoCount=0)
{
	tempHmm = hmm
	tempHmm$transProbs[is.na(hmm$transProbs)]       = 0
	tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
	diff = c()
	for(i in 1:maxIterations)
	{
		# Counts
		vt = viterbiTrainingRecursion(tempHmm, observation)
		T  = vt$TransitionMatrix
		E  = vt$EmissionMatrix
		# Pseudocounts
		T[!is.na(hmm$transProbs)]    = T[!is.na(hmm$transProbs)]    + pseudoCount
		E[!is.na(hmm$emissionProbs)] = E[!is.na(hmm$emissionProbs)] + pseudoCount
		# Relativ Frequencies
		T = (T/apply(T,1,sum))
		E = (E/apply(E,1,sum))
		d = sqrt(sum((tempHmm$transProbs-T)^2)) + sqrt(sum((tempHmm$emissionProbs-E)^2))
		diff = c(diff, d)
		tempHmm$transProbs    = T
		tempHmm$emissionProbs = E
		if(d < delta)
		{
			break
		}
	}
	tempHmm$transProbs[is.na(hmm$transProbs)]       = NA
	tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = NA
	return(list(hmm=tempHmm,difference=diff))
}

viterbiTrainingRecursion = function(hmm, observation)
{
	TransitionMatrix    = hmm$transProbs
	TransitionMatrix[,] = 0
	EmissionMatrix      = hmm$emissionProbs
	EmissionMatrix[,]   = 0
	v = viterbi(hmm,  observation)
	for(i in 1:(length(observation)-1))
	{
		TransitionMatrix[v[i],v[i+1]] = TransitionMatrix[v[i],v[i+1]] + 1
	}
	for(i in 1:length(observation))
	{
		EmissionMatrix[v[i],observation[i]] = EmissionMatrix[v[i],observation[i]] + 1
	}
	return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=EmissionMatrix))
}
