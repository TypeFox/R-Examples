PlotMixturesAndBoundaries <-function(Data, Means, SDs, Weights, IsLogDistribution = rep(FALSE,length(Means)), SingleColor = 'blue', MixtureColor = 'red',DataColor='black', BoundaryColor = 'magenta', xlab, ylab, ...){
#PlotGaussMixesAndBoundaries(Data,Means,SDs,Weights,SingleColor,MixtureColor)
# Plot a Mixture of Gaussian/LogNormal and Bayesian decision boundaries

# INPUT
# Data[1:n]                   data column for which the distribution was modelled
# Means[1:L]                  Means of Gaussians,  L ==  Number of Gaussians
# SDs[1:L]                  estimated Gaussian Kernels = standard deviations
# Weights[1:L]                relative number of points in Gaussians (prior #															probabilities): sum(Weights) ==1
# OPTIONAL
# IsLogDistribution(1:L) 			gibt an ob die jeweilige Verteilung eine 
#															Lognormaverteilung ist,(default ==0*(1:L))
# SingleColor              		PlotSymbol of all the single gaussians, default 
#															magenta
# MixtureColor             		PlotSymbol of the mixture, default black
# BoundaryColor               PlotSymbol of the boundaries default green
# RoundNpower                 Decision Boundaries are rounded by 
#															roundn(DecisionBoundariesRoundNpower)
# ...													other plot arguments, like xlim = c(1,10)

# OUTPUT
# DecisonBoundaries(1:L-1)    Bayes decision boundaries 
# DBY(1:L-1)                  y values at the cross points of the Gaussians
       
# author MT 08/2015

# Labels
if(missing(xlab)){
	xlab = '' # no label for x axis
}
if(missing(ylab)){
	ylab = '' # no label for y axis
}

# Calculate intersections
dec = BayesDecisionBoundaries(Means,SDs,Weights,IsLogDistribution,Ycoor=T)

#MT: Bugifx
#if(is.list(dec)){
  DecisionBoundaries=as.vector(dec$DecisionBoundaries)
  #print('dec was a list, assuming usage of BayesDecisionBoundaries()')
#}
# Plot Gaussians
PlotMixtures(Data,Means,SDs,Weights,IsLogDistribution,SingleColor,MixtureColor,DataColor=DataColor, xlab = xlab, ylab  = ylab,...)

# intersecions
for (i in 1:length(DecisionBoundaries)){
 abline(v = DecisionBoundaries[i], col = BoundaryColor)
}

return (invisible(list(DecisionBoundaries = DecisionBoundaries, DBY = dec$DBY))) 
}

