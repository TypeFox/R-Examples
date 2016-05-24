compute_edgeTList = function(tree)
{	
	nNode = Nnode(tree)
	edgeTList = vector("list", nNode)
	ca = max(branching.times(tree))
									
	for(i in 1:nNode)
  {	
		ntips = Ntip(tree)
		nodeDepths = dist.nodes(tree)[(ntips + 1):(2 * ntips - 1),ntips + 1]        
		nodeDescendants = clade.members.list(tree,tip.labels = TRUE)
			
		sisterNodes = which(sapply(nodeDescendants,length) == 2)
		sisterDepths = nodeDepths[sisterNodes]
	  w = which(sisterDepths==max(sisterDepths))
			
		sisterNodes = which(sapply(nodeDescendants,length) == 2) # find nodes for sisters
		sisterDepths = nodeDepths[sisterNodes] # find age of sister nodes
	  w = which(sisterDepths == max(sisterDepths)) # find youngest sister node
	  toDrop1 = nodeDescendants[[sisterNodes[w]]]
	  toDrop2 = c(which(tree$tip.label == toDrop1[1]),which(tree$tip.label == toDrop1[2]))
	  edgeT = ca - c(dist.nodes(tree)[toDrop2,ntips + 1],sisterDepths[w]) # specify start and finish of edges to be dropped  
	  names(edgeT)[1:2] = toDrop1
	  edgeTList[[i]] = edgeT
	    
	  tree = drop.tip(tree, toDrop1,trim.internal = FALSE)
		tree$tip.label[which(tree$tip.label == "NA")] = paste("p",i - 1,sep = "")
	}  
  return(edgeTList) 	
}			

DAMOCLES_loglik_rhs = function(
   t,
   p,
   pars
   )
{
mu = pars[1]
ga = pars[2]/(1 + pars[3] * (pars[4] - t))

dp = c(ga * (p[2] - p[1]),mu * (p[1] - p[2]))

return(list(dp))
}

DAMOCLES_integrate_ODE = function(
   pars,
   p,
   tt,
   ca
   )
{
   if(pars[3] == 0)
   {
      # TAKE ANALYTICAL SOLUTION
      mu = pars[1]
      ga = pars[2]
      difft = tt[2] - tt[1]
      p0f = p[1] * mu + ga * p[2] + ga * (p[1] - p[2]) * exp(-difft * (ga + mu))
      p1f = p[1] * mu + ga * p[2] - mu * (p[1] - p[2]) * exp(-difft * (ga + mu))
      p = 1/(ga + mu) * c(p0f,p1f)
   } else {
      # SOLVE ODE NUMERICALLY
      y = lsoda(p,tt,DAMOCLES_loglik_rhs,c(pars,ca),rtol = 1E-10,atol = 1E-16)
      p = y[2,2:3]
   }
   return(p)
}

DAMOCLES_loglik = function(
   phy,
   pa,
   pars,
   pchoice = 0,
   edgeTList = compute_edgeTList(phy)
   )
{
#DAMOCLES_loglik computes the likelihood of the DAMOCLES model given a tree and immigration and extinction parameters.
#phy = phylogenetic tree, in nexus-format
#pa = presence-absence table, first column contains tip labels, second column presence (1) or absence (0)
#pars = set of parameters
#- pars[1] = mu = extinction rate
#- pars[2] = gamma_0 = immigration rate parameter in formula gamma(t) = gamma_0/(1+gamma_1 * t)
#- pars[3] = gamma_1 = immigration rate parameter in formula gamma(t) = gamma_0/(1+gamma_1 * t)
#pchoice = choice which probability to optimize:
# pchoice == 0, optimize p_0f + p_1f
# pchoice == 1, optimize p_0f
# pchoice == 2, optimize p_1f
#edgeTList = list of succesive edges to prune; if NULL it is computed using compute_edgeTList(phy)

  patable = cbind(pa[,1],1 - as.numeric(pa[,2]),pa[,2])
	nNode = dim(patable)[1] - 1

	for(i in 1:nNode)
	{
	  ca = max(unlist(edgeTList))
  	branchesTimes = edgeTList[[i]]
    pc = matrix(0,nrow = 2,ncol = 2)
    w2 = c(0,0)
    for(j in 1:2)
    {
       tt = c(branchesTimes[j],branchesTimes[3])
       w2[j] = which(patable[,1] == names(branchesTimes)[j]) 
       p = as.numeric(patable[w2[j],2:3])
       p = DAMOCLES_integrate_ODE(pars,p,tt,ca)
       patable[w2[j],2:3] = p
       pc[j,] = p
    }
    pnew = c(pc[1,1] * pc[2,1], (pc[1,2] * pc[2,1] + pc[2,2] * pc[1,1])/2) 
	  updateTips = as.character(patable[,1])
    updateTips[w2[1]] = as.character(paste("p",i - 1,sep = ""))
    patable[,1] = updateTips
    patable[w2[1],2:3] = pnew
    patable = patable[-w2[2],]
	}	
	if(!is.null(phy$root.edge))
	{
     # SPECIFY START AND FINISH OF BRANCH j
     tt = c(0,phy$root.edge)
     #tt = c(branchesTimes[3],0)
     # READ INITIAL VALUES FOR THIS BRANCH
     p = as.numeric(patable[2:3])
     p = DAMOCLES_integrate_ODE(pars,p,tt,ca)
     patable[2:3] = p
	}
	loglik = log(sum(as.numeric(patable[(2 + (pchoice == 2)):(3 - (pchoice == 1))])))
	#s1 = sprintf('Parameters: %f %f %f',pars[1],pars[2],pars[3])
	#s2 = sprintf(', Loglikelihood: %f',loglik)
	#cat(s1,s2,"\n",sep = "")
	#flush.console()
	return(loglik)
}

