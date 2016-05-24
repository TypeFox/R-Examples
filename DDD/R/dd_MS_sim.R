dd_MS_lamuN = function(ddmodel,pars,N)
{
    laM = pars[1]
    muM = pars[2]
    Kprime = pars[3]
    laS = pars[4]
    muS = pars[5]
    n0 = (ddmodel == 2 | ddmodel == 4)
    if(ddmodel == 1 | ddmodel == 1.3)
    {
        # linear dependence in speciation rate
        laMN = max(0,laM * (1 - N/Kprime))
        muMN = muM
        laSN = max(0,laS * (1 - N/Kprime))
        muSN = muS
    }
    if(ddmodel == 2 | ddmodel == 2.1 | ddmodel == 2.2 | ddmodel == 2.3)
    {
        # exponential dependence in speciation rate
        al = Kprime^(ddmodel != 2.2)
        laMN = laM * (N + n0)^(-al)
        muMN = muM
        laSN = laS * (N + n0)^(-al)
        muSN = muM
    }
    return(c(laMN,muMN,laSN,muSN))
}

dd_MS_sim = function(pars,age,ddmodel = 1.3)
{
# Simulation of diversity-dependent process
#  . start from crown age
#  . no additional species at crown node
#  . no missing species in present
# pars = [laM muM K laS muS tinn]
# - pars1[1] = laM = (initial) speciation rate of main clade
# - pars1[2] = muM = extinction rate of main clade
# - pars1[3] = K' = maximum number of species
# - pars1[4] = laS = (initial) speciation rate of subclade
# - pars1[5] = muS = extinction rate of subclade
# - pars1[6] = tinn = time of key innovation
# age = crown age
# ddmodel = mode of diversity-dependence
#  . ddmodel == 1 : linear dependence in speciation rate with parameter K
#  . ddmodel == 1.3: linear dependence in speciation rate with parameter K'
#  . ddmodel == 2 : exponential dependence in speciation rate
#  . ddmodel == 2.1: variant with offset at infinity
#  . ddmodel == 2.2: 1/n dependence in speciation rate
#  . ddmodel == 2.3: exponential dependence in speciation rate with parameter x

done = 0
if(pars[6] > age)
{
   stop('The key innovation time is before the crown age of the main clade.')
}
if((pars[1] < pars[2]) | (pars[4] < pars[5]))
{
   stop('lambda0 is smaller than mu for one or both clades')
}
if(min(pars) < 0)
{
   stop('One of the parameters is negative')
}
if(!(ddmodel %in% c(1,1.3,2,2.1,2.2,2.3)))
{
   stop('This diversity-dependence model does not exist or is not implemented')
}
while(done == 0)
{
    # number of species N at time t
    # i = index running through t and N
    t = rep(0,1)
    L = matrix(0,2,5)
    i = 1
    t[1] = 0
    NM = 2
    NS = 0
    # L = data structure for lineages,
    # . L[,1] = branching times
    # . L[,2] = index of parent species
    # . L[,3] = index of daughter species
    # . L[,4] = time of extinction
    # . L[,5] = main clade (0) or subclade (1)
    # j = index running through L
    L[1,1:5] = c(0,0,-1,-1,0)
    L[2,1:5] = c(0,-1,2,-1,0)
    linlistM = c(-1,2)
    linlistS = NULL
    newL = 2
    tinn = age - pars[6]
    ff = dd_MS_lamuN(ddmodel,pars,NM[i] + NS[i])
    laMN = ff[1]
    muMN = ff[2]
    laSN = ff[3]
    muSN = ff[4]
    denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
    t[i + 1] = t[i] - log(runif(1)) / denom
    if(t[i + 1] > tinn & t[i] < tinn)
    {
         NM[i] = NM[i] - 1
         NS[i] = NS[i] + 1
         linlistS = sample2(linlistM,1)
         L[abs(linlistS),5] = 1
         linlistM = linlistM[-which(linlistM == linlistS)]
         ff = dd_MS_lamuN(ddmodel,pars,NM[i] + NS[i])
         laMN = ff[1]
         muMN = ff[2]
         laSN = ff[3]
         muSN = ff[4]
         denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
         t[i + 1] = tinn - log(runif(1)) / denom
    }
    while(t[i + 1] <= age)
    {
        event = sample2(x = 1:4,size = 1,prob = c(laMN * NM[i], muMN * NM[i], laSN * NS[i], muSN * NS[i]))
        i = i + 1
        if(event == 1)
        {
            # speciation event in main clade
            ranL = sample2(linlistM,1)
            NM[i] = NM[i - 1] + 1
            NS[i] = NS[i - 1]
            newL = newL + 1
            L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1,0))
            linlistM = c(linlistM,sign(ranL) * newL)
        } else if(event == 3)
        {
            # speciation event in subclade
            ranL = sample2(linlistS,1)
            NM[i] = NM[i - 1]
            NS[i] = NS[i - 1] + 1
            newL = newL + 1
            L = rbind(L,c(t[i],ranL,sign(ranL) * newL,-1,1))
            linlistS = c(linlistS,sign(ranL) * newL)           
        } else if(event == 2)
        {
            # extinction event in main clade
            ranL = sample2(linlistM,1)
            NM[i] = NM[i - 1] - 1
            NS[i] = NS[i - 1]
            L[abs(ranL),4] = t[i]
            w = which(linlistM == ranL)
            linlistM = linlistM[-w]
            linlistM = sort(linlistM)
        } else if(event == 4)
        {
            # extinction event in subclade
            ranL = sample2(linlistS,1)
            NM[i] = NM[i - 1]
            NS[i] = NS[i - 1] - 1
            L[abs(ranL),4] = t[i]
            w = which(linlistS == ranL)
            linlistS = linlistS[-w]
            linlistS = sort(linlistS)        
        }
        if(sum(c(linlistM,linlistS) < 0) == 0 | sum(c(linlistM,linlistS) > 0) == 0)
        {
            t[i + 1] = Inf
        } else {
            ff = dd_MS_lamuN(ddmodel,pars,NM[i] + NS[i])
            laMN = ff[1]
            muMN = ff[2]
            laSN = ff[3]
            muSN = ff[4]
            denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
            t[i + 1] = t[i] - log(runif(1)) / denom
            if(t[i + 1] > tinn & t[i] < tinn)
            {
               NM[i] = NM[i] - 1
               NS[i] = NS[i] + 1
               linlistS = sample2(linlistM,1)
               L[abs(linlistS),5] = 1
               linlistM = linlistM[-which(linlistM == linlistS)]
               ff = dd_MS_lamuN(ddmodel,pars,NM[i] + NS[i])
               laMN = ff[1]
               muMN = ff[2]
               laSN = ff[3]
               muSN = ff[4]
               denom = (laMN + muMN) * NM[i] + (laSN + muSN) * NS[i]
               t[i + 1] = tinn - log(runif(1)) / denom
            }
        }
    }
    if(sum(c(linlistM,linlistS) < 0) == 0 | sum(c(linlistM,linlistS) > 0) == 0)
    {
       done = 0
    } else {
       done = 1
    }
}

L[,1] = age - c(L[,1])
notmin1 = which(L[,4] != -1)
L[notmin1,4] = age - c(L[notmin1,4])
L[which(L[,4] == age + 1),4] = -1
tes = L2phylo(L[,1:4],dropextinct = T)
tas = L2phylo(L[,1:4],dropextinct = F)
tesS = NULL
tes2 = NULL
par(mfrow = c(2,1))
plot(tes)
plot(tas)
cols = c("blue","red")
names(cols) = c(0,1)
if(length(linlistS) > 0)
{
   namesS = paste('t',abs(linlistS), sep = "")
   if(length(linlistS) == 1)
   {
      m = which(tes$tip.label == namesS)
      b2 = 0
   }
   else if(length(linlistS) > 1)
   {
      m = getMRCA(phy = tes,tip = namesS)
      tesS = extract.clade(phy = tes,node = m)
      b2 = age - node.depth.edgelength(tes)[m]
   }  
   m0 = tes$edge[which(tes$edge[,2] == m),1]
   b1 = age - node.depth.edgelength(tes)[m0]
   tes2 = paintSubTree(tes,node = m,state = "1",anc.state = "0",stem = (pars[6] - b2)/(b1 - b2))
   plotSimmap(tes2,cols,lwd = 3,pts = F)
}
tasS = NULL
tas2 = NULL
allS = which(L[,5] == 1)
if(length(allS) > 0)
{
   namesS = paste('t',abs(allS), sep = "")
   if(length(allS) == 1)
   {
      m = which(tas$tip.label == namesS)
      b2 = 0
   }
   else if(length(allS) > 1)
   {
      m = getMRCA(phy = tas,tip = namesS)
      tasS = extract.clade(phy = tas,node = m)
      b2 = age - node.depth.edgelength(tas)[m]
   }
   m0 = tas$edge[which(tas$edge[,2] == m),1]
   b1 = age - node.depth.edgelength(tas)[m0]
   tas2 = paintSubTree(tas,node = m,state = "1",anc.state = "0", stem = (pars[6] - b2)/(b1 - b2))
   plotSimmap(tas2,cols,lwd = 3,pts = F)   
}
out = list(tes = tes,tas = tas,L = L,tesS = tesS,tasS = tasS,tes2 = tes2,tas2 = tas2)
return(out)

}