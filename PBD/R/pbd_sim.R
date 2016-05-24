pbd_sim = function(pars,age,soc = 2,plotit = FALSE)
{
la1 = pars[1]
la2 = pars[2]
la3 = pars[3]
mu1 = pars[4]
mu2 = pars[5]

i = 1
while(i <= soc)
{
   t = 0
   if(i == 1)
   {
      id1 = 0
      id = id1 + 1
      Sid1 = 0
      Sid = 1
      sg = id
      si = NULL
      L = t(as.matrix(c(id,0,-1e-10,t,-1,1)))
   }
   if(i == 2)
   {
      id = id1 + 1
      Sid = Sid1
      sg = NULL
      si = -id
      L = t(as.matrix(c(id,1,t,-1,-1,1)))
      #L = rbind(L1,t(as.matrix(c(id,1,t,-1,-1,1))))
   }

   Ng = length(sg)
   Ni = length(si)  
   probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)   
   denom = sum(probs) 
   probs = probs/denom
   t = t - log(runif(1))/denom

   while(t <= age)
   {
      event = sample2(1:5,size = 1,prob = probs)
      if(event == 1)
      {
          parent = as.numeric(sample2(sg,1))
          id = id + 1
          L = rbind(L,c(id,parent,t,-1,-1,L[abs(parent) - id1,6]))
          si = c(si,-id)
          Ni = Ni + 1
      }
      if(event == 2)
      {
          todie = as.numeric(sample2(sg,1))
          L[todie - id1,5] = t
          sg = sg[-which(sg == todie)]
          Ng = Ng - 1
      }
      if(event == 3)
      {
          tocomplete = abs(as.numeric(sample2(si,1)))
          L[tocomplete - id1,4] = t
          Sid = Sid + 1
          L[tocomplete - id1,6] = Sid      
          sg = c(sg,tocomplete)
          si = si[-which(abs(si) == tocomplete)]
          Ng = Ng + 1
          Ni = Ni - 1
      }
      if(event == 4)
      {
          parent = as.numeric(sample2(si,1))
          id = id + 1
          L = rbind(L,c(id,parent,t,-1,-1,L[abs(parent) - id1,6]))
          si = c(si,-id)
          Ni = Ni + 1
      }
      if(event == 5)
      {
          todie = abs(as.numeric(sample2(si,1)))
          L[todie - id1,5] = t
          si = si[-which(abs(si) == todie)]
          Ni = Ni - 1
      }
      probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)   
      denom = sum(probs) 
      probs = probs/denom
      t = t - log(runif(1))/denom
   }
   if(i == 1)
   {
       if((Ng + Ni) > 0)
       {
           i = i + 1
           L1 = L
           id1 = id
           Sid1 = Sid
           si1 = si
           sg1 = sg
       }       
   } else {
   if(i == 2)
   {
       if(checkgood(L,si,sg) == 1)
       {
           i = i + 1
           L2 = L   
           si2 = si
           sg2 = sg                   
       }       
   }} 
}
L = L1
if(soc == 2)
{
    L = rbind(L1,L2)
}
L0 = L
absL = L
absL[,2] = abs(L[,2])
trans = NULL
igtree.extinct = phytools::read.simmap(text = detphy(absL,age,ig = T,dropextinct = F))
igtree.extant = phytools::read.simmap(text = detphy(absL,age,ig = T,dropextinct = T))
tree = as.phylo(read.tree(text = detphy(absL,age)))

# Random sampling
sL_random = sampletree(absL, age, samplemethod = "random")
stree_random = as.phylo(read.tree(text = detphy(sL_random, age)))
# Sampling the oldest
sL_oldest = sampletree(absL, age, samplemethod = "oldest")
stree_oldest = as.phylo(read.tree(text = detphy(sL_oldest, age)))
# Sampling the youngest
sL_youngest = sampletree(absL, age, samplemethod = "youngest")
stree_youngest = as.phylo(read.tree(text = detphy(sL_youngest, age)))

sL_random[, 3:5][which(sL_random[, 3:5] == -1)] = age + 1
sL_random[, 3:5] = age - sL_random[, 3:5]
sL_random = sL_random[order(sL_random[, 1]), ]
sL_oldest[, 3:5][which(sL_oldest[, 3:5] == -1)] = age + 1
sL_oldest[, 3:5] = age - sL_oldest[, 3:5]
sL_oldest = sL_oldest[order(sL_oldest[, 1]), ] 
sL_youngest[, 3:5][which(sL_youngest[, 3:5] == -1)] = age + 1
sL_youngest[, 3:5] = age - sL_youngest[, 3:5]
sL_youngest = sL_youngest[order(sL_youngest[, 1]), ]

#sL = sampletree(absL,age)
#stree = as.phylo(read.tree(text = detphy(sL,age)))
#sL[,3:5][which(sL[,3:5] == -1)] = age + 1
#sL[,3:5] = age - sL[,3:5]
#sL = sL[order(sL[,1]),]

#rbrts = sort(age - pbd_sim_step2a(L)[[1]])
reconL = pbd_reconstruct(L0)
recontree = as.phylo(read.tree(text = detphy(reconL,age)))
L[,3:5][which(L[,3:5] == -1)] = age + 1
L[,3:5] = age - L[,3:5]
L = L[order(L[,1]),]
#reducL = cbind(L[,3],abs(L[,2:1]),L[,5],L[,4],L[,6])
#fulltree = L2phylo2(reducL,dropextinct = F)
reconL[,3:5][which(reconL[,3:5] == -1)] = age + 1
reconL[,3:5] = age - reconL[,3:5]
reconL = reconL[order(reconL[,1]),]

if(plotit == TRUE)
{
   par(mfrow = c(3,3))
   cols = setNames(c("gray","black"),c("i","g"))
   phytools::plotSimmap(igtree.extinct,colors = cols)
   phytools::plotSimmap(igtree.extant,colors = cols)
   plot(tree, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   plot(stree_random, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)   
   plot(stree_oldest, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)   
   plot(stree_youngest, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)   
   plot(recontree, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
   par(mfrow = c(1,1))
}


Ltreeslist = list(tree = tree,stree_random = stree_random,stree_oldest = stree_oldest,stree_youngest = stree_youngest,L = L,sL_random = sL_random,sL_oldest = sL_oldest,sL_youngest = sL_youngest,igtree.extinct = igtree.extinct,igtree.extant = igtree.extant,recontree = recontree,reconL = reconL,L0 = L0)
return(Ltreeslist)
}