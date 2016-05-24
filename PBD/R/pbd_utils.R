roundn = function(x, digits = 0)
{
    fac = 10^digits
    return(trunc(fac * x + 0.5)/fac)
}

sample2 = function(x,size,replace = FALSE,prob = NULL)
{
    if(length(x) == 1)
    { 
        x = c(x,x)
        prob = c(prob,prob)
    }
    return(sample(x,size,replace,prob))
}

checkgood = function(L,si,sg,id1)
{
    j = 1;
    found = 0
    if(length(sg) > 0)
    {  
       found = 1
    } else {
    if(length(si) == 0) { found = 0 } else {
    while(found == 0 & j <= length(si))
    {           
        rowinc = which(L[,1] == abs(si[j]))
        parent = L[rowinc,2]
        birth = L[rowinc,3]
        #parent = L[si[j] - id1,2]
        #birth = L[si[j] - id1,3]
        while(found == 0 & parent > 1)
        {
            rowpar = which(L[,1] == parent)   
            if(L[rowpar,4] > -1 & L[rowpar,4] < birth)
            #if(L[parent - id1,4] > -1 & L[parent - id1,4] < birth)
            {
                found = 1
            } else
            {
                parent = L[rowpar,2]
                birth = L[rowpar,3]
                #parent = L[parent - id1,2]
                #birth = L[parent - id1,3]
            }
        }
        j = j + 1
    }}}
    invisible(found)
}

detphy = function(L, age, ig = F, dropextinct = T)
{
   #print(L)
   dimL = dim(L)
   if((dimL[1] == 1))
   {
      linlist = paste("(S-1-1-1:",age,");",sep = "")
   } else {
   L = L[order(L[,1]),1:6]
   if(dropextinct == T)
   {
      sall = which(L[,5] == -1)
      tend = age
   } else {
      sall = which(L[,5] >= -1)
      tend = (L[,5] == -1) * age + (L[,5] > -1) * L[,5]
   }   
   
   linlist = matrix(0,nrow = 1,ncol = 8)
   if(length(sall) == 1)
   {
      linlist[1,] = c(L[sall,],paste("S",paste(L[sall,6],L[sall,6],L[sall,1],sep = "-"),sep = ""),tend)
   } else {
      linlist = cbind(L[sall,],paste("S",paste(L[sall,6],L[sall,6],L[sall,1],sep = "-"),sep = ""),tend)
   }
   done = 0
   while(done == 0)
   {
      j = which.max(linlist[,3])      
      daughter = as.numeric(linlist[j,1])
      parent = as.numeric(linlist[j,2])
      parentj = which(linlist[,1] == parent)
      parentinlist = length(parentj)

              #print(linlist)
              #print(paste('j = ',j))
              #print(paste('duaghter = ',daughter))
              #print(paste('parent is ',parent))
              #print(paste('parentj is ',parentj))


      if(parentinlist == 1)
      {
          startedge = as.numeric(linlist[j,3])
          comptime = as.numeric(linlist[j,4])
          endedge = as.numeric(linlist[j,8])
          comptimeparent = as.numeric(linlist[parentj,4])          
          endedgeparent = as.numeric(linlist[parentj,8])
          if(ig == FALSE)
          {
              spec1 = paste(linlist[parentj,7],":",endedgeparent - startedge,sep = "")
              spec2 = paste(linlist[j,7],":",endedge - startedge,sep = "")
          } else {
              if(comptimeparent == -1 | comptimeparent > endedgeparent)
              {
                  comptimeparent = endedgeparent
              }
              itimeparent = max(0,comptimeparent - startedge)
              gtimeparent = endedgeparent - max(comptimeparent,startedge)
              if(itimeparent == 0)
              {
                  spec1 = paste(linlist[parentj,7],":{g,",gtimeparent,":i,",itimeparent,":g,0}",sep = "")
              } else {
                  spec1 = paste(linlist[parentj,7],":{g,",gtimeparent,":i,",itimeparent,"}",sep = "")                  
              }
              if(comptime == -1 | comptime > endedge)
              {
                  comptime = endedge
              }
              itime = comptime - startedge
              gtime = endedge - comptime
              if(itimeparent == 0)
              {
                  spec2 = paste(linlist[j,7],":{g,",gtime,":i,",itime,":g,0}",sep = "")
              } else {
                  spec2 = paste(linlist[j,7],":{g,",gtime,":i,",itime,"}",sep = "")
              }
          }
          linlist[parentj,7] = paste("(",spec1,",",spec2,")",sep = "")
          linlist[parentj,8] = linlist[j,3]
          linlist = linlist[-j,]               
      } else {
          if(as.numeric(parent) != 0)
          {
              parentj2 = which(L[,1] == as.numeric(parent))                            
              comptimeparent2 = L[parentj2,4]
              
              #print(paste('parentj2 is ',parentj2))
              #print(L)
              #print(paste('comptimeparent is ',comptimeparent2))
              
              if(comptimeparent2 > -1 & (comptimeparent2 < as.numeric(linlist[j,3]) | parentj2 <= -1))
              {
                 linlist[j,4] = L[parentj2,4]
              }
              linlist[j,c(1:3,5)] = L[parentj2,c(1:3,5)]
          }
      }
      if(is.null(nrow(linlist)))
      { 
          done = 1
          if(ig == FALSE)
          {
             linlist[7] = paste(linlist[7],":",abs(as.numeric(linlist[3])),";",sep = "")
          } else {
             linlist[7] = paste(linlist[7],";",sep = "")             
          }
      } else {
          if(nrow(linlist) == 1)
          { 
              done = 1
              if(ig == FALSE)
              {
                 linlist[7] = paste("(",linlist[7],":",age,");",sep = "")
              } else {
                 linlist[7] = paste(linlist[7],";",sep = "")              
              }                                                                  
          }
      }
   }
   }
   return(linlist[7])
}

sampletree = function(L,age,samplemethod = "random")
{
   lenL = length(L[,1])
   if(samplemethod == "random")
   {
      neworder = sample2(1:lenL, replace = F) 
   }   
   if(samplemethod == "youngest")
   {
      neworder = order(L[,3],decreasing = T)
   }
   if (samplemethod == "oldest")
   {
      neworder = order(L[,3],decreasing = F)
   }     
   L2 = L[neworder,]
   ss = NULL;
   for(i in 1:lenL)
   {
       if(L2[i,5] == -1)
       {
           if(is.element(L2[i,6],ss) == FALSE)
           {
              ss = c(ss,L2[i,6])
           } else {
              L2[i,5] = age # peudo extinction just before the present
           }
       }
   }
   L2 = L2[rev(order(L2[,3])),]   
   return(L2)
}
 
pbd_reconstruct = function(L)
{
  L2 = L[order(L[,3]),]
  L3 = L2;
  for(i in 1:length(L2[,2]))
  {
      pai = which(abs(L2[,2]) == L2[i,1]);
      L3[pai,2] = sign(L2[pai,2]) * i;
  }
  orglabs = L3[,1];
  numincspec = length(L3[,2])
  id = 1:numincspec;
  L3[,1] = id;
  L = L3;
  L[1,3] = -1E-10;
  if(L[2,3] == 0)
  {
      L[2,3] = 1E-10;
  }
  
  pa = L[,2]; # vector of parent species
  ti = L[,3]; # vector of speciation-initiation times
  tc = L[,4]; # vector of speciation-completion times
  te = L[,5]; # vector of extinction times
  sl = L[,6]; # vector of species labels              
  id2 = id;
  tr = NULL; ###
  ### print(cbind(id,pa,ti,tc,te,sl))
  
  # find the branch that went extinct last
  idx1 = which(te == max(te) & te > 0)
  while(length(idx1) != 0)
  {
      # does this extinct branch have offspring?
      # find the offspring who have the extinct branch as parent
      idx2 = rev(which(abs(pa) == idx1))[1];
      if(is.na(idx2))
      {
          # extinct branch does not have offspring
          # extinct branch can be neglected
          ti[idx1] = 0;
          tc[idx1] = 0;
          te[idx1] = 0;
          pa[idx1] = 0;
          sl[idx1] = 0;
      } else {
          # extinct branch has offspring
          # find the offspring of the offspring of the extinct branch
          idx3 = which(abs(pa) == idx2)
          # was extinct branch good/incipient
          # at time of initiation of offspring?
          if(pa[idx2] > 0)
          {
              # extinct branch was good
              pa[idx3] = idx1;
          } else {
              # extinct branch was incipient
              pa[idx3] = sign(pa[idx3]) * idx1;
              tc[idx1] = tc[idx2];
          }
          te[idx1] = te[idx2];
          sl[idx1] = sl[idx2]; ###
          tr = rbind(tr,c(id2[idx1],id2[idx2])); ##
          id2[idx1] = id2[idx2]; ###
          ti[idx2] = 0;
          tc[idx2] = 0;
          te[idx2] = 0;
          pa[idx2] = 0;
          sl[idx2] = 0;                             
      }
      #idx1 = rev(which(te > 0))[1];
      idx1 = which(te == max(te) & te > 0)
  }
  ### print(cbind(id,pa,ti,tc,te,sl))
  # eliminate zero rows
  idxs = which(ti != 0); 
  diff = (idxs != (1:length(idxs)));
  while(sum(diff) != 0)
  {
      idx1 = (which(diff == 1))[1];
      idx2 = idxs[idx1];
      ti[idx1] = ti[idx2];
      tc[idx1] = tc[idx2];
      te[idx1] = te[idx2];
      pa[idx1] = pa[idx2]; 
      sl[idx1] = sl[idx2];
      id[idx1] = id[idx2]; 
      id2[idx1] = id2[idx2];
      ti[idx2] = 0;
      tc[idx2] = 0;
      te[idx2] = 0;
      pa[idx2] = 0;
      ### pa[abs(pa) == idx2] = sign(pa[abs(pa) == idx2]) * idx1; ###
      sl[idx2] = 0;
      id[idx2] = 0;
      id2[idx2] = 0;
      idxs = which(ti != 0);
      diff = (idxs != (1:length(idxs)));
  }
  ig = rep(0,length(ti)); # good/incipient flags
  ig[te == -1 & tc != -1] = 1;
  ig[te == -1 & tc == -1] = -1;
  if(te[1] == -1)
  { 
     ig[1] = 1;
  }
  zeros = c(which(sl == 0))
  if(length(zeros) > 0)
  {
     id = id[-zeros]
     id2 = id2[-zeros]
     pa = pa[-zeros]
     ti = ti[-zeros]
     te = te[-zeros]
     tc = tc[-zeros]
     sl = sl[-zeros]
     ig = ig[-zeros]
  }
  ### print(cbind(id,pa,ti,tc,te,sl,ig)); ###

  igg = ig; # copy of table of good/incipient flags
  ppa = pa; # copy of table of parent indices
  its = 0; # index that will run through table
  tt = NULL; # table of splitting times
  pp = NULL; # table of parent indices
  dd = NULL; # table of daughter indices
  sls = NULL; # table of species labels
  idxs = which(igg != 0);
  while(idxs[length(idxs)] > 1)
  {
     idx = which.max(ti[idxs]);
     di = idxs[idx]; # daughter index
     parenti = ppa[di]; # parent index (can be negative!)
     pai = which(id == abs(parenti))
     if(igg[pai] == 1 & parenti > 0 & igg[di] == -1)
     {
         #print('1. parent alive, good at event, good at present, daughter inc at present')
         igg[di] = 0;
         #igg[pai] = 1; This was already the case
     } else {
     if(igg[pai] == 1 & parenti > 0 & igg[di] == 1)
     {
         #print('2. parent alive, good at event, good at present, daughter good at present')
         igg[di] = 0;
         #igg[pai] = 1; This was already the case
         its = its + 1;
         tt[its] = ti[di];
         pp[its] = abs(parenti);
         dd[its] = id[di];
         tr = rbind(tr,c(id[di],id2[di])); ##
         sls[its] = sl[di];
     } else {
     if(igg[pai] == 1 & parenti < 0 & igg[di] == -1)
     {
         #print('3. parent alive, inc at event, good at present, daughter inc at present')
         igg[di] = 0;
         #igg[pai] = 1; This was already the case
     } else {
     if(igg[pai] == 1 & parenti < 0 & igg[di] == 1)
     {
         #print('4. parent alive, inc at event, good at present, daughter good at present')
         igg[di] = 0;
         #igg[pai] = 1; This was already the case
         its = its + 1;
         tt[its] = ti[di];
         pp[its] = abs(parenti);
         dd[its] = id[di];
         #dd[its] = id2[di]; ##
         tr = rbind(tr,c(id[di],id2[di])); ##
         sls[its] = sl[di];
     } else {
     if(igg[pai] == -1 & parenti < 0 & igg[di] == -1)
     {
         #print('5. parent alive, inc at event, inc at present, daughter inc at present')
         igg[di] = 0;
         #igg[pai] = -1; This was already the case
     } else {
     if(igg[pai] == -1 & parenti < 0 & igg[di] == 1)
     {
         #print('6. parent alive, inc at event, inc at present, daughter good at present')
         igg[di] = 0;
         igg[pai] = 1;
         pp[which(pp == id[di])] = abs(parenti);
         tr = rbind(tr,c(id2[pai],id2[di])); ###
         #id2[pai] = id2[di]; ##
         sl[pai] = sl[di]; #
     } else {
     if(igg[pai] == 0 & parenti > 0 & igg[di] == -1)
     {
         #print('7. parent dead, good at event, daughter inc at present')
         igg[di] = 0;
         igg[pai] = 1;
         tr = rbind(tr,c(id2[pai],id2[di])); ###
     } else {
     if(igg[pai] == 0 & parenti > 0 & igg[di] == 1)
     {
         #print('8. parent dead, good at event, daughter good at present')
         igg[di] = 0;
         igg[pai] = 1;
         pp[which(pp == id[di])] = abs(parenti);
         sl[pai] = sl[di]
         tr = rbind(tr,c(id2[pai],id2[di])); ###       
         ## The daughter keeps her own species label
     } else {
     if(igg[pai] == 0 & parenti < 0 & igg[di] == -1)
     {
         #print('9. parent dead, inc at event, daughter inc at present')
         igg[di] = 0;
         igg[pai] = -1;
         tr = rbind(tr,c(id2[pai],id2[di])); ###         
     } else {
     if(igg[pai] == 0 & parenti < 0 & igg[di] == 1)
     {
         #print('10. parent dead, inc at event, daughter good at present')
         igg[di] = 0;
         igg[pai] = 1;
         pp[which(pp == id[di])] = abs(parenti);
         tr = rbind(tr,c(id2[pai],id2[di])); ###        
         sl[pai] = sl[di]; 
     }
     }}}}}}}}}
     idxs = which(igg != 0); 
  }
  dd = c(1,dd); ##
  pp = c(0,pp); ##
  tt = c(-1e-10,tt); ##
  tt2 = c(0,tt); ##
  te = c(-1,te); ##
  sls = c(sl[1],sls); ##
  dd2 = dd; ##
  pp2 = pp; ##
  for(i in 1:length(tr[,1]))
  {
      dd2[which(dd2 == tr[i,1])] = tr[i,2]
      pp2[which(pp2 == tr[i,1])] = tr[i,2]
  }
  dd = dd2; ##
  pp = pp2; ##
  reconL = cbind(dd,pp,tt,tt,rep(-1,length(dd)),sls,deparse.level = 0)
  ## reconL = rbind(c(1,0,-1e-10,0,-1,1),cbind(dd,pp,tt,tt,rep(-1,length(dd)),sls,deparse.level = 0))
  ### print(reconL); ###
  L = reconL;
  L[,1] = orglabs[reconL[,1]];
  L[,2] = c(0,orglabs[reconL[,2]]);
  reconL = L;
  ### print(reconL); ###
  return(reconL)
}

L2phylo2 = function(L,dropextinct = T)
# makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
   L = L[order(abs(L[,3])),]  
   age = L[1,1]
   L[,1] = age - L[,1]
   L[1,1] = -1
   notmin1 = which(L[,4] != -1)
   L[notmin1,4] = age - L[notmin1,4]
   if(dropextinct == T)
   {
      sall = which(L[,4] == -1)
      tend = age
   } else {
      sall = which(L[,4] >= -1)
      tend = (L[,4] == -1) * age + (L[,4] > -1) * L[,4]
   }
   specid = L[,6]
   L = L[,-(4:6)]
   linlist = cbind(L[sall,],paste("S",specid,"-",specid,"-",abs(L[sall,3]),sep = ""),tend)
   done = 0
   while(done == 0)
   {
      #print(linlist)
      j = which.max(linlist[,1])
      daughter = linlist[j,3]
      parent = linlist[j,2]
      parentj = which(parent == linlist[,3])
      parentinlist = length(parentj)
      if(parentinlist == 1)
      {
         spec1 = paste(linlist[parentj,4],":",as.numeric(linlist[parentj,5]) - as.numeric(linlist[j,1]),sep = "")
         spec2 = paste(linlist[j,4],":",as.numeric(linlist[j,5]) - as.numeric(linlist[j,1]),sep = "")
         linlist[parentj,4] = paste("(",spec1,",",spec2,")",sep = "")
         linlist[parentj,5] = linlist[j,1]
         linlist = linlist[-j,]
      } else {      
         #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
         linlist[j,1:3] = L[which(L[,3] == as.numeric(parent)),1:3]
      }
      if(is.null(nrow(linlist))) { done = 1 }
   }
   linlist[4] = paste(linlist[4],":",linlist[5],";",sep = "")
   phy = read.tree(text = linlist[4])
   tree = as.phylo(phy)
   return(tree)
}
