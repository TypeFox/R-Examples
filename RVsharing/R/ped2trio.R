ped2trio <- function(ped)
{
# Converts a pedigree object into a recursive list of trio objects
# ped: a pedigree object

dv = kindepth(ped)
# md is the depth of the pedigree, i.e. the number of levels of trios
md = max(dv)

# Extracting the ids of subjects and their parents from ped object
id = as.character(ped$id)
dad.id = mom.id = rep("0",length(id))
dad.id[ped$findex>0] = ped$id[ped$findex]
mom.id[ped$mindex>0] = ped$id[ped$mindex]

trio.flag = logical(length(id))
# The same subject can have multiple trios. This will need to be taken into account
names(trio.flag) = id

trio.list=list()
cumultrio.vec=numeric(md)
foundertrios = logical()

r=1
for (lev in md:2)
  {
  # Create vector of unique moms and dads at the current depth
  currentdads = unique(dad.id[dv==lev])
  currentmoms = unique(mom.id[dv==lev])
  # Create vector of unique dads at the current depth who are non founders
  dad.tokeep = dad.id[id%in%currentdads]!=0
  names(dad.tokeep) = id[id%in%currentdads]
  dad.trio = names(dad.tokeep)[dad.tokeep]
  # Loop over the non-founder dads at the current depth
  if (length(dad.trio)>0)
  {
  for (i in 1:length(dad.trio))
    {
    spousevec = unique(mom.id[dad.id==dad.trio[i]])
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[dad.id==dad.trio[i] & mom.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      # I need to convert the offspring vector into a list to handle sibships that contain both trios and final descendants represented by a character string
      offspring.list = as.list(offspring.vec)
      # For children who already have their own trio, replace their name by the actual trio object in the list of offspring
      for (h in which(trio.flag[offspring.vec]))
        {
#        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[trio.flag[offspring.vec]][h],sep="")]]
#        foundertrios[offspring.vec[trio.flag[offspring.vec]][h]] = FALSE
        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[h],sep="")]]  
        foundertrios[offspring.vec[h]] = FALSE
        }  
      trio.list[[r]]=new("Trio", id=dad.trio[i],spouse=spousevec[j],offspring = offspring.list)
      names(trio.list)[r]=paste("trio",dad.trio[i],sep="")
      foundertrios[r] = TRUE
      names(foundertrios)[r] = dad.trio[i]
      r = r +1
      }
    }
  # We turn the trio.flag of the dad.trio to TRUE to indicate they now have a trio
  trio.flag[dad.trio] = TRUE
  }
  # Create vector of unique moms at the current depth who are non founders
  mom.tokeep = mom.id[id%in%currentmoms]!=0
  names(mom.tokeep) = id[id%in%currentmoms]
  mom.trio = names(mom.tokeep)[mom.tokeep]
  # Loop over the non-founder moms at the current depth
  if (length(mom.trio)>0)
  {
  for (i in 1:length(mom.trio))
    {
    spousevec = unique(dad.id[mom.id==mom.trio[i]])
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[mom.id==mom.trio[i] & dad.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      offspring.list = as.list(offspring.vec)
      # For children who already have their own trio, replace their name by the actual trio object in the list of offspring
      for (h in which(trio.flag[offspring.vec]))
        {
#        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[trio.flag[offspring.vec]][h],sep="")]]  
#        foundertrios[offspring.vec[trio.flag[offspring.vec]][h]] = FALSE
        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[h],sep="")]]  
        foundertrios[offspring.vec[h]] = FALSE
        }  
      trio.list[[r]]=new("Trio",id=mom.trio[i],spouse=spousevec[j],offspring = offspring.list)
      names(trio.list)[r]=paste("trio",mom.trio[i],sep="")
      foundertrios[r] = TRUE
      names(foundertrios)[r] = mom.trio[i]
      r = r +1
      }
    }
  # We turn the trio.flag of the mom.trio to TRUE to indicate they now have a trio
  trio.flag[mom.trio] = TRUE
  }
  # Record number of trios up to now
  cumultrio.vec[md-lev+1] = r - 1
  }
#  print(trio.flag)
# Depth 1
  # Create vector of unique moms and dads at the current depth
  currentdads = unique(dad.id[dv==1])
  currentmoms = unique(mom.id[dv==1])

  # All current dads are founders at depth 1. We give them all a trio
  dad.trio = currentdads
  # Loop over the non-founder dads at the current depth
  if (length(dad.trio)>0)
  {
  for (i in 1:length(dad.trio))
    {
    spousevec = unique(mom.id[dad.id==dad.trio[i]])
#    print (dad.trio[i])
#    print (spousevec)
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[dad.id==dad.trio[i] & mom.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      offspring.list = as.list(offspring.vec)
      # For children who already have their own trio, replace their name by the actual trio object in the list of offspring
      for (h in which(trio.flag[offspring.vec]))
        {
#        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[trio.flag[offspring.vec]][h],sep="")]]  
#        foundertrios[offspring.vec[trio.flag[offspring.vec]][h]] = FALSE
        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[h],sep="")]]  
        foundertrios[offspring.vec[h]] = FALSE
        }  
      trio.list[[r]]=new("Trio",id=dad.trio[i],spouse=spousevec[j],offspring = offspring.list)
      names(trio.list)[r]=paste("trio",dad.trio[i],"_",j,sep="")
      foundertrios[r] = TRUE
      names(foundertrios)[r] = dad.trio[i]
      r = r +1
      }
    # We turn the trio.flag of the spouses to TRUE to indicate they are now members of a trio
    trio.flag[spousevec] = TRUE      
    }
#    print(trio.flag)
  }
  # Create vector of unique moms at the current depth who are not members of a trio
  mom.tokeep = !trio.flag[id%in%currentmoms]
#  print (mom.tokeep)
  mom.trio = names(mom.tokeep)[mom.tokeep]
  # Loop over the non-founder moms at the current depth
  if (length(mom.trio)>0)
  {
  for (i in 1:length(mom.trio))
    {
    spousevec = unique(dad.id[mom.id==mom.trio[i]])
    for (j in 1:length(spousevec))
      {
      # Create vector of children of the couple
      offspring.vec = id[mom.id==mom.trio[i] & dad.id==spousevec[j]]
      # Creating the elements of the trio object in trio.list
      offspring.list = as.list(offspring.vec)
      # For children who already have their own trio, replace their name by the actual trio object in the list of offspring
      for (h in which(trio.flag[offspring.vec]))
        {
#        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[trio.flag[offspring.vec]][h],sep="")]]  
#        foundertrios[offspring.vec[trio.flag[offspring.vec]][h]] = FALSE
        offspring.list[[h]] = trio.list[[paste("trio",offspring.vec[h],sep="")]]  
        foundertrios[offspring.vec[h]] = FALSE
        }  
      trio.list[[r]]=list(id=mom.trio[i],spouse=spousevec[j],offspring = offspring.list)
      names(trio.list)[r]=paste("trio",mom.trio[i],"_",j,sep="")
      foundertrios[r] = TRUE
      names(foundertrios)[r] = mom.trio[i]
      r = r +1
      }
    }
  # We turn the trio.flag of the mom.trio to TRUE to indicate they now have a trio
  trio.flag[mom.trio] = TRUE
  }
# Record the total number of trios
cumultrio.vec[md] = r - 1 
# Compute number of trios per level
#ped.vec = diff(c(0,cumultrio.vec))

#trio.obj = trio.list[[r-1]]
  trio.obj = trio.list[which(foundertrios)]

#return(list(obj.list = trio.list, object = trio.obj, fd.indices=which(dv==0)))
return(list(object = trio.obj, fd.indices=which(dv==0)))
}
