BPEC.TreeEdges = function(MCMCout)
{
  count = MCMCout$countR
  root = which.max(MCMCout$RootProbsR)
  levels = MCMCout$levelsR
  datsiz = MCMCout$NoSamplesR
  clado = MCMCout$cladoR

  SeqLabels =MCMCout$SeqsFileR[MCMCout$SeqLabelsR]

 edgelist = NULL
 plotheight=sqrt(2)
 levels=levels+1

 if(length(SeqLabels)<count)
   {
     SeqLabels=c(SeqLabels,rep(0,count-length(SeqLabels)))
   }
 
 counter=max(SeqLabels)+1
 for(i in 1:count)
   {
     if(SeqLabels[i]==0)
       {
         SeqLabels[i]=counter
         counter=counter+1
       }
   }
 
 NodeOrder=array(0,count)
 newNodeOrder=array(0,count)
 MaxLevel=max(levels)
 NodeOrder[1]=root
 newNodeOrder[1]=root
 descendantcounter=2
 descendants=0
 
 for(i in 1:MaxLevel)
   {
     prevdescendants=descendants
     
     descendants=1
     NodeOrder=newNodeOrder
     for(j in 1:count)
       {
         if(levels[j]==i)
           {
             descendants=descendants+1
           }
       }
     prevord=descendantcounter-1
    # print(descendantcounter)
    # print(prevord)
     descendantcounter=1
     previousone=-1
     
     if(prevord>0)
       {
         for(j in 1:prevord)
           {             
             for(l in 1:count)
               {
                 if((clado[(NodeOrder[j]-1)*count+l]==1||clado[(l-1)*count+NodeOrder[j]]==1)&&levels[l]==i)
                   {
                      # print(l)
                                        # edgelist output
                     edgelist = rbind(edgelist, data.frame(vert.from=SeqLabels[NodeOrder[j]],vert.to=SeqLabels[l],level=i,ssize=datsiz[l]))                   
                     
                     previousone=l
                     newNodeOrder[descendantcounter]=l
                     descendantcounter=descendantcounter+1
                     next
                   }
               }
           }
       }
   }
 return(edgelist)
}

