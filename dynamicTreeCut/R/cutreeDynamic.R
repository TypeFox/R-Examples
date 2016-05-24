#------------------------------------------------------------------------------------------
#
# cutreeDynamic
#
#------------------------------------------------------------------------------------------
# Modification(s) by Peter Langfelder: returns numerical labels (not colors) in a vector (not a factor).
# Several rendundant blocks of code removed; duplicate function definitions removed; unused
# functions removed.
# maxTreeHeight now checked for being too large.

#.minAttachModuleSize = 100;

cutreeDynamicTree = function(dendro, maxTreeHeight=1, deepSplit=TRUE, minModuleSize=50)
{
    if (is.null(maxTreeHeight)) maxTreeHeight = 0.99 * max(dendro$height);
    if (maxTreeHeight > max(dendro$height)) maxTreeHeight = 0.99 * max(dendro$height);
    staticCutCluster = .cutTreeStatic(dendro=dendro, heightcutoff=maxTreeHeight, minsize1=minModuleSize)

    #get tree height for every singleton
    #node_index   tree_height
    demdroHeiAll= rbind( cbind(dendro$merge[,1], dendro$height), cbind(dendro$merge[,2], dendro$height) )

    #singletons will stand at the front of the list
    myorder = order(demdroHeiAll[,1])

    #get # of singletons
    no.singletons = length(dendro$order)
    demdroHeiAll.sort = demdroHeiAll[myorder, ]
    demdroHei.sort = demdroHeiAll.sort[c(1:no.singletons), ]
    demdroHei = demdroHei.sort[seq(no.singletons, 1, by=-1), ]
    demdroHei[,1]     = -demdroHei[,1]

    # combine with prelimilary cluster-cutoff results
    demdroHei         = cbind(demdroHei, as.integer(staticCutCluster))

    # re-order the order based on the dendrogram order dendro$order
    demdroHei.order = demdroHei[dendro$order, ]

    static.clupos = .locateCluster(demdroHei.order[, 3])

    if (is.null(static.clupos) ){
        module.assign   = rep(0, no.singletons)
        return ( module.assign )
    }

    static.no     = dim(static.clupos)[1]
    static.clupos2 =     static.clupos
    static.no2     =     static.no

    #split individual cluster if there are sub clusters embedded
    mcycle=1
    while(1==1){
        clupos = NULL
        for (i in c(1:static.no)){
           mydemdroHei.order = demdroHei.order[ c(static.clupos[i,1]:static.clupos[i,2]), ] #index to [1, clusterSize]
           mydemdroHei.order[, 1] = mydemdroHei.order[, 1] - static.clupos[i, 1] + 1

           #cat("Cycle ", as.character(mcycle), "cluster (", static.clupos[i,1], static.clupos[i,2], ")\n")
           #cat("i=", as.character(i), "\n")

           iclupos = .processIndividualCluster(mydemdroHei.order, 
                                              cminModuleSize       = minModuleSize, 
                                              cminAttachModuleSize = 2*minModuleSize)

           iclupos[,1] = iclupos[,1] + static.clupos[i, 1] -1 #recover the original index
           iclupos[,2] = iclupos[,2] + static.clupos[i, 1] -1

           clupos  = rbind(clupos, iclupos) #put in the final output buffer
        }

        if(deepSplit==FALSE){
          break
        }

        if(dim(clupos)[1] != static.no) {
           static.clupos = clupos
           static.no     = dim(static.clupos)[1]
        }else{
          break
        }
       mcycle = mcycle + 1
       #static.clupos  
    }
    final.cnt = dim(clupos)[1]
    #assign colors for modules
    module.assign = rep(0, no.singletons)
    module.cnt=1
    for (i in c(1:final.cnt ))
    {
       sdx = clupos[i, 1] #module start point
       edx = clupos[i, 2] #module end point
       module.size = edx - sdx +1 
       if(module.size <minModuleSize){
         next
       }
       #assign module lable
       module.assign[sdx:edx] = rep(module.cnt, module.size)
       #update module label for the next module
       module.cnt = module.cnt + 1
    }
    colcode.reduced.order = .assignModuleNumber(module.assign, minsize1=minModuleSize)
   recov.order = order( demdroHei.order[,1])
    colcode.reduced = colcode.reduced.order[recov.order]
    colcode.reduced
}

#leftOrright >0 : running length (with same sign) to right, otherwise to the left
#mysign = -1: negative value, mysign = -1: positive value
.runlengthSign = function(mysequence, leftOrright=-1, mysign=-1){
    seqlen = length(mysequence)   
    if(leftOrright<0){
        pseq = rev(mysequence)
    }else{
        pseq = mysequence
    }

    if(mysign<0){ #see where the first POSITIVE number occurs
       nonezero.bool = (pseq > 0)
    }else{ #see where the first NEGATIVE number occur
       nonezero.bool = (pseq < 0)
    }
    if( sum(nonezero.bool) > 0){
      runlength = min( c(1:seqlen)[nonezero.bool] ) - 1
    }else{
      runlength = 0
    }
}

#"0" is for grey module
#.assignModuleColor = function(labelpred, minsize1=50, anameallmodules=FALSE, auseblackwhite=FALSE) {
    # here we define modules by using a height cut-off for the branches
    #labelpred= cutree(dendro,h=heightcutoff)
    #cat(labelpred)

    #"0", grey module doesn't participate color assignment, directly assigned as "grey"
    #labelpredNoZero = labelpred[ labelpred >0 ]
    #sort1=-sort(-table(labelpredNoZero))
    ## sort1
    #modulename= as.numeric(names(sort1))
    #modulebranch= sort1 >= minsize1
    #no.modules = sum(modulebranch)
#
    #colorcode=GlobalStandardColors;
#
    ##"grey" means not in any module;
    #colorhelp=rep("grey",length(labelpred))
    #if ( no.modules==0){
        #print("No module detected.")
    #}
    #else{
        #if ( no.modules > length(colorcode)  ){
            #print( paste("Too many modules:", as.character(no.modules)) )
        #}
#
        #if ( (anameallmodules==FALSE) || (no.modules <=length(colorcode)) ){
            #labeledModules = min(no.modules, length(colorcode) )
            #for (i in c(1:labeledModules)) {
               #colorhelp=ifelse(labelpred==modulename[i],colorcode[i],colorhelp)
            #}
            #colorhelp=factor(colorhelp,levels=c(colorcode[1:labeledModules],"grey"))
        #}else{#nameallmodules==TRUE and no.modules >length(colorcode)
            #maxcolors=length(colorcode)
            #labeledModules = no.modules
            #extracolors=NULL
            #blackwhite=c("red", "black")
            #for(i in c((maxcolors+1):no.modules)){
              #if(auseblackwhite==FALSE){
                  #icolor=paste("module", as.character(i), sep="")
              #}else{#use balck white alternatively represent extra colors, for display only
                  ##here we use the ordered label to avoid put the same color for two neighboring clusters
                  #icolor=blackwhite[1+(as.integer(modulename[i])%%2) ]
              #}
              #extracolors=c(extracolors, icolor)
            #}
#
            ##combine the true-color code and the extra colorcode into a uniform colorcode for 
            ##color assignment
            #allcolorcode=c(colorcode, extracolors)
#
            #for (i in c(1:labeledModules)) {
               #colorhelp=ifelse(labelpred==modulename[i],allcolorcode[i],colorhelp)
            #}
            #colorhelp=factor(colorhelp,levels=c(allcolorcode[1:labeledModules],"grey"))
        #}
    #}
#

    #colorhelp
#}

# This function written by Peter Langfelder, based on .assignModuleColor above but simplified.
# Assigns module numbers, not colors. All modules are labeled.
#"0" is for grey module
.assignModuleNumber = function(labelpred, minsize1=50) {
    #"0", grey module doesn't participate color assignment, directly assigned as "grey"
    labelpredNoZero = labelpred[ labelpred >0 ]
    sort1=-sort(-table(labelpredNoZero))
    # sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 >= minsize1
    no.modules = sum(modulebranch)

    #"grey" means not in any module;
    colorhelp=rep(0,length(labelpred))

    for (i in c(1:no.modules)) {
        colorhelp=ifelse(labelpred==modulename[i],i,colorhelp)
    }

    colorhelp
}


#locate the start/end positions of each cluster in the ordered cluster label sequence 
#where "-1" indicating no cluster
#3-1 -1 1 1 1 1 2 2 2
#3 3 -1-1 1 1 1 1 2 2 2  (shift)
#---------------------------------
#0-4  0 2 0 0 0 1 0 0 0   (difference)
#       *     * @
.locateCluster = function(clusterlabels)
{
 no.nodes = length(clusterlabels)
 clusterlabels.shift = c(clusterlabels[1], c(clusterlabels[1:(no.nodes-1)]) )
 
 #a non-zero point is the start point of a cluster and it previous point is the end point of the previous
 #cluster 
 label.diff = abs(clusterlabels - clusterlabels.shift)

 #process the first and last positions as start/end points if they belong to a cluster instead of no
 # cluster "-1" 
 if(clusterlabels[1]       >0) {label.diff[1]=1} 
 if(clusterlabels[no.nodes]>0) {label.diff[no.nodes]=1} 

 flagpoints.bool = label.diff > 0
 if( sum(flagpoints.bool) ==0){
     return(NULL)
 }

 flagpoints = c(1:no.nodes)[flagpoints.bool]
 no.points  = length(flagpoints)

 myclupos=NULL
 for(i in c(1:(no.points-1)) ){
   idx = flagpoints[i]
   if(clusterlabels[idx]>0){
      if(flagpoints[i+1]==no.nodes) {#boundary effect
         myclupos = rbind(myclupos, c(idx, flagpoints[i+1]) )
         break
      }else{
         myclupos = rbind(myclupos, c(idx, flagpoints[i+1]-1) )
      }
   }
 }
 myclupos
}


#input is the cluster demdrogram of an individual cluster, we want to find its embbeded subclusters
#execution order: mean-height ==> (mean+max)/2 ==> (mean+min)/2
#useMean: =0 ~ use mean-height   as calibation line
#         =1 ~ use (mean+max)/2  as calibation line to detect relatively a small cluster sitting on the head of a bigger one,
#                      so mean-height is too low to detect the two modules.
#         =-1~ use (mean+min)/2  as calibation line to detect relatively a small cluster sitting on the tail of a bigger one,
#                      so mean-height & (mean+max)/2 are too high to detect the two modules

.processIndividualCluster = function(clusterDemdroHei, cminModuleSize=50, 
     cminAttachModuleSize = 2* cminModuleSize, 
     minTailRunlength= as.integer(cminModuleSize/3)+1, useMean=0)
{
    #for debug: use all genes
    #clusterDemdroHei =demdroHei.order

    no.cnodes = dim(clusterDemdroHei)[1]
    
    cmaxhei   = max(clusterDemdroHei[, 2])
    cminhei   = min(clusterDemdroHei[, 2])
    
    cmeanhei  = mean(clusterDemdroHei[, 2])
    cmidhei = (cmeanhei + cmaxhei)/2.0
    cdwnhei = (cmeanhei + cminhei)/2.0

    if (useMean==1){
        comphei = cmidhei
    }else if (useMean==-1){
        comphei = cdwnhei
    }else{ #normal case
        comphei = cmeanhei
    }
        
    # compute height diffrence with mean height
    heidiff       = clusterDemdroHei[,2] - comphei
    heidiff.shift = .shiftSequence(heidiff, -1)

    # get cut positions
    # detect the end point of a cluster, whose height should be less than meanhei 
    #  and the node behind it is the start point of the next cluster which has a height above meanhei
    cuts.bool = (heidiff<0) & (heidiff.shift > 0)
    cuts.bool[1]         = TRUE
    cuts.bool[no.cnodes] = TRUE
    
    if(sum(cuts.bool)==2){
          if (useMean==0){
             new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }else{
             new.clupos = rbind(c(1, no.cnodes))
          }
          return (new.clupos)
    }

    #a good candidate cluster-end point should have significant # of ahead nodes with head < meanHei
    cutindex =c(1:no.cnodes)[cuts.bool]
    no.cutps = length(cutindex)
    runlens  = rep(999, no.cutps)
    cuts.bool2 = cuts.bool
    for(i in c(2:(no.cutps-1)) ){
       seq = c( (cutindex[i-1]+1):cutindex[i] )
       runlens[i] = .runlengthSign(heidiff[seq], leftOrright=-1, mysign=-1)

       if(runlens[i] < minTailRunlength){
          #cat("run length=", runlens[i], "\n")
          cuts.bool2[ cutindex[i] ] = FALSE
       }
    }

    #attach SMALL cluster to the left-side BIG cluster if the small one has smaller mean height
    cuts.bool3=cuts.bool2
    if(sum(cuts.bool2) > 3) {
       curj = 2
       while (1==1){
           cutindex2 =c(1:no.cnodes)[cuts.bool2]
           no.clus = length(cutindex2) -1
           if (curj>no.clus){
              break
           }
           pre.sdx = cutindex2[ curj-1 ]+1 #previous module start point
           pre.edx = cutindex2[ curj ] #previous module end   point
           pre.module.size = pre.edx - pre.sdx +1 
           pre.module.hei  = mean(clusterDemdroHei[c(pre.sdx:pre.edx) , 2])
         
           cur.sdx = cutindex2[ curj ]+1 #previous module start point
           cur.edx = cutindex2[ curj+1 ] #previous module end   point
           cur.module.size = cur.edx - cur.sdx +1 
           cur.module.hei  = mean(clusterDemdroHei[c(cur.sdx:cur.edx) , 2])

           #merge to the leftside major module, don't change the current index "curj"
           #if( (pre.module.size >minAttachModuleSize)&(cur.module.hei<pre.module.hei)&(cur.module.size<minAttachModuleSize) ){
           if( (cur.module.hei<pre.module.hei)&(cur.module.size<cminAttachModuleSize) ){
                cuts.bool2[ cutindex2[curj] ] = FALSE
           }else{ #consider next cluster
                curj = curj + 1
           }
       }#while
   }#if 

   cutindex2 =c(1:no.cnodes)[cuts.bool2]
   no.cutps = length(cutindex2)

    #we don't want to lose the small cluster at the tail, attch it to the previous big cluster
    #cat("Lclu= ", cutindex2[no.cutps]-cutindex2[no.cutps-1]+1, "\n")
    if(no.cutps > 2){
      if( (cutindex2[no.cutps] - cutindex2[no.cutps-1]+1) < cminModuleSize ){
        cuts.bool2[ cutindex2[no.cutps-1] ] =FALSE  
      }
    }

   cutindex2  = c(1:no.cnodes)[cuts.bool2]
   cutindex2[1]=cutindex2[1]-1 #the first 
   no.cutps2  = length(cutindex2)

   if(no.cutps2 > 2){
     new.clupos = cbind( cutindex2[c(1:(no.cutps2-1))]+1, cutindex2[c(2:no.cutps2)] )
   }else{
     new.clupos = cbind( 1, no.cnodes)
   }

   if ( dim(new.clupos)[1] == 1 ){   
          if (useMean==0){
             new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=1)
          }else if(useMean==1){
             new.clupos=.processIndividualCluster(clusterDemdroHei=clusterDemdroHei, cminModuleSize=cminModuleSize, 
                                         cminAttachModuleSize=cminAttachModuleSize,
                                         useMean=-1)          
          }   
   }
   new.clupos
}


#delta >0 : shift to right, otherwise to the left
.shiftSequence = function(mysequence, delta){
    seqlen = length(mysequence)
    if(delta>0){
        finalseq=c(mysequence[1:delta], mysequence[1:(seqlen-delta)])
    }else{
        posdelta = -delta 
        finalseq=c(mysequence[(posdelta+1):seqlen], mysequence[(seqlen-posdelta+1):seqlen])
    }
    finalseq
}

#use height cutoff to remove
.cutTreeStatic = function(dendro,heightcutoff=0.99, minsize1=50) {

    # here we define modules by using a height cut-off for the branches
    labelpred= cutree(dendro,h=heightcutoff)
    sort1=-sort(-table(labelpred))
    sort1
    modulename= as.numeric(names(sort1))
    modulebranch= sort1 >= minsize1
    no.modules=sum(modulebranch)

    colorhelp = rep(-1, length(labelpred) )
    if ( no.modules==0){
        print("No module detected")
    }
    else{
        for (i in c(1:no.modules)) {
            colorhelp=ifelse(labelpred==modulename[i],i ,colorhelp)
        }
    }
    colorhelp
}

