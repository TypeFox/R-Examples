Gibbs = function(pij, iter, species, abund,  hyperParam, fast.rmultinom.weight, readWeights)
{
  time1 <- Sys.time() 

####### Matrices to record the mixing weights and species assignments through iterations
  abundances<-matrix(0, ncol=(length(species)+1), nrow=iter)
  abundances[,1]=1:iter
  assignments<-matrix(0, ncol=(length(species)+1), nrow=iter)
  assignments[,1]=1:iter
  
#### select  the species you want and at the same time rearrange columns of pij.temp to be same with species vector, so multiplication work

  pij.temp<-pij[,species]
  
#### make matrix that will hold the assignments
  assignm<-matrix(0, nrow=(dim(pij.temp)[1]), ncol=(dim(pij.temp)[2]))
  colnames(assignm)<-colnames(pij.temp)
  rownames(assignm)<-rownames(pij.temp)

  timeA <- Sys.time() 

####### Begin iterations
  for (i in 1:iter) {
    zij <- as.matrix(pij.temp)
    w<-abund
 

#####Do at the same time STEP1 and STEP2
######## STEP1.Calculate  zij {zij=pijwj/sum_j(pijwj)}

######## STEP2. Sample from multinomial ~ Mult(draw=1,objects = 1, prob = zij)
 
    assignm<-fast.rmultinom.weight(proba.matrix=zij, z.matrix=assignm, weights=as.vector(w))

    assignmWeighted<- assignm[rownames(readWeights),] * readWeights[,"weight"] ### multiply assignment (0,1) by read weight


######## STEP3. Calculate nj=sum_i(zij)
    nj <- colSums(assignmWeighted)


######## STEP4. Generate w^(t) from \pi (w|z^(t)). Sample w from Dirichlet (a1+n1, ..., ak+nk)
####new parameters for Dirichlet
    alpha<-nj + hyperParam
    abund<-rdirichlet(1,alpha)
    colnames(abund)<-colnames(assignmWeighted)



##### Combine assignment info with read names    
    assignedReads<-cbind.data.frame(read=rownames(pij.temp), assignmWeighted,  stringsAsFactors=FALSE)


###### record output
    assignments[i,2:(length(species)+1)]<- nj

    abundances[i,2:(length(species)+1)]<-w
    

#### update progress bar
    #setTxtProgressBar(pb, i)

    
} ###end of iterat

  
### Processing time
  time2<-Sys.time()
  timeDiff<-time2-time1
  iterTimeDiff<-time2-timeA

  abundances=data.frame(abundances)
  names(abundances)<-c("Iter", colnames(zij))
  assignments=data.frame(assignments)
  names(assignments)<-c("Iter", colnames(zij))


  
  result<-list("abundances"=abundances, "assignments"=assignments, "RunningTime"=timeDiff, "RunningTimeIterations"=iterTimeDiff, "pijs"= pij.temp, "assignedReads"=assignedReads)
  return(result)
}
