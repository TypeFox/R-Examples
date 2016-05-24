llra.datprep<-function(X, mpoints, groups, baseline=NULL){
 
    Xwide <- X
    if (ncol(Xwide) %% mpoints > 0)
       stop("Number of items must be the same for each timepoint.")
    nitems <- dim(Xwide)[2]/mpoints

   if(missing(groups)) groups <- rep("CG",dim(Xwide)[1])

    covs.prep<-function(groups,baseline){
      groups<-as.matrix(groups)
      grstr<-apply(groups,1,paste,collapse=":")
      grstr <- factor(grstr)    
      if(!is.null(baseline))    
        {
          basel <- paste(baseline,collapse=":")
          grstr <- relevel(grstr,basel) 
        }
      cov.groupvec<-as.numeric(grstr)  
      names(cov.groupvec)<-grstr
      cov.groupvec
    }
    
    # sort data according to cov.groupvec
    cov.groupvec<-covs.prep(groups,baseline)
    Xwide<-Xwide[order(cov.groupvec,decreasing=TRUE),] 
    cov.groupvec<-sort(cov.groupvec)

    # number of people per group
    grp_n<-table(cov.groupvec)
    names(grp_n)<-unique(names(cov.groupvec))
    
    # convert to long format
    Xlong <- matrix(unlist(Xwide), ncol = mpoints)
    
    # assignment vector item x treatment
    assign.vec <- as.vector(sapply(1:nitems, function(i)
                    cov.groupvec + (i-1)*max(cov.groupvec)))
    assign.vec <- rev(assign.vec)
    assign.vec <- abs(assign.vec-max(assign.vec))+1 
    names(assign.vec)<-rev(rep(names(cov.groupvec),nitems)) 
    list(X=Xlong, assign.vec=assign.vec, grp_n=grp_n, nitems=nitems)
}
