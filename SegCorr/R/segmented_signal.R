#################################################
###### Running Multi-seg  for SNP signal ########
#################################################
### Summary: Using the cghseg package, we perform multivariate segmentation on the SNP signal
###          If there the SNP signal is heterogeneous among the different patients, we perform segmentation 
###          on the different variance groups 

### Input: SNP.Chr:   the SNP signal for a given chromosome CHR -- NA's removed
###        group:     the variance classification


### Output: mu.SNP :  mean segmented signal
options(warn=1)
segmented_signal = function(SNP.Chr ,group){
  
  if(missing(SNP.Chr)){
    stop("segmented_signal: SNP matrix missing")
  }else{
    
    
    p = dim(SNP.Chr)[1] # no of probes
    no = dim(SNP.Chr)[2]# no of patients
    
      if(missing(group)){
      warning('segmented_signal: group missing - Default value used ')
      group = rep(1,no)
     }else if(length(group)!=no){ 
      stop('segmented_signal: Check length of group')
      }
    
    mu.SNP = matrix(NA,p,no)
  
    k = length(unique(group))# no of variance groups
    
    for(j in 1:k){
      cat(paste('Performing SNP mean segmentation for group =',j,sep=" "),sep="\n")
      
      # pat.gp = pat.SNP[group.classification==j]
      loc.gp = which(group==j)#which(colnames(SNP.Chr)%in%pat.gp)
      
      datax = SNP.Chr[,loc.gp]
      
      CGHd = new("CGHdata",Y=datax)
      CGHo = new("CGHoptions")
      CGHr = multiseg(CGHd,CGHo)
      mu.SNP[,loc.gp] = as.matrix(getsegprofiles(CGHr))
    }
    
    return(mu.SNP)
    #cat('Done!')
  }
}
