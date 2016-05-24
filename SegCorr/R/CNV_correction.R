#########################################################################
####################### CNV correction step   ###########################
#########################################################################
### Summary: Correcting the gene expression signal for CNV. 
###          Since the two signals are not aligned we need to fix this first.
##           We create a matrix that tells us how which SNP probes belong to a certain  gene.
##           There are many genes for which we dont have a smoothed SNP signal. 
##           Therefore we need to interpolate it. Finally, we perform the regression EXP vs. smoothed mean segmented SNP.


### Input: s.Position.EXP:       start Position of gene
###        e.Position.EXP:       end Position of gene
###        Position.SNP:         SNP position
###        mu.SNP:               mean segmented SNP signal
###        EXP:                  Expression signal

### Output: correction:          CNV corrected signal

options(warn=1)
CNV_correction = function(s.Position.EXP, e.Position.EXP, Position.SNP, mu.SNP, EXP){
  
  if(missing(s.Position.EXP) | missing(e.Position.EXP) ){
    stop('CNV_correction : start/end position for gene expression missing')
  }else if( length(s.Position.EXP)!=length(e.Position.EXP)){ 
    stop('CNV_correction : Dimension mismatch between start and end position for gene expression')
  }
  if(missing(Position.SNP)){
    stop('CNV_correction : SNP position missing')   
  }
  
  if(missing(mu.SNP)){
    stop('CNV_correction : mu.SNP missing')   
  }
  
  if(missing(EXP)){
    stop('CNV_correction : EXP missing')   
  }
  
  n = dim(EXP)[2]
  
  if(dim(mu.SNP)[2]!=dim(EXP)[2]){
    stop("CNV_correction : Dimension mismatch between mu.SNP and end EXP")
  }
  
  Position.common = .aligning.positions(s.Position.EXP, e.Position.EXP, Position.SNP)
  
  mu.SNP.interpolated = .SNP.interpolated(Position.common,Position.SNP,n,mu.SNP)
  correction = .corrected.EXP(Position.common,s.Position.EXP,n, mu.SNP.interpolated, EXP)
 # beta = correction$coef.matrix
  correction = correction$residuals
  
  return(correction)
}



















######################################
########## Aligning Positions ########
######################################
## Summary: 

### Input: s.Position.EXP:       start Position of gene
###        s.Position.EXP:       start Position of gene
###        Position.SNP:         SNP position



.aligning.positions = function(s.Position.EXP, e.Position.EXP, Position.SNP){
  
  pp = rep(NA,length(s.Position.EXP))
    Position.common = c()
    cat('Performing signal alignment')
    for(i in 1:length(s.Position.EXP)){
      
      loc.SNP = which((s.Position.EXP[i]<=Position.SNP) *  (Position.SNP<=e.Position.EXP[i]) == 1)## checking which SNPs belong to that interval
      #loc.SNP = intersect(which(s.Position.EXP[i]<=Position.SNP),which(Position.SNP<=e.Position.EXP[i]))
     pp[i] = p.SNP = length(loc.SNP)
      
      if(p.SNP==0){## if there is no correspondance between the three signals
        xx = cbind(s.Position.EXP[i],e.Position.EXP[i],NA) 
        Position.common = rbind(Position.common,xx)
      }else if( p.SNP!=0){## having just SNP
        xx = cbind(rep(s.Position.EXP[i],p.SNP),rep(e.Position.EXP[i],p.SNP),Position.SNP[loc.SNP]) 
        Position.common = rbind(Position.common,xx)
      }
      
    }
    colnames(Position.common) = c('Expression Start Position','Expression End Position','SNP Position')
    return(Position.common)
  
  
}

###############################################################
######## Interpolating the mus-status for SNPs ################
###############################################################
## Summary: There are many genes for which we dont have a smoothed SNP signal. 
##          Therefore we need to interpolate it.
### Input: Position.common:      matrix indicating the Expression and SNP positions
###        Position.SNP:         SNP position
###        no          :         number of patients in common
###        mu.SNP      :         segmented SNP signal


.SNP.interpolated = function(Position.common,Position.SNP,no,mu.SNP){
  cat("Performing mean SNP interpolation",sep="\n")  
  mu.SNP.interpolated = matrix(NA,dim(Position.common)[1],no)
  Position.EXP.interpolated = apply(Position.common[,1:2],1,mean)
  
  if(sum(is.na(Position.common[,3]))!=0){
     ### checking which positions need to be interpolated
    x       = Position.SNP
    loc.na  = which(is.na(Position.common[,3]))
    x.inter = Position.EXP.interpolated[loc.na]
    
    loc.SNP = rep(NA,length(Position.common[-loc.na,3]))
    for(i in 1:length(loc.SNP)){
      loc.SNP[i] = which(Position.common[-loc.na,3][i]==Position.SNP)
    }
    ## performing the interpolation patient by patient
    
    for(i in 1:no){
      
      
      mu.SNP.interpolated[-loc.na,i]  = mu.SNP[loc.SNP,i]
      y = mu.SNP[,i]
      mu.SNP.interpolated[loc.na,i]  = approx(x,y,xout=x.inter,yright=mu.SNP[dim(mu.SNP)[1],i],yleft=mu.SNP[1,i])$y
      ## we need to define to the code what happens if we need to interpolate a value, for which we have no information
      ##  that is why we define yright(yleft)
      
    }
    
   # colnames(mu.SNP.interpolated) = colnames(mu.SNP)
    
  }else{
    
    for(i in 1:dim(Position.common)[1]){
      loc = which(Position.common[i,3]==Position.SNP)[1]
      mu.SNP.interpolated[i,] = mu.SNP[loc,]
    }
    
  #  colnames(mu.SNP.interpolated) = colnames(mu.SNP)
  }
  mu.SNP.interpolated = sapply(1:dim(mu.SNP.interpolated)[2],function(i,mu.SNP.interpolated,Position.common) tapply(mu.SNP.interpolated[,i],INDEX = as.factor(Position.common),mean),mu.SNP.interpolated = mu.SNP.interpolated, Position.common = Position.EXP.interpolated)
  
  
  
  return(mu.SNP.interpolated)
}


#########################################################################
####################### Regression EXP vs mus ###########################
#########################################################################
### Summary: Performing the regression
###          but first we need to check the values of the SNPs
### Input: Position.common:      matrix indicating the Expression and SNP positions
###        Position.SNP:         SNP position
###        no          :         number of patients in common
###        mu.SNP.interpolated:  segmented SNP signal (interpolated)
###        data.EXP:             Expression signal


.corrected.EXP  = function(Position.common, s.Position.EXP,no, mu.SNP.interpolated, data.EXP){
  
  n.gene = length(s.Position.EXP)
  pvalues = matrix(NA,n.gene,2)
  rgene = c()
  coef.matrix = matrix(NA,n.gene,2)
  residuals.matrix.diff= matrix(NA,n.gene,no)
  
  cat("Performing CNV correction for gene expression",sep="\n")
  for (i in 1:n.gene){
     
    Data = matrix(NA,no,2) 
    Data[,2] = as.numeric(mu.SNP.interpolated[i,])
    
    Data[,1] = as.numeric(data.EXP[i,])
    Data = data.frame(Data)
    
    names(Data)=c("EXP","SNP")
    if(is.na(Data[1,2])==T){
      
      rgene   =c(rgene ,i)
      
    }else{## performing the regression
      
      fit = lm(EXP~SNP,data = Data)
      pvalues[i,] = summary.lm(fit)[[4]][,4] ## keeping the pvalues for the coefficients (beta values of the regression)
      coef.matrix[i,]=fit$coefficients ## beta values
      residuals.matrix.diff[i,]=residuals(fit)## residuals (corrected signal) ## problem with the variance
      #       ### that is why we prefer to use one of the following residuals
      #       for(k in 1:no){
      #         residuals.matrix.diff[i,k] = residuals(fit)[k]/sd(residuals(fit)[-k])## this is the one I have been using so far
      #       }
    }
    
  }
  
  list(p.values.coef=pvalues,problematic.genes = rgene,coef.matrix = coef.matrix, residuals = residuals.matrix.diff )
}



