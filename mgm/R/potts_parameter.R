
potts_parameter <- function(
  graph, #(weighted) adjacency matrix ;matrix
  type,  # vector indicating the type of the variables ;vector
  lev, #how many lev do variables have (>1 only for categorical) ;vector
  thresh #threshold for each variable; for each category for categoricals ;list
)
{
  
  ## create empty matrix with correct dimensions
  model.par.matrix <- matrix(0, sum(lev), sum(lev))
  nVar <- nrow(graph)
  
  r <- 1 #current row; we need to specify that extra, as we use two indices (variables and cat) in one loop
  
  for(v in 1:nVar) { #loop variables
    
    # type[v] == CATEGORICAL
    if(type[v]=="c") {
      
      subm.mpm <- matrix(0,lev[v], sum(lev[-v])) #fill a dummy matrix; then include slip in thresholds and merge with full matrix
      
      
      for(k in 1:lev[v]) { #loop categories
        
        c <- 1 #current column index
        
        for(v2 in (1:nVar)[-v]) #loop other variables
        {
          
          if(graph[v,v2]!=0) { #check whether edge present between variables; if not, we just set all parameters to 0 (means: we dont to anything)
            
            #checking type of interaction (type of 2nd variable)
            if(type[v2]=="c") {
              
              #potts as far as categories overlap
              subm.mpm[k,c:(c+lev[v2]-1)] <- (k==(1:lev[v2])) * graph[v,v2] 
              
              #fill last category up
              if(k>lev[v2]) { 
                subm.mpm[k,c:(c+lev[v2]-1)] <- ((max(1:lev[v2]))==(1:lev[v2])) * graph[v,v2] 
                }
              
              
            } else { #same for all continuous, no further specification
              
              subm.mpm[k,c] <-  (k <= lev[v]/2) * graph[v,v2] # smaller "half" has effect, other not; MISSING: multiply actual edge weight
              
            }
          }#end if: edge present?
          
          c <- c + lev[v2]
        }#end: other variables
        
      }#end: categories
      
      #fill in thresholds
      p1 <- sum(lev[0:(v-1)]) #cols before
      if(v+1 > length(lev)) { p2 <- 0 } else {
        p2 <- sum(lev[(v+1):length(lev)]) } #cols after
      
      dummy <- 1 #necessary in the case where the submatrix is at the right edge of the subm.mpm.t
      if((ncol(subm.mpm)-p2+1)>ncol(subm.mpm)) { dummy <- 0 }
      
      subm.mpm.t <- cbind(subm.mpm[,0:p1],
                          matrix(rep(0,lev[v]^2),lev[v],lev[v]),
                          subm.mpm[,((ncol(subm.mpm)-p2+1)*dummy):(ncol(subm.mpm)*dummy)])
      
      
      #merge with
      model.par.matrix[r:(r+lev[v]-1),] <- subm.mpm.t #include weights from (weighted) adj matrix
      
      # type[v] == CONTINUOUS
    } else {
      
      # here, we only need to get the continuous interactions; the rest we get by symmetry
      
      subv.mpm <- numeric(sum(lev))
      c <- 1 #current column index
      
      for(v4 in (1:nVar)){ #loop other variables
        # only continous
        if(type[v4]!="c") { subv.mpm[c] <- graph[v,v4] }
        c <- c + 1
      }#end: other variables
      
      #merge
      model.par.matrix[r,] <- subv.mpm
      
    }
    
    r <- r + lev[v] #move down the target matrix
    
  } #end for variables

  ## "mirror" matrix
  pars <- ncol(model.par.matrix)
  for(rows in 1:pars)
  {
    for(cols in 1:pars) {
      if(model.par.matrix[rows,cols]!=0) {
        model.par.matrix[cols,rows] <- model.par.matrix[rows,cols]
      }
    }
  }

  # insert thresholds in diagonal
  diag(model.par.matrix) <- unlist(thresh)

  # add labels indicating which row/col corresponds to which variable; this helps to feed the parameters in the sampler
  dummy_matrix <- cbind(lev, 1:length(lev))
  ind.e <- unlist(apply(dummy_matrix,1,function(x) { rep(x[2],x[1])}))
  colnames(model.par.matrix) <- rownames(model.par.matrix) <- ind.e

  return(model.par.matrix)
  
  ### end of FUNCTION
}
