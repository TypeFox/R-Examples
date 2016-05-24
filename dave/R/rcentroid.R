rcentroid <-
function(veg,membership.r,y,...){
# computing a matrix of releve group centroids (rows) with all species(columns)
# available as frequency table (species counts) or as probability (relative
# frequency) 21. 3. 2012
  nrel <- length(veg[,1])
  nspec <- length(veg[1,])
  nele <- nrel*nspec                                 # no. of elements in the entire table
  xall <- as.numeric(as.matrix(t(veg),nrow=nele))    # entire matrix in one dimension
  rowgroups<- rep(seq(1,nspec,1),nrel)
  col <- c(rep(membership.r,nspec)) 
  col <- matrix(col,nrow=nrel,ncol=nspec)
  col <- t(col)
  colgroups <- matrix(col,nrow=nele)
  F <- tapply(sign(xall),list(rowgroups, colgroups),sum)      # contingency table, counts
  F<- t(F)
  colnames(F)<- colnames(veg)
# relativized centroid table in P
  tt<-table(membership.r)
  S<-matrix(rep(tt,nspec),nrow=length(tt))
  P<- F/S
  colnames(P)<- colnames(veg)
 
  dgr<- (1 - cor(t(P^y)))/2
  dgr<- as.dist(dgr,upper=T,diag=T)

# list output
  rcentroid<- list(nrelgroups=length(tt),nspec=nspec,freq.table=F,prob.table=P,dist.mat=dgr)
  }
