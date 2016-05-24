#' @title Multigroup Partial Least Squares Regression
#' 
#' @description 
#' Multigroup PLS regression
#' 
#' @param DataX a numeric matrix or data frame associated with independent dataset
#' @param DataY a numeric matrix or data frame associated with dependent dataset
#' @param Group a vector of factors associated with group structure
#' @param ncomp number of components, if NULL number of components is equal to 2
#' @param Scale scaling variables, by defalt is FALSE. By default data are centered within groups
#' @param Gcenter global variables centering, by defalt is FALSE.
#' @param Gscale global variables scaling, by defalt is FALSE.
#' @return list with the following results:
#' @return \item{DataXm}{Group X  data}
#' @return \item{DataYm}{Group Y  data}
#' @return \item{Concat.X}{Concatenated X data}
#' @return \item{Concat.Y}{Concatenated Y data}
#' @return \item{coefficients}{Coefficients associated with  X data}
#' @return \item{coefficients.Y}{Coefficients associated with regressing  Y on Global components X}
#' @return \item{Components.Global}{Conctenated Components for X and Y}
#' @return \item{Components.Group}{Components associated with groups in X}
#' @return \item{loadings.common}{Common vector of loadings for X and Y}
#' @return \item{loadings.Group}{Group vector of loadings for X and Y}
#' @return \item{expvar}{Explained variance associated with global components X}
#' @return \item{cum.expvar.Group}{Cumulative explained varaince in groups of X and Y}
#' @return \item{Similarity.Common.Group.load}{Cumulative similarity between group and common loadings}
#' @return \item{Similarity.noncum.Common.Group.load}{ NonCumulative  similarity between group and common loadings}
#' @return \item{eigenValue}{EigenValue to calculate the percentage of inertia}
#' @return \item{Percentage.inertia.per.GroupX}{Percentages of inertia per group X}
#' @seealso \code{\link{mgPCA}}, \code{\link{mbmgPCA}}
#' @export
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2013). Multi-group PLS
#' regressMathematics and Statistics, Springer Proceedings (ed), \emph{New Perspectives
#' in Partial Least Squares and Related Methods}, 56, 243-255.
#' 
#' @references A. Eslami, E. M. Qannari, A. Kohler and S. Bougeard (2014). Algorithms
#'  for multi-group PLS. \emph{Journal of Chemometrics}, 28(3), 192-201.
#' 
#' @examples
#' data(oliveoil)
#' DataX = oliveoil[,2:6]
#' DataY = oliveoil[,7:12]
#' Group = as.factor(oliveoil[,1])
#' res.mgPLS = mgPLS (DataX, DataY, Group)
#' barplot(res.mgPLS$noncumper.inertiglobal)
#' #----- Regression coefficients 
#' #res.mgPLS$coefficients[[2]]
#' #----- Similarity index: group loadings are compared to the common structure (in  X and Y spaces)
#' XX1= res.mgPLS$Similarity.noncum.Common.Group.load$X[[1]][-1, 1, drop=FALSE]
#' XX2=res.mgPLS$Similarity.noncum.Common.Group.load$X[[2]][-1, 1, drop=FALSE]
#' simX <- cbind(XX1, XX2)
#' YY1=res.mgPLS$Similarity.noncum.Common.Group.load$Y[[1]][-1, 1, drop=FALSE]
#' YY2=res.mgPLS$Similarity.noncum.Common.Group.load$Y[[2]][-1, 1, drop=FALSE]
#' simY <- cbind(YY1,YY2)
#' XLAB = paste("Dim1, %",res.mgPLS$noncumper.inertiglobal[1])
#' YLAB = paste("Dim1, %",res.mgPLS$noncumper.inertiglobal[2])
#' plot(simX[, 1], simX[, 2], pch=15, xlim=c(0, 1), ylim=c(0, 1), main="Similarity indices in X space", 
#'      xlab=XLAB, ylab=YLAB)
#' abline(h=seq(0, 1, by=0.2), col="black", lty=3)
#' text(simX[, 1], simX[, 2], labels=rownames(simX), pos=2)
#' plot(simY[, 1], simY[, 2], pch=15, xlim=c(0, 1), ylim=c(0, 1), main="Similarity indices in Y space",
#'      xlab=XLAB, ylab=YLAB)
#' abline(h=seq(0, 1, by=0.2), col="black", lty=3)
#' text(simY[, 1], simY[, 2], labels=rownames(simY), pos=2)
mgPLS = function(DataX, DataY, Group, ncomp=NULL, Scale=FALSE, Gcenter=FALSE, Gscale=FALSE){
  #=========================================================================
  #                               Checking the inputs
  #=========================================================================
  check(DataX, Group)
  check(DataY, Group)
  
  if ((nrow(DataX) != nrow(DataY))) 
    stop("Oops unequal number of rows in 'DataX' and 'DataY'.")
  #=========================================================================
  #                               Preparing Data
  #=========================================================================
  if (class(DataX) == 'data.frame') {
    DataX = as.matrix(DataX)
  }
  if (class(DataY) == 'data.frame') {
    DataY = as.matrix(DataY)
  }
  if(is.null(ncomp)) {ncomp = 2}
  if(is.null(colnames(DataX))) {
    colnames(DataX) = paste('VX', 1:ncol(DataX), sep='')
  }
  if(is.null(colnames(DataY))) {
    colnames(DataY) = paste('VY', 1:ncol(DataY), sep='')
  }
  #=========================================================================
  #                        Notations and Presentation of Data
  #=========================================================================
  DataX = as.data.frame(DataX, row.names = NULL)
  DataY = as.data.frame(DataY, row.names = NULL)
  Group = as.factor(Group)
  DataX = as.matrix(DataX)     	 #---- data as matrix form
  rownames(DataX) = Group        #---- rownames of data=groups
  DataY = as.matrix(DataY) 			 #---- data as matrix form
  rownames(DataY) = Group        #---- rownames of data=groups
  M = length(levels(Group))      #---- number of groups: M
  P = dim(DataX)[2]              #---- number of variables: X
  Q = dim(DataY)[2]              #---- number of variables: Y
  n = as.vector(table(Group))    #---- number of individuals in each group
  N = sum(n)                     #---- number of individuals
  H = ncomp
  #=========================================================================
  #					                       Pre-processing 
  #=========================================================================
  #------------------------ Global cenetring and scaling -------------------
  DataX = scale(DataX, center = Gcenter, scale = Gscale) 
  DataY = scale(DataY, center = Gcenter, scale = Gscale)
  #-----------------------------------  Split Data and save them in function
  DataXm = split(DataX, Group)         #----split Data to M parts 
  DataYm = split(DataY, Group)         #----split Data to M parts 
  Concat.X = Concat.Y = NULL
  
  # Centering and scaling if TRUE
  for(m in 1:M){  
    DataXm[[m]] = matrix(DataXm[[m]], nrow=n[m])
    DataXm[[m]] = scale(DataXm[[m]], center=TRUE, scale=Scale)
    
    DataYm[[m]] = matrix(DataYm[[m]], nrow=n[m])
    DataYm[[m]] = scale(DataYm[[m]], center=TRUE, scale=Scale)
    
    Concat.X = rbind(Concat.X, DataXm[[m]])
    Concat.Y = rbind(Concat.Y, DataYm[[m]])
  }
  colnames(Concat.X) = colnames(DataX)
  colnames(Concat.Y) = colnames(DataY)
  rownames(Concat.X) = rownames(DataX)
  rownames(Concat.Y) = rownames(DataY)
  #=========================================================================
  #					                         Output
  #=========================================================================
  res =  list()
  res$DataXm =  DataXm
  res$DataYm =  DataYm
  res$Concat.X = Concat.X
  res$Concat.Y = Concat.Y
  #------------------------ 
  aCommonX = matrix(0, ncol=H, nrow=P)
  TGlobalX = matrix(0, ncol=H, nrow=N)
  UGlobalY = matrix(0, ncol=H, nrow=N)
  bCommonY = matrix(0, ncol=H, nrow=Q)
  aGroupX = vector("list",M)
  bGroupY = vector("list",M)
  TGroupX = vector("list",M)
  UGroupY = vector("list",M)
  for(m in 1:M){
    aGroupX[[m]] = matrix(0,nrow=P, ncol=H)
    bGroupY[[m]] = matrix(0,nrow=Q, ncol=H)
    rownames(aGroupX[[m]]) = colnames(DataX)
    colnames(aGroupX[[m]]) = paste("Dim", 1:H, sep="")
    rownames(bGroupY[[m]]) = colnames(DataY)
    colnames(bGroupY[[m]]) = paste("Dim", 1:H, sep="")
    
    
    TGroupX[[m]] = matrix(0,nrow=n[m],ncol=H)
    UGroupY[[m]] = matrix(0,nrow=n[m],ncol=H)
  }
  
  

  #------------------------ 
  aCommonX_STAR = matrix(0,nrow=P,ncol=H)
  COV  = vector("list",H)
  COVm = vector("list",H)   # sum cov2(um, tm)
  #------------------------ 
  ppbeta = matrix(0,nrow=P,ncol=H)
  betay  = vector("list",H)
  #------------------------ 
  cum.expvar.Xm = matrix(0, ncol=H, nrow=M)
  cum.expvar.Ym = matrix(0, ncol=H, nrow=M)
  #------------------------ to calculate norm 
  nor2x = sum((Concat.X)^2)
  exp_var = c(0)
  nor2y= sum((Concat.Y)^2)
  exp_vary = c(0)
  #------------------------
  res$noncumper.inertiglobal = matrix(0, ncol=ncomp, nrow=1)
  colnames(res$noncumper.inertiglobal) = paste("Dim", 1:ncomp, sep="")
  #---------------
  #eigen value
  res$eigenValue = NULL
  #res$WeigenValue = matrix(0,ncol=ncomp,nrow=ncol(DataX))
  #---------- Percentage.inertia.per.GroupX
  res$Percentage.inertia.per.GroupX = matrix(0, nrow=M, ncol= ncomp)
  rownames(res$Percentage.inertia.per.GroupX) = levels(Group)
  colnames(res$Percentage.inertia.per.GroupX) = paste("Dim", 1:ncomp, sep="")
  
  #=========================================================================
  #                      NIPALS algorithm for multigroup PLS data
  #=========================================================================
  for(h in 1:H){
    eps = 1e-6     # for iteration loops
    #w= svd(t(Concat.X) %*% Concat.Y) $u[,1] # initial value of w
    w = rnorm(P)
    w = w/as.numeric(sqrt(t(w)%*%w))  # normalize w
    #-------------------------------
    x=1.0;
    tTold=0;
    wTold=w;
    uTold=0;
    vTold=0;
    iter =0;
    covv  = 0   # criterion for each iteration
    covvm = 0   # criterion for each iteration
    tt = matrix(0,nrow=N)
    tmm = list()
    umm = list()
    #-------------------------------   
    while (x>eps){
      # ------- compute loading/ components vector associated to Y -----------
      # loading associated to Y
      sumyx=0
      
      for (m in 1:M){ 
        ccc = t(DataYm[[m]]) %*% DataXm[[m]] %*% w
        sumyx = sumyx+ccc
      }
      v= sumyx /  as.numeric(sqrt(sum((sumyx)^2)))
      #---------------------------- Y dataset      
      # component associated to Y
      # and the group component umm
      for (m in 1:M){    
        um = DataYm[[m]] %*% v        # score of group m
        umm[[m]]=um
      }
      u    = (Concat.Y) %*% v  
      # ------------ components and loading vectors asosciated to X -----------
      # update the loading vector associates to X
      #---------------------------- w
      sumxy=0
      for (m in 1:M){ 
        ccc = t(DataXm[[m]]) %*% DataYm[[m]] %*% v
        aGroupX[[m]][,h]=ccc/as.numeric(sqrt(sum((ccc)^2)))
        
        sumxy = sumxy+ccc
      }
      w = sumxy / as.numeric(sqrt(sum((sumxy)^2)))
      # -----  compute t
      # global component t associated to the X data set
      # and the group component tm
      for (m in 1:M){ 
        tm = DataXm[[m]] %*% w       # score of group m
        TGroupX[[m]][,h] = tm
        ddd = t(DataYm[[m]]) %*% tm
        bGroupY[[m]][,h] = ddd/as.numeric(sqrt(sum((ddd)^2)))
        
        tmm[[m]] = tm
      }
      t = (Concat.X) %*% w 
      tt = cbind(tt,t)
      #------------# test accuracy (on the normalized) T and P 
      #-------------check convergence -------------- 
      iter=iter+1
      if (iter > 1){
        xT = crossprod(t-tTold) 
        yw = crossprod(v-vTold)  
        x  = max(xT,yw);
      }
      covv[iter]=(cov(u, t))
      tTold= t;
      wTold= w;
      uTold= u;
      vTold= v;
      #---------------------------- #
      sum_cov=0
      for(m in 1:M){
        sum_cov= sum_cov+( (n[m]-1) * cov(tmm[[m]],umm[[m]]) )
      }
      covvm[iter] = sum_cov
      #---------------------------- #
    }# end of iteration
   
    #--------------------------------
    #eigen value
    VVeigen=0   #--- v=sum (Xm' * Ym)
    for(m in 1: M){
      VVeigen = VVeigen + (t(DataXm[[m]])%*% DataYm[[m]])
    }
    WWeigen=VVeigen%*%t(VVeigen)
    #--- 2. common loadings
   ## res$WeigenValue[,h] =  eigen(WWeigen)$vectors[,1]
    res$eigenValue[h]   =  eigen(WWeigen)$values[1]/(N-1)^2 #eigen(WWeigen)$values[1]/(N-1)^2 
    #---------------------------- #   
    aCommonX[,h] = wTold
    TGlobalX[,h] = t 
    UGlobalY[,h] = u
    bCommonY[,h] = v
    COV[[h]] = covv
    COVm[[h]]=covvm
    project.global=TGlobalX[,h] %*% ginv( t(TGlobalX[,h]) %*% TGlobalX[,h]) %*%   t(TGlobalX[,h])
    project.data.global=project.global %*% res$Concat.X 
    res$noncumper.inertiglobal[h] = round(100*sum(diag(t(project.data.global) %*% project.data.global))/ sum(diag(t(res$Concat.X) %*% res$Concat.X)),1)
    #=========================================================================
    #                      Computation of regression coefficients 
    #=========================================================================
    ppbeta[,h] = t(res$Concat.X)%*% TGlobalX[,h]/(c(( t(TGlobalX[,h])%*% TGlobalX[,h])))
    if(h==1){aCommonX_STAR[, h]=aCommonX[,h]}
    if(h>=2){
      Mat                <- diag(P) 
      for (k in 2 : h) {
        mat                <- t(TGlobalX[, k-1]) %*% res$Concat.X/ as.numeric(t(TGlobalX[, k-1])%*%TGlobalX[, k-1])
        Mat                <- Mat %*% (diag(P) - aCommonX[, k-1] %*% mat)
      }
      aCommonX_STAR[, h] <- Mat %*% aCommonX[, h]
    }
    dd = coefficients(lm(res$Concat.Y ~TGlobalX[,1:h]))
    if(Q==1){betay[[h]]=(dd[-1])}
    if(Q>=2){betay[[h]]=(dd[-1,])}
    #res$BETA[[h]]=as.matrix(aCommonX_STAR[,1:h]) %*% (betay[[h]])
    res$coefficients[[h]] = as.matrix(aCommonX_STAR[,1:h]) %*% (betay[[h]])
    colnames(res$coefficients[[h]]) = colnames(Concat.Y)
    rownames(res$coefficients[[h]]) = colnames(Concat.X)
    #=========================================================================
    #                                   Deflation 
    #=========================================================================
    wp = crossprod(Concat.X, t) / drop(crossprod(t))
    wy = crossprod(Concat.Y, t) / drop(crossprod(t))
    Concat.X = Concat.X -  t%*% t(wp)
    Concat.Y = Concat.Y -  t %*% t(wy)
    
    #------------------------ explained variance 
    exp_var_new= 100* as.numeric(t(t) %*% res$Concat.X %*% t(res$Concat.X) %*% t/ as.vector((t(t)%*%t)) )/nor2x
    exp_var=append(exp_var, exp_var_new)
  
    exp_var_newy= 100* as.numeric(t(t) %*% res$Concat.Y %*% t(res$Concat.Y) %*% t/ as.vector((t(t)%*%t)) )/nor2y
    exp_vary=append(exp_vary, exp_var_newy)
    #---------------------------------- #Split Data
    DataXm = split(Concat.X, Group)  
    DataYm = split(Concat.Y, Group)
    
    for(m in 1:M){  # dataset in form of matrix
      DataXm[[m]] = matrix(DataXm[[m]], ncol=P)
      colnames(DataXm[[m]]) = colnames(Concat.X)
      
      DataYm[[m]] = matrix(DataYm[[m]], ncol=Q)
      colnames(DataYm[[m]]) = colnames(Concat.Y)
    }
    #-----------------------
  }  # END OF DIMENSION
  #=========================================================================
  #                                  Saving outputs
  #=========================================================================
  expvarX = exp_var[-1]
  cumexpvarX = cumsum(expvarX)

  EVX = matrix(c(expvarX, cumexpvarX), ncol=2)
  rownames(EVX) = paste("Dim", 1:H, sep="")
  colnames(EVX) = c("Explained.Var", "Cumulative")  
  
  expvarY = exp_vary[-1]
  cumexpvarY = cumsum(expvarY)
  
  EVY = matrix(c(expvarY, cumexpvarY), ncol=2)
  rownames(EVY) = paste("Dim", 1:H, sep="")
  colnames(EVY) = c("Explained.Var", "Cumulative") 
  #------------------------
  for(m in 1:M){ 
    for(h in 1:H){ 
      
      proj=  TGroupX[[m]][,1:h] %*%  ginv(t(TGroupX[[m]][,1:h]) %*% TGroupX[[m]][,1:h])  %*% t(TGroupX[[m]][,1:h]) 
      xm_hat= proj %*%  DataXm[[m]]
      explain_varx=sum(diag( xm_hat %*% t(xm_hat) ))/sum(diag( DataXm[[m]] %*% t(DataXm[[m]]) ))
      ym_hat= proj %*%  DataYm[[m]]
      explain_vary = sum(diag( ym_hat %*% t(ym_hat) ))/sum(diag( DataYm[[m]] %*% t(DataYm[[m]]) ))
      cum.expvar.Xm[m,h] = explain_varx
      cum.expvar.Ym[m,h] = explain_vary
    }
  }
  #------------------------
  # Y coefficients
  res$coefficients.Y <- t(coefficients(lm(Concat.Y ~ TGlobalX - 1)))
  rownames(res$coefficients.Y) = colnames(res$Concat.Y)
  colnames(res$coefficients.Y) = paste("Dim", 1:H, sep="")
  #=========================================================================
  #			       Similarity between group and common loadings 
  #=========================================================================
  loadings_matricesW = list()
  loadings_matricesW[[1]] = aCommonX
  for(m in 2:(M+1)){
    loadings_matricesW[[m]] = aGroupX[[(m-1)]]
  }
  
  loadings_matricesV = list()
  loadings_matricesV[[1]] = bCommonY
  for(m in 2:(M+1)){
    loadings_matricesV[[m]] = bGroupY[[(m-1)]]
  }
  NAMES = c("Commonload", levels(Group))
  similarityA =similarity_function(loadings_matrices=loadings_matricesW, NAMES)
  similarityB =similarity_function(loadings_matrices=loadings_matricesV, NAMES)
  similarityA_noncum = similarity_noncum(loadings_matrices=loadings_matricesW, NAMES)
  similarityB_noncum = similarity_noncum(loadings_matrices=loadings_matricesV, NAMES)
 
  #-------------------------------------------------------------------
  Pct_vartm = matrix(0, nrow=M, ncol= ncomp)
  for (m in 1:M){
    Pct_vartm[m, ] <- apply(TGroupX[[m]], 2, var) / sum(apply(TGroupX[[m]], 2, var))*100
  }
  
  rownames(Pct_vartm) = levels(Group)
  colnames(Pct_vartm) = paste("Dim", 1:H, sep="")
  #------------------------------------------------------------------
  rownames(aCommonX) = colnames(res$Concat.X)
  rownames(bCommonY) = colnames(res$Concat.Y)
  rownames(TGlobalX) = rownames(res$Concat.X)
  rownames(UGlobalY) = rownames(res$Concat.Y)
  colnames(aCommonX) = paste("Dim", 1:H, sep="")
  colnames(bCommonY) = paste("Dim", 1:H, sep="")
  colnames(TGlobalX) = paste("Dim", 1:H, sep="")
  colnames(UGlobalY) = paste("Dim", 1:H, sep="")
  
  for(m in 1: M){
    colnames(TGroupX[[m]]) = paste("Dim", 1:H, sep="")
  }
  res$Components.Global  = list(X = TGlobalX, Y = UGlobalY)
  res$Components.Group   = TGroupX
  res$loadings.common    = list(X = aCommonX, Y = bCommonY)
  res$loadings.Group     = list(X = aGroupX, Y = bGroupY)
  colnames(cum.expvar.Xm) = paste("Dim", 1:H, sep="")
  colnames(cum.expvar.Ym) = paste("Dim", 1:H, sep="")
  rownames(cum.expvar.Xm) = levels(Group)
  rownames(cum.expvar.Ym) = levels(Group)
  
  res$expvar  = list(X = EVX, Y = EVY)
  res$cum.expvar.Group  = list(X = cum.expvar.Xm, Y = cum.expvar.Ym)
  
  res$Similarity.Common.Group.load          = list(X = similarityA, Y = similarityB)
  res$Similarity.noncum.Common.Group.load   = list(X = similarityA_noncum , Y = similarityB_noncum)
  res$Percentage.inertia.per.GroupX = Pct_vartm
  
  #----------------------------------
  class(res) = c("mgPLS", "mg")
  return(res)
}


#' @S3method print mgPLS
print.mgPLS <- function(x, ...)
{
  cat("\nMultigroup Partial Least Squares Regression\n")
  cat(rep("-",43), sep="")
  cat("\n$loadings.common      ", "Common loadings")
  cat("\n")
  invisible(x)
}


