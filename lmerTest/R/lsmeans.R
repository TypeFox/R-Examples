################################################################################
## LSMEANS and DIFFLSMEANS related functions
################################################################################



################################################################################
## functions for popMatrix for LSMEANS (from doBy package)
################################################################################
.get_xlevels <- function(obj){
  UseMethod(".get_xlevels")
}

.get_xlevels.default <- function(obj){
  obj$xlevels
}


.covariateAve <- function(object, at=NULL, tt=terms(object)){
  tt  <- delete.response(tt)
  att <- attributes(tt)
  rhs.terms <- rownames(att$factors)[rowSums(att$factors)>0]
  rhs.class <- att$dataClass[match(rhs.terms, names(att$dataClass))]
  nums      <- rhs.terms[rhs.class=="numeric"]
  
  ans  <- lapply(model.frame(object)[,nums, drop=FALSE], mean) 
  
  nn <- match(names(ans), names(at))
  nn <- nn[!is.na(nn)]
  at.num <- at[nn]
  ans[names(at[nn])] <- at.num
  attr(ans, "at.num") <- at.num
  ans
}


.get_vartypes <- function(object){
  tt <- terms(object)
  tt  <- delete.response(tt)
  att <- attributes(tt)
  rhs.terms <- rownames(att$factors)[rowSums(att$factors)>0]
  rhs.class <- att$dataClass[match(rhs.terms, names(att$dataClass))]
  nums      <- rhs.terms[rhs.class=="numeric"]
  fact      <- rhs.terms[rhs.class=="factor"]
  list(numeric=nums, factor=fact)
}


.set_xlevels <- function(xlev, at){
  nam    <- names(xlev)
  nn <- match(nam, names(at))
  nn <- nn[!is.na(nn)]
  at.fact <- at[nn]
  xlev[names(at[nn])]  <- at.fact
  attr(xlev, "at.fact") <- at.fact
  xlev
}



.getX <- function(object, newdata){
  tt <- terms(object)
  Terms  <- delete.response(tt)
  mf  <- model.frame(Terms, newdata, xlev = .get_xlevels(object))
  X   <- model.matrix(Terms, mf, contrasts.arg = .get_contrasts(object))
  attr(X,"assign")<-NULL
  attr(X, "contrasts") <- NULL
  X
}


.get_contrasts <- function(obj){
  UseMethod(".get_contrasts")
}

.get_contrasts.default <- function(obj){
  obj$contrasts
}


popMatrix <- function(object, effect=NULL, at=NULL, only.at=TRUE){
  tt <- terms(object)
  Terms   <- delete.response(tt)
  xlev    <- .get_xlevels(object)
  ccc     <- .covariateAve(object,at)
  vartype <- .get_vartypes(object)
  
  
  ##   cat("INPUT: effect:\n"); str(effect)
  ##   cat("INPUT: at:\n"); str(at)
  ##   cat("---------------------------\n")
  xlev   <- .get_xlevels(object)
  
  if (is.null(effect)){
    at.factor <- at[intersect(vartype$factor, names(at))]
    xxx       <- if(length(at.factor)>0)
      at.factor
  } else {
    xlev   <- .set_xlevels(xlev, at=at)
    at.fact <- names(attr(xlev, "at.fact"))
    effect <- setdiff(effect, at.fact)
    xxx    <- xlev[c(effect,at.fact)]
  }
  
  #  print(ccc)
  #  print(xxx)
  
  
  #print(xxx)
  if (is.null(xxx)){
    ## No 'effect' and no 'at'; just to a global average.
    newdata <- expand.grid(xlev)
    newdata[,names(ccc)] <- ccc   
    mf  <- model.frame(Terms, newdata, xlev = .get_xlevels(object))
    X   <- model.matrix(Terms, mf, contrasts.arg = .get_contrasts(object))
    res <- apply(X,2,mean)
    res <- do.call(rbind, list(res))
    attr(res,"at") <- at[intersect(vartype$numeric, names(at))]
  } else {
    eff.grid  <- expand.grid(xxx)
    eff.grid  <- as.data.frame(lapply(eff.grid, as.character),stringsAsFactors=FALSE)
    #cat("eff.grid:\n"); print(eff.grid)
    res <- list()
    for (ii in 1:nrow(eff.grid)){
      conf  <- eff.grid[ii,,drop=FALSE]
      xlev2 <- .set_xlevels(xlev,  at=conf)
      #cat("xlev2 (which defines the grid):\n"); str(xlev2)
      newdata <- expand.grid(xlev2)
      newdata[,names(ccc)] <- ccc   
      
      #print(newdata)
      mm   <- .getX(object, newdata)
      X    <- apply(mm,2,mean)
      res[[ii]] <- X
    }
    
    res <- do.call(rbind, res)
    #    print(eff.grid)
    uuu <- at[intersect(vartype$numeric, names(at))]
    #    print(uuu)
    #    print(vartype)
    #    print(at)
    #    print(ccc)
    #eff.grid[,names(ccc)] <- at[intersect(vartype$numeric, names(at))]
    eff.grid[,names(ccc)] <- ccc
    attr(res,"grid") <- eff.grid
    attr(res,"at") <- at
  }
  class(res) <- c("popMatrix", "conMatrix", "matrix")
  res 
}

################################################################################
################################################################################


###################################################################
#get the combinatoion of the fixed factors for the lsmeans
###################################################################
getFacCombForLSMEANS <- function(split.eff, data)
{
  if(length(split.eff)==1)
    data.merge <- as.data.frame(levels(data[,split.eff]))
  if(length(split.eff)>=2)
    data.merge <- merge(levels(data[,split.eff[1]]),levels(data[,split.eff[2]]))
  if(length(split.eff)>=3)
  {
    for(i in 3:length(split.eff))
    {
      d.split.eff_i <- as.data.frame(levels(data[,split.eff[i]]))
      names(d.split.eff_i) <- paste("l",i)
      data.merge <- merge(data.merge,d.split.eff_i)
    }
    
  }
  names(data.merge) <- split.eff
  return(as.matrix(data.merge))
}

###################################################################
#checks if all the terms in interaction are covariates
###################################################################
checkAllCov <- function(split.eff, data)
{
  for(spleff in split.eff)
  {
    if(!is.factor(data[,spleff]))
    {
      return(TRUE)  
    }
  }
  return(FALSE)
}

###################################################################
#concatenate levels of the effects to form the rownames
###################################################################
concatLevs <- function(matr, row.names)
{
  
  if(is.vector(matr))
    levs <- paste(names(matr),matr)
  else
  {
    levs <- paste(rownames(matr),matr[,1])
    for(i in 2:ncol(matr))
    {
      levs <- paste(levs,matr[,i])
    }    
  }
  
  
  return(levs)
}

#convert facs into numeric
convertFacsToNum <- function(data, begin, end)
{
  
  #convert vars to numeric
  for(i in begin:end)
    data[,i] <- as.numeric(levels(data[,i])[as.integer(data[,i])]) 
  
  return(data)
}


#convert numeric to facs
convertNumsToFac <- function(data, begin, end)
{
  
  #convert vars to numeric
  for(i in begin:end)
    data[,i] <- as.factor(data[,i]) 
  
  return(data)
}


###################################################################
#fills the LSMEANS and DIFF summary matrices
###################################################################
fillLSMEANStab <- function(mat, rho, summ.eff, nfacs, alpha)
{

  newcln <- colnames(mat)[colnames(mat) %in% names(rho$fixEffs)]
  mat <- matrix(mat[,colnames(mat) %in% names(rho$fixEffs)], nrow=nrow(mat), ncol=length(newcln), dimnames=list(rownames(mat),  newcln))
  estim.lsmeans <- mat %*% rho$fixEffs
  summ.eff[,nfacs+1] <- estim.lsmeans  
  ttest.res <- calculateTtestJSS(rho, t(mat), nrow(mat))

  summ.eff[,nfacs+2] <- ttest.res[,4]#stdErrLSMEANS(rho, std.rand, mat)
  #df
  summ.eff[,(nfacs+3)] <- ttest.res[,1]
  #t values
  summ.eff[,(nfacs+4)] <- ttest.res[,2]
  #p values
  summ.eff[,(nfacs+7)] <- ttest.res[,3]
  # CIs
  summ.eff[,nfacs+5] <- estim.lsmeans-abs(qt(alpha/2,ttest.res[,1]))*ttest.res[,4]
  summ.eff[,nfacs+6] <- estim.lsmeans+abs(qt(alpha/2,ttest.res[,1]))*ttest.res[,4]
  return(summ.eff)
}

###################################################################
#round the columns of LSMEANS or DIFFS tables
###################################################################
roundLSMEANStab <- function(summ.eff, nfacs)
{
  summ.eff[,nfacs+1] <- round(summ.eff[,nfacs+1],4)  
  summ.eff[,nfacs+2] <- round(summ.eff[,nfacs+2],4)
  #df
  summ.eff[,(nfacs+3)] <- round(summ.eff[,(nfacs+3)],1)
  #t values
  summ.eff[,(nfacs+4)] <- round(summ.eff[,(nfacs+4)],2)
  #p values
  summ.eff[,(nfacs+7)] <- round(summ.eff[,(nfacs+7)],4)
  # CIs
  summ.eff[,nfacs+5] <- round(summ.eff[,nfacs+5],4)
  summ.eff[,nfacs+6] <- round(summ.eff[,nfacs+6],4)
  return(summ.eff)
}

############################################################################
#function to identify the colors of bar according to significance of effects
############################################################################
calc.cols <- function(x)
{
  if(x<0.001) 
    return("red") 
  if(x<0.01) 
    return("orange") 
  if(x<0.05) 
    return("yellow") 
  return("grey")
}

calc.cols2 <- function(x)
{
  if(x<0.001) 
    return("p-value < 0.001")#return("red")# 
  if(x<0.01) 
    return("p-value < 0.01")#return("orange")# 
  if(x<0.05) 
    return("p-value < 0.05")#return("yellow")# 
  return("NS")#return("grey")#
}

#get names for ploting barplots for the effects
getNamesForPlot <- function(names, ind)
{
  namesForPlot <- unlist(lapply(names, 
                                function(y) 
                                  substring2(y, 1, 
                                             substring.location(y, " ")$first[1]-1)))
  namesForLevels <- unlist(lapply(names,
                                  function(y)  
                                    substring2(y, 
                                               substring.location(y, " ")$first[1]+ind, nchar(y))))
  return(list(namesForPlot=namesForPlot, namesForLevels=namesForLevels))
}

#plots for LSMEANS or DIFF of LSMEANS
plotLSMEANS <- function(table, response, 
                        which.plot=c("LSMEANS", "DIFF of LSMEANS"), 
                        main = NULL, cex = 1.4, effs = NULL, mult = TRUE)
{
  
  if(!is.null(effs)){
    rnames <- rownames(table)
    diffs.facs <- sapply(rnames, 
                         function(x) substring(x, 1, 
                                               substring.location(x, " ")$first[1]-1), 
                         USE.NAMES = FALSE)    
    find.fac <- diffs.facs %in% effs
    table <- table[find.fac,]
  }
  
  if(which.plot=="LSMEANS")
    names <- getNamesForPlot(rownames(table),2)
  else
    names <- getNamesForPlot(rownames(table),1)   
  
  
    
  namesForPlot <- names$namesForPlot
  namesForLevels <- names$namesForLevels
  un.names <- unique(namesForPlot)
  
  
  ### changed code to transfer to ggplot
  ttplot <- table
  ttplot$namesforplots <- namesForPlot
  ttplot$levels <- as.factor(namesForLevels)
  colnames(ttplot)[which(colnames(ttplot)=="p-value")] <- "pvalue"
  colnames(ttplot)[which(colnames(ttplot)=="Lower CI")] <- "lci"
  colnames(ttplot)[which(colnames(ttplot)=="Upper CI")] <- "uci"
  ttplot$col.bars <-  unlist(lapply(ttplot[,"pvalue"], calc.cols2))
  ttplot <- ttplot[,c("levels", "Estimate", "col.bars", "lci", "uci", 
                      "namesforplots")]
  uci <- lci <- col.bars <- Estimate <- NULL
  
  if(mult)
    ggplot(ttplot, aes(x=levels, y = Estimate, fill = col.bars)) + 
    geom_bar(position = "dodge", stat = "identity") +  
    geom_errorbar(aes(ymin = lci, ymax = uci ), colour="black", width=.1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4), 
          axis.title.y = element_text(size = rel(1.4)), 
          axis.text = element_text(size = rel(1)), 
          legend.text = element_text(size = rel(1)), 
          legend.title = element_text(size = rel(1)))  + 
    scale_fill_manual(values  = 
                        c(  "NS" = "grey", "p-value < 0.01" = "orange", 
                            "p-value < 0.05" = "yellow", 
                            "p-value < 0.001" = "red"), name="Significance")  +
    ylab(response) + facet_wrap( ~ namesforplots, scales = "free")
  else{
    for(i in 1:length(un.names)){
      names.plot <- un.names[i]
      subplot <- ttplot[ttplot$namesforplots == names.plot,]
      ggplot(subplot, aes(x=levels, y = Estimate, fill = col.bars)) + 
        geom_bar(position = "dodge", stat = "identity") +  
        geom_errorbar(aes(ymin = lci, ymax = uci ), colour="black", width=.1) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4), 
              axis.title.y = element_text(size = rel(1.4)), 
              axis.text = element_text(size = rel(1)), 
              legend.text = element_text(size = rel(1)), 
              legend.title = element_text(size = rel(1)))  + 
        scale_fill_manual(values  = 
                            c(  "NS" = "grey", "p-value < 0.01" = "orange", 
                                "p-value < 0.05" = "yellow", 
                                "p-value < 0.001" = "red"), name="Significance")  +
        ylab(response) + xlab(names.plot)
    }
  }
}



#calculate DIFFERENCES OF LSMEANS and STDERR for effect
calcDiffsForEff <- function(facs, fac.comb, split.eff, eff, effs, data, rho, 
                            alpha, mat)
{
  ###calculating diffs for 2 way interaction
  if(length(split.eff)>=1 && length(split.eff)<=2)
  {   
    if(length(split.eff)==2)
    {            
      fac.comb.names <- concatLevs(fac.comb)
      main.eff <- effs[effs %in% split.eff]
      
      mat.names.diffs <- combn(fac.comb.names,2)
      mat.nums.diffs <- apply(mat.names.diffs, c(1,2), 
                              function(x) which(fac.comb.names==x))
           
    }
    else
    {
      mat.names.diffs <- combn(fac.comb,2)
      mat.nums.diffs <- apply(mat.names.diffs, c(1,2), 
                              function(x) which(fac.comb==x))
    }
    
    mat.diffs <- matrix(0, nrow=ncol(mat.nums.diffs), ncol=ncol(mat))
    colnames(mat.diffs) <- colnames(mat)
    for(ind.diffs in 1:ncol(mat.nums.diffs))
    {
      mat.diffs[ind.diffs,] <- mat[mat.nums.diffs[1,ind.diffs],]- mat[mat.nums.diffs[2,ind.diffs],]
    }
    names.combn <- apply(mat.names.diffs, 2, 
                         function(x) paste(x[1], x[2], sep=" - "))
    rownames(mat.diffs) <- paste(eff, names.combn)
    
    diffs.summ <-  matrix(NA, ncol=7, nrow=nrow(mat.diffs))
    colnames(diffs.summ) <- c("Estimate","Standard Error", "DF", "t-value", 
                              "Lower CI", "Upper CI", "p-value")
    rownames(diffs.summ) <- rownames(mat.diffs)
    
    diffs.summ <- as.data.frame(fillLSMEANStab(mat.diffs, rho, diffs.summ, 0, 
                                               alpha))
    return(roundLSMEANStab(diffs.summ, 0))
  }
  
}



#calculate LSMEANS and STDERR for effect
calcLsmeansForEff <- function(lsmeans.summ, fac.comb, eff, split.eff, alpha, mat, 
                              rho, facs)
{
  
  summ.eff <- matrix(NA, ncol=ncol(lsmeans.summ), nrow=nrow(fac.comb))
  colnames(summ.eff) <- colnames(lsmeans.summ)
  #rownames(summ.eff) <- rep(eff, nrow(fac.comb))
  summ.eff[,split.eff] <- fac.comb
  names.arg <- concatLevs(summ.eff[,split.eff])
  summ.eff <- as.data.frame(fillLSMEANStab(mat, rho, summ.eff, length(facs), 
                                           alpha))
  summ.eff <- convertFacsToNum(summ.eff, length(facs)+1, ncol(summ.eff))
  
  summ.eff <- roundLSMEANStab(summ.eff, length(facs))
  
   
  rownames(summ.eff) <- paste(rep(eff, nrow(fac.comb)), names.arg)
  return(summ.eff) 
}



###################################################################
#calculate LSMEANS DIFFS and CI for all effects
###################################################################
calcLSMEANS <- function(model, data, rho, alpha, test.effs = NULL, 
                        lsmeansORdiff = TRUE, 
                        l.lmerTest.private.contrast)
{  
  
 
  #####################################
  m <- refitLM(model, l.lmerTest.private.contrast)
 
  effs <- attr(terms(m),"term.labels")
  if(!is.null(test.effs))
    effs <- effs[effs %in% test.effs]
  dclass <- attr(terms(m),"dataClasses")
  facs <- names(dclass[which(dclass=="factor")])
  #Get standard deviation of random parameters from model
  std.rand <- c(unlist(lapply(VarCorr(model), function(x) attr(x,"stddev"))), 
                attr(VarCorr(model), "sc"))^2 #as.numeric(rho$s@REmat[,3])
  
  
  #init lsmeans summary
  if(lsmeansORdiff)
  {
    lsmeans.summ <-  matrix(ncol=length(facs)+7,nrow=0)
    colnames(lsmeans.summ) <- c(facs,"Estimate","Standard Error", "DF", "t-value",
                                "Lower CI", "Upper CI", "p-value")
    summ.data <- as.data.frame(lsmeans.summ)
  }
  else
  {
    #init diff summary
    diff.summ <-  matrix(ncol=7,nrow=0)
    colnames(diff.summ) <- c("Estimate","Standard Error", "DF", "t-value", 
                             "Lower CI", "Upper CI", "p-value")
    summ.data <- as.data.frame(diff.summ)
  }
  
  
  for(eff in effs)
  {
    
    split.eff  <-  unlist(strsplit(eff,":"))
    if(checkAllCov(split.eff, data))
      next
    mat  <-  popMatrix(m, split.eff)
    fac.comb <- getFacCombForLSMEANS(split.eff, data)  
    
    if(!lsmeansORdiff)
      summ.data <- rbind(summ.data,   calcDiffsForEff(facs, fac.comb, split.eff,
                                                      eff, effs, data, rho, alpha,
                                                      mat))
    else
      summ.data <- rbind(summ.data,   calcLsmeansForEff(lsmeans.summ, fac.comb, 
                                                        eff, split.eff, alpha, 
                                                        mat, rho, facs))
  }
  return(list(summ.data = summ.data))
}
