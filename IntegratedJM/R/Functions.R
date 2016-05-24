## THE IntegratedJM PACKAGE ##


## IMPORTS ##

#' @import ggplot2
#' @import nlme
#' @import grid
#' @importFrom Biobase exprs
#' @importFrom Biobase featureNames
    

#####  FUNCTIONS #####


#' @title fitJM
#' 
#' @description The fitJM function fits the model for all the genes for a specific bio-activity vector and a particular fingerprint feature. 
#' @export
#' @param dat Contains the gene expression data matrix for all the genes - can be a matrix or an expression set.
#' @param covariate Vector of 0's and 1's, containing data about the fingerprint feature. 
#' @param responseVector Vector containing the bio-activity data.
#' @param methodMultTest Character string to specify the multiple testing method. Default is the BH-FDR method.
#' @return A data frame, containing the results of the model, to be used later for plots or to identify the top genes.
#' @details The default for the covariate parameter is NULL and if no covariate is specified it returns a data frame containing 5 variables, named as "Pearson","Spearman","p", "adj-p","logratio" 
#' and the data frame is ordered based on the column "p" which is the p-value obtained from the Log-Ratio Test. If there is a covariate, then the output is a dataframe containing 13 variables for all the genes,named as 
#' "adjPearson","adjSpearman","pPearson","Pearson", "Spearman", "pAdjR", "CovEffect1", "adjPeffect1", "CovEffect2", "adjPeffect2", "rawP1", "rawP2","logratio" and sorted based on "rawP1"  and "pPearson" which are
#' p-value corresponding to the effect of the fingerprint feature on the gene expression data as obtained from the t-table after fitting the model using gls and the p-value obtained from the Log-Ratio Test, respectively. 
#' In the first case without any covariate it calls the nullcov function inside it, otherwise the non_nullcov function is called to do the analysis.                    
#' @examples
#' \dontrun{
#' jmRes <- fitJM(dat=gene_eset,responseVector=activity,methodMultTest='fdr')
#' jmRes <- fitJM(dat=gene_eset,responseVector=activity,covariate = fp,methodMultTest='fdr')
#' }

fitJM <- function(dat, responseVector, covariate = NULL, methodMultTest) {
  
  data_error <- try(exprs(dat),silent=TRUE)
  
  if(class(data_error) == "try-error" && is.matrix(dat) == "FALSE" )
  {
    stop("The data type should be ExpressionSet or Matrix..!!")
  } 
  
  data_type = 0 ## 0 - eSet, 1 - matrix
  
  if(class(data_error) == "try-error")
  {
    data_type = 1
  }
  
    
  if(is.null(covariate)){
    
    JMResults <- nullcov(dat,responseVector,methodMultTest,data_type)
    
  } else {
    
    JMResults <- non_nullcov(dat,responseVector,covariate,methodMultTest,data_type)
    
  }
    
  return(JMResults)
}


#' @title nullcov
#' 
#' @description The nullcov function is called while fitting the model when there no covariate is specified in the fitJM function. It returns a data.frame containing the results after fitting the model. The output of this function is also the output of the fitJM function.
#' 
#' @export
#' @param dat Contains the gene expression data matrix for all the genes - can be a matrix or an expression set.
#' @param responseVector Vector containing the bio-activity data.
#' @param methodMultTest Character string to specify the multiple testing method. Default is the BH-FDR method.
#' @param data_type Binary, specifying the type of the parameter dat: 0 - expressionSet, 1 - matrix.
#' @return A data frame, containing the results of the model - same as the output of the fitJM function.
#' @details Fits the model using gls, calculates the correlation, p-values, adjusted p-values (based on the multiple testing method) and logratio from LRT and returns the required results.
#' @examples
#' \dontrun{
#' nullcov(dat=gene_eset,responseVector=activity,methodMultTest='fdr',data_type=0)
#' }


nullcov <- function(dat,responseVector,methodMultTest,data_type){
  
  myResp2 <- as.vector(responseVector)
  
  len <- as.numeric(nrow(dat))
  pearson <- spearman <- pValuePearson <- logratio <- c(1:len) 
  
  if(data_type == 0){
    
    gNames <- featureNames(dat)  #  featureData(fObject)$SYMBOL
        
    for (nGenes in 1:length(gNames)) {
      #  		print(nGenes)
      
      myResp1 <-as.vector(exprs(dat)[nGenes,])
      
      sampleID <- c(1:length(myResp1))							
      
      data <- data.frame(rbind(cbind(myResp1,sampleID), 
                               cbind(myResp2,  sampleID)),
                         respIndex = c(rep(1,length(myResp1)),rep(2,length(myResp2))))
      
      colnames(data) <- c("Responses","sampleID","respIndex")							# respIndex is response indicator
      
      #fit compound symmetry covariance structure
      f1 <- try( gls(Responses ~ as.factor(respIndex),
                     correlation = corSymm(form = ~ respIndex| sampleID),
                     weights = varIdent(form=~ 1|respIndex),
                     data = data, method = "ML"),silent=TRUE)
      
      #fit independence covariance structure
      f2 <- try(gls(Responses ~ as.factor(respIndex),
                    weights = varIdent(form=~ 1|respIndex),
                    data = data, method="ML"),silent=TRUE)
      
      #Likelihood Ratio Test
      LRT <- try(anova(f1,f2),silent=TRUE)
      logratio[nGenes] <- LRT[2,8]
      
      pValuePearson[nGenes] <- try(as.numeric(LRT$"p-value"[2]), 
                                   silent=TRUE)                                                                                                                                                                         
      
      #compute for correlation coeffcients 
      pearson[nGenes] <- try(getVarCov(f1)[1,2]/(sqrt(getVarCov(f1)[1,1])*
                                                   sqrt(getVarCov(f1)[2,2]) ),
                             silent=TRUE)
      
      spearman[nGenes] <- cor(residuals(f1)[1:length(responseVector)],
                              residuals(f1)[(length(responseVector)+1):nrow(data)], 
                              method = "spearman")
    }
    
    #adjust pValues
    pAdjPearson<- p.adjust(pValuePearson, method = methodMultTest, n = length(as.numeric(pValuePearson)))
                           
    
    #results
    result <- data.frame(pearson, spearman,  pValuePearson, pAdjPearson,logratio)
    
    colnames(result) <- c("Pearson","Spearman","p", "adj-p","logratio")
    rownames(result) <- gNames
    
    output <- result[order(result[,"p"]),]
    return(output)
    
} else {
  
  gNames <- rownames(dat)
  
  for (nGenes in 1:length(gNames)) {
    #    	print(nGenes)
    
    myResp1 <-as.vector(dat[nGenes,])
    
    sampleID <- c(1:length(myResp1))							
    
    data <- data.frame(rbind(cbind(myResp1,sampleID), 
                             cbind(myResp2,  sampleID)),
                       respIndex = c(rep(1,length(myResp1)),rep(2,length(myResp2))))
    
    colnames(data) <- c("Responses","sampleID","respIndex")							# respIndex is response indicator
    
    #fit compound symmetry covariance structure
    f1 <- try( gls(Responses ~ as.factor(respIndex),
                   correlation = corSymm(form = ~ respIndex| sampleID),
                   weights = varIdent(form=~ 1|respIndex),
                   data = data, method = "ML"),silent=TRUE)
    
    #fit independence covariance structure
    f2 <- try(gls(Responses ~ as.factor(respIndex),
                  weights = varIdent(form=~ 1|respIndex),
                  data = data, method="ML"),silent=TRUE)
    
    #Likelihood Ratio Test
    LRT <- try(anova(f1,f2),silent=TRUE)
    logratio[nGenes] <- LRT[2,8]
    
    pValuePearson[nGenes] <- try(as.numeric(LRT$"p-value"[2]), 
                                 silent=TRUE)                                                                                                                                                                         
    
    #compute for correlation coeffcients 
    pearson[nGenes] <- try(getVarCov(f1)[1,2]/(sqrt(getVarCov(f1)[1,1])*
                                                 sqrt(getVarCov(f1)[2,2]) ),
                           silent=TRUE)
    
    spearman[nGenes] <- cor(residuals(f1)[1:length(responseVector)],
                            residuals(f1)[(length(responseVector)+1):nrow(data)], 
                            method = "spearman")
  }
  
  #adjust pValues
  pAdjPearson<- p.adjust(pValuePearson, method = methodMultTest, 
                         n = length(as.numeric(pValuePearson)))
  
  #results
  result <- data.frame(pearson, spearman,  pValuePearson, pAdjPearson,logratio)
  
  colnames(result) <- c("Pearson","Spearman","p", "adj-p","logratio")
  rownames(result) <- gNames
  
  output <- result[order(result[,"p"]),]
  return(output)
  
}

}

#' @title non_nullcov
#' 
#' @description The non_nullcov function is called while fitting the model when the covariate is specified in the fitJM function. It returns a data.frame containing the results after fitting the model. The output of this function is also the output of the fitJM function.
#' 
#' @export
#' @param dat Contains the gene expression data matrix for all the genes - can be a matrix or an expression set.
#' @param covariate Vector of 0's and 1's, containing data about the fingerprint feature. 
#' @param responseVector Vector containing the bio-activity data.
#' @param methodMultTest Character string to specify the multiple testing method.
#' @param data_type Binary, specifying the type of the parameter dat: 0 - expressionSet, 1 - matrix.
#' @return A data frame, containing the results of the model - same as the output of the fitJM function.
#' @details Fits the model, adjusting for the covariate effect, using gls, calculates the correlation, p-values, adjusted p-values (based on the multiple testing method) and logratio from LRT and returns the required results.
#' @examples
#' \dontrun{
#' non_nullcov(dat=gene_eset,responseVector=activity,covariate=fp,methodMultTest='fdr',data_type=0)
#' }

non_nullcov <- function(dat,responseVector,covariate,methodMultTest,data_type){

  myResp2 <- as.vector(responseVector)
  
  len <- as.numeric(nrow(dat))
  pValuesEffects <- matrix(0,len,2)
  effects <- pValuePearson <- pearsonAdj <- spearmanAdj <- pearsonUnAdj <- spearmanUnAdj<- logratio <- c(1:len)
  
  
  myCov <- as.vector(covariate)
  sampleID <- c(1:length(myCov))
  
  beta_temp <- lm(responseVector~covariate)
  beta <- beta_temp$coefficients[2]
  pValuesEffects[,2] <- rep(summary(beta_temp)$coefficients["covariate","Pr(>|t|)"],len)
  
  if(data_type == 0){
  
    gNames <- featureNames(dat)  #	featureData(fObject)$SYMBOL
    
    for (nGenes in 1:len) {
      
      myResp1 <-as.vector(exprs(dat)[nGenes,])
      
      data <- data.frame(rbind(cbind(myResp1,myCov,sampleID), 
                               cbind(myResp2, myCov, sampleID)),
                         respIndex = c(rep(1,length(myCov)),rep(2,length(myCov))))
      colnames(data) <- c("Responses","myCov","sampleID","respIndex")				# respIndex is response indicator
      
      #fit compound symmetry covariance structure
      f1 <- try( gls(Responses ~ myCov*as.factor(respIndex),
                     correlation = corSymm(form = ~ respIndex| sampleID),
                     weights = varIdent(form=~ 1|respIndex),
                     data = data, method = "ML"),silent=TRUE)
      
      #extract covariate effects and pvalues
      effects[nGenes] <- summary(f1)$tTable['myCov',"Value"]
      pValuesEffects[nGenes,1] <- summary(f1)$tTable['myCov',"p-value"]
      
      #fit independence covariance structure
      f2 <- try(gls(Responses ~ myCov*as.factor(respIndex),
                    weights = varIdent(form=~ 1|respIndex),
                    data = data, method="ML"),silent=TRUE)
      
      #Likelihood ratio test
      LRT <- try(anova(f1,f2),silent=TRUE)		
      logratio[nGenes] <- LRT[2,8]
      ##LRTpval is same as pValuePearson
      
      #correlation coefficients (adjusted and Unadjusted)
      pearsonAdj[nGenes] <- try(getVarCov(f1)[1,2]/(
        sqrt(getVarCov(f1)[1,1])*
          sqrt(getVarCov(f1)[2,2]) ),
        silent=TRUE)
      pearsonUnAdj[nGenes] <- cor(myResp1,myResp2)
      spearmanUnAdj[nGenes] <- cor(myResp1,myResp2, method = "spearman")
      
      spearmanAdj[nGenes] <- cor(residuals(f1)[1:length(myResp2)],
                                 residuals(f1)[(length(myResp2)+1):nrow(data)], 
                                 method = "spearman")
      
      #pvalue of adjusted Pearson correlation
      pValuePearson[nGenes] <- try(as.numeric(LRT$"p-value"[2]),
                                   silent=TRUE)                                                                                                                                                                         
    }
    
    # adjusted pvalues
    pAdj <- apply(pValuesEffects,2, function(x) p.adjust(x,method = methodMultTest, n = length(x)))
                                                         
    pAdjPearson <- p.adjust(pValuePearson, method = methodMultTest, n = length(as.numeric(pValuePearson)))
                            
    
    #result
    
    result <- data.frame(pearsonAdj, spearmanAdj,  pValuePearson,
                         pearsonUnAdj, spearmanUnAdj, pAdjPearson, 
                         effects, pAdj[,1], rep(beta,len), pAdj[,2], 
                         pValuesEffects,logratio)
    
    colnames(result) <- c("adjPearson","adjSpearman","pPearson",
                          "Pearson", "Spearman", "pAdjR", 
                          "CovEffect1", "adjPeffect1", "CovEffect2", "adjPeffect2",
                          "rawP1", "rawP2","logratio")
    rownames(result) <- gNames
    
    output <- result[order(result[,"rawP1"],result[,"pPearson"]),]
    return(output)  

  } else {
    
    for (nGenes in 1:len) {
      
      gNames <- rownames(dat) 
      
      myResp1 <-as.vector(dat[nGenes,])
      
      data <- data.frame(rbind(cbind(myResp1,myCov,sampleID), 
                               cbind(myResp2, myCov, sampleID)),
                         respIndex = c(rep(1,length(myCov)),rep(2,length(myCov))))
      colnames(data) <- c("Responses","myCov","sampleID","respIndex")  			# respIndex is response indicator
      
      #fit compound symmetry covariance structure
      f1 <- try( gls(Responses ~ myCov*as.factor(respIndex),
                     correlation = corSymm(form = ~ respIndex| sampleID),
                     weights = varIdent(form=~ 1|respIndex),
                     data = data, method = "ML"),silent=TRUE)
      
      #extract covariate effects and pvalues
      effects[nGenes] <- summary(f1)$tTable['myCov',"Value"]
      
      pValuesEffects[nGenes,] <- summary(f1)$tTable[c('myCov', 
                                                      'myCov:as.factor(respIndex)2'),"p-value"]
      
      #fit independence covariance structure
      f2 <- try(gls(Responses ~ myCov*as.factor(respIndex),
                    weights = varIdent(form=~ 1|respIndex),
                    data = data, method="ML"),silent=TRUE)
      
      #Likelihood ratio test
      LRT <- try(anova(f1,f2),silent=TRUE)
      logratio[nGenes] <- LRT[2,8]
      
      #correlation coefficients (adjusted and Unadjusted)
      pearsonAdj[nGenes] <- try(getVarCov(f1)[1,2]/(
        sqrt(getVarCov(f1)[1,1])*
          sqrt(getVarCov(f1)[2,2]) ),
        silent=TRUE)
      pearsonUnAdj[nGenes] <- cor(myResp1,myResp2)
      spearmanUnAdj[nGenes] <- cor(myResp1,myResp2, method = "spearman")
      
      spearmanAdj[nGenes] <- cor(residuals(f1)[1:length(myResp2)],
                                 residuals(f1)[(length(myResp2)+1):nrow(data)], 
                                 method = "spearman")
      
      #pvalue of adjusted Pearson correlation
      pValuePearson[nGenes] <- try(as.numeric(LRT$"p-value"[2]),
                                   silent=TRUE)                                                                                                                                                                         
    }
    
    # adjusted pvalues
    pAdj <- apply(pValuesEffects,2, function(x) p.adjust(x,method = methodMultTest, n = length(x)))
    
    pAdjPearson <- p.adjust(pValuePearson, method = methodMultTest, n = length(as.numeric(pValuePearson)))
    
    
    #result
    
    result <- data.frame(pearsonAdj, spearmanAdj,  pValuePearson,
                         pearsonUnAdj, spearmanUnAdj, pAdjPearson, 
                         effects, pAdj[,1], rep(beta,len), pAdj[,2], 
                         pValuesEffects,logratio)
    
    colnames(result) <- c("adjPearson","adjSpearman","pPearson",
                          "Pearson", "Spearman", "pAdjR", 
                          "CovEffect1", "adjPeffect1", "CovEffect2", "adjPeffect2",
                          "rawP1", "rawP2","logratio")
    rownames(result) <- gNames
    
    output <- result[order(result[,"rawP1"],result[,"pPearson"]),]
    return(output) 
    
  }


}

#' @title topkGenes
#' 
#' @description The topkGenes function is to identify the top genes based on different criteria. 
#' @export
#' @param jointModelResult Data frame, containing the results from the fitJM function.
#' @param subset_type Character string to specify the set of genes. It can have four values: "Effect" for only differentially expressed genes, "Correlation" for only correlated genes, "Effect and Correlation" for genes which are both differentially expressed & correlated and "Other" for the genes which are neither differentially expressed nor correlated.
#' @param ranking Character string, specifying one of the columns of the jointModelResult data frame, based on the genes will be ranked within the selected subset.
#' @param k Integer, specifying the number of genes, to be returned from the list of top genes. Default is 10.
#' @param sigLevel Numeric between 0 and 1, specifying the level of significance, used to select the subset of genes.
#' @return A data frame containing top k genes according to the specified criteria from the specified set of genes.
#' @details Returned data frame contains 6 columns, named as "Genes","FP-Effect", "p-adj(Effect)", "Unadj.Asso.","Adj.Asso.", "p-adj(Adj.Asso.)".
#' @examples
#' \dontrun{
#' topkGenes(jointModelResult=jmRes,subset_type="Effect",ranking="Pearson",k=10,sigLevel = 0.05)
#' }

topkGenes <- function(jointModelResult, subset_type, ranking, k=10, sigLevel = 0.01){
    
    sigLevel <- sigLevel
    
    if(subset_type == "Effect"){
    
      output <- subset(jointModelResult,adjPeffect1 < sigLevel & pAdjR  > sigLevel)
      
    } else if (subset_type == "Correlation" ) {
      
      output <- subset(jointModelResult,pAdjR  < sigLevel & adjPeffect1 > sigLevel )
      
    } else if (subset_type == "Effect and Correlation"){
      
      output <- subset(jointModelResult,pAdjR  < sigLevel & adjPeffect1 < sigLevel )
      
    } else if(subset_type == "Other") {
      
      output <- subset(jointModelResult,pAdjR  > sigLevel & adjPeffect1 > sigLevel )
    } else {
      
      stop("Please specify the type")
    }
    
    ranked_output <- output[order(abs(output[,ranking]),decreasing=TRUE),]
    
    topGenes <- data.frame(rownames(ranked_output), ranked_output$CovEffect1,ranked_output$adjPeffect1, ranked_output$Pearson, ranked_output$adjPearson,ranked_output$pAdjR)
    
    colnames(topGenes) <- c("Genes","FP-Effect", "p-adj(Effect)", "Unadj.Asso.","Adj.Asso.", "p-adj(Adj.Asso.)")
    
    return(topGenes[1:k,])
  
}


#' @title plotAsso
#' 
#' @description The plotAsso function is used to plot the unadjusted association vs the adjusted association for all the genes.
#' @export
#' @param jointModelResult Data frame, containing the results from the fitJM function.
#' @param type Character string, specifying the type of association - Pearson or Spearman. 
#' @return Creates a plot
#' @details Plots the unadjusted association vs the adjusted association for all the genes.
#' @examples
#' \dontrun{
#' plotAsso(jointModelResult=jmRes,type="Pearson")
#' }
                     
plotAsso <- function(jointModelResult,type){
  
   
  if(type == "Pearson"){
    
    adj <- jointModelResult$adjPearson
    unadj <- jointModelResult$Pearson
    
    qplot(unadj, adj, ylab = "Adjusted Association", xlab = "Unadjusted Association",
          main = "Pearson")
    
  } else if(type == "Spearman") {
    
    adj <- jointModelResult$adjSpearman
    unadj <- jointModelResult$Spearman
    
    qplot(unadj, adj, ylab = "Adjusted Association", xlab = "Unadjusted Association",
          main = "Spearman")
  } else {
    stop("Please specify Pearson or Spearman")
  }
}


#' @title plotEff
#' 
#' @description The plotEff function is used to plot the fingerprint effect on gene expression vs the adjusted association for all the genes.
#' @export
#' @param jointModelResult Data frame, containing the results from the fitJM function.
#' @param type Character string, specifying the type of association - Pearson or Spearman. 
#' @return Creates a plot
#' @details Plots the fingerprint effect on gene expression vs the specified type of adjusted association for all the genes.
#' @examples
#' \dontrun{
#' plotEff(jointModelResult=jmRes,type="Pearson")
#' }

plotEff <- function(jointModelResult,type){

  
  if(type == "Pearson"){
    
    adj <- jointModelResult$adjPearson
    alpha <- jointModelResult$CovEffect1
    
    qplot(alpha, adj, ylab = "Adjusted Association(Pearson)", xlab = "FP Effect on Gene Expression")
    
  } else if(type == "Spearman") {
    
    adj <- jointModelResult$adjSpearman
    alpha <- jointModelResult$CovEffect1
    
    qplot(alpha, adj, ylab = "Adjusted Association(Spearman)", xlab = "FP Effect on Gene Expression")
  } else {
    stop("Please specify Pearson or Spearman")
  }
}


#' @title multiplot
#' 
#' @description The multiplot function plots multiple ggplots in the same window.
#' @export
#' @param ... ggplot2 objects, separated by comma.
#' @param cols Integer, specifying the number of plots in one row in the layout.
#' @return Creates multiple ggplots in same window 
#' @details Plots multiple ggplots in the same window - multiplot(p1,p2,p3,p4, cols=2) is similar to the standard R notation par(mfrow=c(2,2)).
#' @examples
#' \dontrun{
#' multiplot(p1,p2,p3,cols=3)
#' }

multiplot <- function(..., cols=1) {
  
  
  plots <- list(...)
  
  numPlots = length(plots)
  
  layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                   ncol = cols, nrow = ceiling(numPlots/cols))
  
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




#' @title plot1gene
#' 
#' @description The plot1gene function plots the data for a single gene.
#' @export
#' @param geneName Character string, specifying the name of the gene.
#' @param fp Vector containing 0's and 1's - the data about the fingerprint feature.
#' @param fpName Character string, used to make the title of the plots. If not specified, the plot title will be blank.
#' @param responseVector Vector containing the bio-activity data.
#' @param dat Contains the gene expression data matrix for all the genes - can be a matrix or an expression set.
#' @param resPlot Logical. If TRUE, also plots the residual from the gls fit. Default is TRUE.
#' @param colP Character string, specifying the colour for the 1's in the fp parameter. Default is blue.
#' @param colA Character string, specifying the colour for the 0's in the fp parameter. Default is white.
#' @return Creates a plot
#' @details Calls the getCorrUnad function and creates the plot(s) accordingly.
#' @examples
#' \dontrun{
#' plot1gene(geneName="Gene21",fp=fp,fpName="Fingerprint",responseVector=activity,dat=gene_eset)
#' }


plot1gene <- function(geneName,fp,fpName = "",responseVector,dat,resPlot=TRUE,colP = "blue",colA = "white"){ 
  
  
  
  if(is.character(geneName) == TRUE && length(fp) == length(responseVector)) {
    
    if(resPlot){
      
      
      data1gene_1 <- getCorrUnad(geneName,fp,fpName,responseVector,dat,resPlot)[[1]]
      data1gene_2 <- getCorrUnad(geneName,fp,fpName,responseVector,dat,resPlot)[[2]]
      UnadjPearson <- getCorrUnad(geneName,fp,fpName,responseVector,dat,resPlot)[[3]]
      AdjPearson <- getCorrUnad(geneName,fp,fpName,responseVector,dat,resPlot)[[4]]
      
      plot1 <- ggplot(data1gene_1,aes(x=as.numeric(as.vector(gene)), y=as.numeric(as.vector(activity)))) +
        xlab("Gene Expression") +  
        ylab ("Activity")+
        facet_grid(.~fpId, scales="free")+
        geom_point(aes(fill=factor(fp)),size=8, shape=21) +
        scale_fill_manual("FP:",values = c("0" = colA,"1" = colP),labels=c("0 - absent", "1 - present"))+
        ggtitle(paste("Unadj. Asso.",UnadjPearson))+
        geom_smooth(method = lm)+
        theme(strip.text.x = element_text(face="bold",size=14))
      
      plot2 <- ggplot(data1gene_2,aes(x=as.numeric(as.vector(gene)), y=as.numeric(as.vector(activity)))) +
        xlab("Gene Expression") +  
        ylab ("Activity")+
        facet_grid(.~fpId, scales="free")+
        geom_point(aes(fill=factor(fp)),size=8, shape=21) +
        scale_fill_manual("FP:",values = c("0" = colA,"1" = colP),labels=c("0 - absent", "1 - present"))+
        ggtitle(paste("Adj. Asso.",AdjPearson))+
        geom_smooth(method = lm)+
        theme(strip.text.x = element_text(face="bold",size=14))
      
      
      multiplot(plot1,plot2,cols=1)
      
    } else {
      
      
      data1gene <- getCorrUnad(geneName,fp,fpName,responseVector,dat,resPlot)[[1]]
      UnadjPearson <- getCorrUnad(geneName,fp,fpName,responseVector,dat,resPlot)[[2]]
      
      ggplot(data1gene,aes(x=as.numeric(as.vector(gene)), y=as.numeric(as.vector(activity)))) +
        xlab("Gene Expression") +  
        ylab ("Activity")+
        facet_grid(.~fpId, scales="free")+
        geom_point(aes(fill=factor(fp)),size=8, shape=21) +
        scale_fill_manual("FP:",values = c("0" = colA,"1" = colP),labels=c("0 - absent", "1 - present"))+
        ggtitle(paste("Unadj. Asso.",UnadjPearson))+
        geom_smooth(method = lm)+
        theme(strip.text.x = element_text(face="bold",size=14))
      
    }
    
  }   else {
    
    stop("missing or wrong data type of parameter:geneName should be character and fp should be numeric")
  }
  
}


#' @title getCorrUnad
#' 
#' @description The getCorrUnad function is a support function for the function plot1gene.
#' @export
#' @param geneName Character string, specifying the name of the gene.
#' @param fp Vector containing 0's and 1's - the data about the fingerprint feature.
#' @param fpName Character string, used to make the title of the plots. If not specified, the plot title will be blank.
#' @param responseVector Vector containing the bio-activity data.
#' @param dat Contains the gene expression data matrix for all the genes - can be a matrix or an expression set.
#' @param resPlot Logical. If TRUE, creates the plot data for the residual plot
#' @return A list containing the data to create the respective plots and the unadjusted association between the gene expression and bio-activity data.
#' @details Works as a support function for plot1gene.
#' @examples
#' \dontrun{
#' getCorrUnad(geneName="Gene21",fp=fp,fpName="Fingerprint",
#' responseVector=activity,dat=gene_eset,resPlot=TRUE)
#' }

getCorrUnad <- function(geneName,fp,fpName,responseVector,dat,resPlot){
  
  data_error <- try(exprs(dat),silent=TRUE)
  
  if(class(data_error) == "try-error" && is.matrix(dat) == "FALSE" )
  {
    stop("The data type should be ExpressionSet or Matrix..!!")
  } 
  
  data_type = 0 ## 0 - eSet, 1 - matrix
  
  if(class(data_error) == "try-error")
  {
    data_type = 1
  }
  
  if(data_type == 0){
    myResp1 <-as.vector(exprs(dat)[geneName,])
  } else {
    myResp1 <- as.vector(dat[geneName,])
  }
  
  
  myResp2 <- as.vector(responseVector)
  pearson <- round(cor(myResp1, myResp2),4)
  
  if(resPlot){
    
    ## calculate residuals and prepare plot data accordingly
    
    myCov <- as.vector(fp)
    sampleID <- c(1:length(myCov))  						
    
    data <- data.frame(rbind(cbind(myResp1,myCov,sampleID), 
                             cbind(myResp2, myCov, sampleID)),
                       respIndex = c(rep(1,length(myCov)),rep(2,length(myCov))))
    colnames(data) <- c("Responses","myCov","sampleID","respIndex")				# respIndex is response indicator
    
    f1 <- try( gls(Responses ~ myCov*as.factor(respIndex),
                   correlation = corSymm(form = ~ respIndex| sampleID),
                   weights = varIdent(form=~ 1|respIndex),
                   data = data, method = "ML"),silent=TRUE)
    
    gcorr<-round(try(getVarCov(f1)[1,2]/(sqrt(getVarCov(f1)[1,1])* sqrt(getVarCov(f1)[2,2]) ),silent=TRUE),4)						
    
    plot_dat1 <- data.frame(cbind(gene=myResp1,activity=responseVector, fp=fp,fpId=c(rep(paste0(geneName,",Observed:",fpName),length(fp)))))
    plot_dat2 <- data.frame(cbind(gene=residuals(f1)[1:length(myResp1)],activity=residuals(f1)[(length(myResp1)+1):length(data$Responses)], 
                                  fp=fp,fpId=c(rep(paste0(geneName,",Residuals:",fpName),length(fp)))))
    dataoneGene <- list(plot_dat1,plot_dat2,pearson,gcorr) 
    return(dataoneGene)
    
    
  } else {
    
    
    plot_dat <- data.frame(cbind(gene=myResp1,activity=responseVector, fp=fp,fpId=c(rep(paste0(geneName,",FP:",fpName),length(fp)))))
    dataoneGene <- list(plot_dat,pearson) 
    return(dataoneGene)
    
    
    
  }
  
  
  
  
}


#' @title volcano
#' 
#' @description The volcano function produces the volcano plot for logratio / fp-effect vs corresponding p-values.
#' @export
#' @param x Numeric vector of logratios or covariate effect values to be plotted.
#' @param pValue Numeric vector of corresponding p-values obtained from some statistical test.
#' @param pointLabels Character vector providing the texts for the points to be labelled in the plot.
#' @param topPValues Number of top p-values to be labelled. Default value is 10.
#' @param topXvalues Number of top logratios or covariate effect values to be labelled. Default value is 10.
#' @param smoothScatter Logical parameter to decide if a smooth plot is expected or not. Default is TRUE.
#' @param xlab Text for the x-axis of the plot. Default is NULL.
#' @param ylab Text for the y-axis of the plot. Default is NULL.
#' @param main Text for the main title of the plot. Default is NULL.
#' @param newpage Logical parameter
#' @param additionalPointsToLabel Set of points other than the top values to be labelled in the plot. Default is NULL.
#' @param additionalLabelColor Colour of the additionally labelled points. Default colour is red.
#' @param dir Logical parameter deciding if the top values should be in decreasing (= TRUE) or increasing (= FALSE) order. Default is TRUE.
#' @return A plot which looks like a volcano.
#' @details Creates a plot which looks like a volcano with the interesting points labelled within the plot.
#' @examples
#' \dontrun{
#' volcano(x=jmRes$CovEffect1,pValue=jmRes$rawP1,pointLabels=rownames(jmRes),
#' topPValues = 10, topXvalues = 10,xlab="FP Effect (alpha)",ylab="-log(p-values)")
#' }



volcano <- function (x, pValue, pointLabels, topPValues = 10, topXvalues = 10, 
                smoothScatter = TRUE, xlab = NULL, ylab = NULL, main = NULL,
                newpage = TRUE, additionalPointsToLabel = NULL, 
                additionalLabelColor = "red", dir=TRUE) {
  logRatio <- x
  pVals <- -log10(pValue)
  topLR <- order(abs(logRatio), decreasing = dir)[seq(length.out = topXvalues)]
  topP <- order(abs(pVals), decreasing = dir)[seq(length.out = topPValues)]
  pointsToLabel <- union(topP, topLR)
  pointsToLabel <- union(pointsToLabel, which(names(logRatio) %in% additionalPointsToLabel))
                                                
  if (is.null(additionalPointsToLabel)) {
    colPointsToLabel <- rep("black", length(pointsToLabel))
  } else {
    if (is.null(names(logRatio))) {
      stop("labeling additional points required a named vector for the logRatios")
    } else {
      colPointsToLabel <- ifelse(names(logRatio)[pointsToLabel] %in% 
                                   additionalPointsToLabel, additionalLabelColor, 
                                 "black")
    }
  }
  
  if (newpage) 
  
  grid.newpage()
  
  pvp <- plotViewport(c(5, 6, 5, 3))
  pushViewport(pvp)
  tg <- textGrob(label = pointLabels[pointsToLabel], x = unit(logRatio[pointsToLabel], 
                                                              "native"), y = unit(pVals[pointsToLabel], "native"), 
                 gp = gpar(cex = 0.65, col = colPointsToLabel))
  maxLabelWidth <- max(grobWidth(tg))
  nMaxLabelWidth <- convertHeight(maxLabelWidth, "native", valueOnly = TRUE)
  
  xr <- max(abs(min(logRatio)),abs(max(logRatio)))
  xr_ep <- c(-xr,xr)
  
  dvp <- dataViewport(xscale = xr_ep + c(-nMaxLabelWidth/8, 
                                                   nMaxLabelWidth/8), yscale = range(pVals, na.rm = TRUE))
  pushViewport(dvp)
  atPositionsY <- seq(0,ceiling(max(current.viewport()$yscale)),2)
  grid.yaxis(name = "ya", at = atPositionsY, label = atPositionsY)
  xa <- xaxisGrob(name = "xa")
  moveUnit <- unit(-0.5, "char")
  xa <- editGrob(xa, edits = gEditList(gEdit("major", y = moveUnit),gEdit("ticks", y0 = moveUnit), gEdit("ticks", y1 = unit(-0.5, "lines") + moveUnit),
                                         gEdit("labels", y = unit(-1.5, "lines") + moveUnit)))
                                                                                               
                                                                                                                                              
  grid.draw(xa)
  
  dotColors <- if (smoothScatter) {
    densCols(x = logRatio[-pointsToLabel], y = pVals[-pointsToLabel],colramp = colorRampPalette(blues9[-(1:3)]))
  } else {
    "#9ECAE1"
  }
  
  grid.points(x = unit(logRatio[-pointsToLabel], "native"), 
              y = unit(pVals[-pointsToLabel], "native"), pch = 20, 
              gp = gpar(col = dotColors))
  grid.draw(tg)
  if (!is.null(main)) {
    grid.text(label = main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))
                                                      
  }

  if (!is.null(xlab)) {
    grid.text(label = xlab, y = unit(-3, "lines") + 0.5 * 
                moveUnit)
  }
  if (!is.null(ylab)) {
    grid.text(label = ylab, x = unit(-4.5, "lines"), rot = 90)
  }
}






if(getRversion() >= "3.0.0"){
	globalVariables(c("gene", "activity", "adjPeffect1", "pAdjR"))
}


