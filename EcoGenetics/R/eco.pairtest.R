#' Kruskall - Wallis + Wilcoxon (Mann-Whitney U) and aov + Tukey-HSD tests 
#' for an ecogen object
#' 
#' @details This program returns the Wilcoxon (Mann-Whitney U) or Tukey-HSD 
#' statistics and p values for the multiple comparisons of the variables contained
#' in the selected data frame, among the levels of a factor of the slot "S".
#' 
#' @param eco Object of class "ecogen".
#' @param df The data frame for the analysis. Could be "P", "E" or "C".
#' @param x The name of the S slot column with the groups for the analysis.
#' @param test Test to perform ("wilcoxon", "tukey").
#' @param  only.p  Should be only returned a matrix with P-values? 
#' Default TRUE.
#' @param adjust P-values correction method for multiple tests 
#' passed to \code{\link[stats]{p.adjust}}. Defalut is "fdr".
#' @param ... Additional arguments passed to \code{\link{wilcox.test}} 
#' or  \code{\link{TukeyHSD}}.
#' 
#' @seealso \code{\link{wilcox.test}} \code{\link{TukeyHSD}}
#' 
#' @examples 
#' \dontrun{
#' data(eco3)
#' wil <- eco.pairtest(eco = eco3, df = "P", x = "structure")
#' wil
#' wil <- eco.pairtest(eco = eco3,df = "E", x = "structure")
#' wil
#' wil <- eco.pairtest(eco = eco3, df = "P", x = "structure", only.p = FALSE)
#' wil
#' wil <- eco.pairtest(eco = eco3,df = "P", x = "structure", test = "tukey")
#' wil
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.pairtest", 
					 
					 function(eco, df = c("P", "E","A", "C"),
					 				 x, test = c("wilcoxon", "tukey"),
					 				 adjust = "fdr",
					 				 only.p = TRUE, ...) {
					 	
	test <- match.arg(test)

	grupo <- eco@S
	fact <- match(x, colnames(eco@S), nomatch = 0)
	fact <- fact[!fact == 0]
	if(length(fact) == 0) {
		stop("incorrect factor name")
	}
	
	grupos <- as.factor(eco@S[, fact])
	
	P <- match.arg(df)
	
	if(P == "P") {
		P<-eco@P
	} else if(P == "E") {
		P<-eco@E
	} else if(P == "A") {
		P<- eco@A
	} else if(P == "C") {
		P<-eco@C
	}
  
  if(test == "wilcoxon") {
    
    
    niveles <- as.numeric(max(levels(grupos)))
    lev <- list()
    
    tabla <- table(1:length(grupos), grupos)
    
    a <- matrix(0, nrow = niveles, ncol = niveles)
    a <- upper.tri(a)
    index <- which(a == TRUE, arr.ind = TRUE)
    index <- data.frame(index, rep(0, nrow(index)), rep(0, nrow(index)))
    colnames(index) <- c("pop1", "pop2", "est", "P")
    res <- list()
    krus <- list()
    for(i in 1:ncol(P)) {
      tabla <- index
      temp <- kruskal.test(P[,i]~grupos)
      krus[[i]] <- c(temp$statistic, temp$p.value)
      for(j in 1:nrow(index)) {
        wil <- wilcox.test(P[grupos == index[j, 1], i], 
                           P[grupos == index[j, 2],i], ...)
        tabla[j, 3] <- wil$statistic
        tabla[j, 4] <- round(p.adjust(wil$p.value, method = adjust), 4)
      }
      res[[i]] <- tabla
    }
    krus <- sapply(krus, c)
    krus <- round(krus, 4)
    colnames(krus) <- colnames(P)
    rownames(krus) <- c("statistic", "P")
    names(res) <- names(P)
    
    if(only.p == TRUE) {
    pvalues <- sapply(res, function(y) {y[, 4]})
    colnames(pvalues) <- colnames(P)
    rownames(pvalues) <- paste(index[, 1], index[, 2], sep = "-")
    return(list(kruskall.test = krus, wilcoxon.test = pvalues, 
    						p.correction = adjust))
    } else {
      return(list(kruskall.test = krus, wilcoxon.test = res, 
      						p.correction = adjust))
    } 
    
  } else if(test == "tukey") {
    
    m.tukey <- function(x, y) {
    general <- listafac <- tukeyfac <- list()
      for(i in 1:ncol(x)) {
        listafac[[i]] <- aov(x[, i] ~ grupos)
        temporal <- anova(listafac[[i]])
        general[[i]] <- c(temporal[,4][1], temporal[, 5][1])
        tukeyfac[[i]] <- post <- TukeyHSD(listafac[[i]], ...)
      }
      general <- sapply(general, c)
      rownames(general) <- c("F", "P")
      colnames(general) <- colnames(P)
      names(tukeyfac) <- names(P)
      res <- list(general, tukeyfac)
    }
    resultados <- m.tukey(P, grupos)
    tuk <- resultados[[2]]
    general <- round(resultados[[1]], 4)
    p.values <- sapply(tuk, function(y) round(y$grupos[, 4], 4))
    rownames(p.values) <- rownames(tuk[[1]]$grupos)
    colnames(p.values) <- colnames(P)
    
    if(only.p == FALSE) {
      return(list(aov = general, pairtest = tuk))
    } else  {
    return(list(aov = general, pairtest = p.values))
    }
    
  }
})
