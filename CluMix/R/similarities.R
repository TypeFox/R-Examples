
## fast calculation of Goodman&Kruskal's gamma
myGKgamma <- function(x, y){
  
  Hmisc::rcorr.cens(as.numeric(x), as.numeric(y), outx=TRUE)["Dxy"]  # outx=T is important to leave out tied observations
}
# -> "fast" for continuous vs categorical
# -> faster for categorical vs categorical: GoodmanKruskalGamma from package DescTools


## calculate chi^2 test p-value; relevant lines taken from chisq.test()
chisq <- function(x, group){
  OK <- complete.cases(x, group)
  x <- factor(x[OK])
  y <- factor(group[OK])
  x <- table(x, y)
  n <- sum(x)
  nr <- as.integer(nrow(x))
  nc <- as.integer(ncol(x))
  sr <- rowSums(x)
  sc <- colSums(x)
  E <- outer(sr, sc, "*")/n
  
  STATISTIC <- sum((abs(x - E))^2/E)
  PARAMETER <- (nr - 1L) * (nc - 1L)
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  return(PVAL)
}



## calculate new measure for association btw. ranked (continuous or ordinal) variable A and categorical (K>2) variable B
assoc.rank.cat <- function(x, a, alpha=0.05){
  # x: continuous or ordinal vector
  # a: factor (not ordinal)
  # alpha: significance level for Kruskal-Wallis test
  
  # test if there might be any association
  p <- kruskal.test(x, a)$p.value 
  
  # only go on with optimizing ordering of a in case of a significant result
  if(p < alpha){ 
    
    # order levels of a variable according to mean ranks of x
    a <- factor(as.numeric(a))
    r <- rank(x, na.last="keep")
    m <- tapply(r, a, median, na.rm=TRUE) 
    a <- ordered(a, levels=order(m))
  }
  
  # continuous x -> Spearman
  if(data.class(x) == "numeric")
    res <- abs(cor(x, as.numeric(a), method="spearman", use="complete.obs"))
  
  # ordinal x -> Goodman & Kruskal's gamma
  if(data.class(x) == "ordered"){
    res <- abs(DescTools::GoodmanKruskalGamma(x, a))
  }
  
  return(res)
}



## function to "impose an order" in categorical (K>2) variables by "diagonalization" of cross table
assoc.cat.cat <- function(a, b, alpha=0.05){
  # a, b: factors (not ordinal) with >2 categories
  # alpha: significance level for Kruskal-Wallis test
  
  # when both variables are binary, of course no reordering is necessary
  ka <- length(table(a))
  kb <- length(table(b))
  if(ka > 2 | kb > 2){
    
    # test if there might be any association
    p <- chisq(a, b)
    
    # only go on with diagonalization of cross table in case of a significant result
    if(p < alpha){ 
      tab <- table(a, b)
      
      # get "optimal" ordering of categories (large frequencies in the diagonal)
      tab.new <- extracat::optile(tab)
      
      # reorder categories correspondingly
      a <- ordered(a, levels=rownames(tab.new))
      b <- ordered(b, levels=colnames(tab.new))  
    }
  }
  
  # calculate Goodman & Kruskal's gamma on "optimized" cross table
  abs(DescTools::GoodmanKruskalGamma(a, b))  
}



