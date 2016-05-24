#' Profile Analysis by Group: Testing Parallelism, Equal Levels, and Flatness
#'
#' The \code{pbg} function implements three hypothesis tests. These tests are whether the profiles are parallel, have equal levels, and are flat across groups defined by the grouping variable. If parallelism is rejected, the other two tests are not necessary. In that case, flatness may be assessed within each group, and various within- and between-group contrasts may be analyzed. 
#'
#' @importFrom graphics par matplot axis legend
#' @importFrom stats pf var aov manova na.omit
#' @param data A matrix or data frame with multiple scores; rows represent individuals, columns represent subscores. Missing subscores have to be inserted as NA.
#' @param group A vector or data frame that indicates a grouping variable. It can be either numeric or character (e.g., male-female, A-B-C, 0-1-2). The grouping variable must have the same length of x. Missing values are not allowed in y.
#' @param original.names Use original column names in x. If FALSE, variables are renamed using v1, v2, ..., vn for subscores and "group" for the grouping variable. Default is FALSE.
#' @param profile.plot Print a profile plot of scores for the groups. Default is FALSE.
#'
#' @return An object of class \code{profg} is returned, listing the following components:
#' \itemize{
#' \item \code{data.summary} - Means of observed variables by the grouping variable	
#' \item \code{corr.table} - A matrix of correlations among observed variables splitted by the grouping variable
#' \item \code{profile.test} - Results of F-tests for testing parallel, coincidential, and level profiles across two groups.
#' }
#' @examples
#' \dontrun{
#'data(spouse)
#'mod <- pbg(data=spouse[,1:4], group=spouse[,5], original.names=TRUE, profile.plot=TRUE)
#'print(mod) #prints average scores in the profile across two groups
#'summary(mod) #prints the results of three profile by group hypothesis tests
#' }
#' @seealso \code{\link{pr}}, \code{\link{profileplot}}
#'@export

pbg <- function(data, group, original.names=FALSE, profile.plot=FALSE) {
  
  x <- as.data.frame(data)
  y <- as.factor(group)
  n=nrow(x)
  m=length(group)
  k=ncol(x)
  g=length(unique(y))
  
  if(!(n == m)){stop("The number of rows are not equal between the data and the group.")}
  else {
    z <- as.data.frame(cbind(x,y))
    cor.table=by(z[,1:k],z[,(1+k)],cor)
    
    if(original.names) {colnames(z) <- c(labs <- colnames(x),"group")}
    else {colnames(z) <- c(labs <- paste("v",1:k,sep=""),"group")}
    z$group <- as.factor(z$group)
    
    average <- matrix(NA,k,g) #average scores
    colnames(average) <- levels(z[,(k+1)])
    rownames(average) <- labs
    
    for(groups in colnames(average)) # for each column of grouping variable
      average[,groups] <- colMeans(z[z$group==groups,1:k],na.rm = TRUE)
    rm(groups)
    
    #The following part creates a profile plot of observed variables across two groups
    if(profile.plot) {
      par(mar=c(4.1,4.1,0.5,0.5))
      p.plot <- matplot(1:k, average, type="b", pch=21:(21+g-1),
                        xaxt="n", ylab="Mean Score", xlab="Observed Variables",cex=1.5)
      axis(1,at=1:k,labels=labs)
      legend(x="topleft", legend=colnames(average), lty=1:g, pch=21:(21+g-1)) 
      }
    
    #Testing parallel profiles
    #1) Create a constrast matrix
    cnt <- c(1,-1,rep(0,k-2))
    
    for (i in 2:(k-1)) {
      cnt <- c(cnt, c(rep(0,(i-1)),1,-1,rep(0,(k-i-1))))
    }
    
    cont1 <- matrix(cnt, nrow=k, ncol=(k-1))
    
    #2) Multiple contrast matrix with the data
    xm <- as.matrix(x)
    deviation1 <- xm%*%cont1
    
    #3) Test parallel profiles using multivariate tests
    parallel <- data.frame(Multivariate.Test=c("Wilks","Pillai","Hotelling-Lawley","Roy"),Statistic=rep(0,4),Approx.F=rep(0,4),
                           num.df=rep(0,4), den.df=rep(0,4),p.value=rep(0,4))
    
    fit1 <- manova(deviation1 ~ y, na.action=na.omit)
    parallel[1,2:6] <- summary(fit1, test="Wilks")$stats[1,2:6]
    parallel[2,2:6] <- summary(fit1, test="Pillai")$stats[1,2:6]
    parallel[3,2:6] <- summary(fit1, test="Hotelling-Lawley")$stats[1,2:6]
    parallel[4,2:6] <- summary(fit1, test="Roy")$stats[1,2:6]

    #Equal levels: Sum all variables, get the average, run univariate anova
    z$sum <- rowMeans(data, na.rm = TRUE)
    fit2 <- aov(sum ~ group, data=z)
    equal.levels <- summary(fit2)
    
    #Flatness
    cont2 <- diag(x = -1, k, k)
    for(i in 1:(k-1)) {cont2[i,(i+1)]=1}
    cont2 <- cont2[1:(k-1),]
    allS <- var(z[,1:k]) #overall variance
    xbar <- apply(z[,1:k],2,mean) #grand mean
    xbar.dif <- cont2 %*% xbar
    
    c.vector <- vector()
    for (i in 1:g) {
      c.vector[i] <- nrow(z[z$group==levels(z$group)[i],])
    }
    
    F.val <- (sum(c.vector, na.rm = TRUE))*sum(xbar.dif*(solve(cont2%*%allS%*%t(cont2)) %*% xbar.dif))
    p.val <- pf(F.val*(sum(c.vector, na.rm = TRUE)-k+1)/(sum(c.vector, na.rm = TRUE)-1)/(k-1), (k-1), (sum(c.vector, na.rm = TRUE)-k+1), 
             lower.tail=F)
    df1 <- (k-1)
    df2 <- (sum(c.vector, na.rm = TRUE)-k+1)
    flatness <- data.frame(F.val,df1,df2,p.val)
    names(flatness) <- c("F","df1","df2","p-value")
    
    #Combine the results from parallelism, equal levels, and flatness
    profile.test <- list(parallel, equal.levels, flatness)
    names(profile.test) <- c("Ho: Profiles are parallel","Ho: Profiles have equal levels","Ho: Profiles are flat")
    
    call <- match.call()
    output <- list(call=call, data.summary=average, corr.table=cor.table, profile.test=profile.test)
    class(output) <- "profg"
    return(output)
  }
}
