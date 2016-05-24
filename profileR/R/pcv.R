#' Cross-Validation for Profile Analysis
#' 
#' Implements the cross-validation described in Davison & Davenport (2002).
#'
#' The \code{pcv} function requires two arguments: criterion and predictor. The criterion corresonds to the dependent variable and the predictor corresponds to the matrix of predictor variables. The function performs the cross-validation technique described in Davison & Davenport (2002) and an object of class \code{critpat} is returned. There the following s3 generic functions are available: \code{summary()},\code{anova()}, \code{print()}, and \code{plot()}. These functions provide a summary of the cross-validation (namely, R2); performs ANOVA of the R2 based on the split for the level, pattern, and overall; provide output similar to \code{lm()}; and plot the estimated parameters for the random split.  Missing data are presently handled by specifying \code{na.action = "na.omit"}, which performs listwise deletion and \code{na.action = "na.fail"}, the default, which causes the function to fail. A seed may also be set for reproducibility by setting the \code{seed}. 
#'
#' @export
#' @importFrom stats anova glm coef cor fitted pf
#' @param formula An object of class \code{\link{formula}} of the form \code{response ~ terms}.
#' @param data An optional data frame, list or environment containing the variables in the model.
#' @param seed Should a seed be set? Function defaults to a random seed.
#' @param na.action How should missing data be handled? Function defaults to failing if missing data are present.
#' @param family A description of the error distribution and link function to be used in the model. See \code{\link{family}}.
#' @param weights An option vector of weights to be used in the fitting process.
#' @return An object of class \code{critpat} is returned, listing the f ollowing components: 
#' \itemize{
#' \item \code{R2.full}, test of the null hypothesis that R2 = 0  
#' \item \code{R2.pat}, test that the R2_pattern = 0
#' \item \code{R2.level}, test that the R2_level = 0
#' \item \code{R2.full.lvl}, test that the R2_full = R2_level = 0
#' \item \code{R2.full.pat}, test that the R2_full = R2_pattern = 0}
#' @references Davison, M., & Davenport, E. (2002). Identifying criterion-related patterns of predictor scores using multiple regression. \emph{ Psychological Methods, 7}(4), 468-484. DOI: 10.1037/1082-989X.7.4.468.
#' @seealso \code{\link{cpa}},\code{\link{print.critpat}},\code{\link{summary.critpat}},\code{\link{anova.critpat}},\code{\link{plot.critpat}} 
#' @keywords methods


pcv <- function(formula, data, seed=NULL, na.action = "na.fail", family = "gaussian", weights = NULL){
  k = 100
  if(is.numeric(seed))
  	set.seed(seed)
  
  index <- 1:nrow(data)
  index.samp <- sample(index,nrow(data)/2)
  
  if(is.null(weights)){
  	regweg.X1<- glm(formula=formula,data=data[index.samp,],family = family,na.action = na.action)
  	regweg.X2<- glm(formula=formula,data=data[-(index.samp),],family = family,na.action = na.action)
  	}
  else {
  	regweg.X1<- glm(formula=formula,data=data[index.samp,],family = family,na.action = na.action,weights=weights)
  	regweg.X2<- glm(formula=formula,data=data[-(index.samp),],family = family,na.action = na.action,,weights=weights)
  }
  X1.b <- coef(regweg.X1)[-1]
  X2.b <- coef(regweg.X2)[-1]
  X1.bstar <- X1.b - mean(X1.b)
  X2.bstar <- X2.b - mean(X2.b)
  X1.xc <- k*X1.bstar # criterion-pattern
  X2.xc <- k*X2.bstar
  
  if(is.null(weights)){
    X1 <- regweg.X1$model[,-1]
    X2 <- regweg.X2$model[,-1]
    }
    else{
    	X1 <- regweg.X1$model[,c(-1,-ncol(regweg.X1$model))]
    	X2 <- regweg.X2$model[,c(-1,-ncol(regweg.X2$model))]
    	}

  k <- 1

  Y1 <- regweg.X1$model[,1]
  Y2 <- regweg.X2$model[,1]
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  V1 <- 1/ncol(X1)
  V2 <- 1/ncol(X2)
  
  
  Xp.X1 <- apply(X1,1,mean)
  Xp.X2 <- apply(X2,1,mean)
  pat.compX1 <- X1 - apply(X1,1,mean)
  pat.compX2 <- X2 - apply(X2,1,mean)
  Covpc.X1 <- V1*(as.matrix(pat.compX1)%*%as.matrix(X2.xc))
  Covpc.X2 <- V2*(as.matrix(pat.compX2)%*%as.matrix(X1.xc))
  R2.lvl.rs1 <- cor(Xp.X1,Y1)^2
  R2.lvl.rs2 <- cor(Xp.X2,Y2)^2
  R2.pat.rs1 <- cor(Covpc.X1,Y1)^2
  R2.pat.rs2 <- cor(Covpc.X2,Y2)^2
  ypredrs1 <- fitted(lm(Y1 ~ 1 + Covpc.X1 + Xp.X1))
  ypredrs2 <- fitted(lm(Y2 ~ 1 + Covpc.X2 + Xp.X2))
  
  R2.f.rs1 <- cor(ypredrs1,Y1)^2
  R2.f.rs2 <- cor(ypredrs2,Y2)^2
  
  df.f1 <- c(1,N1-3)
  df.f2 <- c(1,N2-3)
  
  df.rs1 <- c(1,N1-2)
  df.rs2 <- c(1,N2-2)
  
  
  F.R2.fullrs1 <- (R2.f.rs1*df.f1[2])/((1-R2.f.rs1)*df.f1[1])
  F.R2.fullrs2 <- (R2.f.rs2*df.f2[2])/((1-R2.f.rs2)*df.f2[1])
  p.value.F.R2.fullrs1 <- pf(F.R2.fullrs1,df.f1[1],df.f1[2],lower.tail = FALSE)
  p.value.F.R2.fullrs2 <- pf(F.R2.fullrs2,df.f2[1],df.f2[2],lower.tail = FALSE)
  r2.frs <- rbind(R2.f.rs1,R2.f.rs2)
  r2.frF <- rbind(F.R2.fullrs1,F.R2.fullrs2)
  r2.fr.p <- rbind(p.value.F.R2.fullrs1,p.value.F.R2.fullrs2)
  df.full <- rbind(df.f1,df.f2)
  r2.full <- cbind(r2.frs,df.full, r2.frF,r2.fr.p)
  rownames(r2.full) <- c("Random Sample 1","Random Sample 2")
  colnames(r2.full) <- c("R2","df1","df2","F value","Pr(>F)")
  
  F.R2.patOrs1 <- (R2.pat.rs1*df.rs1[2])/((1-R2.pat.rs1)*df.rs1[1])
  F.R2.patOrs2 <- (R2.pat.rs2*df.rs2[2])/((1-R2.pat.rs2)*df.rs2[1])
  p.value.F.R2.patOrs1 <- pf(F.R2.patOrs1,df.rs1[1],df.rs1[2],lower.tail = FALSE)
  p.value.F.R2.patOrs2 <- pf(F.R2.patOrs2,df.rs2[1],df.rs2[2],lower.tail = FALSE)
  r2.patOrs <- rbind(R2.pat.rs1,R2.pat.rs2)
  r2.pOrF <- rbind(F.R2.patOrs1,F.R2.patOrs2)
  r2.pOr.p <- rbind(p.value.F.R2.patOrs1,p.value.F.R2.patOrs2)
  df.pat <- rbind(df.rs1,df.rs2)
  R2.patO <- cbind(r2.patOrs,df.pat,r2.pOrF,r2.pOr.p)
  rownames(R2.patO) <- c("Random Sample 1","Random Sample 2")
  colnames(R2.patO) <- c("R2","df1","df2","F value","Pr(>F)")
  
  F.R2.lvlOrs1 <- (R2.lvl.rs1*df.rs1[2])/((1-R2.lvl.rs1)*df.rs1[1])
  F.R2.lvlOrs2 <- (R2.lvl.rs2*df.rs2[2])/((1-R2.lvl.rs2)*df.rs2[1])
  p.value.F.R2.lvlOrs1 <- pf(F.R2.lvlOrs1,df.rs1[1],df.rs1[2],lower.tail = FALSE)
  p.value.F.R2.lvlOrs2 <- pf(F.R2.lvlOrs2,df.rs2[1],df.rs2[2],lower.tail = FALSE)
  r2.lvlOrs <- rbind(R2.lvl.rs1,R2.lvl.rs2)
  r2.lOrF <- rbind(F.R2.lvlOrs1,F.R2.lvlOrs2)
  r2.lOr.p <- rbind(p.value.F.R2.lvlOrs1,p.value.F.R2.lvlOrs2)
  df.lvl <- rbind(df.rs1,df.rs2)
  R2.lvlO <- cbind(r2.lvlOrs,df.lvl, r2.lOrF,r2.lOr.p)
  rownames(R2.lvlO) <- c("Random Sample 1","Random Sample 2")
  colnames(R2.lvlO) <- c("R2","df1","df2", "F value","Pr(>F)")
  
  F.R2.patrs1 <- ((R2.f.rs1 - R2.lvl.rs1)*df.rs1[2])/((1-R2.f.rs1)*df.rs1[1])
  F.R2.patrs2 <- ((R2.f.rs2 - R2.lvl.rs2)*df.rs2[2])/((1-R2.f.rs2)*df.rs2[1])
  p.value.F.R2.patrs1 <- pf(F.R2.patrs1,df.rs1[1],df.rs1[2],lower.tail = FALSE)
  p.value.F.R2.patrs2 <- pf(F.R2.patrs2,df.rs2[1],df.rs2[2],lower.tail = FALSE)
  r2.prs <- rbind(R2.pat.rs1,R2.pat.rs2)
  r2.prF <- rbind(F.R2.patrs1,F.R2.patrs2)
  r2.pr.p <- rbind(p.value.F.R2.patrs1,p.value.F.R2.patrs2)
  r2.pattern <- cbind(r2.prs,r2.prF,df.pat,r2.pr.p)
  rownames(r2.pattern) <- c("Random Sample 1","Random Sample 2")
  colnames(r2.pattern) <- c("R2","F value","df1","df2","Pr(>F)")
  
  F.R2.lvl.rs1 <- ((R2.f.rs1 - R2.pat.rs1)*df.rs1[2])/((1-R2.f.rs1))
  F.R2.lvl.rs2 <- ((R2.f.rs2 - R2.pat.rs2)*df.rs2[2])/((1-R2.f.rs2))
  p.value.F.R2.lvl.rs1 <- pf(F.R2.lvl.rs1,df.rs1[1],df.rs1[2],lower.tail=FALSE)
  p.value.F.R2.lvl.rs2 <- pf(F.R2.lvl.rs2,df.rs2[1],df.rs2[2],lower.tail=FALSE)
  r2.lrs <- rbind(R2.lvl.rs1,R2.lvl.rs2)
  r2.lrF <- rbind(F.R2.lvl.rs1,F.R2.lvl.rs2)
  df.lvl <- rbind(df.rs1,df.rs2)
  r2.lr.p <- rbind(p.value.F.R2.lvl.rs1,p.value.F.R2.lvl.rs2)
  r2.level <- cbind(r2.lrs,df.lvl,r2.lrF,r2.lr.p)
  rownames(r2.level) <- c("Random Sample 1","Random Sample 2")
  colnames(r2.level) <- c("R2","df1","df2","F value","Pr(>F)")
  r2.full <- round(r2.full,digits=6)
  R2.patO <- round(R2.patO,digits=6)
  R2.lvlO <- round(R2.lvlO,digits=6)
  r2.pattern <- round(r2.pattern,digits=6)
  r2.level <- round(r2.level,digits=6)
  b = list(X1.b,X2.b)
  names(b) = c("Random Sample 1", "Random Sample 2")
  call<- match.call()
  ftable <- list(r2.full, R2.patO,R2.lvlO,r2.pattern,r2.level)
  r2 <- list(r2.full[,1],R2.patO[,1],R2.lvlO[,1],r2.pattern[,1],r2.level[,1])
  names(ftable) <- c("R2.full = 0", "R2.pat = 0", "R2.lvl = 0", "R2.full = R2.lvl", "R2.full = R2.pat")
  names(r2) <- names(ftable)
  
  output = list(call=call,b = b, ftable = ftable,r2 = r2)
  class(output) <- "critpat"
  return(output)
}




