MDR.sing.mod <-
function(x,y,order=NULL,trace=NULL,...){
  if(is.null(trace)) trace <- FALSE
  if(is.null(order)) order <- 3
  res <- MDR.sing.mod1(x=x,y=y,order=order,trace=trace)
  return(res)
}


MDR.sing.mod1 <-
function(x,y,order,trace,...){	
  n <- ncol(x)
  index <- t(combn(n,order))

  n <- nrow(index)
  c <- ncol(index)
  #MDR.R2 <- numeric()
  MDR.AdjR2 <- numeric()
  MDR.Pvalue <- numeric() 
  #MLR.R2 <- numeric()
  MLR.AdjR2 <- numeric()
  MLR.Pvalue <- numeric() 
  Pvalue <- numeric()

  for(i in 1:n){
	if(trace)cat("Stage 1, Step:",i,"/",n,"\n")

	X1 <- x[,index[i,]]
	ok <- complete.cases(X1) ## Missing 3/3/2014	
	colnames(X1)= NULL
	X1 <- data.frame(X=X1)
	X1 <- X1[ok,]
	y1 <- matrix(y,,1)
	y1 <- y1[ok,]
        
	
	#### Without Interaction
        reg1 <- lm(y1~.,data=X1) 
        a <- summary(reg1)
        MLR.R2 <- a$r.squared
        nR <- nrow(X1)
        pR <- ncol(X1)
        MLR.AdjR2[i] <- 1-(1-MLR.R2)*(nR-1)/(nR-pR-1) # Adjusted R square 

        ANOVA <- aov(y1~.,data=X1)
        A <- summary(ANOVA)
        P <- A[[1]][[5]]
        s <- length(P)
        P <- P[-s]
        MLR.Pvalue[i] <- max(P)
	  

	#### With Interaction
	if(c>1){
	   X2 <- X1[,1]
           for(j in 1:(c-1)){X2<-paste(X2,X1[,j+1],sep=":")}    
	   X2 <- factor(X2)

	   reg2 <- lm(y1~X2)
	   B <- summary(reg2)
           MDR.R2 <- B$r.squared
           nR <- nrow(X1)
           pR <- 1
           MDR.AdjR2[i] <- 1-(1-MDR.R2)*(nR-1)/(nR-pR-1) # Adjusted R square
   
	   ANOVA <- aov(y1~X2)
	   A <- summary(ANOVA)
	   P <- A[[1]][[5]]
	   s <- length(P)
	   P <- P[-s]
	   MDR.Pvalue[i] <- max(P)		

	   #### Compare These Tow Models and Extract the P-value
	   Comp <- anova(reg1,reg2)
	   Pvalue[i] <- Comp$"Pr(>F)"[2]
	} 
  }

  #res = cbind(MDR_R2 = MDR.R2,MDR_Adj_R2 = MDR.AdjR2, MDR_Pvalue = MDR.Pvalue,
  #            MLR_R2 = MLR.R2,MLR_Adj_R2 = MLR.AdjR2, MLR_Pvalue = MLR.Pvalue,
  #            Model_Pvalue = Pvalue)
  res = cbind(index,MDR_R2 = MDR.AdjR2, MDR_Pvalue = MDR.Pvalue,
              MLR_R2 = MLR.AdjR2, MLR_Pvalue = MLR.Pvalue,
              Model_Pvalue = Pvalue)
  return(res)

}
