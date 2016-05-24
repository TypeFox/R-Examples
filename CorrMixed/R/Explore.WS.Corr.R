Explore.WS.Corr <- function(OLS.Model=" ", Dataset, Id, Time, Alpha=0.05, 
                            Smoother.Span=.2, Number.Bootstrap=100, Seed=1){
  
  F_val <- Smoother.Span  
  Max.Time <- length(unique(Dataset[,Time]))
  
  options(warn=-1)
  
  Id_val <- Dataset[, Id]  
  
  Res <- lm(formula = OLS.Model, data=Dataset)$residuals
  
  Dataset.Short <-
    data.frame(cbind(Dataset[,colnames(Dataset)%in%Id], Res, Dataset[,colnames(Dataset)%in%Time]))
  names(Dataset.Short) <- c("Id", "Res", "Time")
  
  Data_Wide <- reshape(data = Dataset.Short, timevar= "Time", idvar= "Id", direction="wide")
  Data_Wide <- Data_Wide[, 2:(1+Max.Time)]
  
  cor <- cor(Data_Wide, use = "pairwise.complete.obs")
  cor <- matrix(cor, nrow = nrow(cor))
  
  alles <- NULL
  for (i in 1: max(Max.Time)){
    cor_hier <- cor[row(cor) == (col(cor) - i)]
    try(erbij <- cbind(i, cor_hier), silent=TRUE)    
    try(alles <- rbind(alles, erbij), silent=TRUE)
    rm(erbij)
  } 
  
  Est.Corr <- na.exclude(unique(lowess(x = alles[,1], y=alles[,2], f=F_val)$y))
  
  ### Bootstrap CI
  all_y <- Est.Corr
  Dataset <- cbind(Dataset, Dataset[, Id]) 
  colnames(Dataset[dim(Dataset)[2]]) <- c("Id_name")
  
  for (i in 1: Number.Bootstrap){
    set.seed(Seed+i)
    boot <- 
      sample(unique(Dataset[,colnames(Dataset)%in%Id]), 
             size = length(unique(Dataset[,colnames(Dataset)%in%Id])), replace = T)
    
    sample.boot <- NULL
    for (j in 1: length(boot)){
      Dataset <- data.frame(Dataset)
      samen <- 
        Dataset[Dataset[dim(Dataset)[2]]==boot[j], ]
      unit <- rep(j, dim(samen)[1])
      samen <- cbind(samen, unit)
      sample.boot <- rbind(sample.boot, samen) 
    }
    
    Res_hier <-  lm(formula = OLS.Model, data=sample.boot)$residuals
    
    Dataset.Short <-
      data.frame(cbind(sample.boot[,dim(sample.boot)[2]], Res_hier, 
                       sample.boot[,colnames(sample.boot)%in%Time]))
    names(Dataset.Short) <- c("Unit", "Res", "Time")
    
    Data_Wide <- reshape(data = Dataset.Short, timevar= "Time", 
                         idvar= "Unit", direction="wide")
    
    cor <- cor(Data_Wide[,2: dim(Data_Wide)[2]], use = "pairwise.complete.obs")
    cor <- matrix(cor, nrow = nrow(cor))
    
    alles <- NULL
    for (i in 1: Max.Time){
      cor_hier <- cor[row(cor) == (col(cor) - i)]
      try(erbij <- cbind(i, cor_hier), silent=TRUE)    
      try(alles <- rbind(alles, erbij), silent=TRUE)
      rm(erbij)
    } 
    
    Y_hier <- unique(lowess(x = alles[,1], y=alles[,2], f=F_val)$y)
    Y_hier[Y_hier>c(1)] <- 1; Y_hier[Y_hier< c(-1)] <- -1
    
    if (length(Y_hier) != length(Est.Corr)){
      diff = length(Est.Corr) - length(Y_hier)
      if (diff > 0){Y_hier <- c(Y_hier, rep(NA, times=diff) )}
      if (diff < 0){Y_hier <- Y_hier <- Y_hier[1: length(Est.Corr)]}
    }
    all_y <- rbind(all_y, Y_hier)
  }
  
  all_y <- all_y[2: dim(all_y)[1],]
  all_y <- (matrix(data = all_y, nrow=nrow(all_y)))
  
  ci <- NULL
  for (i in 1: dim(all_y)[2]){
    ci <- cbind(ci, t(t(quantile(all_y[,i], probs = c(Alpha/2, (1-(Alpha/2))), 
                                 na.rm = T)))) 
  }
  
  CI.Upper <- ci[2,]
  CI.Lower <- ci[1,]
  CI.Upper[CI.Upper > 1] <- 1
  CI.Lower[CI.Lower < -1] <- -1
  colnames(alles) <- c("Time_lag", "R")
  
  fit <- 
    list(Est.Corr=Est.Corr, All.Corrs=alles, Bootstrapped.Corrs=all_y, Alpha=Alpha, 
         CI.Upper=CI.Upper, CI.Lower=CI.Lower, Call=match.call())   
  
  class(fit) <- "Explore.WS.Corr"
  fit
}