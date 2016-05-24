modeval <- function (calculated, measured,stat=c("N","pearson","MBE","RMBE","MAE","RMAE","RMSE","RRMSE","R2","slope","intercept","EF","SD","CRM","MPE","AC","ACu","ACs"),minlength=4) 
{
  stopifnot(length(calculated)==length(measured), length(stat)!=0)
  rval <- list()
  for (st in stat) rval[[st]] <- NA 
  x <- na.omit(as.data.frame(cbind(calculated, measured)))
  if (dim(x)[1]>minlength) {    
    calculated <- x[, 1]
    measured <- x[, 2]
    
    if("N"%in%stat) rval$N <- nrow(x)
    
    dif <- calculated - measured
    mdif <- mean(dif, na.rm = T)
    sdif <- dif^2
    msdif <- mean(sdif, na.rm = T)
    mmeasured <- mean(measured, na.rm = T)
    mcalculated <- mean(calculated, na.rm = T)  
    
    if("pearson"%in%stat) rval$pearson <- cor(calculated, measured)
    if("MBE"%in%stat) rval$MBE <- mdif
    if("RMBE"%in%stat) rval$RMBE <- 100 * mdif/mmeasured    
    if("MAE"%in%stat) rval$MAE <- mean(abs(dif), na.rm = T)
    if("RMAE"%in%stat) rval$RMAE <- 100 * mean(abs(dif), na.rm = T)/mmeasured 
    if("RMSE"%in%stat) rval$RMSE <- sqrt(msdif)
    if("RRMSE"%in%stat) rval$RRMSE <- 100 * sqrt(msdif)/mmeasured
    
    if(TRUE%in%(c("R2","slope","intercept")%in%stat)) m <- lm(calculated ~ measured)
    if("R2"%in%stat) rval$R2 <- summary(m)$r.squared
    if("slope"%in%stat) rval$slope <- coef(m)[2]
    if("intercept"%in%stat) rval$intercept <- coef(m)[1]
    
    if("EF"%in%stat)rval$EF <- 1 - (sum(sdif)/sum((measured - mmeasured)^2))
    if("SD"%in%stat)rval$SD <- sd(dif, na.rm = T)
    if("CRM"%in%stat)rval$CRM <- (mean(measured - calculated, na.rm = T))/mmeasured
    if("MPE"%in%stat)rval$MPE <- mean(dif/measured, na.rm = T) * 100
    
    if(TRUE%in%(c("AC","ACs","ACu")%in%stat)){
      SSD <- sum(sdif)
      SPOD <- sum((abs(mcalculated - mmeasured) + abs(calculated - mcalculated)) * (abs(mcalculated - mmeasured) + abs(measured - mmeasured)))    
      b <- sqrt((sum((measured-mean(measured))^2))/(sum((calculated-mean(calculated))^2)))
      if (!is.na(cor(calculated,measured))  & cor(calculated,measured) < 0) b <- -b   
      a <-  mmeasured - b * mcalculated
      dY <- a + b*calculated 
      dX <- -a/b + (1/b)*measured   
      SPDu <- sum((abs(calculated-dX))*(abs(measured-dY)))
      SPDs <- SSD-SPDu }
    if("AC"%in%stat)rval$AC <- 1 - (SSD/SPOD)
    if("ACu"%in%stat)rval$ACu <- 1 - (SPDu/SPOD)
    if("ACs"%in%stat)rval$ACs <- 1 - (SPDs/SPOD)
  }
  rval
}

