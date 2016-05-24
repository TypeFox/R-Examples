assign("dblm",
       function(data, y, sc=0.003, ev.min=0.007, ...){
         Delta <- gowdis(data, ...)
         mds <- cmdscale(sqrt(Delta), k = nrow(data)-1, eig = TRUE)
         m <- sum(mds$eig > ev.min)
         mds <- cmdscale(sqrt(Delta), k = m, eig = TRUE)
         X <- mds$points
         eigenvalues <- mds$eig
         SquareCor <- as.vector(cor(y,X)^2)                   
         Percent.Iner <- eigenvalues/sum(eigenvalues) 
         o<-data.frame(1:length(SquareCor),round(eigenvalues[1:length(SquareCor)],10),round(SquareCor,10),round(Percent.Iner[1:length(SquareCor)],10))
         names(o)<-c("ID", "Eigenvalues.cp", "Sq.Cor", "PercentIner")
         o1<-o[o$Sq.Cor>sc,]                                          
         Xr <- X[,o1$ID]
         rdb <- lm(y ~ Xr)                                     
         model.db <- summary(rdb)
         model.db1 <- model.db["coefficients"]$coefficients[-1,][model.db["coefficients"]$coefficients[-1,4]<0.05,]
         lv <- length(rownames(model.db["coefficients"]$coefficients[-1,][model.db["coefficients"]$coefficients[-1,4]<0.05,]))
         an <- as.numeric(unlist(strsplit(rownames(model.db["coefficients"]$coefficients[-1,][model.db["coefficients"]$coefficients[-1,4]<0.05,]), "[r]"))[seq(2,lv*2,2)])
         Xr.0 <- Xr[,an]
         db.lm <- list(table=o1,ev=eigenvalues,cp=Xr.0,dbmodel=model.db)
         db.lm
       }
)
