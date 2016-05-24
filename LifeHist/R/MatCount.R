MatCount <-
function(matdat,fem.key,mal.key,stage.key,season.key) 
             {
              MatDat                 <- data.frame(pred=round(matdat[,1]),
                                                   sex=matdat[,2],
                                                   season=ifelse(matdat[,3] < season.key[1] | matdat[,3] > season.key[2], 0, 1),                                          
                                                   mat=ifelse(matdat[,4] <= stage.key,0,1));
              MatDat.mal.nrep        <- MatDat[MatDat$sex==mal.key & MatDat$season==0,c(1,4)];
              MatDat.mal.rep         <- MatDat[MatDat$sex==mal.key & MatDat$season==1,c(1,4)];
              MatDat.fem.nrep        <- MatDat[MatDat$sex==fem.key & MatDat$season==0,c(1,4)];
              MatDat.fem.rep         <- MatDat[MatDat$sex==fem.key & MatDat$season==1,c(1,4)];
              #Males Non-Reproductive season
              juv                    <- aggregate(MatDat.mal.nrep$mat[MatDat.mal.nrep$mat==0],list(MatDat.mal.nrep$pred[MatDat.mal.nrep$mat==0]),length);
              names(juv)             <- c("Pred","Count"); 
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(juv$Pred),unique(pred.series$Pred)),2] <- juv$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              juv                    <- pred.series;
              juv$matkey             <- 0;
              rm(pred.series);
              mat                    <- aggregate(MatDat.mal.nrep$mat[MatDat.mal.nrep$mat==1],list(MatDat.mal.nrep$pred[MatDat.mal.nrep$mat==1]),length);
              names(mat)             <- c("Pred","Count");
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(mat$Pred),unique(pred.series$Pred)),2] <- mat$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              mat                    <- pred.series;
              mat$matkey             <- 1; 
              rm(pred.series);
              MatDat.mal.nrep         <- rbind(juv,mat);
              rm(juv,mat);
              MatDat.mal.nrep         <- MatDat.mal.nrep[do.call(order,list(MatDat.mal.nrep$Pred,MatDat.mal.nrep$matkey)),];
              MatDat.mal.nrep         <- MatDat.mal.nrep[,-3];
              MatDat.mal.nrep         <- data.frame(pred=MatDat.mal.nrep$Pred,
                                                   mat=rep.int(x=c(0,1),times=dim(MatDat.mal.nrep)[1]/2),
                                                   juv=rep.int(x=c(1,0),times=dim(MatDat.mal.nrep)[1]/2),
                                                   count=MatDat.mal.nrep$Count);
              class(MatDat.mal.nrep)  <- "MatData";
              #Males Reproductive season
              juv                    <- aggregate(MatDat.mal.rep$mat[MatDat.mal.rep$mat==0],list(MatDat.mal.rep$pred[MatDat.mal.rep$mat==0]),length);
              names(juv)             <- c("Pred","Count"); 
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(juv$Pred),unique(pred.series$Pred)),2] <- juv$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              juv                    <- pred.series;
              juv$matkey             <- 0;
              rm(pred.series);
              mat                    <- aggregate(MatDat.mal.rep$mat[MatDat.mal.rep$mat==1],list(MatDat.mal.rep$pred[MatDat.mal.rep$mat==1]),length);
              names(mat)             <- c("Pred","Count");
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(mat$Pred),unique(pred.series$Pred)),2] <- mat$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              mat                    <- pred.series;
              mat$matkey             <- 1; 
              rm(pred.series);
              MatDat.mal.rep        <- rbind(juv,mat);
              rm(juv,mat);
              MatDat.mal.rep        <- MatDat.mal.rep[do.call(order,list(MatDat.mal.rep$Pred,MatDat.mal.rep$matkey)),];
              MatDat.mal.rep        <- MatDat.mal.rep[,-3];
              MatDat.mal.rep        <- data.frame(pred=MatDat.mal.rep$Pred,
                                                   mat=rep.int(x=c(0,1),times=dim(MatDat.mal.rep)[1]/2),
                                                   juv=rep.int(x=c(1,0),times=dim(MatDat.mal.rep)[1]/2),
                                                   count=MatDat.mal.rep$Count);
              class(MatDat.mal.rep)  <- "MatData";
              #Females Non-Reproductive season
              juv                    <- aggregate(MatDat.fem.nrep$mat[MatDat.fem.nrep$mat==0],list(MatDat.fem.nrep$pred[MatDat.fem.nrep$mat==0]),length);
              names(juv)             <- c("Pred","Count"); 
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(juv$Pred),unique(pred.series$Pred)),2] <- juv$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              juv                    <- pred.series;
              juv$matkey             <- 0;
              rm(pred.series);
              mat                    <- aggregate(MatDat.fem.nrep$mat[MatDat.fem.nrep$mat==1],list(MatDat.fem.nrep$pred[MatDat.fem.nrep$mat==1]),length);
              names(mat)             <- c("Pred","Count");
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(mat$Pred),unique(pred.series$Pred)),2] <- mat$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              mat                    <- pred.series;
              mat$matkey             <- 1; 
              rm(pred.series);
              MatDat.fem.nrep         <- rbind(juv,mat);
              rm(juv,mat);
              MatDat.fem.nrep         <- MatDat.fem.nrep[do.call(order,list(MatDat.fem.nrep$Pred,MatDat.fem.nrep$matkey)),];
              MatDat.fem.nrep         <- MatDat.fem.nrep[,-3];
              MatDat.fem.nrep         <- data.frame(pred=MatDat.fem.nrep$Pred,
                                                   mat=rep.int(x=c(0,1),times=dim(MatDat.fem.nrep)[1]/2),
                                                   juv=rep.int(x=c(1,0),times=dim(MatDat.fem.nrep)[1]/2),
                                                   count=MatDat.fem.nrep$Count);
              class(MatDat.fem.nrep)  <- "MatData";
              #Females Reproductive season
              juv                    <- aggregate(MatDat.fem.rep$mat[MatDat.fem.rep$mat==0],list(MatDat.fem.rep$pred[MatDat.fem.rep$mat==0]),length);
              names(juv)             <- c("Pred","Count"); 
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(juv$Pred),unique(pred.series$Pred)),2] <- juv$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              juv                    <- pred.series;
              juv$matkey             <- 0;
              rm(pred.series);
              mat                    <- aggregate(MatDat.fem.rep$mat[MatDat.fem.rep$mat==1],list(MatDat.fem.rep$pred[MatDat.fem.rep$mat==1]),length);
              names(mat)             <- c("Pred","Count");
              pred.seq               <- min(MatDat$pred):max(MatDat$pred);
              pred.series            <- data.frame(pred.seq, NA);
              names(pred.series)     <- c("Pred","Count");
              pred.series[match(unique(mat$Pred),unique(pred.series$Pred)),2] <- mat$Count;
              pred.series$Count[is.na(pred.series$Count)] <- 0;
              mat                    <- pred.series;
              mat$matkey             <- 1; 
              rm(pred.series);
              MatDat.fem.rep        <- rbind(juv,mat);
              rm(juv,mat);
              MatDat.fem.rep        <- MatDat.fem.rep[do.call(order,list(MatDat.fem.rep$Pred,MatDat.fem.rep$matkey)),];;
              MatDat.fem.rep        <- MatDat.fem.rep[,-3];
              MatDat.fem.rep        <- data.frame(pred=MatDat.fem.rep$Pred,
                                                   mat=rep.int(x=c(0,1),times=dim(MatDat.fem.rep)[1]/2),
                                                   juv=rep.int(x=c(1,0),times=dim(MatDat.fem.rep)[1]/2),
                                                   count=MatDat.fem.rep$Count);
              class(MatDat.fem.rep)  <- "MatData";
              #
              MatDat                 <- list(mal.nrep=MatDat.mal.nrep,mal.rep=MatDat.mal.rep,fem.nrep=MatDat.fem.nrep,fem.rep=MatDat.fem.rep);
              rm(MatDat.mal.nrep,MatDat.mal.rep,MatDat.fem.nrep,MatDat.fem.rep);
              return(MatDat);
              }
