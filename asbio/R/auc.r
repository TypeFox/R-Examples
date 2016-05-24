auc <- function(obs, fit, plot=FALSE){
compf <- outer(fit[obs==0],fit[obs==1],function(x,y)paste(x,y,sep=","))
cf<-as.data.frame(unlist(strsplit(compf,",")))

n <- prod(summary(factor(obs)))

seq1 <- seq(1,2*n,by = 2)
seq2 <- seq(2,2*n,by = 2)

cf1 <- rep(NA,length(seq1)); cf2 <- rep(NA,length(seq1))

for(i in 1:length(seq1)){
cf1[i] <- as.numeric(levels(cf[[1]][seq1[i]])[cf[[1]][seq1[i]]])
cf2[i] <- as.numeric(levels(cf[[1]][seq2[i]])[cf[[1]][seq2[i]]])
}

U <- rep(NA, length(cf2))
for(i in 1:length(cf2)) U[i] <- ifelse(cf2[i] > cf1[i], 1, 0)


if(plot==TRUE){
        
        int <- seq(0.01, 0.99, by =.01)
        sens <- rep(NA, length(int))
        omspec <- rep(NA, length(int))
       
        for(i in 1:length(int)){
        t <- table(ifelse(fit > int[i], 1, 0), obs)
                    
        sens[i] <- t[,2][2]/sum(t[,2])
        omspec[i] <- 1 - t[,1][1]/sum(t[,1])
        }
        
        plot(c(1,omspec,0), c(1,sens,0), xlab = "False positive rate", ylab = "True positive rate", type = "l", xlim=c(0,1), ylim=c(0,1), lwd = 2.5, main = "ROC curve",cex.lab=1.3,cex.main=1.3)
        abline(0, 1, lty = 2)
        }
sum(U)/n
}



