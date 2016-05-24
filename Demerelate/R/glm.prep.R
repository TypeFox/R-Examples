glm.prep <- function(empirical.list, offfull.list, offhalf.list, offnon.list)
{
  

  # Preparing output data.frame

	  non <- data.frame(rep("NON",length(offnon.list[!is.na(offnon.list)])),offnon.list[!is.na(offnon.list)])
    hs <- data.frame(rep("HS",length(offhalf.list[!is.na(offhalf.list)])),offhalf.list[!is.na(offhalf.list)])
    fs <- data.frame(rep("FS",length(offfull.list[!is.na(offfull.list)])),offfull.list[!is.na(offfull.list)])
		rel.Y <- c(as.character(non[,1]),as.character(hs[,1]),as.character(fs[,1]))
		rel.X <- c(as.numeric(non[,2]),as.numeric(hs[,2]),as.numeric(fs[,2]))
		relate.data <- data.frame(as.factor(rel.Y),as.numeric(rel.X))
		names(relate.data) <- c("Sib","Mxy")

	
	# Multinominal logistic regression
	
  	redata <- mlogit.data(relate.data,varying=NULL,choice="Sib",shape="wide")
  	mlogit.model <- mlogit(Sib~1|Mxy, data=redata, reflevel="NON")
  	sumlrm <- summary(mlogit.model)

	# Prepares data for plotting mlr
  	empirical <- as.vector(empirical.list[!is.na(empirical.list)])
  	logits1 <- exp(rep(0,length(empirical)))
  	nonlog2 <- (log(1)-sumlrm[[1]][1])/sumlrm[[1]][3]
  	nonlog3 <- (log(1)-sumlrm[[1]][2])/sumlrm[[1]][4]
  	half <- min(nonlog2,nonlog3)
		
  	
    return(list(half, sumlrm))

}

      
      
