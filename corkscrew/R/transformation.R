transformation <-
function(x,y,data)
{
	Transformation <- NULL
	Variable <- NULL
    correlation <- data.frame()
    out <- data.frame()
    z <- data.frame()
    for (i in seq(1,length(x),1)){
      print(paste("Working on variable #",i,sep=""))

      if((class(data[,y])=="numeric" | class(data[,y])=="integer") & var(data[,y],na.rm=T) > 0)
      {
        if((class(data[,x[i]])=="numeric" | class(data[,x[i]])=="integer") & var(data[,x[i]],na.rm=T) > 0)
        {
          correlation[1,1] <- "linear"
          correlation[2,1] <- "log"
          correlation[3,1] <- "power_0.1"
          correlation[4,1] <- "power_0.3"
          correlation[5,1] <- "power_0.5"
          correlation[6,1] <- "power_0.7"
          correlation[7,1] <- "power_0.9"
          correlation[8,1] <- "arctangent"
          
          if(min(data[,x[i]],na.rm=T) > 0)
          {
            correlation[1,2] <- cor(data[,x[i]],data[,y],use="na.or.complete")
            correlation[2,2] <- cor(log(data[,x[i]]),data[,y],use="na.or.complete")
            correlation[3,2] <- cor((data[,x[i]]^0.1),data[,y],use="na.or.complete")
            correlation[4,2] <- cor((data[,x[i]]^0.3),data[,y],use="na.or.complete")
            correlation[5,2] <- cor((data[,x[i]]^0.5),data[,y],use="na.or.complete")
            correlation[6,2] <- cor((data[,x[i]]^0.7),data[,y],use="na.or.complete")
            correlation[7,2] <- cor((data[,x[i]]^0.9),data[,y],use="na.or.complete")
            correlation[8,2] <- cor((1/(1+(exp(-(data[,x[i]]-((max(data[,x[i]],na.rm=T)+ min(data[,x[i]],na.rm=T))/2)))))),data[,y],use="na.or.complete")
          }else{
            correlation[1,2] <- cor(data[,x[i]],data[,y],use="na.or.complete")
            correlation[2,2] <- cor(log(data[,x[i]] + (abs(min(data[,x[i]],na.rm=T))+1)),data[,y],use="na.or.complete")
            correlation[3,2] <- cor(((data[,x[i]] + (abs(min(data[,x[i]],na.rm=T))+1))^0.1),data[,y],use="na.or.complete")
            correlation[4,2] <- cor(((data[,x[i]] + (abs(min(data[,x[i]],na.rm=T))+1))^0.3),data[,y],use="na.or.complete")
            correlation[5,2] <- cor(((data[,x[i]] + (abs(min(data[,x[i]],na.rm=T))+1))^0.5),data[,y],use="na.or.complete")
            correlation[6,2] <- cor(((data[,x[i]] + (abs(min(data[,x[i]],na.rm=T))+1))^0.7),data[,y],use="na.or.complete")
            correlation[7,2] <- cor(((data[,x[i]] + (abs(min(data[,x[i]],na.rm=T))+1))^0.9),data[,y],use="na.or.complete")
            correlation[8,2] <- cor((1/(1+(exp(-(data[,x[i]]-((max(data[,x[i]],na.rm=T)+ min(data[,x[i]],na.rm=T))/2)))))),data[,y],use="na.or.complete")
            
            warning(sprintf("Variable %s is adjusted with an offset (making the minimum to 1 due to negative values) for computing correlations using log and power 
transformations",x[i]))
            
          }
          correlation[,3] = x[i]
          
          out = rbind(out,correlation)
          
          z[i,"Variables"] <- x[i]
          z[i,"Transformations"] <- correlation[which(abs(correlation$V2) == max(abs(correlation$V2), na.rm = T)),1]
          z[i,"Correlation"] <- ifelse(abs(max(correlation$V2, na.rm = T))>=abs(min(correlation$V2, na.rm = T)),max(correlation$V2, na.rm = T),min(correlation$V2, 
na.rm = T))
          z[i,"Influence"] <- ifelse(z[i,"Correlation"]>0,"Positive","Negative")
          
        }else{
          warning(sprintf("Variable %s is not numeric / integer or has 0 variance",x[i]))
        }
      }
      else{
        warning("y variable is not numeric / integer or has 0 variance")
      }
    }
    out = out[,c(3,1,2)]
    names(out) = c("Variable","Transformation","cor")
    print(qplot(x=Variable, y=Transformation, data=out, fill=cor, geom="tile")+scale_fill_gradient2(limits=c(-1, 1),low="red4",high="seagreen4",mid="white"))
    return(list(z,out))
}
