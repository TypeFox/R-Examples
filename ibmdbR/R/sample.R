# 
# Copyright (c) 2013, 2014, IBM Corp. All rights reserved. 
# 		
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
#

idaSample <- function(bdf, n, stratCol=NULL,stratVals=NULL,stratProbs=NULL,dbPreSamplePercentage=100,fetchFirst=F) {
  
  idaCheckConnection();
  if (!is.ida.data.frame(bdf))
    stop("idaSample can only be applied to objects of type ida.data.frame.")
  
  
  #Check if there are any generated columns or a where condition
  if(((length(idadf.defined.columns(bdf))>0)|nchar(bdf@where))&(dbPreSamplePercentage<100)) {
    stop("Presampling is not supported for a ida.data.frame with constructed columns or selection.")
  }
  
  cols <- paste('"', bdf@cols,'"', sep='');	
  
  if(!is.null(stratCol)) {
    
    stats <- idaQuery("SELECT \"",stratCol,"\" AS GRP,COUNT(*) AS CNT FROM ",idadf.from(bdf)," ",ifelse(nchar(bdf@where),paste(" WHERE ",bdf@where,sep=''),'')," GROUP BY \"",stratCol,"\"");
    stats$CNT <- as.numeric(stats$CNT)
    
    if(!is.null(stratVals)) {
      stats$CNT <- ifelse(stats$GRP %in% stratVals,stats$CNT,0);
      for(i in 1:length(stratVals)) {
        if(nrow(stats[stats$GRP==stratVals[i],])==0)
          warning(paste("Value ",stratVals[i]," not in column."))
      }
    }
    
    if(is.null(stratProbs)) {
      numRows <- sum(stats$CNT)
      fraction <- n/numRows;
      stats$n <- as.integer(ceiling(stats$CNT*fraction));
    } else {
      
      if(is.null(stratVals))
        stop("You need to specify stratVals along with the probabilities");
      
      if(length(stratVals) != length(stratProbs))
        stop("The length of the vector stratVals and stratProbs needs to be the same")
      
      if(abs(sum(stratProbs)-1.0)>0.001)
        warning("The probabilities do not sum up to one, ")
      
      stratMap <- list();
      for(i in 1:length(stratVals)) {
        stratMap[stratVals[i]] <- stratProbs[i];
      }
      
      for(i in 1:nrow(stats)) {			
        stats$n[i] <- as.integer(ifelse(stats$GRP[i] %in% stratVals,ceiling(stratMap[[stats$GRP[i]]]*n),0));
      }
  
      for(i in 1:nrow(stats)) {
        if(stats$CNT[i]<stats$n[i]) {
          warning(paste("Not enough samples for ", stats$GRP[i], ": required ",stats$n[i], " available ", stats$CNT[i], sep=''))
        }	
      }
    }
    
    queries <- c();
    if(nrow(stats)>50) {
      stop("Only a maximum of 50 different strata allowed at the moment.")
    }
    for(i in 1:nrow(stats)) {
      if(stats$n[i]>0) {
        
        if(!fetchFirst) {
          query <- paste("(SELECT ",paste(cols, collapse=","),",RAND() AS P_RAND_INDEX FROM ",idadf.from(bdf)," A ",ifelse(dbPreSamplePercentage<100,paste(" TABLESAMPLE BERNOULLI(",dbPreSamplePercentage,") ",sep=''),"")," WHERE \"",stratCol,"\" = '",stats$GRP[i],"'", 
              ifelse(nchar(bdf@where),paste(" AND ",bdf@where,sep=''),'')," ORDER BY P_RAND_INDEX FETCH FIRST ",format(stats$n[i], scientific = FALSE), " ROWS ONLY)",sep='');
        } else  {
          query <- paste("(SELECT ",paste(cols, collapse=",")," FROM ",idadf.from(bdf)," A  WHERE \"",stratCol,"\" = '",stats$GRP[i],"'", 
              ifelse(nchar(bdf@where),paste(" AND ",bdf@where,sep=''),'')," FETCH FIRST ",format(stats$n[i], scientific = FALSE), " ROWS ONLY)",sep='');
        }
        
        queries <- c(queries,query);
        
      }
    }
    result <- idaQuery(paste(queries,collapse=" UNION "))
    
  } else {
    sampleSize <- as.integer(n);
    
    if(!fetchFirst) {
      result <- idaQuery("SELECT ",paste( cols, collapse=","),",RAND() AS P_RAND_INDEX FROM ",idadf.from(bdf),ifelse(dbPreSamplePercentage<100,paste(" TABLESAMPLE BERNOULLI(",dbPreSamplePercentage,") ",sep=''),""), ifelse(nchar(bdf@where),paste(" WHERE ",bdf@where,sep=''),'')," ORDER BY P_RAND_INDEX FETCH FIRST ",format(sampleSize, scientific = FALSE), " ROWS ONLY");
    } else {
      result <- idaQuery("SELECT ",paste( cols, collapse=",")," FROM ",idadf.from(bdf),ifelse(nchar(bdf@where),paste(" WHERE ",bdf@where,sep=''),'')," FETCH FIRST  ",format(sampleSize, scientific = FALSE), " ROWS ONLY");
    }
  }
  result$P_RAND_INDEX <- NULL;
  result;	
}