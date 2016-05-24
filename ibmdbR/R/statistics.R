# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
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

################ idaTable ############################
idaTable<-function (idadf,max.entries=1000) 
{
  
  
  if (!is.ida.data.frame(idadf)) 
    stop(simpleError("idaTable is valid only for ida.data.frame objects"))
   
  ####### Identifying Categorical Fields #########
  res <- idaTableDef(idadf,F)
  categorical <- as.vector(res[res$valType=='CATEGORICAL','name'])
  args<-categorical
  
  if (!length(args)) 
    stop("no categorical columns to tabulate")       
  
  dims <- integer()
  nm <- NULL
  
  for (a in args) {
    
    levels<-idaLevels(idadf,a) 	 #Select Levels of the Column
    o<-order(levels[1])		 #Order them in Ascending
    nm<-c(nm,list(levels[o,1]))	 #List of Level Names
    nl <- length(unlist(levels))	 #Number of Levels			 
    dims <- c(dims, nl) 
  }
  if (prod(dims) > max.entries) {
    print("Value counts for factors:")
    names(dims)<-args;
    print(dims)
    stop("Attempt to make a table with more than max.entries elements. Either increase max.entries or remove columns with too many levels.")
  }
  
  nmz<-c(nm,list(0))		# An extra column for counting occurences
  com<-expand.grid(nmz)		# Get all possible combinations
  bin<-idaFactors(idadf,args)		# Get actually Existing Combinations	
  names(com)<-names(bin)
  total<-rbind(bin,com)
  mat<-total[!duplicated(total[-length(total)]),]   # Remove duplicated combinations generated previously
  mat<-mat[ do.call(order, mat), ]
  names(nm)<-args
  
  z<-array(mat$COUNT,dim=rev(dims))               # Preparing the result for display
  z<-aperm(z)
  dim(z)<-dims
  dimnames(z)<-nm
  class(z)<-"table"
  return(z)
}


########################## idaLevels #####################################################

idaLevels<-function(idadf,column) {
  query <- paste("SELECT DISTINCT" ,paste("\"", column, "\"", collapse=",",sep=''), " FROM ",idadf.from(idadf)," ",ifelse(nchar(idadf@where),paste(" WHERE ",idadf@where,sep=''),''))
  return(idaQuery(query));
}

########################## idaFactors #####################################################

idaFactors<-function(idadf,args) {
  
  query <- paste("SELECT" ,paste("\"", args, "\"", collapse=",",sep=''),",COUNT(*) as COUNT FROM ",idadf.from(idadf)," ",ifelse(nchar(idadf@where),paste(" WHERE ",idadf@where,sep=''),'')," GROUP BY",paste("\"", args, "\"", collapse=",",sep=''))
  return(idaQuery(query))
}

########################## cor #####################################################
setMethod(f="cor", signature=c(x="ida.data.frame"),
    
    function(x,y=NULL,use = NULL,method = NULL ) {   
      
      # use
      if (!missing(use) && !is.null(use))
        stop(simpleError("use option is not implemented yet"))
      
      # method
      if (!missing(method) && !is.null(method))
        stop(simpleError("method option is not implemented yet"))
      
      if(nrow(x)==0) {
        stop("No rows in ida.data.frame.")
      }
      
      ####### Identifying Numeric Fields #########
      res <- idaTableDef(x,F)
      xCols <- as.vector(res[res$valType=='NUMERIC','name'])
      
      # Check if any columns are left in the end
      if (!length(xCols)) 
        stop("nothing to calculate") 
      
      xMean<-idaMean(x,xCols)
      n<-NROW(x)
      queryList <- c();
      
      
      if (!missing(y) && !is.null(y)){
        
        if (!is.ida.data.frame(y))
          stop("cor is valid only for ida.data.frame objects")
        
        if(idadf.from(x)!= idadf.from(y))
          stop("x and y must be from the same database table")
        
        ####### Identifying Numeric Fields #########
        res <- idaTableDef(y,F)
        yCols<- as.vector(res[res$valType=='NUMERIC','name'])
        
        
        # Check if any columns are left in the end
        if (!length(yCols)) 
          stop("nothing to calculate") 
        
        yMean<-idaMean(y,yCols)
        
        for(i in 1:length(xCols)) {
          for(j in 1:length(yCols)) {
            queryList <- c(queryList,paste("SUM((", paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")*(", paste("\"", yCols[j], "\"", collapse=",",sep=''),"-",yMean[j],"))/(SQRT(SUM((",paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")*(",paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")))*SQRT(SUM((",paste("\"", yCols[j], "\"", collapse=",",sep=''),"-",yMean[j],")*(",paste("\"", yCols[j], "\"", collapse=",",sep=''),"-",yMean[j],"))))",sep=''));
          }	
        }
      }
      
      
      else{
        for(i in 1:length(xCols)) {
          for(j in i:length(xCols)) {
            queryList <- c(queryList,paste("SUM((", paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")*(", paste("\"", xCols[j], "\"", collapse=",",sep=''),"-",xMean[j],"))/(SQRT(SUM((",paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")*(",paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")))*SQRT(SUM((",paste("\"", xCols[j], "\"", collapse=",",sep=''),"-",xMean[j],")*(",paste("\"", xCols[j], "\"", collapse=",",sep=''),"-",xMean[j],"))))",sep=''));
          }	
        }
        yCols<-xCols
      }
      
      queryList<-paste("SELECT ", paste(queryList,sep=',',collapse=',')," FROM ",idadf.from(x)," ",ifelse(nchar(x@where),paste(" WHERE ",x@where,sep=''),''),sep='');
      cor<-idaQuery(queryList)
      mdat <- matrix(1:(length(xCols))*(length(yCols)),nrow=length(xCols),ncol=length(yCols),dimnames = list(c(xCols),c(yCols)),byrow=T)
      
      # Arrange matrix values
      r <- 1;
      c <- 1;
      for(i in 1:ncol(cor)) {
        mdat[r,c] <- cor[[i]][1];
        if (is.null(y))
          mdat[c,r] <- mdat[r,c];
        c <- c + 1;
        if(c>length(yCols)) {
          r <- r+1;
          c <- if (is.null(y)) r else 1
        }
      }
      mdat
    }
)

########################## sd #############################################################
setMethod(f="sd", signature=c(x="ida.data.frame"),
    function(x,na.rm=NULL) { 
      idadf<-x 
      
      # na.rm
      if (!missing(na.rm) && !is.null(na.rm))
        warning("Missing values are not removed before computation")
      
      if(nrow(x)==0) {
        stop("No rows in ida.data.frame.")
      }
      
      col<-idadf@cols
      if(length(col)>1)
        warning("the condition has length > 1 and only the first numeric element will be used")
      
      queryList <- c();
      
      ####### Identifying Numeric Fields #########
      res <- idaTableDef(idadf,F)
      numeric <- as.vector(res[res$valType=='NUMERIC','name'])
      
      col<-numeric[1];
      
      # Check if any column is left
      if (!length(col)) 
        stop("nothing to calculate") 
      
      n<-NROW(idadf)
      mean<-idaMean(idadf,col)
      queryList <-paste("SUM((", paste("\"", col, "\"", collapse=",",sep=''),"-",mean,")*(",paste("\"", col, "\"", collapse=",",sep=''),"-",mean,"))/",n-1,sep='');
      queryList<-paste("SELECT ", queryList," FROM ",idadf.from(idadf)," ",ifelse(nchar(idadf@where),paste(" WHERE ",idadf@where,sep=''),''),sep='');
      sd<-sqrt(idaQuery(queryList))
      return(sd[1,1])
    }
)

########################## cov #############################################################
idacov<-function(x,y=NULL,use = NULL,method = NULL ) {   
  if (!is.ida.data.frame(x))
    stop("cov is valid only for ida.data.frame objects")
  
  if(nrow(x)==0) {
    stop("No rows in ida.data.frame.")
  }
  # use
  if (!missing(use) && !is.null(use))
    stop(simpleError("use option is not implemented yet"))
  
  # method
  if (!missing(method) && !is.null(method))
    stop(simpleError("method option is not yet available"))
  
  ####### Identifying Numeric Fields #########
  res <- idaTableDef(x,F)
  xCols <- as.vector(res[res$valType=='NUMERIC','name'])
  
  # Check if any columns are left in the end
  if (!length(xCols)) 
    stop("nothing to calculate") 
  
  xMean<-idaMean(x,xCols)
  n<-NROW(x)
  queryList <- c();
  
  if (!missing(y) && !is.null(y)){
    
    if (!is.ida.data.frame(y))
      stop("cov is valid only for ida.data.frame objects")
    
    if(idadf.from(x)!= idadf.from(y))
      stop("x and y must be from the same database table")
    
    ####### Identifying Numeric Fields #########
    res <- idaTableDef(y,F)
    yCols<- as.vector(res[res$valType=='NUMERIC','name'])
    
    # Check if any columns are left in the end
    if (!length(yCols)) 
      stop("nothing to calculate") 
    
    yMean<-idaMean(y,yCols)
    
    for(i in 1:length(xCols)) {
      for(j in 1:length(yCols)) {
        queryList <- c(queryList,paste("SUM((", paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")*(", paste("\"", yCols[j], "\"", collapse=",",sep=''),"-",yMean[j],"))/",n-1,sep=''));
      }	
    }
  }
  
  
  else{
    for(i in 1:length(xCols)) {
      for(j in i:length(xCols)) {
        queryList <- c(queryList,paste("SUM((", paste("\"", xCols[i], "\"", collapse=",",sep=''),"-",xMean[i],")*(", paste("\"", xCols[j], "\"", collapse=",",sep=''),"-",xMean[j],"))/",n-1,sep=''));
      }	
    }
    yCols<-xCols
  }
  
  queryList<-paste("SELECT ", paste(queryList,sep=',',collapse=',')," FROM ",idadf.from(x)," ",ifelse(nchar(x@where),paste(" WHERE ",x@where,sep=''),''),sep='');
  
  cov<-idaQuery(queryList)
  mdat <- matrix(1:(length(xCols))*(length(yCols)),nrow=length(xCols),ncol=length(yCols),dimnames = list(c(xCols),c(yCols)),byrow=T)
  
  # Arrange matrix values
  r <- 1;
  c <- 1;
  for(i in 1:ncol(cov)) {
    mdat[r,c] <- cov[[i]][1];
    if (is.null(y))
      mdat[c,r] <- mdat[r,c];
    c <- c + 1;
    if(c>length(yCols)) {
      r <- r+1;
      c <- if (is.null(y)) r else 1
    }
  }
  mdat
}

########################## idaVar #############################################################
setMethod(f="var", signature=c(x="ida.data.frame"),
    function(x,y=NULL,na.rm=NULL,use=NULL) {
      idacov(x,y) 
    }
)

###############################################################################################
## workaround for generic function
cov.ida.data.frame <- function (x, y = NULL, use = NULL, method = NULL) {
  return(idacov(x,y))
}

setMethod("cov", signature(x="ida.data.frame"), cov.ida.data.frame)
setMethod("cov", signature(y="ida.data.frame"), cov.ida.data.frame)
setMethod("cov", signature(x="ida.data.frame", y="ida.data.frame"), cov.ida.data.frame)

############################### summary ###############################################################

setMethod(f="summary", signature=c("ida.data.frame"),
    function (object,digits=max(3L, getOption("digits") -3L), maxsum = 7L, ...) {
      idadf<-object
          
      options(scipen=999)
      
      ####### Identifying Categorical Fields #########
      res <- idaTableDef(idadf,F)
      numeric <- as.vector(res[res$valType=='NUMERIC','name'])
      categorical <- as.vector(res[res$valType=='CATEGORICAL','name'])
      
      ####### Calculating Summary Values #############
      catRes<-idaCategorical(idadf,categorical,maxsum) 
      maxsum<-max(7,maxsum)
      quantiles <- idaQuantiles(idadf,numeric,maxsum,digits)
      
      final<-c(quantiles,catRes)
      final<-final[intersect(idadf@cols,c(numeric,categorical))]
      final<-unlist(final)
      dim(final) <- c(maxsum, length(numeric)+length(categorical))
      
      ####### Naming the summary columns ############
      blanks <- paste(character(max(10, na.rm = TRUE) + 2L),collapse = " ")
      pad <- floor(nchar(final[1,])/2 - nchar(intersect(idadf@cols,c(numeric,categorical)))/2)
      names <- paste0(substring(blanks, 1, pad), intersect(idadf@cols,c(numeric,categorical)))
      dimnames(final)<-list(rep.int("",maxsum),names)
      attr(final, "class") <- c("table")
      final		
    }
)


########################## Categorical ##################################################


idaCategorical<-function(idadf,categorical,maxsum=7L) {
  
  res<-c();
  output<-list();
  
  if(length(categorical)){
    
    for(a in categorical){
      
      query <- paste("SELECT",paste("\"", a, "\"", collapse=",",sep=''),",COUNT(*) FROM",idadf.from(idadf)," ",ifelse(nchar(idadf@where),paste(" WHERE ",idadf@where,sep=''),'')," GROUP BY",paste("\"", a, "\"", collapse=",",sep=''))
      
      cf<-idaQuery(query)
      if (length(cf[,1]) > maxsum) {
        o<-order(cf[,2],decreasing=TRUE)
        cf<-cf[o,]
        drop <- maxsum:length(cf[,1])
        cf[maxsum,1]<-"(Others)"
        cf[maxsum,2]<- sum(cf[drop,2])
        cf<-cf[1:maxsum,]
      }
      
      
      #### Dont we need to Display the NAs ####
      
      nas<-cf[is.na(cf[,1]==c('NA')),]
      cf<-cf[!is.na(cf[,1]==c('NA')),]
      
      if(length(nas)){
        nas[1,1]<-"NA's"
        cf[length(cf[,1])+1,]<-nas}
      ########################################
      
      res<-paste0(format(unlist(cf[1])), ":", format(unlist(cf[2])), "  ")
      length(res)<-maxsum
      length(res)<-max(7,maxsum)
      output<-c(output,list(res))	}
  }
  names(output)<-c(categorical)
  output
}

########################## Quantiles Method 2 (R "default", type=7)########################
idaQuantiles<-function(idadf,numeric,maxsum=7L,digits=4) {
  
  quantiles<-list();
  probs = seq(0, 1, 0.25)
  
  label<-c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")  
  mean<- matrix(idaMean(idadf,numeric), ncol = length(numeric))
  
  for(i in 1:length(numeric)){
    
    col <- paste("\"", numeric[i], "\"", collapse=",",sep='')
    query<-paste("SELECT COUNT(*) FROM", idadf.from(idadf),"WHERE ",col,"IS NULL "," ",ifelse(nchar(idadf@where),paste(" AND ",idadf@where,sep=''),''))
    nas<-as.numeric(idaQuery(query))
    
    n <- NROW(idadf)
    
    if(nas<n) {
      
      n<-n- nas
      index <- 1 + (n - 1) * probs
      lo <- floor(index)
      hi <- ceiling(index)
      unique<-sort(unique(c(lo,hi)))
      j <- which(index > lo)
      h <- (index - lo)[j]
      where <- paste("", unique, collapse=",",sep='')
      
      query <- paste("SELECT ",col,"FROM ","(SELECT ROW_NUMBER() over (ORDER BY",col,") as rn,",col,"FROM (", idadf.query(idadf),")) WHERE rn in(",where,")")
      
      full<-data.matrix(idaQuery(query))
      qs<-full[unique%in%lo,]
      xhi<-full[unique%in%hi,]
      
      if(n<5){
        qs<-qs[lo]
        xhi<-xhi[hi]}
      
      qs[j] <- (1 - h) * qs[j] + h * xhi[j]
      qs<-c(qs[1:3],mean[i],qs[4:5])	
      qs <- paste0(format(label), ":", format(qs,digits=digits,nsmall=3), "  ")
      if(nas>0){
        qs <- c(qs,paste0("NA's   :", format(nas), "  "))
      }
    } else {
      qs <- paste0("NA's   :", format(nas), "  ")
    }
    
    length(qs)<-maxsum;
    quantiles<-c(quantiles,list(c(qs)))
    
  }
  names(quantiles)<-numeric
  quantiles
}
########################## idaMean #######################################################

idaMean<-function(idadf,numeric) {
  query <- paste(paste("AVG(\"", numeric, "\") AS ", numeric, collapse=",",sep='')," FROM ", idadf.from(idadf), " ",ifelse(nchar(idadf@where),paste(" WHERE ",idadf@where,sep=''),''))
  m<-idaQuery ("SELECT ", query)
  m<-matrix(unlist(m),nrow=1)
  return(m)
}