MDR_forward <-
function(Index,dat,alpha=NULL,rsquared=NULL,trace=NULL,ModifyName=NULL,max.step=NULL,...){
	if(is.null(alpha))alpha=0.01 
	if(is.null(rsquared))rsquared=0.02
	if(is.null(trace))trace=FALSE
	if(is.null(ModifyName))ModifyName=TRUE
	if(is.null(max.step))max.step = nrow(Index) 
	res <- MDR_forward1(Index=Index,dat=dat,alpha=alpha,rsquared=rsquared,trace=trace,ModifyName=ModifyName,max.step=max.step)
	return(res)

}

MDR_forward1 <-
function(Index,dat,alpha=NULL,rsquared=NULL,trace=NULL,ModifyName=NULL,max.step=NULL,...){
	
  
	if(ModifyName){
		data=Data_Name(dat) # This fuction spends lots of time
  		dat=data$dat #head(dat)
		dat_Original = dat
  		datColNames=data$datColNames
		data=NULL
	}

	if(!ModifyName){
		dat_Original = dat
  		datColNames=colnames(dat)
	}

	c=ncol(Index) #order number of combination (c=1 means no combination) 
	dat2 = data.frame(dat$y) #dat[,ncol(dat)]
	colnames(dat2)=c("y") #head(dat2)

  	marker =  numeric()
  	R=0
  	Rmax=0 ## save max R square in every step, temporary value
  	row=1
  	R2 = numeric() ## output of R square for every step
	R2_o = numeric() ## save the original R square
	
	DeleteRow = numeric() 
	R2_OneWay = numeric() ## output of R square of single variables for every step
	R2_OneWay_o = numeric()
	R2_singal = numeric() 
	R2_singal_o = numeric() # save the R square
		
  	Ind=NULL
  	k=1
  	IndexNew = matrix(0,,)

  	if (c<1) stop ("'c' must be c >= 1") ##########[1]c<1
	
  	#################[2]c=1#################
  	if (c==1) { ##(1)##

		IndexNew=numeric()
		StepNum=0

      		while(row !=0){ ##(2)##
			StepNum=StepNum+1
	   		row=0
         		r=nrow(Index) 
			colname<- colnames(dat)
	
   	   		for(i in 1:r){ ##(3)##
				if(trace)cat("Step",StepNum," ")
				if(trace)cat("Scan:", i,"\n")
				Name <- c(colnames(dat2),colname[Index[i]])
				ok = complete.cases(factor(dat[,Index[i]]))
				Data1 <- data.frame(dat2,factor(dat[,Index[i]]))
          			colnames(Data1) <- Name
				Data1 = Data1[ok,]
				reg <- lm(y~.,data=Data1) 
				a=summary(reg)
				r2=a$r.squared # R square
				r2_o = r2 # save the R square 1/16/2013
		
				nR = nrow(Data1)
 				pR = ncol(Data1)-1
				r2 = 1-(1-r2)*(nR-1)/(nR-pR-1)

		
	
				if(r2<1 & r2 > Rmax & r2 > (R + rsquared)){
					ANOVA <- aov(y~.,data=Data1)
		 			A=summary(ANOVA)
		  			P = A[[1]][[5]]
		  			s=length(P)
		  			P = P[-s]
		  			Pvalue= max(P)
		  			if(Pvalue < alpha){
						Data = Data1
						Rmax = r2
						Rmax_o = r2_o		
						row = i
						OK = ok		
		  			}
				}
				Data1 = NULL

	   		}##end of(3)##

	   		if(Rmax < (R+rsquared)){row = 0}
	   		if(max.step < StepNum) {row = 0}

	   		if(row!=0){ ##(4)##
				data_singial = data.frame(dat_Original$y,dat_Original[,Index[row]])
				colnames(data_singial) = c("y","x")
				ok1 = complete.cases(data_singial)
				data_singial = data_singial[ok1,]
				reg <- lm(y~.,data=data_singial) 
				a=summary(reg)
				r2=a$r.squared				
				nR = nrow(data_singial)		
				R2_singal[k] = 1-(1-r2)*(nR-1)/(nR-2)
				R2_singal_o[k] = r2 #1/16/2013

				dat2 = Data
				Data = NULL
				R = Rmax
				R2_o[k] = Rmax_o
				R2[k] = R
				IndexNew[k]=Index[row]
				Index=matrix(Index[-row])				
				dat=dat[OK,]
				DeleteRow[k]=length(which(OK=="FALSE"))	
				k=k+1 	
	   		}##end of(4)##

        	}##end of(2)##

	  	IndexNew=matrix(IndexNew)
	  	R2_OneWay = R2 ##2/23
		R2_OneWay_o = R2_o ##1/16/2013

	} ##end of(1)##


  	##################[3]c>1#################
  	if (c > 1) { ##(1)##
		StepNum=0
    		while(row !=0){ ##(2)##
			StepNum=StepNum+1
			row=0
			r=nrow(Index)

			for(i in 1:r){##(3)##
				if(trace)cat("Step",StepNum," ")
				if(trace)cat("Scan:", i,"\n")
				colname<- colnames(dat)
            			ok=complete.cases(dat[,Index[i,]])
				coldata = dat[,Index[i,1]]
				for(j in 1:(c-1)){		
					coldata = paste(coldata,dat[,Index[i,j+1]])			
				colname[Index[i,j+1]]<- paste(colname[Index[i,j]],colname[Index[i,j+1]],sep="")
				}

				Name <- c(colnames(dat2),colname[Index[i,j+1]])		
				Data1 <- data.frame(dat2,factor(coldata))
				colnames(Data1) <- Name
	      			Data1 = Data1[ok,]

				reg <- lm(y~.,data=Data1)
				a=summary(reg)
				r2=a$r.squared		
				r2_o = r2
				nR = nrow(Data1)
 				pR = ncol(Data1)-1
				r2 = 1-(1-r2)*(nR-1)/(nR-pR-1)

				if(r2<1 & r2 > Rmax & r2 > (R + rsquared)){
		  			ANOVA <- aov(y~.,data=Data1)
		  			A=summary(ANOVA)
		  			P = A[[1]][[5]]
		  			s=length(P)
		  			P = P[-s]
		  			Pvalue= max(P)		
		  			if(Pvalue < alpha){
						Data = Data1
						Rmax = r2  
						Rmax_o = r2_o 
						row = i
						OK = ok
		  			}
				}
				Data1=NULL		
			}##end of(3)## 

			if(Rmax < (R+rsquared)){row = 0}
			if(max.step < StepNum) {row = 0} 

			if(row!=0){##(4)##
 				ok1=complete.cases(dat_Original[,Index[row,]])	    
				coldata = dat_Original[,Index[row,1]]
				for(j in 1:(c-1)){		
					coldata = paste(coldata,dat_Original[,Index[row,j+1]])			
				}						
				data_singial = data.frame(dat_Original$y,factor(coldata)) 
				colnames(data_singial) = c("y","x")
				data_singial = data_singial[ok1,]
				reg <- lm(y~.,data=data_singial) 
				a=summary(reg)
				r2=a$r.squared				
				nR = nrow(data_singial)		
				R2_singal[k] = 1-(1-r2)*(nR-1)/(nR-2)
				R2_singal_o[k] = r2 #1/16/2013

				dat2 = Data
				Data=NULL
				R = Rmax

				R2[k] = R
				R2_o[k]=Rmax_o
				IndexNew=data.frame(IndexNew,Index[row,])
				Ind1=Index[row,]
				Index= Index[-row,]		
	
	      			dat = dat[OK,]
            			DeleteRow[k]=length(which(OK=="FALSE"))
			}##end of(4)##
 
			if(row!=0){ ##(5)##	
				col=colnames(dat)
				Ind=union(Ind,Ind1)
				Name0 <- c("y",col[Ind])    
				dat4 <- data.frame(dat$y,dat[,Ind])
	  			colnames(dat4)<-Name0    
	    			reg0 <- lm(y~.,data=dat4)
	    			b=summary(reg0)
				r2_1way = b$r.squared
				R2_OneWay_o[k]=r2_1way
				nR_1way = nrow(dat4)
 				pR_1way = ncol(dat4)-1
				r2_1way = 1-(1-r2_1way)*(nR_1way - 1)/(nR_1way - pR_1way -1) 
	    			R2_OneWay[k]=r2_1way
				k=k+1
				dat4 =NULL
	 		}##end of(5)##

		} ##end of(2)##
	} ##end of(1)##

   	if (c==1){IndexNew=t(t(IndexNew))}
   	if(c>1){IndexNew=t(IndexNew[,-1])}

   	############get marker name##############
   	Num = IndexNew[,1]
   	for(i in 1:c){Num=union(Num,IndexNew[,i])}
   	Num=sort(Num)
   	for(i in 1:length(Num)){marker[i]=datColNames[Num[i]]}
   	Num = t(t(Num))
   	marker = t(t(marker))
   	markers = cbind(Num,marker) 
   	summary = cbind(IndexNew,R2,R2_OneWay,R2_singal,R2_o,R2_OneWay_o,R2_singal_o,DeleteRow)
   	result=list(summary=summary,markers=markers)
   	return(result)

}
