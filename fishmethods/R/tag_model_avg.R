tag_model_avg<-function(...,global=NULL){
             foo.call <- as.character(match.call())[-1] 
             all.x <- list(...) 
             n<-length(all.x)
             names(all.x) <- unique(foo.call)
            if(is.null(global)) chat<-all.x[[1]]$statistics[8] else
		chat<-eval(parse(text = paste("all.x$",global,"$statistics[8]",sep="")))           
            if(chat<1) chat<-1
#Check model type
       check<-NULL
       chyrs<-NULL
	 for(i in 1:n){ 
             check[i]<-all.x[[i]]$type
             chyrs[i]<-length(all.x[[i]]$fishing_mortality$Year)
          }
        if(length(unique(check))>1) stop ("Model outputs are not of the same type (irm_cr or irm_h)")
	  if(length(unique(chyrs))>1) stop("Model year intervals are not the same length")
        
   #Create Statistics Table
        outpt<-array(NA,dim=c(n,11))
        rownames(outpt)<-c(unique(foo.call))
        #rownames(outpt)<-c(dodo[-length(dodo)])
        colnames(outpt)<-c("Likelihood","No. Parms","AIC","AICc","N","QAIC"
                 ,"QAICc","dQAICc","e(-0.5*dQAICc)","QAICc Wgts","Global chat")
         cnt<-0
         for(x in list(...)){  
               cnt<-cnt+1
          outpt[cnt,1]<-x$statistics[1]
          outpt[cnt,2]<-x$statistics[2]
          outpt[cnt,3]<-x$statistics[3]
          outpt[cnt,4]<-x$statistics[4]
          outpt[cnt,5]<-x$statistics[5]
          outpt[cnt,6]<-((-2*outpt[cnt,1])/chat)+2*(outpt[cnt,2]+1)
          outpt[cnt,7]<-outpt[cnt,6]+(2*(outpt[cnt,2]+1)*(outpt[cnt,2]+2))/(outpt[cnt,5]-outpt[cnt,2])
         }
         outpt[,8]<-outpt[,7]-min(outpt[,7])
         outpt[,9]<-exp(-0.5*outpt[,8])
         outpt[,10]<-outpt[,9]/sum(outpt[,9])
         outpt[,11]<-chat

  #################Generate Adjusted SE          
       for(i in 1:n){
         all.x[[i]]$fishing_mortality$SE<-sqrt(all.x[[i]]$fishing_mortality$VAR*chat)
         all.x[[i]]$natural_mortality$SE<-sqrt(all.x[[i]]$natural_mortality$VAR*chat)
         all.x[[i]]$total_mortality$SE<-sqrt(all.x[[i]]$total_mortality$VAR*chat)
         all.x[[i]]$survival$SE<-sqrt(all.x[[i]]$survival$VAR*chat)
         if(unique(check)=="cr") all.x[[i]]$tag_mortality$SE<-sqrt(all.x[[i]]$tag_mortality$VAR*chat)
        }
 #Wgt SE
	  yrlen<-unique(chyrs)
	  F<-data.frame(Year=all.x[[1]]$fishing_mortality$Year,avgF=NA,Wgt_SE=NA,Uncond_SE=NA )
        M<-data.frame(Year=all.x[[1]]$natural_mortality$Year,avgM=NA,Wgt_SE=NA,Uncond_SE=NA )
	  Z<-data.frame(Year=all.x[[1]]$total_mortality$Year,avgZ=NA,Wgt_SE=NA,Uncond_SE=NA )
        S<-data.frame(Year=all.x[[1]]$survival$Year,avgS=NA,Wgt_SE=NA,Uncond_SE=NA )
	  if(unique(check)=="cr") FA<-data.frame(Year=all.x[[1]]$tag_mortality$Year,avgFA=NA,Wgt_SE=NA,Uncond_SE=NA )
        tempF<-NULL;tempFWSE<-NULL;tempFUSE<-NULL
        tempM<-NULL;tempMWSE<-NULL;tempMUSE<-NULL
        tempZ<-NULL;tempZWSE<-NULL;tempZUSE<-NULL
        tempS<-NULL;tempSWSE<-NULL;tempSUSE<-NULL
	   if(unique(check)=="cr"){
          tempFA<-NULL;tempFAWSE<-NULL;tempFAUSE<-NULL
          }
        for(i in 1:n){
		tempF<-cbind(tempF,all.x[[i]]$fishing_mortality$F*outpt[i,10])
		tempFWSE<-cbind(tempFWSE,all.x[[i]]$fishing_mortality$SE*outpt[i,10])
            tempM<-cbind(tempM,all.x[[i]]$natural_mortality$M*outpt[i,10])
		tempMWSE<-cbind(tempMWSE,all.x[[i]]$natural_mortality$SE*outpt[i,10])
            tempZ<-cbind(tempZ,all.x[[i]]$total_mortality$Z*outpt[i,10])
		tempZWSE<-cbind(tempZWSE,all.x[[i]]$total_mortality$SE*outpt[i,10])
            tempS<-cbind(tempS,all.x[[i]]$survival$S*outpt[i,10])
		tempSWSE<-cbind(tempSWSE,all.x[[i]]$survival$SE*outpt[i,10])
		if(unique(check)=="cr"){
            	tempFA<-cbind(tempFA,all.x[[i]]$tag_mortality$FA*outpt[i,10])
			tempFAWSE<-cbind(tempFAWSE,all.x[[i]]$tag_mortality$SE*outpt[i,10])
	        }

        }
       F$avgF<-rowSums(tempF)
       F$Wgt_SE<-rowSums(tempFWSE)
       M$avgM<-rowSums(tempM)
       M$Wgt_SE<-rowSums(tempMWSE)
       Z$avgZ<-rowSums(tempZ)
       Z$Wgt_SE<-rowSums(tempZWSE)
       S$avgS<-rowSums(tempS)
       S$Wgt_SE<-rowSums(tempSWSE)
		if(unique(check)=="cr"){
      	 FA$avgFA<-rowSums(tempFA)
      	 FA$Wgt_SE<-rowSums(tempFAWSE)
		}
#calculate unconditional SE
  	 for(i in 1:n){
	    tempFUSE<-cbind(tempFUSE,outpt[i,10]*(all.x[[i]]$fishing_mortality$SE^2+
                 (all.x[[i]]$fishing_mortality$F-F$avgF)^2))
   	    tempMUSE<-cbind(tempMUSE,outpt[i,10]*(all.x[[i]]$natural_mortality$SE^2+
                 (all.x[[i]]$natural_mortality$M-M$avgM)^2))
   	    tempZUSE<-cbind(tempZUSE,outpt[i,10]*(all.x[[i]]$total_mortality$SE^2+
                 (all.x[[i]]$total_mortality$Z-Z$avgZ)^2))
	    tempSUSE<-cbind(tempSUSE,outpt[i,10]*(all.x[[i]]$survival$SE^2+
                 (all.x[[i]]$survival$S-S$avgS)^2))
          if(unique(check)=="cr"){
			tempFAUSE<-cbind(tempFAUSE,outpt[i,10]*(all.x[[i]]$tag_mortality$SE^2+
                 (all.x[[i]]$tag_mortality$FA-FA$avgFA)^2))
        	   }
       }
 	F$Uncond_SE<-sqrt(rowSums(tempFUSE))
  	M$Uncond_SE<-sqrt(rowSums(tempMUSE))
  	Z$Uncond_SE<-sqrt(rowSums(tempZUSE))
  	S$Uncond_SE<-sqrt(rowSums(tempSUSE))
     if(unique(check)=="cr"){
          FA$Uncond_SE<-sqrt(rowSums(tempFAUSE))
       }
      ans<-NULL
      ans$statistics<-outpt
	ans$model_averaged_F<-F
	ans$model_averaged_M<-M
	ans$model_averaged_Z<-Z
	ans$model_averaged_S<-S
     if(unique(check)=="cr") ans$model_averaged_FA<-FA
          return(ans)  
} 

