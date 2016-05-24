summary.mcposthoc<-function(object,ph.list=NULL,term=NULL,print=TRUE,...){
                if(is.null(ph.list)){
			tmp<-strsplit(names(object$summaries),"_")
			tmp<-sort(unique(unlist(lapply(tmp,FUN=function(x)x[1]))))
			tmp.mat<-as.matrix(tmp)
			if(length(tmp)==1){
				num<-1
			}else{
				cat("display results for which posthoc list?\n")
				print(tmp.mat)
				num<-as.numeric(readline(prompt="Enter a number: "))
			}
			ph.list<-tmp[num]
		}

		tmp<-names(object$summaries)[grep(ph.list,names(object$summaries))]
		for(i in 1:length(tmp)){
			if(i==1){
				terms<-rownames(object$summaries[[tmp[i]]])
			}else{
				terms<-c(terms,rownames(object$summaries[[tmp[i]]]))
			}
			rem<-which(terms=="(Intercept)")
			if(length(rem)>0){
				terms<-terms[-rem]
			}
			terms<-sort(unique(terms))
			terms.mat<-as.matrix(terms)
		}
                if(is.null(term)){
			cat("display results for which reference level?\n")
			print(terms.mat)
			num<-as.numeric(readline(prompt="Enter a number: "))
                        cat("\n")
			term<-terms.mat[num]
                }

		for(i in 1:length(tmp)){
			if(i==1){
				nms<-vector("character")
				#ref.lev<-vector("character")
				res<-object$summaries[[tmp[i]]][term,]
				res<-res[-c(1:nrow(res)),]
			}
			if(term%in%rownames(object$summaries[[names(object$summaries)[i]]])){
				nms<-c(nms,tmp[i])
				res<-rbind(res,object$summaries[[tmp[i]]][term,])
			}#else{
                                #ref.lev<-c(ref.lev,tmp[i])
                        #}
		}
		rownames(res)<-gsub(paste(ph.list,"_",sep=""),"",nms)
		#ref.lev<-gsub(paste(ph.list,"_",sep=""),"",ref.lev)
		res<-as.data.frame(res)
                res$TMP<-rownames(res)

		var<-object$var[[ph.list]]
		levs<-strsplit(rownames(res),"_")
		for(i in var){
			tmp<-vector("character")
			for(j in 1:length(levs)){
				tmp<-c(tmp,levs[[j]][grep(i,var)])
			}
			res[,i]<-tmp
		}
	        res<-res[,c((ncol(res)-length(var)+1):ncol(res),1:(ncol(res)-length(var)))]
		res<-na.omit(res)
                res<-res[order(res[,"TMP"]),]
                res<-res[,-grep("TMP",colnames(res))]
	        rownames(res)<-1:nrow(res)

#                if(length(var)>1){
#                        tmp<-res[,var]
#                        for(v1 in 1:nrow(tmp)){
#                                if(v1==1){
#                                        for(v2 in 1:ncol(tmp)){
#                                                if(v2==1){
#                                                        tmp2<-paste(var[v2],tmp[v1,v2],sep="")
#                                                }else{
#                                                        tmp2<-paste(tmp2,paste(var[v2],tmp[v1,v2],sep=""),sep="_")
#                                                }
#                                        }
#                                        tmp3<-tmp2
#                                }else{
#                                        for(v2 in 1:ncol(tmp)){
#                                                if(v2==1){
#                                                        tmp2<-paste(var[v2],tmp[v1,v2],sep="")
#                                                }else{
#                                                        tmp2<-paste(tmp2,paste(var[v2],tmp[v1,v2],sep=""),sep="_")
#                                                }
#                                        }
#                                        tmp3<-c(tmp3,tmp2)
#                                }
#                        }
#                        res$Contrast<-tmp3
#                        ref.lev<-strsplit(ref.lev,"_")
#                        for(i in nrow(res)){
#                                
#                        }
#                        for(v1 in 1:nrow(tmp)){
#                                if(v1==1){
#                                        for(v2 in 1:ncol(tmp)){
#                                                if(v2==1){
#                                                        tmp2<-paste(var[v2],ref.lev[[v1]][v2],sep="")
#                                                }else{
#                                                        tmp2<-paste(tmp2,paste(var[v2],ref.lev[[v1]][v2],sep=""),sep="_")
#                                                }
#                                        }
#                                        tmp3<-tmp2
#                                }else{
#                                        for(v2 in 1:ncol(tmp)){
#                                                if(v2==1){
#                                                        tmp2<-paste(var[v2],ref.lev[[v1]][v2],sep="")
#                                                }else{
#                                                        tmp2<-paste(tmp2,paste(var[v2],ref.lev[[v1]][v2],sep=""),sep="_")
#                                                }
#                                        }
#                                        tmp3<-c(tmp3,tmp2)
#                                }
#                        }
#                }

		if(print){
	          cat("--------------------\n")
	          cat("POSTHOC LIST:",ph.list,"\n")
	          cat("\n")
	          cat("TERM:",term,"\n")
	          #cat("(first ref.lev of the difference)\n")
	          cat("\n")
	          #cat("VARIABLES:",object$var[[ph.list]],"\n")
	          #cat("\n")
	          cat("SUMMARY:\n")
	          #cat("CONTRAST LEVEL:\n")
	          print(res)
	          #cat("\n")
	          #cat("NOTE: \"Estimate\" is equal to contrast - reference level.\n")
		}
		return(invisible(list(ph.list=ph.list,term=term,summary=res)))
}
