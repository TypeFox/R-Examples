tt.Summary <-
function(tabtabres,trees,focal){
				no.trees<-length(c(trees[focal]))
				treeComp<-c(trees[focal])
				num.edges<-NA
					for(Tr in 1:no.trees){#
					thisTre<-treeComp[[Tr]]
					num.edges[Tr]<-sum(thisTre$edge[,2]>length(thisTre$tip.label))
					}
				n<-sum(num.edges)

	
	colNames<-c("k","Groups")
	totgroups<-(dim(tabtabres)[2]-2)/5


					ratetocalc<-totgroups
					rowtocal<-dim(tabtabres)[1]
					P1s<-seq(from=3,to=dim(tabtabres)[2],by=5)
					P2s<-seq(from=4,to=dim(tabtabres)[2],by=5)
					mIDs<-seq(from=7,to=dim(tabtabres)[2],by=5)###Changed this from from=6 to from=7

					rAtes<-array(dim=c(rowtocal,ratetocalc))
					ratcn<-NA
						for (i in 1:ratetocalc){
							ratcn[i]<-paste("g",i,"_rate",sep="")
							}
					
					for (i in 1:rowtocal){
						for(j in 1:ratetocalc){
							p.1<-tabtabres[i,P1s[j]]
							p.2<-tabtabres[i,P2s[j]]
							m.id<-tabtabres[i,mIDs[j]]
							if(is.na(m.id)==TRUE){rAtes[i,j]<-NA}
							else{
							if(m.id==1){rAtes[i,j]<-p.1}
							if(m.id==2)rAtes[i,j]<-calc.rate.weib(c(p.1,p.2))
							if(m.id==3)rAtes[i,j]<-calc.rate.lnorm(c(p.1,p.2))
							if(m.id==4)rAtes[i,j]<-calc.rate.vrat(c(p.1,p.2))}
							}	
						}
					colnames(rAtes)<-ratcn



	for (i in 1:totgroups){
			P1n<-paste("g",i,"_P1",sep="")
			P2n<-paste("g",i,"_P2",sep="")
			Lin<-paste("g",i,"_LL",sep="")
			Min<-paste("g",i,"_mod.id",sep="")
			npn<-paste("g",i,"_n.param",sep="")
			colNames<-c(colNames,P1n,P2n,Lin,Min,npn)
		}
	
	col.liks<-seq(5,length(colNames),5)
	col.nparm<-seq(7,length(colNames),5)
		lik.sums<-apply(tabtabres[,col.liks],1,sum,na.rm=TRUE)
		nparm.sums<-apply(tabtabres[,col.nparm],1,sum,na.rm=TRUE)
		aic<-(-2*lik.sums)+(2*nparm.sums)
		aicc<-aic+((2*nparm.sums*(nparm.sums+1))/(n-nparm.sums-1)) 
		dAICc<-aicc-min(aicc)
	result<-data.frame(tabtabres,nparm.sums,lik.sums,aic,aicc,dAICc) 
	colnames(result)<-c(colNames,"total.param","LL","AIC","AICc","dAICc")
	result<-cbind(result[,1:(length(result[1,])-5)],rAtes,result[,(length(result[1,])-4):length(result[1,])])
	return(result)
	}

