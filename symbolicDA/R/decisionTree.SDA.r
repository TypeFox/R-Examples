#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Warstwowe drzewo klasyfikacyjne dla danych symbolicznych
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach projektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach li-cencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************
decisionTree.SDA<-function(sdt,formula,testSet,treshMin=0.0001,treshW=-1e10,tNodes=NULL,minSize=2,epsilon=1e-4,useEM=FALSE,multiNominalType="ordinal",rf=FALSE,rf.size=nrow(sdt$variables)%/%2+1,objectSelection=1:nrow(sdt$individuals))
{
  debug<-FALSE
	if (!(.is.symbolic(sdt)))
		stop("First argument has to be symbolic data table")
	ile<-nrow(sdt$individuals)
	if (is.null(tNodes))
		tNodes=ile%/%2

  f<-.parseFormula(formula,sdt)
  #print(f$variableSelection)
  classes<-f$classes

	nodes<-c(0,0,0,0,1)
	dim(nodes)<-c(1,5)
	colnames(nodes)<-c("question_type","variable","value","parent","terminal")

	nodeObjects<-array(0,c(2*tNodes,ile))
	nodeObjects[1,]<-rep(1,ile)
	
	conditionalProbab<-array(0,c(2*tNodes,max(classes)))
	conditionalProbab[1,]<-1/(max(classes))


	minimalSizeExists<-TRUE
	minimalWExists<-TRUE
	##print("conditionalProbab")
	##print(conditionalProbab[1,])
  vs<-f$variableSelection
  if(rf){
    vs<-sample(vs,rf.size)
    print("Random forest variable sample")
    print(vs)
  }
	petla<-0
	if(debug){write("STARTUJE",file="w.txt",append=FALSE)}
	cond_L_W<-cond_R_W<-array(0,c(max(classes)))
	while(sum(nodes[,"terminal"]==1)<tNodes && minimalSizeExists && minimalWExists && petla<3)
	{
		petla<-petla+1
		print(paste("petla",petla,Sys.time()))
		#if(petla==2)
		#	stop("TYMCZASOWE PRZERWANIE")
		W_maximal<--1e100
		minimalSizeExists<-FALSE
		#Psck<-array(0,c(2*tNodes,max(classes)))
		#for(i in 1:petla)
		#for(m in 1:max(classes))
		#{
		#	Psck[i,m]<-prod(conditionalProbab[i,classes==m])
		#}
		for(i in 1:nrow(nodes))
		{
			#print(i)
			lambda<-array(0,c(ile))
			for (k in 1:ile) 
			{
				if(nrow(nodes)==1)		# zastanowic sie
				{
					lambda[k]<-0.05
				}
				else
				{
					for (s in 1:nrow(nodes))
					{
						if (s!=i && nodes[s,"terminal"]==1)
						#if (nodes[s,"terminal"]==1)
						{
							lambda[k]<-lambda[k]+nodeObjects[s,k]*conditionalProbab[s,classes[k]]
						}
					}
				}
			}
			#print("lambda")
			#print(lambda)
			#stop()
			if(nodes[i,"terminal"]==1)
			{
				for (v in vs)
				{
          #print(paste("petla,zmienna",petla,v))
          #print(paste("Typ zmiennej:",sdt$variables[v,"type"],petla,Sys.time()))
					if (sdt$variables[v,"type"]=="IC")
					{
						centers<-unique(sort(c(sdt$indivIC[,v,1],sdt$indivIC[,v,2])))
						if(length(centers)>10){
               fact<-length(centers)%/%10
               sfact<-seq(fact,length(centers),fact)
               centers<-centers[sfact]
						}
						#print(paste("variable",v,"centers before"))
						#print(centers) 
						centers<-apply(cbind(centers[1:length(centers)-1],centers[2:length(centers)]),1,mean)
						#print(paste("variable",v))
						#print(centers) 
						p<-array(0,c(ile,length(centers)))
						dlugoscq<-length(centers)
						question_type<-1
					}
					if (sdt$variables[v,"type"]=="NM")
					{
						#print("NM")
						values<-as.vector((sdt$detailsListNom[sdt$detailsListNom[,"details_no"]==as.integer(sdt$variables[v,"details"]),"num"]))						
						values<-sort(values[1:(length(values)-1)])
						p<-array(0,c(ile,(2^(length(values)-1)-1)))
						dlugoscq<-(2^(length(values)-1)-1)
						question_type<-2
					}
					if ((sdt$variables[v,"type"]=="N" || sdt$variables[v,"type"]=="MN") && multiNominalType!="ordinal")
					{
					
						##print(sdt$variables[v,"details"])
						##print(sdt$detailsListNom[,"details_no"]==sdt$variables[v,"details"])
						values<-as.vector((sdt$detailsListNom[sdt$detailsListNom[,"details_no"]==as.integer(sdt$variables[v,"details"]),"num"]))						
						p<-array(0,c(ile,(2^(length(values)-1)-1)))
						dlugoscq<-(2^(length(values)-1)-1)
						question_type<-3
					}
					if ((sdt$variables[v,"type"]=="N" || sdt$variables[v,"type"]=="MN") && multiNominalType=="ordinal")
					{
					
						##print(sdt$variables[v,"details"])
						##print(sdt$detailsListNom[,"details_no"]==sdt$variables[v,"details"])
						values<-as.vector((sdt$detailsListNom[sdt$detailsListNom[,"details_no"]==as.integer(sdt$variables[v,"details"]),"num"]))						
						p<-array(0,c(ile,length(values)))
						dlugoscq<-length(values)
						question_type<-2
					}
					for (q in 1:dlugoscq)
					{
						leftNode<-NULL
						rightNode<-NULL
						if ((sdt$variables[v,"type"]=="N" || sdt$variables[v,"type"]=="MN") && multiNominalType!="ordinal" )
						{
							t<-q
							lewy<-NULL
							for(j in 1:length(values))
							{
								if((t %% 2)==1) 
								{
									lewy<-c(lewy,values[j])
								}
								t<-t%/%2
							}
						}
						if ((sdt$variables[v,"type"]=="N" || sdt$variables[v,"type"]=="MN") && multiNominalType=="ordinal" )
						{
						  lewy<-values[1:q]
						}
						##print(nodeObjects[nodeObjects[,"node"]==i,"object"])
						for(k in 1:ile)
						{
              #print(paste("k,q",k,q,Sys.time()))
							if (sdt$variables[v,"type"]=="IC")
							{
								przedzial<-sdt$indivIC[k,v,]
								if(przedzial[2]<centers[q])
								{
									p[k,q]<-1
								}
								else 
								{
									if(przedzial[1]>centers[q])
									{
										p[k,q]<-0
									}
									else
									{
										if (przedzial[2]==przedzial[1])
										{
											p[k,q]<-0
										}
										else
										{
											p[k,q]<-(centers[q]-przedzial[1])/(przedzial[2]-przedzial[1])
										}
									}
								}
							}
							if (sdt$variables[v,"type"]=="N" || sdt$variables[v,"type"]=="MN")
							{
								t<-sdt$indivN[sdt$indivN[,"indiv"]==k,]
								t<-t[t[,"variable"]==v,"value"]
								if(!is.null(t)&&length(t)!=0){
								p[k,q]=sum(sapply(as.integer(t),function(x){sum(lewy==x)}))/length(t)
								}
							}
							if (sdt$variables[v,"type"]=="NM")
							{
								t<-sdt$indivNM[sdt$indivNM[,"indiv"]==k,]
								##print(t)
								t1<-t[t[,"variable"]==v,"frequency"]
								t2<-t[t[,"variable"]==v,"value"]
								##print(k)
								##print(v)
								##print(t1)
								##print(t2)
								p[k,q]<-sum(t1[t2<=as.integer(values[q])])
								#TYMCZASOWO !!!!
								#p[k,q]<-t1[q]
								if (p[k,q]>1.02)
									stop("Prawdopodobienstwo >1")
							}
#							tryCatch(if (p[k,q]<0.5){p[k,q]}, error = function(e){print(paste("BLAD   V:", v, "   ; k: ",k,"  ; q: ",q,"  ; dlugoscq: ",dlugoscq,"   ; p[k,q] :",p[k,q],przedzial[1],przedzial[2]))})
#							if(p[k,q]<0.5)
#							{
#								leftNode<-c(leftNode,k)
#							}
#							else
#							{
#								rightNode<-c(rightNode,k)
#							}
						}
						#print("values")
						#try(print(values))
						#print("dlugoscq")
						#print(dlugoscq)
						#print("p")
						#print(p)
						#zastanowic sie
						pkl<-p[,q]*nodeObjects[i,]
						pkr<-(1-p[,q])*nodeObjects[i,]
						pkl[pkl<0]<-0		#  Bledy zaokraglenia
						pkr[pkr<0]<-0		#	Bledy zaokraglenia
						#print("pkl")
						#print(pkl)
						#print("pkr")
						#print(pkr)
						if (sum(pkl)>=treshMin && sum(pkr)>=treshMin)
						{
						  if(useEM){
  							#EM alghoritm
  							Pr<-Pl<-array(0,c(max(classes)))	
  							Er<-El<-array(0,c(ile))	
  							conditional<-array(0,c(max(classes)))
  							increase<-2*abs(epsilon);
  							W_old<--1e10;
  							for(m in 1:max(classes))
  							{
  									Pl[m]<-sum(pkl[classes==m])/sum(pkl)
  									Pr[m]<-sum(pkr[classes==m])/sum(pkr)
  							}
  							#print("pl")
  							#print(Pl)
  							#print("pr")
  							#print(Pr)
  							petlaEM<-0
  							#print(paste("przed em",Sys.time()))
  							while (abs(increase)>=epsilon)
  							{
  								petlaEM<-petlaEM+1
  								for (k in 1:ile)							
  								{
  									# !!!!!!!!!!!!!!!!!! kolejnosc EM , Czy Tu nie ma zadnego sumowania
  									#E Step
  									if (Pl[classes[k]]*pkl[k]+lambda[k]!=0)
  										El[k]<-Pl[classes[k]]*pkl[k]/(Pl[classes[k]]*pkl[k]+lambda[k])
  									else
  										El[k]<-0
  									if (Pr[classes[k]]*pkr[k]+lambda[k]!=0)
  										Er[k]<-Pr[classes[k]]*pkr[k]/(Pr[classes[k]]*pkr[k]+lambda[k])
  									else
  										Er[k]<-0
  								}
  								#print("El")
  								#print(El)
  								#print("Er")
  								#print(Er)
  								for (m in 1:max(classes))							
  								{
  									#M Step
  									if (TRUE) #(sum(El)!=0)
  									{
  										Pl[m]<-sum(El[classes==m])/sum(El)			# ZASTANOWIC SIE, CZY TO OGRANICZYC TYLKO DO OBIEKTOW ZAWARTYCH W WEZLE
  									}
  									else
  									{
  										Pl[m]<-0
  									}
  									if (TRUE) #(sum(Er)!=0)
  									{
  										Pr[m]<-sum(Er[classes==m])/sum(Er)
  									}
  									else
  									{
  										Pr[m]<-0
  									}
  								}
  								W_actual=1
  								#print("lambda")
  								#print(lambda)
  								#print("Pl")
  								#print(Pl)
  								#print("Pr")
  								#print(Pr)
  								if(debug){write(paste("----------------->!!!!!!!!!!!PETLAEM",petlaEM,sep=" : "),file="w.txt",append=TRUE)}
  								for (k in 1:ile)
  								{
  									if(debug){write(paste("----------------->",pkl[k],Pl[classes[k]],pkr[k],Pr[classes[k]],lambda[k],sep=" : "),file="w.txt",append=TRUE)}
  									W_actual<-W_actual*(pkl[k]*Pl[classes[k]]+pkr[k]*Pr[classes[k]]+lambda[k])     # !!!!!!!!!!!!!!!!! czy na pewno?
  								}
  								if(W_actual==1) 
  									stop("Jedynka w W")
  								#print(W_actual)
  								W_actual<-log(W_actual)
  								if(debug){write(paste(petla,v,q,W_actual,sep=" : "),file="w.txt",append=TRUE)}
  								if(debug){write(paste("INCREASE",W_actual-W_old,sep=" : "),file="w.txt",append=TRUE)}
  								#print("W_actual")
  								#print(W_actual)
  								increase<-W_actual-W_old
  								W_old<-W_actual
  								#print(paste("Pytanie : ",q))
  								#print(paste("Zmienna : ",v))
  
  							}
  							#print(paste("po em",Sys.time()))
              }
              if(!useEM){
  							Pr<-Pl<-array(0,c(max(classes)))	
  							Er<-El<-array(0,c(ile))	
  							conditional<-array(0,c(max(classes)))
  							for(m in 1:max(classes))
  							{
  									Pl[m]<-sum(pkl[classes==m])/sum(pkl)
  									Pr[m]<-sum(pkr[classes==m])/sum(pkr)
  							}
                W_actual<-1
		            for (k in 1:ile)
  								{
  								  if(!any(testSet==k)){
    									if(debug){write(paste("----------------->",pkl[k],Pl[classes[k]],pkr[k],Pr[classes[k]],lambda[k],sep=" : "),file="w.txt",append=TRUE)}
    									W_actual<-W_actual*(pkl[k]*Pl[classes[k]]+pkr[k]*Pr[classes[k]]+lambda[k])     # !!!!!!!!!!!!!!!!! czy na pewno?
  									}
  								}
 								W_actual<-log(W_actual)
              }
							#print("Skonczylem While")
							czyByloPytanie<-nodes[nodes[,"variable"]==v,]
							if(is.null(czyByloPytanie))
								czyByloPytanie<-0
							else
							{
								if(is.null(dim(czyByloPytanie)))
									dim(czyByloPytanie)<-c(1,length(czyByloPytanie))
								czyByloPytanie<-sum(czyByloPytanie[,3]==q)
							}
							if (W_actual>W_maximal && czyByloPytanie==0)
							{
								node_W<-i
								W_maximal<-W_actual
								question_type_W<-question_type
								question_variable_W<-v
								if(question_type==1)
								{
									question_value_W<-centers[q]
								}
								else
								if(question_type==2)
								{
									question_value_W<-q
								}
								else
								{
									question_value_W<-q	#przemyslec
								}
								prob_L_W<-pkl
								prob_R_W<-pkr
								for (m in 1:max(classes))
								{
									cond_L_W[m]<-sum(pkl[classes==m])/sum(pkl)
									cond_R_W[m]<-sum(pkr[classes==m])/sum(pkr)
								}
							}
							minimalSizeExists<-TRUE
						}
					}
				}				
			}
		}	
		# Node actualization
		if(W_maximal<treshW)
		{
			minimalWExists<-FALSE
		}
		if(minimalWExists && minimalSizeExists)
		{
			nodes[node_W,1]<-question_type_W
			nodes[node_W,2]<-question_variable_W
			nodes[node_W,3]<-question_value_W
			nodes[node_W,5]<-0
			actualNodes<-nrow(nodes)
			nodes<-rbind(nodes,c(0,0,0,node_W,1))
			nodes<-rbind(nodes,c(0,0,0,node_W,1))
			conditionalProbab[actualNodes+1,]<-cond_L_W
			conditionalProbab[actualNodes+2,]<-cond_R_W
			
			#zmienic
			#print(cond_L_W)
			#print(cond_R_W)
			#print(prob_L_W)
			#print(prob_R_W)
			nodeObjects[actualNodes+1,]<-prob_L_W
			nodeObjects[actualNodes+2,]<-prob_R_W
			#print(nodeObjects)
			
			if (sum(nodeObjects[i,]!=1)!=0)
			{
				#stop (i)
			}

		}
	}
	#print("DEBUG po while");
	#print(nodes);
	#print(conditionalProbab)
	prediction<-NULL
	indivIC<-sdt$indivIC
  for(i in testSet){
  	print(paste("DEBUG test set ",i));
    pr<-rep(0,max(f$classes))
    for(j in 1:nrow(nodes)){
      jj<-j
      if(nodes[j,"terminal"]==1){
        prt<-conditionalProbab[j,]
        while(jj!=1){
        	#print(paste("DEBUG jj ",jj));
        
#        v<-nodes[j,"variable"]
#          if(nodes[j,"question_type"]==1){
#            if(indivIC[i,v,2]<nodes[j,"value"]){
#              pr<-pr+1
#            }
#            else{
#              if(indivIC[i,v,1]>nodes[j,"value"]){
#                pr<-pr+1*conditionalProbab[2*j+1,]
#              }
#              else{
#                pr<-pr+abs(nodes[j,"value"]-indivIC[i,v,2])*conditionalProbab[2*j,]
#                pr<-pr+abs(indivIC[i,v,1]-nodes[j,"value"])*conditionalProbab[2*j+1,]
#              }
#            }
#          }
           prt<-prt*nodeObjects[jj,i]
           jj<-nodes[jj,"parent"]
           #print(c("parent node",jj))
        }
        pr<-pr+prt
	  #DOROBIC DLA INNYCH TYPOW PYTAN !!!!
      }
    }
    #print(paste("obiekt",i))
    #print(pr)
    prediction<-c(prediction,which.max(pr))
  }
	resul<-list(nodes=nodes,nodeObjects=nodeObjects,conditionalProb=conditionalProbab,prediction=prediction,variables=sdt$variables)
	resul
	#print(resul)
}



.nodeNumber<-function(nodes,k){
    if(nodes[k,"parent"]==0){
      resul<-1
    }
    else{
      resul<-.nodeNumber(nodes,nodes[k,"parent"])*2+k%%2
    }
    resul
}

.centerF<-function(i,j,hunit,vunit,levels){
    c(-hunit+j*2*hunit/(2^(i-1)+1),-vunit/2+(levels-i-1)*2*vunit)
}

.plot.decisionTree.SDA<-function(decisionTree.SDA,boxWidth=1,boxHeight=3){
    vunit<-10
    hunit<-15
    dt<-decisionTree.SDA
    levelsA<-NULL
    for(i in 1:nrow(dt$nodes)){
      levels<-1
      if(dt$nodes[i,"terminal"]==1){
        node<-i
        while(dt$nodes[node,"parent"]!=0){
          node<-dt$nodes[node,"parent"]
          levels<-levels+1
        }
      }
      levelsA<-c(levelsA,levels)
      
    }
    levels<-max(levelsA)
    #print(levels)
    plot(NULL,xlim=c(-hunit,hunit),ylim=c(-levels*vunit,levels*vunit),type="n",axes=FALSE,xlab="",ylab="")
    for(k in 1:nrow(dt$nodes)){
        nodeNr<-.nodeNumber(dt$nodes,k)
        i<-as.integer(log(nodeNr,2))+1
        j<-nodeNr%%(2^(i-1))+1
        center<-.centerF(i,j,hunit,vunit,levels)
        rect(center[1]-boxWidth,center[2]-boxHeight,center[1]+boxWidth,center[2]+boxHeight)
        text(center[1],center[2],paste(paste(format(dt$conditionalProb[k,]*100,digits=3),"%",sep=""),collapse="\n"),cex=0.4)  
        if(dt$nodes[k,"terminal"]==0){
          text(center[1],center[2]-vunit,paste(dt$variables[dt$nodes[k,"variable"],"label"],">",dt$nodes[k,"value"]),cex=0.7)  
        }
        if(k!=1){
          nodeNr<-nodeNr%/%2
          i<-as.integer(log(nodeNr,2))+1
          j<-nodeNr%%(2^(i-1))+1
          centerParent<-.centerF(i,j,hunit,vunit,levels)
          lines(c(centerParent[1],center[1]),c(centerParent[2]-boxHeight,center[2]+boxHeight))     
        }
    }
  }


draw.decisionTree.SDA<-function(decisionTree.SDA,boxWidth=1,boxHeight=3){
  .plot.decisionTree.SDA(decisionTree.SDA,boxWidth,boxHeight)
}

boosting.SDA<-function (sdt,formula,testSet, mfinal = 20,...) 
{
    debug=FALSE
    f<-.parseFormula(formula,sdt)
    variableSelection<-f$variableSelection
    classes<-f$classes[testSet]
    formula <- as.formula(formula)
    vardep <- classes
    n <- nrow(sdt$individuals)
    nclases <- length(unique(sdt$variables))
    arboles <- list()
    pesos <- rep(1/n,n)
    replicas <- array(0, c(n, mfinal))
    pred<-NULL
    pond <- rep(0, mfinal)
    for (m in 1:mfinal) {
        boostrap <- sort(sample((1:n)[-testSet],replace=TRUE,prob = pesos[-testSet]))
        currBoostrap<-sort(c(boostrap,testSet))
        btestSet<-NULL
        for(ts in testSet){
          for(i in 1:length(currBoostrap)){
            if(currBoostrap[i]==ts){
              btestSet<-c(btestSet,i)
            }
          }
        }
        currTree<-subsdt.SDA(sdt,currBoostrap)
#        currTestSet<<-btestSet
        if(debug)print("DEBUG przed fit")
        fit<-decisionTree.SDA(currTree,formula,1:length(currBoostrap),...)
        flearn <- fit$prediction
        ind<-rep(0,n)
        ind1 <- as.numeric(f$classes[currBoostrap] != flearn)
#        cind1<<-ind1
        for(t in 1:length(currBoostrap)){
          if(ind1[t]==1){
            ind[currBoostrap[t]]<-1
          }
        }
        err <- sum(ind)/n
        c <- 1/2*log((1 - err)/err)
        pesos <- pesos * exp(c * ind)
        pesos <- pesos/sum(pesos)
        maxerror <- min(1 - max(summary(vardep))/sum(summary(vardep)),0.5)
        if (err >= maxerror) {
            pesos <- rep(1/n, n)
            c <- 0.001
        }
        if (err == 0) {
            pesos <- rep(1/n, n)
            c <- 3
        }
        pond[m] <- c
        fit<-decisionTree.SDA(currTree,formula,btestSet,...)
        arboles[[m]] <- fit
        pred<-cbind(pred,fit$prediction)
    }

#    currPred<<-pred
#    currVardep<<-vardep
    classfinal <- array(0, c(length(testSet), length(unique(f$classes))))
    for (i in sort(unique(f$classes))) {
	  print(i)
        classfinal[, i] <- matrix(as.numeric(pred == sort(unique(f$classes))[i]),
		nrow = length(testSet)) %*% as.vector(pond)

    }
    predclass <- rep("O", length(testSet))
    for (i in 1:length(testSet)) {
        predclass[i] <- as.character(sort(unique(f$classes))[(order(classfinal[i, 
            ], decreasing = TRUE)[1])])
    }
#    cpredClass<<-predclass
#    cvardep<<-vardep
#    ctestSet<<-testSet
    
    error <- 1 - sum(predclass == vardep)/length(testSet)
    ans <- list(formula = formula, trees = arboles, weights = pond, 
        votes = classfinal, class = predclass,error=error)
    ans
}


bagging.SDA<-function (sdt,formula,testSet, mfinal = 20,rf=FALSE,...) 
{
    #sdt<-generate.SO(20,3,4,4,TRUE,separationIndex=20)
    #mfinal<-3
    #rf<-FALSE
    #testSet<-sample(1:60,13)
    #formula<-var_9~.
    debug=FALSE
    f<-.parseFormula(formula,sdt)
    variableSelection<-f$variableSelection
    classes<-f$classes
    formula <- as.formula(formula)
    vardep <- classes
    n <- nrow(sdt$individuals)
#    currN<<-n
#    currTest<<-testSet
    nclases <- length(unique(sdt$variables))
    arboles <- list()
    replicas <- array(0, c(n, mfinal))
    pred<-NULL
    for (m in 1:mfinal) {
        boostrap <- sort(sample((1:n)[-testSet],replace=TRUE))
        currBoostrap<-sort(c(boostrap,testSet))
        btestSet<-NULL
        for(ts in testSet){
          for(i in 1:length(currBoostrap)){
            if(currBoostrap[i]==ts){
              btestSet<-c(btestSet,i)
            }
          }
        }
        currTree<-subsdt.SDA(sdt,currBoostrap)
#        currTestSet<<-btestSet
        if(debug)print("DEBUG przed fit")
        fit<-decisionTree.SDA(currTree,formula,btestSet,rf=rf)
        if(debug)print("DEBUG po fit")
#        currFit<<-fit
#        currReplicas<<-replicas
        arboles[[m]] <- fit
        #flearn <- predict(fit, newdata = data, type = "class")
        #ind <- as.numeric(vardep != flearn)
        #err <- sum(ind)/n
        replicas[, m] <- currBoostrap
        pred<-cbind(pred,fit$prediction)
    }
#    currPred<<-pred
#    currVardep<<-vardep
    classfinal <- array(0, c(length(testSet), length(unique(f$classes))))
    for (i in sort(unique(f$classes))) {
	  print(i)
        classfinal[, i] <- matrix(as.numeric(pred == sort(unique(f$classes))[i]),
		nrow = length(testSet)) %*% rep(1, mfinal)
    }
    predclass <- rep("O", length(testSet))
    for (i in 1:length(testSet)) {
        predclass[i] <- as.character(sort(unique(f$classes))[(order(classfinal[i, 
            ], decreasing = TRUE)[1])])
    }
    tabla <- table(predclass, vardep[testSet], dnn = c("Predicted Class", 
        "Observed Class"))
    error <- 1 - sum(predclass == vardep[testSet])/length(testSet)
    output <- list(class = predclass, confusion = tabla, error = error,pred=pred,classfinal=classfinal)
}

random.forest.SDA<-function (sdt,formula,testSet, mfinal = 100,...) {
  bagging.SDA(sdt,formula,testSet,rf=TRUE)
}

subsdt.SDA<-function(sdt,objectSelection=1:nrow(sdt$individuals)){
  resul<-sdt
  resul$individuals<-resul$individuals[objectSelection,]
  resul$indivIC<-resul$indivIC[objectSelection,,]
  dimnames(resul$indivIC)[[1]]<-objectSelection
  indivN<-NULL
  indivNM<-NULL
  for(os in objectSelection){
    if(!is.null(resul$indivN))
    for(i in 1:nrow(resul$indivN)){
      if(any(resul$indivN[i,"indiv"]==os)){
      indivN<-rbind(indivN,resul$indivN[i,])
      }
    }
    if(!is.null(resul$indivNM))
    for(i in 1:nrow(resul$indivNM)){
      if(any(resul$indivNM[i,"indiv"]==os)){
      indivNM<-rbind(indivNM,resul$indivNM[i,])
      }
    }
  }
  detailsListNom<-NULL
  dn<-0;
  dnm<-0;
  resul$detailsIC<-as.matrix(resul$detailsIC)
  resul$detailsN<-as.matrix(resul$detailsN)
  for(i in 1:nrow(resul$variables)){
    
    if(resul$variables[i,"type"]=="MN" || resul$variables[i,"type"]=="N"){
      t<-unique(as.matrix(indivN)[as.matrix(indivN)[,"variable"]==i,"value"])
      #print(t)
      for(tt in t){
        detailsListNom<-rbind(detailsListNom,c(resul$variables[i,5],tt,paste("x",tt,sep="_"),paste("x",tt,sep="_")))
      }
        resul$detailsN[resul$variables[i,5],3]<-length(t)
      dn<-dn+1;
    }
    if(resul$variables[i,"type"]=="IC"){
      
      #print("DEBUG IC START")
      resul$detailsIC[resul$variables[i,5],3]<-min(resul$indivIC[,i,1])
      #print("DEBUG IC MIDDLE")
      resul$detailsIC[resul$variables[i,5],4]<-max(resul$indivIC[,i,2])
      #print("DEBUG IC STOP")
    }
#    if(resul$variables[i,"type"]=="NM" || resul$variables[i,"type"]=="N"){
#      dnm<-dn+1;
#      detailsN<-rbind(detailsN,c(0,0,0)
#    }
    #if(resul$variables[i,"type"]=="NM" || resul$variables[i,"type"]=="N"){
    #  detailsN<-rbind(detailsN,
      
    #}
  }
  #print("po")
  detailsListNom<-as.data.frame(detailsListNom)
  names(detailsListNom)<-c("details_no","num","name","label")
  names(indivN)<-c("indiv","variable","value")
  resul$detailsIC<-as.data.frame(resul$detailsIC)
  resul$detailsN<-as.data.frame(resul$detailsN)
  #dimnames(resul$individuals)[[1]]<-1:length(boostrap)
  #dimnames(resul$indivN)[[1]]<-1:nrow(sdt$indivN)
  resul$detailsListNom<-detailsListNom
  
  # update of indexes
  if(!is.null(indivN)){
  newIndices<-as.matrix(resul$individuals)
  newIndices[,1]<-1:nrow(resul$individuals)
  
  individuals<-as.matrix(resul$individuals)
  for(i in 1:nrow(individuals)){
    #print(i)
    individuals[i,1]<-  newIndices[newIndices[,"name"]==paste("",individuals[i,"name"],"",sep=""),1]
  }
  for(i in 1:nrow(indivN)){
    indivN[i,1]<-  newIndices[paste("",indivN[i,1],"",sep=""),1]
  }
  dimnames(individuals)[[1]]<-1:nrow(individuals)
  dimnames(indivN)[[1]]<-1:nrow(indivN)
  indivN<-as.matrix(indivN)
  resul$individuals<-as.data.frame(individuals)
  resul$indivN<-as.data.frame(indivN)
  }
  resul
}
