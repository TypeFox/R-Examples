#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Odległości dla danych symbolicznych
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach projektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach li-cencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************

.nomValues<-function(indivN,i,k){
indivN_t1<-indivN[indivN[,"indiv"]==i ,]		
indivN1<-indivN_t1[indivN_t1[,"variable"]==k  ,]		
as.matrix(indivN1)
}

.p<-function(list,l){
  resul<-0
  if(sum(list[,"value"]==l)!=0){
    if(any(dimnames(list)[[2]]=="frequency")){
      resul<-list[list[,"value"]==l,]["frequency"]
    }
    else{
      resul<-1/nrow(list)
    }
  }
  resul
}

.pp<-function(list,indivN,k,l){
  resul<-0
  if(sum(list[,"value"]==l)!=0){
    indivN1<-indivN[indivN[,"variable"]==k  ,]
    resul<-1/sum(indivN1[,"value"]==l)
  }
  resul
}


.gsum<-function(table.Symbolic,i,j,k){
indivIC<-table.Symbolic$indivIC
indivN<-table.Symbolic$indivN
variables<-table.Symbolic$variables
  gsum=0
  if(is.null(k))stop("wrong k null")
  if(length(k)!=1)stop(paste("wrong k",k))
  if(variables[k,"type"]=="IC"){
      gsum<-max(indivIC[i,k,2],indivIC[j,k,2])-min(indivIC[i,k,1],indivIC[j,k,1])
  }
  if (variables[k,"type"]=="N" || variables[k,"type"]=="MN"){
		indivN_t1<-indivN[indivN[,"indiv"]==i ,]		
		indivN_t2<-indivN[indivN[,"indiv"]==j ,]		
		indivN1<-indivN_t1[indivN_t1[,"variable"]==k  ,]		
		indivN2<-indivN_t2[indivN_t2[,"variable"]==k  ,]
		#print(paste("nrow(indivN1)",nrow(indivN1)))
		#print(paste("nrow(indivN2)",nrow(indivN2)))

		gsum<-nrow(indivN1)+nrow(indivN2)-.gprod(table.Symbolic,i,j,k)
  }
  gsum
}


.gprod<-function(table.Symbolic,i,j,k){
indivIC<-table.Symbolic$indivIC
indivN<-table.Symbolic$indivN
variables<-table.Symbolic$variables
  gprod=0
    if(is.null(k))stop("wrong k null")
  if(length(k)!=1)stop(paste("wrong k",k))

  if(variables[k,"type"]=="IC"){
      if(indivIC[i,k,1]>=indivIC[j,k,2] || indivIC[i,k,2]<=indivIC[j,k,1]){
           gprod<-0
      }
      else{
          if (indivIC[i, k, 1] <= indivIC[j, k, 1]) {
      
              if (indivIC[i, k, 2] <= indivIC[j, k, 2]) 
                  gprod <- abs(indivIC[i, k, 2] - indivIC[j, k, 
                    1])
              else gprod <- abs(indivIC[j, k, 2] - indivIC[j, k, 
                  1])
          }
          else {
          #print("mniejszy")
              if (indivIC[j, k, 2] <= indivIC[i, k, 2]) 
                  gprod <- abs(indivIC[j, k, 2] - indivIC[i, k, 
                    1])
              else gprod <- abs(indivIC[i, k, 2] - indivIC[i, k, 
                  1])
          }
      }
  }
  if (variables[k,"type"]=="N" || variables[k,"type"]=="MN"){
		gprod<-0
		indivN_t1<-indivN[indivN[,"indiv"]==i ,]		
		indivN_t2<-indivN[indivN[,"indiv"]==j ,]		
		indivN1<-indivN_t1[indivN_t1[,"variable"]==k  ,]		
		indivN2<-indivN_t2[indivN_t2[,"variable"]==k  ,]
		if (!is.null(indivN1[,"value"]) && !is.null(indivN1[,"value"]))
		{
			for (z  in 1:(nrow(indivN1)))
			{
				if ( sum(indivN2[,"value"]==indivN1[z,"value"])!=0)
				{
					gprod<-gprod+1
				}			
			}
		}
	}
  gprod
}



.gabs<-function(table.Symbolic,i,k){
indivIC<-table.Symbolic$indivIC
indivN<-table.Symbolic$indivN
variables<-table.Symbolic$variables
  gabs<-0
  if(is.null(k))stop("wrong k null")
  if(length(k)!=1)stop(paste("wrong k",k))
  if(variables[k,"type"]=="IC"){
      gabs<-abs(indivIC[i,k,2]-indivIC[i,k,1])
  }
  if (variables[k,"type"]=="N" || variables[k,"type"]=="MN"){
		indivN_t1<-indivN[indivN[,"indiv"]==i ,]		
		indivN1<-indivN_t1[indivN_t1[,"variable"]==k  ,]		
		gabs<-nrow(indivN1)
  }
  gabs
}

.gdomainLength<-function(table.Symbolic,k){
  indivIC<-table.Symbolic$indivIC
  indivN<-table.Symbolic$indivN
  variables<-table.Symbolic$variables
  if(is.null(k))stop("wrong k null")
  if(length(k)!=1)stop(paste("wrong k",k))
  if(variables[k,"type"]=="IC"){
      gdomainLength<-abs(max(indivIC[,k,1])-min(indivIC[,k,1]))
  }
  if (variables[k,"type"]=="N"){
		gdomainLength<-length(unique(indivN[indivN[,"variable"]==k,]))
  }
  if (variables[k,"type"]=="MN"){
		gdomainLength<-length(unique(indivN[indivN[,"variable"]==k,]))
  }
  gdomainLength
}


.dpobject<-function(table.Symbolic,i,variableSelection){
  dp<-1
  for(k in variableSelection){
    dp<-dp*.gabs(table.Symbolic,i,k)
  }
  dp
}

.dpsum<-function(table.Symbolic,i,j,variableSelection){
  dp<-1
  for(k in variableSelection){
    dp<-dp*.gsum(table.Symbolic,i,j,k)
  }
  dp
}


.dpprod<-function(table.Symbolic,i,j,variableSelection){
  dp<-1
  for(k in variableSelection){
    dp<-dp*.gprod(table.Symbolic,i,j,k)
  }
  dp
}

.dpmax<-function(table.Symbolic,variableSelection){
  dpobjects<-NULL
  for(i in 1:length(table.Symbolic$individuals)){
    dpobjects<-c(dpobjects,.dpobject(table.Symbolic,i,variableSelection))
  }
  max(dpobjects)
}

.C1alpha<-function(table.Symbolic,i,j,k){
    .gprod(table.Symbolic,i,j,k)
}
.C1beta<-function(table.Symbolic,i,j,k){
  indivIC<-table.Symbolic$indivIC
  indivN<-table.Symbolic$indivN
  variables<-table.Symbolic$variables
  if(table.Symbolic$variables[k,"type"]=="N" || table.Symbolic$variables[k,"type"]=="MN"){
    beta<-.gsum(table.Symbolic,i,j,k)-.gabs(table.Symbolic,j,k)
  }
  if(table.Symbolic$variables[k,"type"]=="IC"){
  if (indivIC[i,k,1]<=indivIC[j,k,1])
			{
				if (indivIC[i,k,2]<=indivIC[j,k,2])
					beta<-abs(indivIC[i,k,1]-indivIC[j,k,1])
				else
					beta<-abs(indivIC[j,k,2]-indivIC[j,k,2])+abs(indivIC[i,k,2]-indivIC[i,k,2])
			}
			else
			{
				if (indivIC[j,k,2]<=indivIC[i,k,2])
					beta<-abs(indivIC[i,k,2]-indivIC[j,k,2])
				else
					beta<-0
			}

  }  
  beta
}

.C1gamma<-function(table.Symbolic,i,j,k){
  indivIC<-table.Symbolic$indivIC
  indivN<-table.Symbolic$indivN
  variables<-table.Symbolic$variables
  if(table.Symbolic$variables[k,"type"]=="N" || table.Symbolic$variables[k,"type"]=="MN"){
    gamma<-.gsum(table.Symbolic,i,j,k)-.gabs(table.Symbolic,j,k)
  }
  if(table.Symbolic$variables[k,"type"]=="IC"){
  if (indivIC[j,k,1]<=indivIC[j,k,1])
			{
				if (indivIC[j,k,2]<=indivIC[i,k,2])
					gamma<-abs(indivIC[i,k,1]-indivIC[j,k,1])
				else
					gamma<-abs(indivIC[j,k,2]-indivIC[j,k,2])+abs(indivIC[i,k,2]-indivIC[i,k,2])
			}
			else
			{
				if (indivIC[i,k,2]<=indivIC[j,k,2])
					gamma<-abs(indivIC[i,k,2]-indivIC[j,k,2])
				else
					gamma<-0
			}
  }  
  gamma
}

.probDist<-function(indivNM,i,j,k,probType,s,p){
  resul<-0
  list1<-.nomValues(indivNM,i,k)
  #print(c(i,j,k,probType,s,p))
  #print(list1)
  list2<-.nomValues(indivNM,j,k)
  #print(list2)
  if(probType=="J"){
    for(l in 1:max(c(list1[,"value"],list2[,"value"]))){
      if(.p(list1,l)!=0 && .p(list2,l)!=0){
        resul<-resul+(.p(list2,l)*log(.p(list2,l)/.p(list1,l))+.p(list1,l)*log(.p(list1,l)/.p(list2,l)))/2
      }
    }
  }
  if(probType=="CHI"){
    for(l in 1:max(c(list1[,"value"],list2[,"value"]))){
      if(.p(list1,l)!=0){
        resul<-resul+(.p(list1,l)-.p(list2,l))^2/.p(list1,l)
      }
    }
  }
  if(probType=="CHER"){
    for(l in 1:max(c(list1[,"value"],list2[,"value"]))){
      resul<-resul+(.p(list1,l)^(1-s)*.p(list2,l)^s)
    }
    resul<--log(resul)
  }
  if(probType=="REN"){
    for(l in 1:max(c(list1[,"value"],list2[,"value"]))){
      resul<-resul+(.p(list1,l)^(1-s)*.p(list2,l)^s)
    }
    resul<--log(resul)/(s-1)
  }
  if(probType=="LP"){
    for(l in 1:max(c(list1[,"value"],list2[,"value"]))){
      resul<-resul+abs(.p(list1,l)-.p(list2,l))^p
    }
    resul<-resul^(1/p)
  }
  resul
}


dist.SDA<-function(table.Symbolic,type="U_2",subType=NULL,gamma=0.5,power=2,probType="J",probAggregation="P_1",s=0.5,p=2,variableSelection=NULL,weights=NULL){
  types<-c(paste("U",2:4,sep="_"),"C_1",paste("SO",1:5,sep="_"),"H","L_1","L_2")
  subTypes<-c(paste("D",1:5,sep="_"))
  probTypes=c("J","CHI2","REN","CHER","LP")
  #probAggregation=
  if(sum(types==type)==0){
    stop(paste("Dissimilarity type should be one of the following:",paste(types,collapse=",")))
  }
  if(!is.null(subType) && sum(subTypes==subType)==0){
    stop(paste("Dissimilarity subType for(C_1) and (SO_1) distance should be one of the following:",paste(subTypes,collapse=",")))
  }
  if(is.null(subType) && (type=="C_1" || type=="SO_1")){
    stop(paste("Dissimilarity subType for(C_1) and (SO_1) distance should be one of the following:",paste(subTypes,collapse=",")," but not null"))
  }

  if(sum(types==type)==0){
    stop(paste("Dissimilarity type should be one of the following:",paste(types,collapse=",")))
  }
  if (!.is.symbolic(table.Symbolic))
    stop ("Cannot count dissimilarity measures on non-symbolic data")
  #print("step 0")
  indivIC<-table.Symbolic$indivIC
  indivN<-table.Symbolic$indivN
  indivNM<-table.Symbolic$indivNM
  variables<-table.Symbolic$variables
  individuals<-table.Symbolic$individuals
  individualsNo<-nrow(individuals)
  variablesNo<-nrow(variables)
  #print("step 1")
  if(is.null(weights)) weights<-rep(1,nrow(variables))
  if(length(weights)==1) weights<-rep(weights,nrow(variables))
  if(length(weights)!=nrow(variables)){
    stop ("Weights vector should have the same length as the number of variables")
  }
  #print("step 2")
  if (is.null(variableSelection))
    variableSelection=(1:variablesNo)
  distS<-array(0,c(individualsNo,individualsNo))
  for (i in 1:(individualsNo-1))
  for (j in (i+1):individualsNo){
    if(type=="U_2"){
      D<-0
      for(k in variableSelection){
        d<-.gsum(table.Symbolic,i,j,k)-.gprod(table.Symbolic,i,j,k)+gamma*(2*.gprod(table.Symbolic,i,j,k)-.gabs(table.Symbolic,i,k)-.gabs(table.Symbolic,j,k))
        D<-D+d^power
        if(variables[k,"type"]=="NM"){
          D<-D+.probDist(indivNM,i,j,k,probType,s,p)^power
        }
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="U_3"){
      D<-0
      for(k in variableSelection){
        #print(paste("krok",i,j,k))
        #print(paste("gsum",.gsum(table.Symbolic,i,j,k)))
        #print(paste("gprod",.gprod(table.Symbolic,i,j,k)))
        #print(paste("gdomain",.gdomainLength(table.Symbolic,k)))
        d<-.gsum(table.Symbolic,i,j,k)-.gprod(table.Symbolic,i,j,k)+gamma*(2*.gprod(table.Symbolic,i,j,k)-.gabs(table.Symbolic,i,k)-.gabs(table.Symbolic,j,k))
        d<-d/.gdomainLength(table.Symbolic,k)
        #if(is.inf(d))stop("niesk")
        #print(d)
        D<-D+d^power
        if(variables[k,"type"]=="NM"){
          D<-D+.probDist(indivNM,i,j,k,probType,s,p)^power
        }
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="U_4"){
      D<-0
      for(k in variableSelection){
      print(paste("k=",k))
        d<-.gsum(table.Symbolic,i,j,k)-.gprod(table.Symbolic,i,j,k)+gamma*(2*.gprod(table.Symbolic,i,j,k)-.gabs(table.Symbolic,i,k)-.gabs(table.Symbolic,j,k))
        d<-d/.gdomainLength(table.Symbolic,i)
        D<-D+weights[k]*d^power
        if(variables[k,"type"]=="NM"){
          D<-D+weights[k]*.probDist(indivNM,i,j,k,probType,s,p)^power
        }
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="C_1"){
      D<-0
      for(k in variableSelection){
        alpha=.C1alpha(table.Symbolic,i,j,k)
        beta=.C1gamma(table.Symbolic,i,j,k)
        gamma=.C1gamma(table.Symbolic,i,j,k)
        if(subType=="D_1"){
          d<-alpha/(alpha+beta+gamma)
        }
        if(subType=="D_2"){
          d<-2*alpha/(2*alpha+beta+gamma)
        }
        if(subType=="D_3"){
          d<-alpha/(alpha+2*beta+2*gamma)
        }
        if(subType=="D_4"){
          d<-alpha/(2*alpha+2*beta)+alpha(2*alpha+2*gamma)
        }
        if(subType=="D_5"){
          d<-alpha/(sqrt(alpha+beta)*sqrt(alpha+gamma))
        }
        if(!is.nan(d))
        D<-D+(d^power)/length(variables)
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="SO_1"){
      D<-0
      for(k in variableSelection){
        alpha=as.numeric(.C1alpha(table.Symbolic,i,j,k))
        bbeta=as.numeric(.C1beta(table.Symbolic,i,j,k))
        gamma=as.numeric(.C1gamma(table.Symbolic,i,j,k))
        if(subType=="D_1"){
          d<-alpha/(alpha+bbeta+gamma)
        }
        if(subType=="D_2"){
          d<-2*alpha/(2*alpha+bbeta+gamma)
        }
        if(subType=="D_3"){
          d<-alpha/(alpha+2*bbeta+2*gamma)
        }
        if(subType=="D_4"){
          d<-alpha/(2*alpha+2*bbeta)+alpha(2*alpha+2*gamma)
        }
        if(subType=="D_5"){
          d<-alpha/(sqrt(alpha+bbeta)*sqrt(alpha+gamma))
        }
        if(!is.nan(d))
        D<-D+weights[k]*(d^power)
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="SO_2"){
      D<-0
      for(k in variableSelection){
        d<-.gsum(table.Symbolic,i,j,k)-.gprod(table.Symbolic,i,j,k)+gamma*(2*.gprod(table.Symbolic,i,j,k)-.gabs(table.Symbolic,i,k)-.gabs(table.Symbolic,j,k))
        d<-d/.gsum(table.Symbolic,i,j,k)
        D<-D+(d^power)/length(variables)
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="SO_3"){
      D<-0
      for(k in variableSelection){
        d<-.dpsum(table.Symbolic,i,j,variableSelection)-.dpprod(table.Symbolic,i,j,variableSelection)+gamma*(2*.dpprod(table.Symbolic,i,j,variableSelection)-.dpobject(table.Symbolic,i,variableSelection)-.dpobject(table.Symbolic,j,variableSelection))
        D<-d
      }
      distS[i,j]<-D		
      distS[j,i]<-distS[i,j]
    }
    if(type=="SO_4"){
      D<-0
      for(k in variableSelection){
        d<-.dpsum(table.Symbolic,i,j,variableSelection)-.dpprod(table.Symbolic,i,j,variableSelection)+gamma*(2*.dpprod(table.Symbolic,i,j,variableSelection)-.dpobject(table.Symbolic,i,variableSelection)-.dpobject(table.Symbolic,j,variableSelection))
        D<-d/.dpmax(table.Symbolic,k)
      }
      distS[i,j]<-D/.dpmax(table.Symbolic,variableSelection)		
      distS[j,i]<-distS[i,j]
    }
    if(type=="SO_5"){
      D<-0
      for(k in variableSelection){
        d<-.dpsum(table.Symbolic,i,j,variableSelection)-.dpprod(table.Symbolic,i,j,variableSelection)+gamma*(2*.dpprod(table.Symbolic,i,j,variableSelection)-.dpobject(table.Symbolic,i,variableSelection)-.dpobject(table.Symbolic,j,variableSelection))
        D<-d/.dpsum(table.Symbolic,i,j,variableSelection)
      }
      distS[i,j]<-D		
      distS[j,i]<-distS[i,j]
    }
    if(type=="L_1" || type=="L_2"){
      if(type=="L_1"){
        L_i<-1
      }
      else{
        L_i<-2
      }
      D<-0
      for(k in variableSelection){
        if(variables[k,"type"]=="IC"){
          D<-D+abs(indivIC[i,k,1]-indivIC[j,k,1])^2++(indivIC[i,k,2]-indivIC[j,k,2])^L_i
        }
        if(variables[k,"type"]=="MN" || variables[k,"type"]=="N" || variables[k,"type"]=="NM"){
          if(variables[k,"type"]=="NM"){
            list1<-.nomValues(indivNM,i,k)
            list2<-.nomValues(indivNM,j,k)
          }
          else{
            list1<-.nomValues(indivN,i,k)
            list2<-.nomValues(indivN,j,k)
          }
          for(l in 1:max(c(list1[,"value"],list2[,"value"]))){
            if(variables[k,"type"]=="NM"){
              D<-D+abs(.p(list1,l)-.p(list2,l))^L_i
            }
            else{
              D<-D+abs(.pp(list1,indivN,k,l)-.pp(list2,indivN,k,l))^L_i
            }
          }
        }
        #print(paste(variables[k,"type"],D))
      }
      distS[i,j]<-D^(1/power)		
      distS[j,i]<-distS[i,j]
    }
  }
	as.dist(distS)
}




