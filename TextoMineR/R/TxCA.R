TxCA <-
function(DocVar, DocTerm, num.agg=NULL, idiom="en", Fmin=5, Dmin=5, Fmax=NULL,
 equivalence=NULL, stop.word.user=NULL, stop.word.tm=FALSE,
  lmd=3, lmk=3,  ncp=5, row.sup = NULL, col.sup = NULL, quanti.sup=NULL,
   quali.sup = NULL, graph = TRUE, axes = c(1,2), row.w = NULL){
    
  if(nrow(DocTerm)<3 & ncol(DocTerm)<3 | nrow(DocTerm)>3 & ncol(DocTerm)<2)
  stop("DocTerm must have at least three rows and three columns")
  if(nrow(DocTerm)!= nrow(DocVar))
  stop("the number of  rows of DocTerm  is different to the number of  rows of DocVar")
  if (!is.null(num.agg)&!is.null(quanti.sup))   
    stop("When there is aggregation quanti.sup  should be NULL")
  
    if (!is.null(equivalence)){
	DocTermB<-DocTerm
      equivalence<-equivalence[which(equivalence[,1]%in%colnames(DocTermB)),]
      if(nrow(equivalence)<1 | ncol(equivalence)<2 )
      stop("equivalence must have at least one rows and two columns")
	for(i in 1:nrow(equivalence)){
      if (equivalence[i,1]%in%colnames(DocTermB) & equivalence[i,2]%in%colnames(DocTermB) ){
        DocTermB[,which(colnames(DocTermB)==equivalence[i,2])]<-(DocTermB[,which(colnames(DocTerm)==equivalence[i,1])]
							  + DocTermB[,which(colnames(DocTerm)==equivalence[i,2])])
         }else{
            pos<-which(colnames(DocTermB)==equivalence[i,1]) 
            colnames(DocTermB)[pos]<-as.character(equivalence[i,2])
       }
      }
       DocTerm<-DocTermB[,which(!colnames(DocTermB)%in%equivalence[,1])]
     }

	dtmM<-as.matrix(DocTerm)
      FreqWord<-apply(dtmM,MARGIN=2,FUN=sum)
      dtmA<-dtmM
      dtmA[dtmM[,]>0]<-1
      NumDoc<-apply(dtmA,MARGIN=2,FUN=sum)
	DocTermR<-DocTerm[,which(FreqWord>=Fmin&NumDoc>=Dmin)]
      if(nrow(DocTermR)<3 & ncol(DocTermR)<3 | nrow(DocTermR)>3 & ncol(DocTermR)<2)
      stop(cat("\n After selection the dimension of data frame is ",dim(DocTermR),"
      It must have at least three rows and three columns.\n"))
  
     if(!is.null(Fmax)){
	  FrecMax<-apply(DocTermR,2,sum)
        PalSel<-which(FrecMax>=Fmax)
        DocTermR<-DocTermR[,-PalSel]
   }

  if (!is.null(stop.word.user)){
	DocTermR<-DocTermR[,which(!colnames(DocTermR)%in%stop.word.user)]
   }

  if(stop.word.tm){
  	stopword<-stopwords(idiom)
	DocTermR<-DocTermR[,which(!colnames(DocTermR)%in%stopword)]
   }
 
  if(nrow(DocTermR)<3 & ncol(DocTermR)<3 | nrow(DocTermR)>3 & ncol(DocTermR)<2)
  stop(cat("\n After selection the dimension of data frame is ",dim(DocTermR),"
  It must have at least three rows and three columns.\n"))

      TableSummary<-matrix(c(sum(DocTerm),sum(DocTermR),nrow(DocTerm),nrow(DocTermR),
       ncol(DocTerm),ncol(DocTermR),round(sum(DocTerm)/nrow(DocTerm),2),
       round(sum(DocTermR)/nrow(DocTermR),2)), nrow =2, ncol = 4, byrow = FALSE,
       dimnames = list(c("Before-selection", "After-selection"),
                    c("Occurrences", "Documents", "Words", "Mean-length")))     
      dtmMR<-as.matrix(DocTermR)
      Nfreqword<-apply(dtmMR,MARGIN=2,FUN=sum)
      dtmAR<-dtmMR
      dtmAR[dtmMR[,]>0]<-1
      Ndocword<-apply(dtmAR,MARGIN=2,FUN=sum)
      Nfreqwords<-as.data.frame(Nfreqword)
      Ndocwords<-as.data.frame(Ndocword)
      Table<-cbind(Nfreqwords,Ndocwords)
	colnames(Table)<-c("Frequency", "N.Documents")
   	Glossary<-with(Table,Table[order(Nfreqwords,Ndocwords,decreasing=TRUE),])

      if (!is.null(col.sup)){
         if(is.character(col.sup))
           col.sup<-which(colnames(DocTermR)%in%col.sup)
          if(is.numeric(col.sup))
	     col.sup<-col.sup
       }

       base<-DocVar
  if(!is.null(num.agg)){ 
         if(is.character(num.agg))
           num.agg<-which(colnames(base)%in%num.agg)
          if(is.numeric(num.agg))
	     num.agg<-num.agg
	 agg<-as.factor(base[,num.agg])
       DocTermRA<-DocTermR
     	 dis.X<-tab.disjonctif(agg)
	 Tagreg<-t(DocTermRA)%*%dis.X
	 Tagreg<-t(Tagreg)
	 Torig<-t(DocTerm)%*%dis.X
	 Torig<-t(Torig)
         if (!is.null(row.sup)){
         if(is.character(row.sup))
          row.sup<-which(rownames(Tagreg)%in%row.sup)
          if(is.numeric(row.sup))
	     row.sup<-row.sup
    } 

    Tagreg1<-Tagreg    
    if (!is.null(quali.sup)){
	  if(is.character(quali.sup))
           quali.sup<-which(colnames(base)%in%quali.sup)
         if(is.numeric(quali.sup))
	    quali.sup<-quali.sup 
       for(i in 1:length(quali.sup)){
       nquali<-quali.sup[i]
       qag<-as.factor(base[,nquali])
       DocTermRA<-DocTermR
     	 dis.X<-tab.disjonctif(qag)
	 Tquali<-t(DocTermRA)%*%dis.X
	 Tquali<-t(Tquali)
       if(i==1){
        Tagquali<-Tquali
        }else{
        Tagquali<-rbind(Tagquali,Tquali) 
        }
       }
       TagregEx<-rbind(Tagreg,Tagquali)
       if (!is.null(row.sup)&!is.null(quali.sup)){
       row.sup=c(row.sup,(nrow(Tagreg)+1):nrow(TagregEx))
       }else{
         row.sup=c((nrow(Tagreg)+1):nrow(TagregEx))
        }
       	Tagreg1<-TagregEx
    }
		ncp <- min(ncp, (nrow(Tagreg1) - 1), (ncol(Tagreg1) - 1))
      	res.ca<-CA(Tagreg1[apply(Tagreg1,1,sum)>0,],ncp, row.sup, col.sup,
           quanti.sup=NULL, quali.sup=NULL,graph=FALSE) 

     Categ<- as.data.frame(base[,num.agg])
     colnames(Categ)<-"GrupoCat"
     DocTermRA<- cbind(DocTermRA,Categ)
     DocTermRA$GrupoCat<-as.factor(DocTermRA$GrupoCat)
     DocTermOrg<- cbind(DocTerm,Categ)
     DocTermOrg$GrupoCat<-as.factor(DocTermOrg$GrupoCat)
     NumClas<-as.matrix(summary(DocTermRA[,ncol(DocTermRA)]))
     NumForms<- as.matrix(apply(Torig,1,sum))
     Ptot<- round((NumForms/sum(DocTerm))*100,2)
     NumFormsRet<- as.matrix(apply(Tagreg,1,sum))
     SumIndGen<-apply(DocTerm[,1:ncol(DocTermR)],1,sum) 
     IndNoRes<-which(SumIndGen==0)
   
   listG<-split(DocTermRA,DocTermRA[,ncol(DocTermRA)])        
   listOrg<-split(DocTermOrg,DocTermOrg[,ncol(DocTermOrg)])
    NoRes<-vector()
    NDistF<-vector()
    NoResOr<-vector()
    NDistFOr<-vector()
     for(i in 1:length(listG)){ 
            Gp1<- as.data.frame(listG[[i]])
            Gp<-Gp1[,1:(ncol(Gp1)-1)]
            SumInd<-length(which(apply(Gp,1,sum)==0)) 
              NoRes[i]<-SumInd
               NDistG<-ncol(Gp[,which(apply(Gp,2,sum)>0)])
               NDistF[i]<-NDistG
             GpOr1<- as.data.frame(listOrg[[i]])
             GpOr<-GpOr1[,1:ncol( GpOr1)-1]
		 SumIndOr<-length(which(apply(GpOr,1,sum)==0))  
              NoResOr[i]<-SumIndOr
              NDistGOr<-ncol(GpOr[,which(apply(GpOr,2,sum)>0)])
              NDistFOr[i]<-NDistGOr
              }
         NoRes<-as.matrix(NoRes)
	   NoResOr<-as.matrix(NoResOr)
	   NumResOr<- NumClas-NoResOr 
         NumRes<- NumClas-NoRes
          MedRes<-round((NumForms/NumResOr),2)
         TableLexCat<-cbind(NumClas, NumResOr)
	   TableLexCat<-cbind(TableLexCat, NumRes)
         overall<-apply(TableLexCat,2,sum)
         TableLexCat<-rbind(TableLexCat, overall)
         colnames(TableLexCat)<- c("Documents", "Non-empty-before", "Non-empty-after")
         TableDistForm1<-cbind(NumForms, Ptot)
         TableDistForm2<- cbind(MedRes,NDistFOr)
         TableDistForm3<- cbind(TableDistForm1, TableDistForm2)
         TableDistForm4<-cbind(TableDistForm3,NDistF )
         TableDistForm<-cbind(TableDistForm4, NumFormsRet )
	   overall<-apply(TableDistForm,2,sum)
         TableDistForm<-rbind(TableDistForm, overall)
         TableDistForm[nrow(TableDistForm),5]<-TableSummary[2,3]
         TableDistForm[nrow(TableDistForm),4]<-TableSummary[1,3]
         TableDistForm[nrow(TableDistForm),3]<- TableSummary[1,4]
         colnames( TableDistForm)<- c("Occurrences-before", "Proportion", "Mean-length", 
                 "Words-before", "Words-after","Occurrences-after" )
         TableDistForm<-TableDistForm[,c(1,2,6,3:5)]
        res.agg<-list(Tagreg=Tagreg,TableLexCat=TableLexCat,TableDistForm=TableDistForm)
	}else{
       if (!is.null(quali.sup)|!is.null(quanti.sup)){     
       if (!is.null(quali.sup)){
	  if(is.character(quali.sup))
           quali.sup<-which(colnames(base)%in%quali.sup)
         if(is.numeric(quali.sup))
	    quali.sup<-quali.sup
          }
       if (!is.null(quanti.sup)){
	  if(is.character(quanti.sup))
           quanti.sup<-which(colnames(base)%in%quanti.sup)
         if(is.numeric(quanti.sup))
	    quanti.sup<-quanti.sup 
           }
		DocTermRAC<-DocTermR[apply(DocTermR,1,sum)>0,]
            DocTermRAC<-DocTermRAC[,apply(DocTermRAC,2,sum)>0]
 		base<-as.data.frame(base[rownames(DocTermRAC),])
         if (!is.null(quali.sup)&!is.null(quanti.sup))
		VarAgreg<-base[,c(quali.sup,quanti.sup)]
	   if (!is.null(quali.sup)& is.null(quanti.sup))
             VarAgreg<-base[,quali.sup]
         if (is.null(quali.sup)& !is.null(quanti.sup))
 		  VarAgreg<-base[,quanti.sup]
		 VarAgreg<-as.data.frame(VarAgreg)
            DocTermRA<-cbind.data.frame(DocTermRAC, VarAgreg) 
        if (!is.null(quali.sup)&!is.null(quanti.sup)){
 		q<-length(quali.sup)
            quali.sup<-(ncol(DocTermR)+1):(ncol(DocTermRA)-length(quanti.sup))
            quanti.sup<-(ncol(DocTermR)+1+q):ncol(DocTermRA)
	  }
	if (!is.null(quali.sup)& is.null(quanti.sup))
          quali.sup<-(ncol(DocTermR)+1):ncol(DocTermRA)
        if (is.null(quali.sup)& !is.null(quanti.sup))
 	   quanti.sup<-(ncol(DocTermR)+1):ncol(DocTermRA)

          if (!is.null(row.sup)){
             if(is.character(row.sup))
             row.sup<-which(rownames(DocTermRA)%in%row.sup)
             if(is.numeric(row.sup))
	        row.sup<-row.sup
            }
 		ncp <- min(ncp, (nrow(DocTermRA) - 1), (ncol(DocTermRA) - 1))
            res.ca<-CA(DocTermRA, ncp, row.sup, col.sup,
	                quanti.sup, quali.sup,graph)
		   res.agg=NULL
               Tagreg=NULL
         }else{
      if (!is.null(row.sup)){
         if(is.character(row.sup))
          row.sup<-which(rownames(DocTermR)%in%row.sup)
          if(is.numeric(row.sup))
	     row.sup<-row.sup
       }
	 ncp <- min(ncp, (nrow(DocTermR) - 1), (ncol(DocTermR) - 1))
	 res.ca<-CA(DocTermR[apply(DocTermR,1,sum)>0,],ncp, row.sup, col.sup,
	    quanti.sup, quali.sup,graph)
        res.agg=NULL
        Tagreg=NULL
      }
	 }
      if(!is.null(num.agg))
         Table<-Tagreg[apply(Tagreg,1,sum)>0,]
       else Table<-DocTermR[apply(DocTermR,1,sum)>0,]	
	   
      Inertia<-round(sum(res.ca$eig[,1]),4)
	if(graph){
            barplot(res.ca$eig[,1],main="Eigenvalues",names.arg=paste("dim",1:nrow(res.ca$eig)))
	      VCr<-round(sqrt(sum(res.ca$eig[,1])/min((nrow(Table)-1),(ncol(Table)-1))),4)
 
            if((!is.null(col.sup))&(is.null(num.agg))){
		dev.new()
              plot(res.ca,invisible=c("col","row","row.sup"),selectCol = "cos2 30",unselect=1,axes=axes, title="Supplementary" )
            }
		
          	 dev.new()
      		res.meta<-META.CA(res.ca,naxes=ncp,axe.x=axes[1],axe.y=axes[2],lmd=-Inf,lmk, main="CA documents/words")

		dev.new() 
		plot(res.ca,axes=axes, title="CA documents/words" )

		dev.new()
		res.meta1<-META.CA(res.ca,naxes=ncp,axe.x=axes[1],axe.y=axes[2],lmd=Inf,lmk, main="CA words")

		dev.new() 
		plot(res.ca,invisible=c("col","col.sup"),axes=axes, title="CA documents" )
           		
 	}else{
		res.meta<-META.CA(res.ca,naxes=ncp,axe.x=axes[1],axe.y=axes[2],lmd,lmk, main="CA documents/words",graph = FALSE)
		   VCr<-round(sqrt(sum(res.ca$eig[,1])/min((nrow(Table)-1),(ncol(Table)-1))),4)
            
	}
       InVc<-cbind(Inertia,VCr)
       colnames(InVc)<-c("Inertia", "Cramer V")
       rownames(InVc)<-"Total"
  	 res<-list(TableSummary=TableSummary,DocTermR=DocTermR,Tagreg=Tagreg,Table=Table,Nfreqword=Nfreqword,Ndocword=Ndocword,
                 Glossary=Glossary, res.ca=res.ca, VCr=VCr, res.meta=res.meta,Inertia=Inertia, res.agg=res.agg, Inertia.VCr=InVc)	
  	 class(res)<-c("TxCA","list")

	return(res)
}
