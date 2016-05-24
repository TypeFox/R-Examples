mrMLM<-function(){
  mrenv <- new.env()
  gnewtable<-function (items, multiple = FALSE, chosencol = 1, icon.FUN = NULL, 
                       filter.column = NULL, filter.labels = NULL, filter.FUN = NULL, 
                       handler = NULL, action = NULL, container = NULL, ..., toolkit = guiToolkit()) 
  {
    if (!missing(items)) {
      if (is.vector(items)) 
        items <- data.frame(.= items, stringsAsFactors = FALSE)
      if (is.matrix(items)) 
        items <- data.frame(items, stringsAsFactors = FALSE)
    }
    widget <- .gtable(toolkit, items = items, multiple = multiple, 
                      chosencol = chosencol, icon.FUN = icon.FUN, filter.column = filter.column, 
                      filter.labels = filter.labels, filter.FUN = filter.FUN, 
                      handler = handler, action = action, container = container, 
                      ...)
    obj <- new("gTable", widget = widget, toolkit = toolkit)
    return(obj)
  }
  
  window<-gwindow(title="Multilocus Random-SNP-effect Mixed Linear Model (mrMLM)",visible=TRUE,width=1260,height=730,expand=TRUE)
  importwin<-gwindow("Input Dataset",visible=FALSE,width=250,height=420)
  gimpwin<-ggroup(container=importwin,expand=FALSE)
  
  plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
  gpw<-ggroup(container=plotwin)
  ggpw<-ggraphics(container=gpw)
  plotwin1<-gwindow("Q-Q Plot",visible=FALSE,width=600,height=360)
  gpw1<-ggroup(container=plotwin1)
  ggpw1<-ggraphics(container=gpw1)
  choicekk<-gwindow("Choose Kinship option",visible=FALSE,width=320,height=150)
  gkk<-ggroup(container=choicekk,expand=FALSE)
  includeps<-gwindow("Include population structure?",visible=FALSE,width=320,height=150)
  gps<-ggroup(container=includeps,expand=FALSE)
  parsetwin<-gwindow("Parameter Setting",visible=FALSE,width=260,height=280)
  gpar<-ggroup(container=parsetwin,expand=FALSE)
  choicesave<-gwindow("Save as ...",visible=FALSE,width=250,height=150)
  gcsave<-ggroup(container=choicesave,expand=FALSE)
  
  lyt<-glayout(container=window,spacing=13)
  
  importdata<-gbutton("Input Dataset",container=lyt)
  parset<-gbutton("Parameter Setting",container=lyt)
  manhattan<-gbutton("Manhattan Plot",container=lyt)
  qqplot<-gbutton("QQ Plot",container=lyt)
  
  savefile<-gbutton(" Save ",container=lyt)
  run<-gbutton("Run",container=lyt)
  exit<-gbutton("Exit",container=lyt)
  gwline<-glabel("Critical value for Manhattan Plot",container=lyt)
  gwedit<-gedit("3",width=20,coerce.with=as.numeric,container=lyt)
  svgwline<-svalue(gwedit)
  gwstandp<-glabel("Critical P-value for QQ plot",container=lyt)
  gwedit1<-gedit("0.992",width=20,coerce.with=as.numeric,container=lyt)
  svgwstandp<-svalue(gwedit1)
  
  lyt[1,1]<-importdata
  lyt[4,1]<-parset
  lyt[5,1]<-run
  lyt[6,1]<-savefile
  lyt[9,1]<-gwline
  lyt[10,1]<-gwedit
  lyt[11,1]<-manhattan
  lyt[12,1]<-gwstandp
  lyt[13,1]<-gwedit1
  lyt[14,1]<-qqplot
  lyt[17,1]<-exit
  
  nb1<-gnotebook(tab.pos=3,closebuttons=TRUE,dontCloseThese=TRUE,container=lyt,expand=TRUE)
  size(nb1)<-c(680,540)
  tb<-gnewtable("     
                1. mrMLM is a R software package for genome-wide association studies based on a multi-locus random-SNP-effect mixed linear model.
                
                2. Please cite: Wang Shi-Bo, Feng Jian-Ying, Ren Wen-Long, Huang Bo, Zhou Ling, Wen Yang-Jun, Zhang Jin, Jim M. Dunwell, Xu Shizhong (*), Zhang Yuan-Ming (*). 2016. 
                Improving power and accuracy of genome-wide association studies via a multi-locus mixed linear model methodology.Scientific Reports 6: 19444. 
                
                3. The software package is developed by Wen-Long Ren, Shi-Bo Wang, Bo Huang & Yuan-Ming Zhang.
                
                
                Version 1.2, Realeased April 2016",multiple=TRUE,container=nb1,expand=TRUE,label="About the software")
  
  
  lyt[1:20,2,expand=TRUE]<-nb1
  
  staprogress<-gtkButton()
  lyt[21,2,expand=TRUE]<-staprogress
  
  
  addHandlerClicked(importdata,handler=function(h,...){
    if(isExtant(importwin)==FALSE)
    {
      importwin<-gwindow("Import Dataset",visible=FALSE,width=250,height=420)
      gimpwin<-ggroup(container=importwin,expand=FALSE)
    }
    lytimp<-glayout(container=gimpwin,spacing=13)
    impchoose<-glabel("1. Choose dataset format",container=lytimp)
    impfile1<-glabel("2. Input Genotypic and Phenotypic files",container=lytimp)
    impprepare<-glabel("3. Sort & Transform for dataset",container=lytimp)
    impfile2<-glabel("4. Input Kinship and Population-structure files",container=lytimp)
    radioimp<-gradio(c("mrMLM numeric format","mrMLM character format","Hapmap (TASSEL) format"),selected=3,horizontal=FALSE,container=lytimp)
    genotype<-gbutton("Genotype",container=lytimp)
    phenotype<-gbutton("Phenotype",container=lytimp)
    kinship<-gbutton("Kinship",container=lytimp)
    population<-gbutton("Population Structure",container=lytimp)
    preimp<-gbutton("Do",container=lytimp)
    lytimp[1,2:5]<-impchoose
    lytimp[2:4,2:5]<-radioimp
    lytimp[5,2:5]<-impfile1
    lytimp[6,2:4]<-genotype
    lytimp[7,2:4]<-phenotype
    lytimp[8,2:5]<-impprepare
    lytimp[9,2:4]<-preimp
    lytimp[10,2:5]<-impfile2
    lytimp[11,2:4]<-kinship
    lytimp[12,2:4]<-population
    visible(importwin)<-TRUE
    
    addHandlerClicked(genotype,handler=function(h,...){
      mrenv$flagps<-1
      input1<-gfile(text="Select a file...",type="open",
                    filter=list("All files"=list(patterns=c("*")),
                                "CSV files"=list(patterns=c("*.csv"))))
        
      if(is.na(input1))
      {
        gmessage("Please input correct genotype data !","Warning",icon="warning")
        return
       }else{
         mrenv$genRaw<-as.matrix(read.csv(input1,header=FALSE))
         showgenRaw<-mrenv$genRaw[-1,]
         colnames(showgenRaw)<-mrenv$genRaw[1,]
         showgenRaw<-as.data.frame(showgenRaw)
         tbdfe1<-gdfedit(showgenRaw,container=nb1,expand=TRUE,label="Raw_Genotype") 
      }
    })
    
    addHandlerClicked(phenotype,handler=function(h,...){
      input2<-gfile(text="Select a file...",type="open",
                    filter=list("All files"=list(patterns=c("*")),
                                "CSV files"=list(patterns=c("*.csv"))))
      if(is.na(input2))
      {
        gmessage("Please input correct phenotype data !","Warning",icon="warning")
        return
      }else{
        mrenv$pheRaw<-as.matrix(read.csv(input2,header=FALSE)) 
        showpheRaw<-mrenv$pheRaw[-1,]
        colnames(showpheRaw)<-c(mrenv$pheRaw[1,1],"   ")
        showpheRaw<-as.data.frame(showpheRaw)
        tbdfe2<-gdfedit(showpheRaw,container=nb1,expand=TRUE,label="Raw_Phenotype")
      }
    })
    
    addHandlerClicked(preimp,handler=function(h,...){
      if(svalue(radioimp)=="mrMLM numeric format"){
        mrenv$inputform<-1
        nameGen <- as.matrix(mrenv$genRaw[1,],1,)
        namePhe <- as.matrix(mrenv$pheRaw[,1],,1)
        mrenv$sameName <- intersect(nameGen,namePhe)
        ##########To find the location of the same name 
        locGen <- match(mrenv$sameName,nameGen)
        locPhe <- match(mrenv$sameName,namePhe)
        ##########Produce new genotype matrix and phenotype matrix
        hapName <- matrix(c("rs#","chrom","pos","genotype for code 1"),1,)
        hapHave <- intersect(nameGen,hapName)
        locHap <- match(hapHave,nameGen)
        newGenloc <- c(locHap,locGen)
        newPheloc <- locPhe
        newGen <- as.matrix(mrenv$genRaw[-1,newGenloc])
        newPhe <- as.matrix(mrenv$pheRaw[newPheloc,])
        nnhap <- length(hapHave)
        rownewGen <- dim(newGen)[1]
        colnewGen <- dim(newGen)[2]
        rownewPhe <- dim(newPhe)[1]
        ###########To show on the table ----newGen
        mrenv$newGen <-rbind(mrenv$genRaw[1,newGenloc],newGen)
        ###########To be computed ----gen
        locChr <- as.numeric(which(mrenv$newGen[1,]=="chrom"))
        locPos <- as.numeric(which(mrenv$newGen[1,]=="pos"))
        mrenv$needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
        mrenv$needGen <- mrenv$newGen[,mrenv$needloc]
        mrenv$gen<-as.matrix(mrenv$needGen[-1,])
        mrenv$gen<-matrix(as.numeric(mrenv$gen),nrow=nrow(mrenv$gen))
        ###########To show on the table ----newPhe
        mrenv$pheRaw[1,2]<-"  "
        mrenv$newPhe<-rbind(mrenv$pheRaw[1,],newPhe)
        ###########To be computed ----phe
        mrenv$phe<-as.matrix(mrenv$newPhe[-1,-1])
        mrenv$phe<-matrix(as.numeric(mrenv$phe),nrow=nrow(mrenv$phe))
        shownewGen<-mrenv$newGen[-1,]
        colnames(shownewGen)<-mrenv$newGen[1,]
        shownewGen<-as.data.frame(shownewGen)
        shownewPhe<-mrenv$newPhe[-1,]
        colnames(shownewPhe)<-c(mrenv$newPhe[1,1],"   ")
        shownewPhe<-as.data.frame(shownewPhe)
        tbdfe3<-gdfedit(shownewGen,container=nb1,expand=TRUE,label="Genotype")
        tbdfe4<-gdfedit(shownewPhe,container=nb1,expand=TRUE,label="Phenotype")
      }else if(svalue(radioimp)=="mrMLM character format"){
        mrenv$inputform<-2
        ##########To find the same name between genotype and phenotype
        nameGen <- as.matrix(mrenv$genRaw[1,],1,)
        namePhe <- as.matrix(mrenv$pheRaw[,1],,1)
        mrenv$sameName <- intersect(nameGen,namePhe)
        ##########To find the location of the same name 
        locGen <- match(mrenv$sameName,nameGen)
        locPhe <- match(mrenv$sameName,namePhe)
        ##########Produce new genotype matrix and phenotype matrix
        hapName <- matrix(c("rs#","chrom","pos"),1,)
        hapHave <- intersect(nameGen,hapName)
        locHap <- match(hapHave,nameGen)
        newGenloc <- c(locHap,locGen)
        newPheloc <- locPhe
        newGen <- as.matrix(mrenv$genRaw[-1,newGenloc])
        newPhe <- as.matrix(mrenv$pheRaw[newPheloc,])
        ##########Transfer ATCG to numeric
        nnhap <- length(hapHave)
        rownewGen <- dim(newGen)[1]
        colnewGen <- dim(newGen)[2]
        rownewPhe <- dim(newPhe)[1]
        computeGen <- newGen[,(nnhap+1):colnewGen]
        colComGen <- ncol(computeGen) 
        referSam <- as.vector(computeGen[,1])
        ATCGloc <- c(which(computeGen[,1]=="A"),which(computeGen[,1]=="T"),which(computeGen[,1]=="C"),which(computeGen[,1]=="G"))
        NNRRloc <- setdiff(c(1:rownewGen),ATCGloc)
        for(i in 2:colComGen)
        {
          if(length(NNRRloc)>0){
            referSam[NNRRloc] <- as.vector(computeGen[NNRRloc,i])
            ATCGlocLoop <- c(which(computeGen[NNRRloc,i]=="A"),which(computeGen[NNRRloc,i]=="T"),which(computeGen[NNRRloc,i]=="C"),which(computeGen[NNRRloc,i]=="G"))
            NNRRloc <- setdiff(NNRRloc,NNRRloc[ATCGlocLoop]) 
          }else{
            break
          }
        }
        for(i in 1:rownewGen)
        {
          tempSel1 <- as.vector(c(which(computeGen[i,]=="A"),which(computeGen[i,]=="T"),which(computeGen[i,]=="C"),which(computeGen[i,]=="G")))
          tempSel2 <- as.vector(c(which(computeGen[i,]==referSam[i])))
          notRef <- setdiff(tempSel1,tempSel2)
          notATCG <- setdiff(c(1:colComGen),tempSel1)
          computeGen[i,tempSel2] <- as.numeric(1)
          computeGen[i,notRef] <- as.numeric(-1)
          computeGen[i,notATCG] <- as.numeric(0)
        }
        mrenv$outATCG<-referSam
        ###########To show on the table ----newGen
        newGen <- cbind(newGen[,1:nnhap],computeGen)
        mrenv$newGen <-rbind(mrenv$genRaw[1,newGenloc],newGen)
        ###########To be computed ----gen
        locChr <- as.numeric(which(mrenv$newGen[1,]=="chrom"))
        locPos <- as.numeric(which(mrenv$newGen[1,]=="pos"))
        mrenv$needloc <- c(locChr,locPos,(nnhap+1):colnewGen)
        mrenv$needGen<-mrenv$newGen[,mrenv$needloc]
        mrenv$gen<-as.matrix(mrenv$needGen[-1,])
        mrenv$gen<-matrix(as.numeric(mrenv$gen),nrow=nrow(mrenv$gen))
        ###########To show on the table ----newPhe
        mrenv$pheRaw[1,2]<-"  "
        mrenv$newPhe<-rbind(mrenv$pheRaw[1,],newPhe)
        ###########To be computed ----phe
        mrenv$phe<-as.matrix(mrenv$newPhe[-1,-1])
        mrenv$phe<-matrix(as.numeric(mrenv$phe),nrow=nrow(mrenv$phe))
        shownewGen<-mrenv$newGen[-1,]
        colnames(shownewGen)<-mrenv$newGen[1,]
        shownewGen<-as.data.frame(shownewGen)
        shownewPhe<-mrenv$newPhe[-1,]
        colnames(shownewPhe)<-c(mrenv$newPhe[1,1],"   ")
        shownewPhe<-as.data.frame(shownewPhe)
        tbdfe3<-gdfedit(shownewGen,container=nb1,expand=TRUE,label="Genotype")
        tbdfe4<-gdfedit(shownewPhe,container=nb1,expand=TRUE,label="Phenotype")
      }else if(svalue(radioimp)=="Hapmap (TASSEL) format"){
        mrenv$inputform<-3
        ##########To find the same name between genotype and phenotype
        nameGen<-as.matrix(mrenv$genRaw[1,],1,)
        namePhe<-as.matrix(mrenv$pheRaw[,1],,1)
        mrenv$sameName<-intersect(nameGen,namePhe)
        ##########To find the location of the same name 
        locGen<-match(mrenv$sameName,nameGen)
        locPhe<-match(mrenv$sameName,namePhe)
        ##########Produce new genotype matrix and phenotype matrix
        hapName<-matrix(c("rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panel","QCcode"),1,)
        hapHave<-intersect(nameGen,hapName)
        locHap<-match(hapHave,nameGen)
        newGenloc<-c(locHap,locGen)
        newPheloc<-locPhe
        newGen<-as.matrix(mrenv$genRaw[-1,newGenloc])
        newPhe<-as.matrix(mrenv$pheRaw[newPheloc,])   
        ##########Transfer ATCG to numeric
        nnhap<-length(hapHave)
        rownewGen<-dim(newGen)[1]
        colnewGen<-dim(newGen)[2]
        rownewPhe<-dim(newPhe)[1]
        computeGen<-newGen[,(nnhap+1):colnewGen]
        colComGen<-ncol(computeGen) 
        referSam<-as.vector(computeGen[,1])
        ATCGloc<-c(which(computeGen[,1]=="AA"),which(computeGen[,1]=="TT"),which(computeGen[,1]=="CC"),which(computeGen[,1]=="GG"))
        NNRRloc<-setdiff(c(1:rownewGen),ATCGloc)
        for(i in 2:colComGen)
        {
          if(length(NNRRloc)>0){
            referSam[NNRRloc]<-as.vector(computeGen[NNRRloc,i])
            ATCGlocLoop<-c(which(computeGen[NNRRloc,i]=="AA"),which(computeGen[NNRRloc,i]=="TT"),which(computeGen[NNRRloc,i]=="CC"),which(computeGen[NNRRloc,i]=="GG"))
            NNRRloc<-setdiff(NNRRloc,NNRRloc[ATCGlocLoop]) 
          }else{
            break
          }
        }
        for(i in 1:rownewGen)
        {
          tempSel1<-as.vector(c(which(computeGen[i,]=="AA"),which(computeGen[i,]=="TT"),which(computeGen[i,]=="CC"),which(computeGen[i,]=="GG")))
          tempSel2<-as.vector(c(which(computeGen[i,]==referSam[i])))
          notRef<-setdiff(tempSel1,tempSel2)
          notATCG<-setdiff(c(1:colComGen),tempSel1)
          computeGen[i,tempSel2]<-as.numeric(1)
          computeGen[i,notRef]<-as.numeric(-1)
          computeGen[i,notATCG]<-as.numeric(0)
        }
        mrenv$outATCG<-referSam
        ###########To show on the table ----mrenv$newGen
        newGen<-cbind(newGen[,1:nnhap],computeGen)
        mrenv$newGen<-rbind(mrenv$genRaw[1,newGenloc],newGen)
        ###########To be computed ----mrenv$gen
        locChr<-as.numeric(which(mrenv$newGen[1,]=="chrom"))
        locPos<-as.numeric(which(mrenv$newGen[1,]=="pos"))
        mrenv$needloc<-c(locChr,locPos,(nnhap+1):colnewGen)
        mrenv$needGen<-mrenv$newGen[,mrenv$needloc]
        mrenv$gen<-as.matrix(mrenv$needGen[-1,])
        mrenv$gen<-matrix(as.numeric(mrenv$gen),nrow=nrow(mrenv$gen))
        ###########To show on the table ----mrenv$newPhe
        mrenv$pheRaw[1,2]<-"  "
        mrenv$newPhe<-rbind(mrenv$pheRaw[1,],newPhe)
        ###########To be computed ----mrenv$phe
        mrenv$phe<-as.matrix(mrenv$newPhe[-1,-1])
        mrenv$phe<-matrix(as.numeric(mrenv$phe),nrow=nrow(mrenv$phe))
        shownewGen<-mrenv$newGen[-1,]
        colnames(shownewGen)<-mrenv$newGen[1,]
        shownewGen<-as.data.frame(shownewGen)
        shownewPhe<-mrenv$newPhe[-1,]
        colnames(shownewPhe)<-c(mrenv$newPhe[1,1],"   ")
        shownewPhe<-as.data.frame(shownewPhe)
        tbdfe3<-gdfedit(shownewGen,container=nb1,expand=TRUE,label="Genotype")
        tbdfe4<-gdfedit(shownewPhe,container=nb1,expand=TRUE,label="Phenotype")
      }
    })
    
    addHandlerClicked(kinship,handler=function(h,...){
      if(isExtant(choicekk)==FALSE)
      {
        choicekk<-gwindow("Choose Kinship option",visible=FALSE,width=320,height=150)
        gkk<-ggroup(container=choicekk,expand=FALSE)
      }
      lytkk<-glayout(container=gkk,spacing=13)
      mrenv$okkk<-gbutton("     OK    ",container=lytkk)
      mrenv$cancelkk<-gbutton(" Cancel ",container=lytkk)
      mrenv$radiokk<-gradio(c("Input the Kinship matrix file","Calculate the Kinship matrix by this software"),selected=1,horizontal=FALSE,container=lytkk)
      lytkk[2:3,2:5]<-mrenv$radiokk
      lytkk[5,2]<-mrenv$okkk
      lytkk[5,5]<-mrenv$cancelkk
      visible(choicekk)<-TRUE
      addHandlerClicked(mrenv$okkk,handler=function(h,...){
      if(svalue(mrenv$radiokk)=="Input the Kinship matrix file"){
          input3<-gfile(text="Select a file...",type="open",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          if(is.na(input3))
          {
            gmessage("Please input correct kinship data !","Warning",icon="warning")
            return
          }else{
            mrenv$kkRaw<-read.csv(input3,header=FALSE)
            nnkk<-dim(mrenv$kkRaw)[1]
            mrenv$kkRaw[1,2:nnkk]<-"  "
            tbdfe5<-gdfedit(mrenv$kkRaw,container=nb1,expand=TRUE,label="Kinship")
            kkPre<-as.matrix(mrenv$kkRaw[-1,-1])
            nameKin<-as.matrix(mrenv$kkRaw[-1,1])
            sameGenKin<-intersect(mrenv$sameName,nameKin)
            locKin<-match(sameGenKin,nameKin)
            mrenv$kk<-kkPre[locKin,locKin]
            mrenv$kk<-matrix(as.numeric(mrenv$kk),nrow=nrow(mrenv$kk))
            dispose(choicekk)
          }
        }else{
          envgen <- mrenv$gen
          if(exists("envgen")==FALSE)
          {
            gmessage("Please input correct genotype data !","Warning",icon="warning")
            return
          }else{
            envgen<-envgen[,3:(ncol(envgen))]
            envgen<-t(envgen)
            m<-ncol(envgen)
            n<-nrow(envgen)
            kk1<-matrix(0,n,n)
            for(k in 1:m){
              z<-as.matrix(envgen[,k])
              kk1<-kk1+z%*%t(z)
            }
            cc<-mean(diag(kk1))
            kk1<-kk1/cc
            mrenv$kk<-as.matrix(kk1)
            rowsize<-dim(mrenv$kk)[1]
            aa<-as.character()
            for(i in 1:(rowsize+1))
            {
              a<-paste("V",i,sep="")
              aa<-c(aa,a)
            }
            mrenv$kkShow<-cbind(matrix(mrenv$newPhe[-1,1],,1),round(mrenv$kk,5))
            tempFirst<-rep("  ",rowsize)
            tempFirst<-c(as.character(rowsize),tempFirst)
            mrenv$kkShow<-as.matrix(rbind(tempFirst,mrenv$kkShow))
            colnames(mrenv$kkShow)<-aa
            rownames(mrenv$kkShow)<-c(1:(rowsize+1))
            tbdfe5<-gdfedit(mrenv$kkShow,container=nb1,expand=TRUE,label="Kinship")
            dispose(choicekk)
          }     
        }
      })
      addHandlerClicked(mrenv$cancelkk,handler=function(h,...){
        dispose(choicekk)
      })
    })
    
      
      addHandlerClicked(population,handler=function(h,...){
        if(isExtant(includeps)==FALSE)
        {
          includeps<-gwindow("Include population structure?",visible=FALSE,width=320,height=150)
          gps<-ggroup(container=includeps,expand=FALSE)
        }
        lytps<-glayout(container=gps,spacing=13)
        okps<-gbutton("     OK    ",container=lytps)
        cancelps<-gbutton(" Cancel ",container=lytps)
        radiops<-gradio(c("Population structure has no effect on GWAS","Input Population Structure file"),selected=1,horizontal=FALSE,container=lytps)
        lytps[2:3,2:5]<-radiops
        lytps[5,2]<-okps
        lytps[5,5]<-cancelps
        visible(includeps)<-TRUE
        addHandlerClicked(okps,handler=function(h,...){
          if(svalue(radiops)=="Input Population Structure file"){
            mrenv$flagps<-0
            input4<-gfile(text="Select a file...",type="open",
                          filter=list("All files"=list(patterns=c("*")),
                                      "CSV files"=list(patterns=c("*.csv"))))
            if(is.na(input4))
            {
              gmessage("Please input correct population data !","Warning",icon="warning")
              return
            }else{
              mrenv$psmatrixRaw<-as.matrix(read.csv(input4,header=FALSE))
              nnpprow<-dim(mrenv$psmatrixRaw)[1]
              nnppcol<-dim(mrenv$psmatrixRaw)[2]
              mrenv$psmatrixRaw[1,2:nnppcol]<-"  "
              psmatrixPre<-mrenv$psmatrixRaw[3:nnpprow,]
              namePop<-as.matrix(psmatrixPre[,1])
              sameGenPop<-intersect(mrenv$sameName,namePop)
              locPop<-match(sameGenPop,namePop)
              mrenv$psmatrix<-psmatrixPre[locPop,-1:-2]
              mrenv$psmatrix<-matrix(as.numeric(mrenv$psmatrix),nrow=nrow(mrenv$psmatrix))
              tbdfe6<-gdfedit(mrenv$psmatrixRaw,container=nb1,expand=TRUE,label="Population Structure")
              dispose(includeps)
              dispose(importwin)
            }
          }else{
            mrenv$flagps<-1
            enabled(population)<-FALSE
            dispose(includeps)
            dispose(importwin)
          }
        })
        addHandlerClicked(cancelps,handler=function(h,...){
          dispose(includeps)
        })
      })
    
  })
  
  addHandlerClicked(parset,handler=function(h,...){
    if(isExtant(parsetwin)==FALSE)
    {
      parsetwin<-gwindow("Parameter Setting",visible=FALSE,width=260,height=280)
      gpar<-ggroup(container=parsetwin,expand=FALSE)
    }
    lytpar<-glayout(container=gpar,spacing=13)
    mrenv$pvallabel<-glabel("1. Critical P-value in rMLM:",container=lytpar)
    mrenv$pvaledit<-gedit("0.01",width=20,coerce.with=as.numeric,container=lytpar)
    mrenv$radlabel<-glabel("2. Search radius of candidate gene (kb):",container=lytpar)
    mrenv$radedit<-gedit("20",width=20,coerce.with=as.numeric,container=lytpar)
    mrenv$mlodlabel<-glabel("3. Critical LOD score in mrMLM:",container=lytpar)
    mrenv$mlodedit<-gedit("3",width=20,coerce.with=as.numeric,container=lytpar)
    mrenv$okpar<-gbutton("     OK    ",container=lytpar)
    mrenv$cancelpar<-gbutton(" Cancel ",container=lytpar)
    lytpar[1,1:5]<-mrenv$pvallabel
    lytpar[2,1:5]<-mrenv$pvaledit
    lytpar[3,1:5]<-mrenv$radlabel
    lytpar[4,1:5]<-mrenv$radedit
    lytpar[5,1:5]<-mrenv$mlodlabel
    lytpar[6,1:5]<-mrenv$mlodedit
    lytpar[7,1]<-mrenv$okpar
    lytpar[7,4]<-mrenv$cancelpar
    visible(parsetwin)<-TRUE
    addHandlerClicked(mrenv$okpar,handler=function(h,...){
      mrenv$svpal<-svalue(mrenv$pvaledit)
      mrenv$svrad<-svalue(mrenv$radedit)
      mrenv$svmlod<-svalue(mrenv$mlodedit)
      if((mrenv$svpal<0)||(mrenv$svpal>1))
      {
        gmessage("Please input critical P-value more than 0 and less than 1!","Warning",icon="warning")
        return
      }
      if(mrenv$svrad<0)
      {
        gmessage("Please input search radius of candidate gene more than 0!","Warning",icon="warning")
        return
      }
      if(mrenv$svmlod<0)
      {
        gmessage("Please input critical LOD score more than 0!","Warning",icon="warning")
        return
      }
      if((mrenv$svpal>0)&&(mrenv$svpal<1)&&(mrenv$svrad>=0)&&(mrenv$svmlod>=0))
      {
        dispose(parsetwin)
      }
    })
    
    addHandlerClicked(mrenv$cancelpar,handler=function(h,...){
      mrenv$svpal<-svalue(mrenv$pvaledit)
      mrenv$svrad<-svalue(mrenv$radedit)
      mrenv$svmlod<-svalue(mrenv$mlodedit)
      dispose(parsetwin)
    })
  })
    
  addHandlerClicked(exit,handler=function(h,...){
    gconfirm("Yes or no?",handler=function(h,...){dispose(window)})
  })
  
  addHandlerClicked(run,handler=function(h,...){
    gen<-mrenv$gen
    phe<-mrenv$phe
    kk<-mrenv$kk
    flagps<-mrenv$flagps
    psmatrix<-mrenv$psmatrix
    if(exists("gen")==FALSE)
    {
      gmessage("Please input correct genotype data !","Warning",icon="warning")
      return
    }
    if(exists("phe")==FALSE)
    {
      gmessage("Please input correct phenotype data !","Warning",icon="warning")
      return
    }
    if(exists("kk")==FALSE)
    {
      gmessage("Please input correct kinship data !","Warning",icon="warning")
      return
    }
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
    {
      gmessage("Sample size between genotype and phenotype is inconsistent!","Error",icon="error")
      return
    }
    
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2))))
    {
      progress_bar <- gtkProgressBar ( )
      staprogress$add(progress_bar)
      progress_bar$setText ( "Please be patient ..." )
      progress_bar$setFraction(2/100)
      
      mixed<-function(x,y,kk){
        
        loglike<-function(theta){
          lambda<-exp(theta)
          logdt<-sum(log(lambda*delta+1))
          h<-1/(lambda*delta+1)
          yy<-sum(yu*h*yu)
          yx<-matrix(0,q,1)
          xx<-matrix(0,q,q)
          for(i in 1:q){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:q){
              xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
          }
          loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
          return(-loglike)
        }
        
        fixed<-function(lambda){
          h<-1/(lambda*delta+1)
          yy<-sum(yu*h*yu)
          yx<-matrix(0,q,1)
          xx<-matrix(0,q,q)
          for(i in 1:q){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:q){
              xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
          } 
          beta<-solve(xx,yx)
          sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
          sigma2<-drop(sigma2)
          var<-diag(solve(xx)*sigma2)
          stderr<-sqrt(var)
          return(c(beta,stderr,sigma2))
        }
        
        qq<-eigen(kk)
        delta<-qq[[1]]
        uu<-qq[[2]]
        q<-ncol(x)
        n<-ncol(kk)
        vp<-var(y)
        yu<-t(uu)%*%y
        xu<-t(uu)%*%x
        theta<-0
        parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
        lambda<-exp(parm$par)
        conv<-parm$convergence
        fn1<-parm$value
        fn0<-loglike(-Inf)
        lrt<-2*(fn0-fn1)
        hess<-parm$hessian
        parmfix<-fixed(lambda)
        beta<-parmfix[1:q]
        stderr<-parmfix[(q+1):(2*q)]
        sigma2<-parmfix[2*q+1]
        lod<-lrt/4.61
        p_value<-1-pchisq(lrt,1)
        sigma2g<-lambda*sigma2
        goodness<-(vp-sigma2)/vp
        par<-data.frame(lrt,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
        return(par)
      }
      
      
      loglike<-function(theta){
        xi<-exp(theta)
        tmp0<-zz*xi+1
        tmp<-xi*solve(tmp0)
        yHy<-yy-t(zy)%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        logdt2<-log(det(tmp0))
        loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
        return(-loglike)
      }
      
      fixed<-function(xi){
        tmp0<-zz*xi+diag(1)
        tmp<-xi*solve(tmp0)
        yHy<-yy-t(zy)%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        zHy<-zy-zz%*%tmp%*%zy
        zHx<-zx-zx%*%tmp%*%zz
        zHz<-zz-zz%*%tmp%*%zz
        beta<-solve(xHx,yHx)
        tmp2<-solve(xHx)
        sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
        gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
        var<-abs((xi*diag(1)-xi*zHz*xi)*as.numeric(sigma2))
        stderr<-sqrt(diag(var))
        result<-list(gamma,stderr,beta,sigma2)
        return(result)
      }
      
      name<-gen[,1:2]
      gen<-gen[,3:(ncol(gen))]
      gen<-t(gen)
      n<-nrow(gen)
      m<-ncol(gen)
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        x<-matrix(1,n,1)
      }else if(flagps==0)
      {
        x<-cbind(matrix(1,n,1),psmatrix)
      }
      ll<-numeric()
      s<-ncol(x)
      kk<-as.matrix(kk)
      qq<-eigen(kk)
      delta<-qq[[1]]
      uu<-qq[[2]]
      xu<-t(uu)%*%x
      for(ii in 1:1)
      {yy<-phe[,1]
       y<-as.matrix(yy)
       parm<-mixed(x=x,y=y,kk=kk)
       lambda<-parm$lambda[1]
       h<-1/(delta*lambda+1)
       yu<-t(uu)%*%y
       xx<-matrix(0,s,s)
       for(i in 1:s){
         for(j in 1:s){
           xx[i,j]<-sum(xu[,i]*h*xu[,j])
         }
       }
       yy<-sum(yu*h*yu)
       yx<-matrix(0,s,1)
       for(i in 1:s){
         yx[i]<-sum(yu*h*xu[,i])
       }
       
       qq<-numeric()
       for(k in 1:m){
         progress_bar$setFraction((2+(90/m)*k)/100)
         z<-as.matrix(gen[,k])
         zu<-t(uu)%*%z
         zy<-as.matrix(sum(yu*h*zu))
         zz<-as.matrix(sum(zu*h*zu))
         zx<-matrix(0,s,1)
         for(i in 1:s){
           zx[i]<-sum(xu[,i]*h*zu)
         }
         theta<-c(0)
         par<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
         xi<-exp(par$par)
         conv<-par$convergence
         fn1<-par$value
         hess<-par$hessian
         parmfix<-fixed(xi)
         gamma<-parmfix[[1]]
         stderr<-parmfix[[2]]
         if((flagps==1)||(exists("psmatrix")==FALSE))
         {
           beta<-parmfix[[3]]
         }else if(flagps==0)
         {
           beta<-parmfix[[3]][1]
         }
         sigma2<-parmfix[[4]]
         lambda<-xi
         sigma2g<-lambda*sigma2
         fn0<-loglike(-Inf)
         lrt<-2*(fn0-fn1)
         p_lrt<-1-pchisq(lrt,1)
         wald<-(gamma/stderr)^2
         p_wald<-1-pchisq(wald,1)
         parm0<-c(ii,name[k,1],name[k,2],beta,sigma2,sigma2g,gamma,stderr,wald,p_wald)
         qq<-rbind(qq,parm0)
       }
       ll<-rbind(ll,qq)
      }
      mrenv$parms<-ll
      mrenv$parms<-matrix(mrenv$parms,,10)
      
      multinormal<-function(y,mean,sigma)
      {
        pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
        return (pdf_value)
      }
      #LOD value test
      likelihood<-function(xxn,xxx,yn)
      {
        nq<-ncol(xxx)
        ns<-nrow(yn)
        at1<-0
        ww1<-1:ncol(xxx)
        ww1<-as.matrix(ww1)
        at1<-dim(ww1)[1]
        lod<-matrix(rep(0,nq),nq,1)
        if(at1>0.5)
          ad<-cbind(xxn,xxx[,ww1])
        else
          ad<-xxn
        #if(abs(det(crossprod(ad,ad)))<1e-6)
        if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
          bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
        else
          bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
        vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
        ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
        
        sub<-1:ncol(ad);
        if(at1>0.5)
        {
          for(i in 1:at1)
          {
            #ij<-which(sub!=sub[i+1])
            ij<-which(sub!=sub[i+ncol(xxn)])
            ad1<-ad[,ij]
            #if(abs(det(crossprod(ad1,ad1)))<1e-6)
            if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
              bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
            else
              bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
            vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
            ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
            lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
          }
        }
        return (lod)
      }
      
      #2010 EM_Bayes
      ebayes_EM<-function(x,z,y)
      {
        n<-nrow(z);k<-ncol(z)
        if(abs(min(eigen(crossprod(x,x))$values))<1e-6)
          b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
        else
          b<-solve(crossprod(x,x))%*%crossprod(x,y)
        v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
        u<-matrix(rep(0,k),k,1)
        v<-matrix(rep(0,k),k,1)
        s<-matrix(rep(0,k),k,1)
        for(i in 1:k)
        {
          zz<-z[,i]
          s[i]<-((crossprod(zz,zz))^(-1))*v0
          u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
          v[i]<-u[i]^2+s[i]
        }
        vv<-matrix(rep(0,n*n),n,n);
        for(i in 1:k)
        {
          zz<-z[,i]
          vv=vv+tcrossprod(zz,zz)*v[i]
        }
        vv<-vv+diag(n)*v0
        iter<-0;err<-1000;iter_max<-100;err_max<-1e-8
        tau<-0;omega<-0
        while((iter<iter_max)&&(err>err_max))
        {
          iter<-iter+1
          v01<-v0
          v1<-v
          b1<-b
          vi<-solve(vv)
          xtv<-crossprod(x,vi)
          if(ncol(x)==1)
          {
            b<-((xtv%*%x)^(-1))*(xtv%*%y)
          }else
          {
            if(abs(min(eigen(xtv%*%x)$values))<1e-6){
              b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
            }
            else{
              b<-solve(xtv%*%x)%*%(xtv%*%y)
            }
          }
          r<-y-x%*%b
          ss<-matrix(rep(0,n),n,1)
          for(i in 1:k)
          {
            zz<-z[,i]
            zztvi<-crossprod(zz,vi)
            u[i]<-v[i]*zztvi%*%r
            s[i]<-v[i]*(1-zztvi%*%zz*v[i])
            v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
            ss<-ss+zz*u[i]
          }
          v0<-as.numeric(crossprod(r,(r-ss))/n)
          vv<-matrix(rep(0,n*n),n,n)
          for(i in 1:k)
          {
            zz<-z[,i]
            vv<-vv+tcrossprod(zz,zz)*v[i]
          }
          vv<-vv+diag(n)*v0
          err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
          beta<-t(b)
          sigma2<-v0
        }
        wang<-matrix(rep(0,k),k,1)
        for (i in 1:k){
          stderr<-sqrt(s[i]+1e-20)
          t<-abs(u[i])/stderr
          f<-t*t
          p<-1-pchisq(f,1)
          wang[i]<-p
        }
        return (wang)
      }
           
      gen<-t(gen) 
      chr_pos<-mrenv$parms[,2:3]
      pfit<-which(mrenv$parms[,10]<=(mrenv$svpal))
      pfit<-as.matrix(pfit)
      pfitrow<-nrow(pfit)
      no_p<-cbind((1:(nrow(mrenv$parms))),mrenv$parms[,10])
      no_porder<-order(no_p[,2])
      no_p<-no_p[no_porder,]
      choose_orderp<-no_p[1:pfitrow,]
      orderno<-no_p[1:pfitrow,1]
      orderno<-as.matrix(orderno)
      sigma2g_SNPerr<-cbind(mrenv$parms[,6],mrenv$parms[,8])
      correct_each<-matrix(1,(nrow(sigma2g_SNPerr)),1)-sigma2g_SNPerr[,2]*sigma2g_SNPerr[,2]/sigma2g_SNPerr[,1]
      k0<-which(correct_each<0)
      k0<-as.matrix(k0)
      if(nrow(k0)>0){
        correct_each[k0,1]<-matrix(0,(nrow(k0)),1)
      }
      correct_sum<-sum(correct_each)
      newp<-0.05/correct_sum
      mrenv$mannewp<-newp
      mrenv$manstandchoice<-1
      no_porder<-which(no_p[,2]<=newp)
      no_porder<-as.matrix(no_porder)
      no_porderrow<-nrow(no_porder)
      gg<-orderno
      for (ii in 1:(nrow(orderno)-1)){
        for (jj in (ii+1):(nrow(orderno))){
          ci<- chr_pos[orderno[ii],1]
          cj<- chr_pos[orderno[jj],1]
          if (ci==cj){
            ye<-abs(chr_pos[orderno[ii],2]-chr_pos[orderno[jj],2])
            if (ye<=((mrenv$svrad)*1000)){
              gg[jj,1]<-0
            }
          }
        }
      }
      progress_bar$setFraction(95/100)
      if(mrenv$inputform==1){
        #output result1 using mrMLM numeric format
        mrenv$parmsShow<-mrenv$parms[,-1]
        mrenv$parmsShow<-cbind(mrenv$genRaw[-1,1],mrenv$parms[,2:3],round(mrenv$parms[,4:10],5))
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,10)
        mrenv$parmsShow<-cbind(mrenv$parmsShow,mrenv$genRaw[-1,4])
        meadd<-matrix(1,nrow(mrenv$parms),1)
        meadd[which(mrenv$parms[,10]<newp),1]<-round(newp,8)
        meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
        mrenv$parmsShow<-cbind(mrenv$parmsShow,meadd)
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,12)
        colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype for code 1","Significance")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      if(mrenv$inputform==2){
        #output result1 using mrMLM character format
        mrenv$parmsShow<-mrenv$parms[,-1]
        mrenv$parmsShow<-cbind(mrenv$genRaw[-1,1],mrenv$parms[,2:3],round(mrenv$parms[,4:10],5))
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,10)
        mrenv$outATCG<-matrix(mrenv$outATCG,,1)
        mrenv$parmsShow<-cbind(mrenv$parmsShow,mrenv$outATCG)
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,11)
        meadd<-matrix(1,nrow(mrenv$parms),1)
        meadd[which(mrenv$parms[,10]<newp),1]<-round(newp,8)
        meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
        mrenv$parmsShow<-cbind(mrenv$parmsShow,meadd)
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,12)
        colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype  for code 1","Significance")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      if(mrenv$inputform==3){
        #output result1 using TASSEL format
        mrenv$parmsShow<-mrenv$parms[,-1]
        mrenv$parmsShow<-cbind(mrenv$genRaw[-1,1],mrenv$parms[,2:3],round(mrenv$parms[,4:10],5))
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,10)
        mrenv$outATCG<-matrix(mrenv$outATCG,,1)
        mrenv$outATCG<-unlist(strsplit(mrenv$outATCG,""))
        mrenv$outATCG<-matrix(mrenv$outATCG[c(TRUE,FALSE)],,1)
        mrenv$parmsShow<-cbind(mrenv$parmsShow,mrenv$outATCG)
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,11)
        meadd<-matrix(1,nrow(mrenv$parms),1)
        meadd[which(mrenv$parms[,10]<newp),1]<-round(newp,8)
        meadd[which(mrenv$parms[,10]>=newp),1]<-"  "
        mrenv$parmsShow<-cbind(mrenv$parmsShow,meadd)
        mrenv$parmsShow<-matrix(mrenv$parmsShow,,12)
        colnames(mrenv$parmsShow)<-c("RS#","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald","Genotype  for code 1","Significance")
        tbdfe7<-gdfedit(mrenv$parmsShow,container=nb1,expand=TRUE,label="Result1")
      }
      
      gg<-as.matrix(gg)
      misfit<-numeric()
      kk<- numeric()
      kk0<- numeric()
      l0<- numeric()
      bong<-no_porderrow
      if (bong>0){
        g0<-gg[1:no_porderrow,1]
        g0<-as.matrix(g0)
        kk0<-no_porderrow
        no_porderrow<-which(g0>0)
        no_porderrow<-as.matrix(no_porderrow)
        g0<-g0[no_porderrow,1]
        g0<-as.matrix(g0)
        xxx0<-gen[g0,]
        if(dim(g0)[1]==1){
          xxx0<-as.matrix(xxx0)
        }
        if(dim(g0)[1]>1)
        {
          xxx0<-as.matrix(xxx0)
          xxx0<-t(xxx0)
        }
        phe<-as.matrix(phe)
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          par<-likelihood(matrix(1,(nrow(xxx0)),1),xxx0,phe)
          lod<-par
        }else if(flagps==0)
        {
          temp<-cbind(matrix(1,(nrow(xxx0)),1),psmatrix)
          par<-likelihood(temp,xxx0,phe)
          lod<-par
        }
        kk<-which(lod>=1.5)
        kk<-as.matrix(kk)
        kk1<-which(lod<1.5)
        kk1<-as.matrix(kk1)
        if ((nrow(kk1))>0){
          misfit<-g0[kk1,1]
          misfit<-as.matrix(misfit)
        }
        if ((nrow(kk))>0){
          g0<-as.matrix(g0)
          g0<-g0[kk,1]
          xx0<-xxx0[,kk]
          lo<-lod[kk,1]
        }
        if ((nrow(kk))==0){kk<-0}
      }
      if (bong==0){
        kk0<-0
        kk<-0
      }
      nleft<-as.matrix(gg[(kk0+1):(nrow(gg)),1])
      if ((length(misfit))>0){gg<-rbind(nleft,misfit)}
      if ((length(misfit))==0){gg<-nleft}
      a1<-which(gg>0)
      a1<-as.matrix(a1)
      a2<-gg[a1,1]
      a2<-as.matrix(a2)
      xx<-t(gen[a2,])
      xx<-as.matrix(xx)
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        if (length(kk)>1){xin<-cbind(matrix(1,(nrow(xx)),1),xx0)}
        if (length(kk)==1){
          if(kk==0){
            xin<- matrix(1,(nrow(xx)),1)
          }
          if(kk>0){
            xin<-cbind(matrix(1,(nrow(xx)),1),xx0)
          }
        }
      }else if(flagps==0)
      {
        temp<-cbind(matrix(1,(nrow(xx)),1),psmatrix)
        if (length(kk)>1){xin<-cbind(temp,xx0)}
        if (length(kk)==1){
          if(kk==0){
            xin<-temp
          }
          if(kk>0){
            xin<-cbind(temp,xx0)
          }
        }
      }
      xin<-as.matrix(xin)
      par<-ebayes_EM(xin,xx,phe)
      w2<-which(par[,1]<=0.01)
      w2<-as.matrix(w2)
      ww<- numeric()
      if ((nrow(w2))>0){
        orderno<-a2[w2,1]
        orderno<-as.matrix(orderno)
        x3<-cbind(xin,xx[,w2])
        x3<-as.matrix(x3)
        lodfix<-matrix(x3[,1],nrow(x3),)
        lodrand<-matrix(x3[,2:(ncol(x3))],nrow(x3),)
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          lod<-likelihood(lodfix,lodrand,phe)
        }else if(flagps==0)
        {
          temp<-cbind(psmatrix,lodfix)
          lod<-likelihood(temp,lodrand,phe)
        }
        w3<-which(lod[,1]>=(mrenv$svmlod))
        w3<-as.matrix(w3)
        if ((kk[1])>0){
          g0<-as.matrix(g0)
          orderno<-rbind(g0,orderno)
          orderno<-as.matrix(orderno)
        }
        if ((w3[1])>0){
          if((flagps==1)||(exists("psmatrix")==FALSE))
          {
            lo<-lod[w3,1]
            ww<-orderno[w3,]
          }else if(flagps==0)
          {
            lo<-lod[w3,1]
            no_loc<-w3-ncol(psmatrix)
            ww<-orderno[no_loc,]
          }
        }
        if ((nrow(w3))==0){ww<-0}
        #if (length(kk)==1){ww<-0}
      }
      if ((nrow(w2))==0){
        g0<-as.matrix(g0)
        lo<-as.matrix(lo)
        yang<-which(lo>=(mrenv$svmlod))
        yang<-as.matrix(yang)
        if ((nrow(yang))>0){
          ww<-g0[yang,1]
          lo<-lo[yang,1]
        }
        if ((nrow(yang))==0){ww<-0}
      }
      ww<-as.matrix(ww)
      mrenv$needww<-ww
      if (length(ww)>1){
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          ex<-cbind(matrix(1,(nrow(xx)),1),t(gen[ww,]))
        }else if(flagps==0)
        {
          ex<-cbind(cbind(matrix(1,(nrow(xx)),1),psmatrix),t(gen[ww,]))
        }
        ex<-as.matrix(ex)
        cui<-det(t(ex)%*%ex)
        p1<-rep(1,ncol(ex))
        p2<-diag(p1)
        if (cui<1e-6){bbbb<-solve(t(ex)%*%ex+p2*0.01)%*%t(ex)%*%phe}
        if (cui>=1e-6){ bbbb<-solve(t(ex)%*%ex)%*%t(ex)%*%phe }
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          eeff<-bbbb[2:(nrow(bbbb)),1]
        }else if(flagps==0)
        {
          eeff<-bbbb[(2+ncol(psmatrix)):(nrow(bbbb)),1]
        }
        
        eeff<-as.matrix(eeff)
        er<-as.numeric()
        her<-as.numeric()
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          excol<-ncol(ex)
          for(i in 1:(excol-1))
          {
            em<-ex[,(1+i)]
            as1<-length(which(em==1))/nrow(ex)
            as2<-1-as1
            er<-rbind(er,(1-(as1-as2)*(as1-as2))*eeff[i]*eeff[i])
          }
          v0<-(1/(nrow(ex)-1))*(t(phe-ex%*%bbbb)%*%(phe-ex%*%bbbb))
          her<-er/as.numeric(sum(er)+v0)
        }else if(flagps==0)
        {
          excol<-ncol(ex)
          for(i in 1:(excol-1-ncol(psmatrix)))
          {
            em<-ex[,(1+ncol(psmatrix)+i)]
            as1<-length(which(em==1))/nrow(ex)
            as2<-1-as1
            er<-rbind(er,(1-(as1-as2)*(as1-as2))*eeff[i]*eeff[i])
          }
          v0<-(1/(nrow(ex)-1))*(t(phe-ex%*%bbbb)%*%(phe-ex%*%bbbb))
          her<-er/as.numeric(sum(er)+v0)
        }
        
        mrenv$wan<-cbind(chr_pos[ww,],round(eeff,5),round(lo,5),round(her,5))
        mrenv$wan<-matrix(mrenv$wan,dim(mrenv$wan)[1],)

        mrenv$wan<-cbind(matrix(mrenv$parmsShow[mrenv$needww,1],,1),mrenv$wan,matrix(mrenv$parmsShow[mrenv$needww,11],,1))
        mrenv$wan<-matrix(mrenv$wan,dim(mrenv$wan)[1],)
        colnames(mrenv$wan)<-c("RS#","Chromosome","Position","QTN effect","LOD score","R2","Genotype  for code 1")
        
      }
      wan<-mrenv$wan
      if(exists("wan")==FALSE||is.null(wan)==TRUE)
      {
        gmessage("There is no result meets the requirements in the second step!","Info",icon="info")
      }else{
        tbdfe8<-gdfedit(wan,container=nb1,expand=TRUE,label="Result2")
      }
      progress_bar$setFraction(100/100)
      progress_bar$setText("All done.")
    }
    return
  })
  
  addHandlerClicked(savefile,handler=function(h,...){
    if(isExtant(choicesave)==FALSE)
    {
      choicesave<-gwindow("Save as ...",visible=FALSE,width=250,height=150)
      gcsave<-ggroup(container=choicesave,expand=FALSE)
    }
    lytsavere<-glayout(container=gcsave,spacing=13)
    mrenv$oksa<-gbutton("     OK    ",container=lytsavere)
    mrenv$cancelsa<-gbutton(" Cancel ",container=lytsavere)
    mrenv$radiosa<-gradio(c("Result1","Result2"),selected=1,horizontal=FALSE,container=lytsavere)
    lytsavere[2:3,2:5]<-mrenv$radiosa
    lytsavere[5,3]<-mrenv$oksa
    lytsavere[5,5]<-mrenv$cancelsa
    visible(choicesave)<-TRUE
    
    addHandlerClicked(mrenv$oksa,handler=function(h,...){
      if(svalue(mrenv$radiosa)=="Result1"){
        parms<-mrenv$parmsShow
        if(exists("parms")==FALSE||is.null(parms)==TRUE)
        {
          gmessage("There is something wrong in the first step!","Info",icon="info")
          return
        }else{
          output<-gfile(text="Save a file...",type="save",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          write.table(parms,output,sep = ",",row.names=FALSE,col.names = TRUE) 
        }
      }else{
        wan<-mrenv$wan
        if(exists("wan")==FALSE||is.null(wan)==TRUE)
        {
          gmessage("There is no result meets the requirements in the second step!","Info",icon="info")
          return
        }else{
          output<-gfile(text="Save a file...",type="save",
                        filter=list("All files"=list(patterns=c("*")),
                                    "CSV files"=list(patterns=c("*.csv"))))
          write.table(wan,output,sep = ",",row.names=FALSE,col.names = TRUE) 
        }
      }
    })
    
    addHandlerClicked(mrenv$cancelsa,handler=function(h,...){
      dispose(choicesave)
    })
  })  

  addHandlerClicked(manhattan,handler=function(h,...){
    if((exists("svgwline")==FALSE)||(svalue(gwedit)<=0))
    {
      gmessage("Please input correct genomewideline value!","Warning",icon="warning")
      return
    }else{
      svgwline<-svalue(gwedit)
      if(mrenv$manstandchoice==1)
      {
        mrenv$standline<--log10(mrenv$mannewp)
        mrenv$manstandchoice<-mrenv$manstandchoice+1
      }else{
        mrenv$standline<-svgwline
      }
      if(isExtant(plotwin)==FALSE)
      {
        plotwin<-gwindow("Manhattan Plot",visible=FALSE,width=600,height=360)
        gpw<-ggroup(container=plotwin)
        ggpw<-ggraphics(container=gpw)
      }
      addHandlerChanged(ggpw, handler=function(h,...) {
        colnames(mrenv$parms)<-c("Trait","Chromosome","Position","Mean","Sigma2","Sigma2_k","SNP effect","Sigma2_k_posteriori","Wald","P_wald")
        parms<-as.data.frame(mrenv$parms)
        plotman<-manhattan(parms,chr = "Chromosome",bp ="Position",p ="P_wald",snp="Sigma2",suggestiveline=FALSE,genomewideline = mrenv$standline)
      })
      visible(plotwin)<-TRUE
    }
  })
  
  addHandlerClicked(qqplot,handler=function(h,...){
    if((exists("svgwstandp")==FALSE)||(svalue(gwedit1)<=0))
    {
      gmessage("Please input correct standard P-value!","Warning",icon="warning")
      return
    }else{
      svgwstandp<-svalue(gwedit1)
      mrenv$standp<-svgwstandp
      if(isExtant(plotwin1)==FALSE)
      {
        plotwin1<-gwindow("Q-Q Plot",visible=FALSE,width=600,height=360)
        gpw1<-ggroup(container=plotwin1)
        ggpw1<-ggraphics(container=gpw1)
      }
      addHandlerChanged(ggpw1, handler=function(h,...) {
        pvalue<-matrix(mrenv$parms[,10],,1)
        observed<-sort(pvalue[,1])
        newobserved<-observed[which(observed<mrenv$standp)]
        lobs<--(log10(newobserved))
        expected<-c(1:length(newobserved))
        lexp<--(log10(expected/(length(expected)+1)))
        plot(lexp,lobs,xlab=expression('Expected -log'[10]*'(P)'),ylab=expression('Observed -log'[10]*'(P)'))
        abline(0,1,col="red")
      })
      visible(plotwin1)<-TRUE
    }
  })
}