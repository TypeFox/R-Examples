# Several functions to calculate GO semantic similarities (TCSSGetAncestors, TCSSCompute_ICA, GetOntology, GetGOParents, GetLatestCommonAncestor) are re-used from R package GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
# Reference: G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
# Modification time: 2011.11

.initial<-function(pos = 1,envir = as.environment(pos)){
  if(!exists("ppiPreEnv") || length(ppiPreEnv)<1) {
    packageStartupMessage("initializing ppiPre...")
    assign("ppiPreEnv",new.env(),envir=envir)  
    assign("ppiPreCache", new.env(),envir=envir)
    assign("ICEnv", new.env(),envir=envir)
    
    packageStartupMessage("done")
  }
}
################
KEGGSim <- function(protein1, protein2)    # KEGG-based similarity of two proteins
{
  
  if(!requireNamespace("KEGG.db")){ stop("package KEGG.db is needed.")}
  Pathway1 <- KEGG.db::KEGGEXTID2PATHID[[protein1]]
  Pathway2 <- KEGG.db::KEGGEXTID2PATHID[[protein2]]
  intersec <- length(na.omit(match(Pathway1, Pathway2)))
  if(intersec==0)
    sim<-0
  else
    sim<-intersec/(length(Pathway1)+length(Pathway2)-intersec)
  return(sim) 
}

GOKEGGSims <- function(gene1, gene2, organism="yeast", drop ="IEA")  #KEGG- and GO-based similarity of two proteins
{
  wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
  
  dropcodes <- drop
  Sims <- data.frame(protein1=gene1,protein2=gene2,BPWang=0,MFWang=0,CCWang=0,BPTCSS=0,MFTCSS=0,CCTCSS=0,BPIG=0,MFIG=0,CCIG=0,KEGGSim=0)
  
  Sims[[3]][1]<-GOSemSim::geneSim(gene1,gene2,ont="BP",organism=wh_organism,measure="Wang",drop = dropcodes)$geneSim
  Sims[[4]][1]<-GOSemSim::geneSim(gene1,gene2,ont="MF",organism=wh_organism,measure="Wang",drop = dropcodes)$geneSim
  Sims[[5]][1]<-GOSemSim::geneSim(gene1,gene2,ont="CC",organism=wh_organism,measure="Wang",drop = dropcodes)$geneSim
  
  Sims[[6]][1]<-TCSSGeneSim(gene1,gene2,ont="BP",organism=wh_organism,drop = dropcodes)$geneSim
  Sims[[7]][1]<-TCSSGeneSim(gene1,gene2,ont="MF",organism=wh_organism,drop = dropcodes)$geneSim
  Sims[[8]][1]<-TCSSGeneSim(gene1,gene2,ont="CC",organism=wh_organism,drop = dropcodes)$geneSim
  
  Sims[[9]][1]<-IntelliGOGeneSim(gene1,gene2,ont="BP",organism=wh_organism,drop = dropcodes)$geneSim
  Sims[[10]][1]<-IntelliGOGeneSim(gene1,gene2,ont="MF",organism=wh_organism,drop = dropcodes)$geneSim
  Sims[[11]][1]<-IntelliGOGeneSim(gene1,gene2,ont="CC",organism=wh_organism,drop = dropcodes)$geneSim
  
  Sims[[12]][1]<-KEGGSim(gene1,gene2)
  return(Sims)
}

GOKEGGSimsFromFile <- function(input,output="GOKEGGSims-ppiPre.csv",header=TRUE,sep=",", organism="yeast", drop ="IEA") ##KEGG- and GO-based similarity of protein pairs in an input file
{	
  cache<-read.csv(file=input,header=header,sep=sep)
  wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
  
  dropcodes <- drop
  SimsFromFile<-data.frame(protein1=cache[1],protein2=cache[2],BPWang=0,MFWang=0,CCWang=0,BPTCSS=0,MFTCSS=0,CCTCSS=0,BPIG=0,MFIG=0,CCIG=0,KEGGSim=0)
  i<-1
  for(n in 1:length(cache[[1]]))
  {
    message(paste("Computing GO- & KEGG-based similarities of",as.character(cache[[1]][i]),"and", as.character(cache[[2]][i])))
    SimsFromFile[[3]][i]<-GOSemSim::geneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="BP",organism=wh_organism,measure="Wang",drop = dropcodes )$geneSim
    SimsFromFile[[4]][i]<-GOSemSim::geneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="MF",organism=wh_organism,measure="Wang",drop = dropcodes )$geneSim
    SimsFromFile[[5]][i]<-GOSemSim::geneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="CC",organism=wh_organism,measure="Wang",drop = dropcodes )$geneSim
    SimsFromFile[[6]][i]<-TCSSGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="BP",organism=wh_organism,drop = dropcodes )$geneSim
    SimsFromFile[[7]][i]<-TCSSGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="MF",organism=wh_organism,drop = dropcodes )$geneSim
    SimsFromFile[[8]][i]<-TCSSGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="CC",organism=wh_organism,drop = dropcodes )$geneSim
    SimsFromFile[[9]][i]<-IntelliGOGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="BP",organism=wh_organism,drop = dropcodes )$geneSim
    SimsFromFile[[10]][i]<-IntelliGOGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="MF",organism=wh_organism,drop = dropcodes )$geneSim
    SimsFromFile[[11]][i]<-IntelliGOGeneSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]),ont="CC",organism=wh_organism,drop = dropcodes )$geneSim
    SimsFromFile[[12]][i]<-KEGGSim(as.character(cache[[1]][i]),as.character(cache[[2]][i]))	
    i <- i+1
  }
  write.csv(SimsFromFile,file=output,row.names=FALSE)
}
################

TCSSGetChildren <- function(ont="MF") {
  if(!exists("ppiPreEnv")) .initial()
  
  wh_Children <- switch(ont,
                        MF = "MFChildren",
                        BP = "BPChildren",
                        CC = "CCChildren"
  )	
  
  Children <- switch(ont,
                     MF = AnnotationDbi::as.list(GOMFCHILDREN) ,
                     BP = AnnotationDbi::as.list(GOBPCHILDREN) , 
                     CC = AnnotationDbi::as.list(GOCCCHILDREN)	
  )
  assign(eval(wh_Children), Children, envir=ppiPreEnv)
}

# Re-used from GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
# Reference: G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
# Modification time: 2011.11
TCSSGetAncestors <- function(ont="MF") { 
  if(!exists("ppiPreEnv")) .initial()
  wh_Ancestors <- switch(ont,
                         MF = "MFAncestors",
                         BP = "BPAncestors",
                         CC = "CCAncestors"	
  )		
  Ancestors <- switch(ont,
                      MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
                      BP = AnnotationDbi::as.list(GOBPANCESTOR) , 
                      CC = AnnotationDbi::as.list(GOCCANCESTOR)	
  )
  assign(eval(wh_Ancestors), Ancestors, envir=ppiPreEnv)
}

CheckAnnotationPackage <- function(species){
  if (species == "human")
    if(!requireNamespace("org.Hs.eg.db"))
      stop("The package org.Hs.eg.db is needed.")
  if (species == "yeast")
    if(!requireNamespace("org.Sc.sgd.db"))
      stop("The package org.Sc.sgd.db is needed.")
  if (species == "fly")
    if(!requireNamespace("org.Dm.eg.db"))
      stop("The package org.Dm.eg.db is needed.")
  if (species == "mouse")
    if(!requireNamespace("org.Mm.eg.db"))
      stop("The package org.Mm.eg.db is needed.")
  if (species == "rat")
    if(!requireNamespace("org.Rn.eg.db"))
      stop("The package org.Rn.eg.db is needed.")
  if (species == "zebrafish")
    if(!requireNamespace("org.Dr.eg.db"))
      stop("The package org.Dr.eg.db is needed.")
  if (species == "worm")
    if(!requireNamespace("org.Ce.eg.db"))
      stop("The package org.Ce.eg.db is needed.")
  if (species == "arabidopsis")
    if(!requireNamespace("org.At.tair.db"))
      stop("The package org.At.tair.db is needed.")
  if (species == "ecolik12")
    if(!requireNamespace("org.EcK12.eg.db"))
      stop("The package org.EcK12.eg.db is needed.")
  if (species == "bovine")
    if(!requireNamespace("org.Bt.eg.db"))
      stop("The package org.Bt.eg.db is needed.")
  if (species == "canine")
    if(!requireNamespace("org.Cf.eg.db"))
      stop("The package org.Cf.eg.db is needed.")
  if (species == "anopheles")
    if(!requireNamespace("org.Ag.eg.db"))
      stop("The package org.Ag.eg.db is needed.")
  if (species == "ecsakai")
    if(!requireNamespace("org.EcSakai.eg.db"))
      stop("The package org.EcSakai.eg.db is needed.")
  if (species == "chicken")
    if(!requireNamespace("org.Gg.eg.db"))
      stop("The package org.Gg.eg.db is needed.")
  if (species == "chimp")
    if(!requireNamespace("org.Pt.eg.db"))
      stop("The package org.Pt.eg.db is needed.")
  if (species == "malaria")
    if(!requireNamespace("org.Pf.plasmo.db"))
      stop("The package org.Pf.plasmo.db is needed.")
  if (species == "rhesus")
    if(!requireNamespace("org.Mmu.eg.db"))
      stop("The package org.Mmu.eg.db is needed.")
  if (species == "pig")
    if(!requireNamespace("org.Ss.eg.db"))
      stop("The package org.Ss.eg.db is needed.")
  if (species == "xenopus")
    if(!requireNamespace("org.Xl.eg.db"))
      stop("The package org.Xl.eg.db is needed.")
}


# Re-used from GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
# Reference: G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
# Modification time: 2011.11
`GetOntology` <-  function(gene, organism, ontology, dropCodes) {
  .initial() 
  species <- switch(organism,
                    human = "Hs",
                    fly = "Dm",
                    mouse = "Mm",
                    rat = "Rn",
                    yeast = "Sc",
                    zebrafish = "Dr",
                    worm = "Ce",
                    arabidopsis = "At",
                    ecolik12 = "EcK12",
                    bovine	= "Bt",
                    canine	= "Cf", 
                    anopheles	=	"Ag", 
                    ecsakai	=	"EcSakai", 
                    chicken	=	"Gg", 
                    chimp	=	"Pt", 
                    malaria	=	"Pf", 
                    rhesus	=	"Mmu", 
                    pig	= "Ss", 
                    xenopus	=	"Xl"
  )
  
  if (!exists(species, envir=ppiPreEnv)) {
    
    CheckAnnotationPackage(organism) #download and install the packages
    
    gomap <- switch(organism,
                    human = org.Hs.eg.db::org.Hs.egGO,
                    fly = org.Dm.eg.db::org.Dm.egGO,
                    mouse = org.Mm.eg.db::org.Mm.egGO,
                    rat = org.Rn.eg.db::org.Rn.egGO,
                    yeast = org.Sc.sgd.db::org.Sc.sgdGO,
                    zebrafish = org.Dr.eg.db::org.Dr.egGO,
                    worm = org.Ce.eg.db::org.Ce.egGO,
                    arabidopsis = org.At.tair.db::org.At.tairGO,
                    ecoli = org.EcK12.eg.db::org.EcK12.egGO,
                    bovine	= org.Bt.eg.db::org.Bt.egGO,
                    canine	= org.Cf.eg.db::org.Cf.egGO, 
                    anopheles	=	org.Ag.eg.db::org.Ag.egGO, 
                    ecsakai	=	org.EcSakai.eg.db::org.EcSakai.egGO, 
                    chicken	=	org.Gg.eg.db::org.Gg.egGO, 
                    chimp	=	org.Pt.eg.db::org.Pt.egGO, 
                    malaria	=	org.Pf.plasmo.db::org.Pf.plasmoGO, 
                    rhesus	=	org.Mmu.eg.db::org.Mmu.egGO, 
                    pig	= org.Ss.eg.db::org.Ss.egGO, 
                    xenopus	=	org.Xl.eg.db::org.Xl.egGO
    )
    assign(eval(species), gomap, envir=ppiPreEnv) 
  }
  gomap <- get(species, envir=ppiPreEnv) 
  
  allGO <- gomap[[gene]] 
  if (is.null(allGO)) {
    return (NA)
  }
  if (sum(!is.na(allGO)) == 0) {
    return (NA)
  }
  if(!is.null(dropCodes)) { 
    evidence<-sapply(allGO, function(x) x$Evidence) 
    drop<-evidence %in% dropCodes 
    allGO<-allGO[!drop] 
  }
  
  category<-sapply(allGO, function(x) x$Ontology) 
  allGO<-allGO[category %in% ontology] 
  
  if(length(allGO)==0) return (NA)
  return (unlist(unique(names(allGO)))) #return the GOIDs
}

# Re-used from GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
# Reference: G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
# Modification time: 2011.11
GetLatestCommonAncestor<-function(GOID1, GOID2, ont, organism){
  #message("Calulating Latest Common Ancestor...")
  if(!exists("ppiPreEnv")) .initial()
  fname <- paste("Info_Contents",organism, ont, sep="_")
  tryCatch(utils::data(list=fname, package="GOSemSim", envir=ICEnv))  
  InfoContents <- get("IC", envir=ICEnv)
  
  rootCount <- max(InfoContents[InfoContents != Inf])
  InfoContents["all"] = 0
  p1 <- InfoContents[GOID1]/rootCount
  p2 <- InfoContents[GOID2]/rootCount    
  if(is.na(p1) || is.na(p2)) return (NA)
  if (p1 == 0 || p2 == 0) return (NA)
  Ancestor.name <- switch(ont,MF = "MFAncestors",BP = "BPAncestors",CC = "CCAncestors")	
  if (!exists(Ancestor.name, envir=ppiPreEnv)) {
    TCSSGetAncestors(ont)
  }
  
  Ancestor <- get(Ancestor.name, envir=ppiPreEnv)					
  ancestor1 <- unlist(Ancestor[GOID1])
  ancestor2 <- unlist(Ancestor[GOID2])
  if (GOID1 == GOID2) { 
    commonAncestor <- GOID1
  } else if (GOID1 %in% ancestor2) { 
    commonAncestor <- GOID1
  } else if (GOID2 %in% ancestor1) {
    commonAncestor <- GOID2
  } else { 
    commonAncestor <- intersect(ancestor1, ancestor2)
  }
  if (length(commonAncestor) == 0)
    LCA<-NULL
  max<- -100
  LCA<-NULL
  for(a in commonAncestor){
    if(!is.na(InfoContents[a])) {
      if(InfoContents[a]>max){
        max<-InfoContents[a]
        LCA<-a
      }
    }
  }
  #message("done...")
  return (LCA)
  
}

# Re-used from GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
# Reference: G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
# Modification time: 2011.11
TCSSCompute_ICA<- function(dropCodes="IEA", ont, organism) {
  
  wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
  wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
  
  ICA.name <- paste(wh_organism,wh_ont,"ICA",sep="_")
  if (exists(ICA.name, envir=ppiPreEnv)) {
    return(get(ICA.name, envir=ppiPreEnv))
  }
  CheckAnnotationPackage(wh_organism)
  #message("Calulating ICA...")
  gomap <- switch(wh_organism,
                  human = org.Hs.eg.db::org.Hs.egGO,
                  fly = org.Dm.eg.db::org.Dm.egGO,
                  mouse = org.Mm.eg.db::org.Mm.egGO,
                  rat = org.Rn.eg.db::org.Rn.egGO,
                  yeast = org.Sc.sgd.db::org.Sc.sgdGO,
                  zebrafish = org.Dr.eg.db::org.Dr.egGO,
                  worm = org.Ce.eg.db::org.Ce.egGO,
                  arabidopsis = org.At.tair.db::org.At.tairGO,
                  ecoli = org.EcK12.eg.db::org.EcK12.egGO,
                  bovine	= org.Bt.eg.db::org.Bt.egGO,
                  canine	= org.Cf.eg.db::org.Cf.egGO, 
                  anopheles	=	org.Ag.eg.db::org.Ag.egGO, 
                  ecsakai	=	org.EcSakai.eg.db::org.EcSakai.egGO, 
                  chicken	=	org.Gg.eg.db::org.Gg.egGO, 
                  chimp	=	org.Pt.eg.db::org.Pt.egGO, 
                  malaria	=	org.Pf.plasmo.db::org.Pf.plasmoGO, 
                  rhesus	=	org.Mmu.eg.db::org.Mmu.egGO, 
                  pig	= org.Ss.eg.db::org.Ss.egGO, 
                  xenopus	=	org.Xl.eg.db::org.Xl.egGO
  )
  mapped_genes <- mappedkeys(gomap)
  gomap = AnnotationDbi::as.list(gomap[mapped_genes])
  if (!is.null(dropCodes)){
    gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
    gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
    gomap<-gomap[sapply(gomap,length) >0]		
  }else {
    gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
  }
  
  goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE)	
  goids <- toTable(GOTERM)
  
  goids <- unique(goids[goids[,"Ontology"] == wh_ont, "go_id"])  	
  gocount <- table(goterms)
  goname <- names(gocount) 
  
  go.diff <- setdiff(goids, goname)
  m <- double(length(go.diff)) 
  names(m) <- go.diff
  gocount <- as.vector(gocount)
  names(gocount) <- goname
  gocount <- c(gocount, m)
  Children.name <- switch(wh_ont,
                          MF = "MFChildren",
                          BP = "BPChildren",
                          CC = "CCChildren"	
  )	
  if (!exists(Children.name, envir=ppiPreEnv)) {
    TCSSGetChildren(wh_ont)
  }
  Children<- get(Children.name, envir=ppiPreEnv)	
  cnt <- sapply(goids,function(x){ c=gocount[unlist(Children[x])]; gocount[x]+sum(c[!is.na(c)])})
  cnt<-cnt+1;
  names(cnt)<-goids;
  ICA<- -log(cnt/sum(gocount))
  #message("done...")
  assign(eval(ICA.name), ICA, envir=ppiPreEnv)
  return (ICA)
}	

`TCSSGoSim` <- function(GOID1, GOID2, ont, organism, ICA) {
  if(!exists("ppiPreEnv")) .initial()
  AnnotationIC <- ICA
  
  AnnotationIC["all"] = 0
  lca<-GetLatestCommonAncestor(GOID1, GOID2, ont, organism)
  if(length(lca)==0)
    sim<-0 
  ica_max<-max(AnnotationIC)
  ica_lca<-max(AnnotationIC[lca])
  sim<-ica_lca/ica_max
  return(round(sim, digits=3))
}

`TCSSGeneSim` <-
  function(gene1, gene2, ont="MF", organism="yeast", drop="IEA"){
    
    wh_ont <- match.arg(ont, c("MF", "BP", "CC")) 
    wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
    
    go1 <- GetOntology(gene1, organism= wh_organism, ontology= wh_ont, dropCodes=drop) 
    go2 <- GetOntology(gene2, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
        
    if (is.na(go1)|| is.na(go2)) {
      return (list(geneSim=NA, GO1=go1, GO2=go2)) 
    }
    go1<-unlist(go1)
    go2<-unlist(go2)
    m<-length(go1)
    n<-length(go2)
    scores <- matrix(nrow=m, ncol=n)
    rownames(scores) <- go1
    colnames(scores) <- go2
    ICA<-TCSSCompute_ICA(drop,wh_ont,wh_organism) 
    
    for( i in 1:m) {
      for (j in 1:n) {
        scores[i,j] <-TCSSGoSim(go1[i], go2[j], wh_ont, wh_organism, ICA)
      }
    }
    if (!sum(!is.na(scores))) return (list(geneSim=NA, GO1=go1, GO2=go2)) 
    sim<-max(scores)
    sim <- round(sim, digits=3)
    
    return (list(geneSim=sim, GO1=go1, GO2=go2)) 
  }	

# Re-used from GOSemSim authored by Guangchuang Yu <guangchuangyu@gmail.com>. 
# Reference: G. Yu, F. Li, Y. Qin, X. Bo, Y. Wu, and S. Wang, "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products", Bioinformatics, vol. 26, no. 7, pp. 976-978, Apr. 2010.
# Modification time: 2011.11
GetGOParents <- function(ont="MF") {
  if(!exists("ppiPreEnv")) .initial()
  
  wh_Parents <- switch(ont,
                       MF = "MFParents",
                       BP = "BPParents",
                       CC = "CCParents"	
  )
  
  Parents <- switch(ont,
                    MF = AnnotationDbi::as.list(GOMFPARENTS) ,
                    BP = AnnotationDbi::as.list(GOBPPARENTS) , 
                    CC = AnnotationDbi::as.list(GOCCPARENTS)	
  )
  assign(eval(wh_Parents), Parents, envir=ppiPreEnv)
}

IntelliGOInverseAnnotationFrequency<-function(dropCodes="IEA", goid,ont,organism){ 
  wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
  wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
  
  CheckAnnotationPackage(wh_organism)
  gomap <- switch(wh_organism,
                  human = org.Hs.eg.db::org.Hs.egGO,
                  fly = org.Dm.eg.db::org.Dm.egGO,
                  mouse = org.Mm.eg.db::org.Mm.egGO,
                  rat = org.Rn.eg.db::org.Rn.egGO,
                  yeast = org.Sc.sgd.db::org.Sc.sgdGO,
                  zebrafish = org.Dr.eg.db::org.Dr.egGO,
                  worm = org.Ce.eg.db::org.Ce.egGO,
                  arabidopsis = org.At.tair.db::org.At.tairGO,
                  ecoli = org.EcK12.eg.db::org.EcK12.egGO,
                  bovine	= org.Bt.eg.db::org.Bt.egGO,
                  canine	= org.Cf.eg.db::org.Cf.egGO, 
                  anopheles	=	org.Ag.eg.db::org.Ag.egGO, 
                  ecsakai	=	org.EcSakai.eg.db::org.EcSakai.egGO, 
                  chicken	=	org.Gg.eg.db::org.Gg.egGO, 
                  chimp	=	org.Pt.eg.db::org.Pt.egGO, 
                  malaria	=	org.Pf.plasmo.db::org.Pf.plasmoGO, 
                  rhesus	=	org.Mmu.eg.db::org.Mmu.egGO, 
                  pig	= org.Ss.eg.db::org.Ss.egGO, 
                  xenopus	=	org.Xl.eg.db::org.Xl.egGO
  )
  
  if (!is.null(dropCodes)){
    gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
    gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
    gomap<-gomap[sapply(gomap,length) >0]		
  }else {
    gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
  }
  
  Gti <- switch(wh_organism,
                human = length(org.Hs.eg.db::org.Hs.egGO2EG[[goid]])	,
                fly = length(org.Dm.eg.db::org.Dm.egGO2EG[[goid]]),
                yeast = length(org.Sc.sgd.db::org.Sc.sgdGO2ORF[[goid]]),
                worm = length(org.Ce.eg.db::org.Ce.egGO2EG[[goid]]),	
                mouse = length(org.Mm.eg.db::org.Mm.egGO2EG[[goid]]),
                rat = length(org.Rn.eg.db::org.Rn.egGO2EG[[goid]]),
                zebrafish = length(org.Dr.eg.db::org.Dr.egGO2EG[[goid]]),
                arabidopsis = length(org.At.tair.db::org.At.tairGO2TAIR[[goid]]),
                ecoli = length(org.EcK12.eg.db::org.EcK12.egGO2EG[[goid]]),
                bovine	= length(org.Bt.eg.db::org.Bt.egGO2EG[[goid]]),
                canine	= length(org.Cf.eg.db::org.Cf.egGO2EG[[goid]]), 
                anopheles	=	length(org.Ag.eg.db::org.Ag.egGO2EG[[goid]]), 
                ecsakai	=	length(org.EcSakai.eg.db::org.EcSakai.egGO2EG[[goid]]), 
                chicken	=	length(org.Gg.eg.db::org.Gg.egGO2EG[[goid]]), 
                chimp	=	length(org.Pt.eg.db::org.Pt.egGO2EG[[goid]]), 
                malaria	=	length(org.Pf.plasmo.db::org.Pf.plasmoGO2ORF[[goid]]), 
                rhesus	=	length(org.Mmu.eg.db::org.Mmu.egGO2EG[[goid]]), 
                pig	= length(org.Ss.eg.db::org.Ss.egGO2EG[[goid]]), 
                xenopus	=	length(org.Xl.eg.db::org.Xl.egGO2EG[[goid]])
  )
  Gtot <- length(mappedkeys(gomap))
  IAF <- log(Gtot/Gti)	
  return (IAF)
}

IntelliGOGetParents<-function(goid,ont,organsim){
  if(!exists("ppiPreEnv")) .initial()
  Parents.name <- switch(ont,
                         MF = "MFParents",
                         BP = "BPParents",
                         CC = "CCParents"	
  )
  if (!exists(Parents.name, envir=ppiPreEnv)){ 
    GetGOParents(ont)
  }
  Parents <- get(Parents.name, envir=ppiPreEnv)
  p<- Parents[goid]
  par<-p[[1]][[1]]
  return (par)
}

IntelliGOGetDepth<-function(goid,ont,organsim){
  wh_ont<-ont
  depth<- 0
  parent<- IntelliGOGetParents(goid,wh_ont,organsim)
  if(length(parent)==0)
    return (depth+1)
  while(length(parent)!=0){
    depth<- depth+1
    parent<- IntelliGOGetParents(parent,wh_ont,organsim)
  }
  return (depth+1)
}
IntelliGOComputeEiEj<-function(GOID1,GOID2,ont,organsim){
  EiEj <- 0
  lca <- GetLatestCommonAncestor(GOID1,GOID2,ont,organsim)
  depth_ei <- IntelliGOGetDepth(GOID1,ont,organsim)
  depth_ej <- IntelliGOGetDepth(GOID2,ont,organsim)
  depth_lca <- IntelliGOGetDepth(lca,ont,organsim)
  EiEj <- 2*depth_lca/(depth_ei+depth_ej)
  return (EiEj)
}
IntelliGOGOSim<-function(GOID1,GOID2,w1,w2,ont,organsim){
  sim<-0
  go1<-unlist(GOID1)
  go2<-unlist(GOID2)
  if(length(go1)==0||length(go2)==0)
    return (0)
  a<-IntelliGOComputeEiEj(go1,go2,ont,organsim)
  x<-sqrt(IntelliGOComputeEiEj(go1,go1,ont,organsim))
  y<-sqrt(IntelliGOComputeEiEj(go2,go2,ont,organsim))
  sim<- a/(x*y)
  return(round(sim, digits=3))
}

`IntelliGOGeneSim` <-
  function(gene1, gene2, w1=1,w2=1,ont="MF", organism="yeast", drop="IEA"){
    wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
    wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine","canine","anopheles","ecsakai","chicken","chimp","malaria","rhesus","pig","xenopus"))
    
    go1 <- GetOntology(gene1, organism= wh_organism, ontology= wh_ont, dropCodes=drop) 
    go2 <- GetOntology(gene2, organism= wh_organism, ontology= wh_ont, dropCodes=drop)
    if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
      return (list(geneSim=NA, GO1=go1, GO2=go2)) 
    }
    weight_go1<-w1
    weight_go2<-w2
    go1<-unlist(go1)
    go2<-unlist(go2)
    m<-length(go1)
    n<-length(go2)
    scores <- matrix(nrow=m, ncol=n)
    rownames(scores) <- go1
    colnames(scores) <- go2
    
    for( i in 1:m ) {
      for (j in 1:n) {
        scores[i,j] <-IntelliGOGOSim(go1[i], go2[j],weight_go1,weight_go2,wh_ont, wh_organism)
      }
    }
    if (!sum(!is.na(scores))) return(list(geneSim=NA, GO1=go1, GO2=go2)) 
    if (n ==1 || m == 1) {
      sim<-max(scores)
    }
    
    sim <- (sum(sapply(1:m, function(x) {max(scores[x,], na.rm=TRUE)})) + sum(sapply(1:n, function(x) {max(scores[,x], na.rm=TRUE)})))/(m+n)
    sim <- round(sim, digits=3)
    return (list(geneSim=sim, GO1=go1, GO2=go2))
  }