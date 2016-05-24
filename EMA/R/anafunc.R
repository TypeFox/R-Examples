runHyperGO <- function(list, pack.annot, categorySize = 1, verbose = TRUE, name = "hyperGO", htmlreport = TRUE, txtreport = TRUE, tabResult = FALSE, pvalue = 0.05) {
    if (require(GOstats))
    {
        require(pack.annot, character.only=TRUE)
        
        filehtmlGO <- paste(name, ".GO.html", sep="")
        filetxtCC <- paste(name, ".GO.CC.txt", sep="")
        filetxtBP <- paste(name, ".GO.BP.txt", sep="")
        filetxtMF <- paste(name, ".GO.MF.txt", sep="")
        
        pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep = "")))
        pack.annot.ACC <- eval(as.name(paste(gsub(".db", "", pack.annot), "ACCNUM", sep = "")))
        pack.annot.GO <- eval(as.name(paste(gsub(".db", "", pack.annot), "GO", sep = "")))
        
        listALL <- list    
        
        ## ENTREZ
        designALL <- names(unlist(as.list(pack.annot.ACC)))
        entrezIds <- mget(designALL, envir = pack.annot.EID, ifnotfound=NA)
        haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
        print ("ENTREZ ID done.")
        
        ## GO
        haveGo <- sapply(mget(haveEntrezId, pack.annot.GO, ifnotfound=NA),
                         function(x) {
                             if (length(x) == 1 && is.na(x)) {
                                 FALSE
                             }else {
                                 TRUE
                             }
                         })
        numNoGO <- sum(!haveGo)
        if (verbose)
            print(paste(numNoGO, "pbsets in Universe have no GO ids"))
        
        designGO <- haveEntrezId[haveGo]
        
        ## UNIQUE UNIVERSE
        entrezUniverse <- na.omit(unique(unlist(mget(designALL, pack.annot.EID, ifnotfound=NA))))
        if (any(duplicated(entrezUniverse))){
            stop("error in gene universe: can't have duplicate Ent")
        }
        
        ## SELECTED GENE LIST
        listGO <- intersect(listALL, designGO)
        numNoGO <- length(listALL) - length(listGO)
        if (verbose)
            print (paste(numNoGO, " pbsets in the list have no GO ids"))
        
        selectedEntrezIds <- na.omit(unique(unlist(mget(listGO, pack.annot.EID, ifnotfound=NA))))
        
        ## GO ANALYSIS
        
        htmlheader(paste(date(), "<br>Gene Ontology HyperGeometric Analysis<br>","<br>List length:", length(listALL),"<br><br>"), filename = filehtmlGO)
        
        params.BP.over <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = entrezUniverse, annotation = pack.annot, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
        params.MF.over <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = entrezUniverse, annotation = pack.annot, ontology = "MF", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
        params.CC.over <- new("GOHyperGParams", geneIds = selectedEntrezIds, universeGeneIds = entrezUniverse, annotation = pack.annot, ontology = "CC", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")
        
        hgOver.BP <- hyperGTest(params.BP.over)
        hgOver.MF <- hyperGTest(params.MF.over)
        hgOver.CC <- hyperGTest(params.CC.over)
        
        if (htmlreport){
            htmlresult(hgOver.BP, filename = filehtmlGO, app = TRUE, categorySize = categorySize, pvalue = pvalue)
            htmlresult(hgOver.MF, filename = filehtmlGO, app = TRUE, categorySize = categorySize, pvalue = pvalue)
            htmlresult(hgOver.CC, filename = filehtmlGO, app = TRUE, categorySize = categorySize, pvalue = pvalue)
        }
        
        if (txtreport){
            goReport(hgOver.BP, fileout = filetxtBP, type = "BP", pack.annot = pack.annot, categorySize = categorySize, pvalue = pvalue)
            goReport(hgOver.MF, fileout = filetxtMF, type = "MF", pack.annot = pack.annot, categorySize = categorySize, pvalue = pvalue)
            goReport(hgOver.CC, fileout = filetxtCC, type = "CC", pack.annot = pack.annot, categorySize = categorySize, pvalue = pvalue)
        }
        
        if (tabResult){
            res1 <- summary(hgOver.BP, pvalue = pvalue, categorySize = categorySize)
            res2 <- summary(hgOver.MF, pvalue = pvalue, categorySize = categorySize)
            res3 <- summary(hgOver.CC, pvalue = pvalue, categorySize = categorySize)
            res <- list(BP = res1, MF = res2, CC = res3)
            return(res)
        } else{
            return(list(BP = hgOver.BP, MF = hgOver.MF, CC = hgOver.CC))
        }
    } else stop("Failed to load required package GOstats.")
    
}

runHyperKEGG <- function(list, pack.annot, categorySize = 1, name = "hyperKEGG", htmlreport = TRUE, txtreport = TRUE, tabResult = FALSE, pvalue = 0.05) {
    if (require(GOstats)){
        require(pack.annot, character.only=TRUE)
        
        pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep="")))
        pack.annot.ACC <- eval(as.name(paste(gsub(".db", "", pack.annot), "ACCNUM", sep="")))
        pack.annot.KEGG <- gsub(".db", "",  pack.annot)

        listALL <- list
        
        filehtmlKEGG <- paste(name, ".KEGG.html", sep="")
        filetxtKEGG <- paste(name, ".KEGG.txt", sep="") 
        
        ## ENTREZ
        designALL <- names(unlist(as.list(pack.annot.ACC)))
        entrezIds <- mget(designALL, envir = pack.annot.EID, ifnotfound=NA)
        haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
        designKEGG<-haveEntrezId
        
        ## UNIQUE UNIVERSE
        entrezUniverse <- unique(unlist(mget(designKEGG, pack.annot.EID, ifnotfound=NA)))
        
        listKEGG <- intersect(listALL, designKEGG)
        selectedEntrezIds <- unlist(mget(listKEGG, pack.annot.EID, ifnotfound=NA))
        
        ## UNIQUE LISTE
        selectedEntrezIds <- na.omit(unique(unlist(selectedEntrezIds)))
        params.KEGG.over <-new("KEGGHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse, annotation = pack.annot.KEGG, pvalueCutoff=0.05, testDirection="over")
        
        hgOver.KEGG <- hyperGTest(params.KEGG.over)
        
        if (htmlreport){
            htmlheader(paste(date(), "<br>KEGG HyperGeometric Analysis<br>", "<br>List length:", length(listALL),"<br><br>"), filename = filehtmlKEGG)
            htmlresult(hgOver.KEGG, filename = filehtmlKEGG, app = TRUE, categorySize, pvalue = pvalue)
        }
        
        if (txtreport)
            keggReport(hgOver.KEGG, fileout = filetxtKEGG, pack.annot = pack.annot, pvalue = pvalue)
        
        if (tabResult) {
            res <- summary(hgOver.KEGG, pvalue = pvalue)
            return(res)
        }
    } else stop("Failed to load required package GOstats.")
}
####
##
## HTML reports
##
####

htmlresult <- function(hgOver, filename, app = FALSE, categorySize = 1, pvalue = 0.05) {
  write(paste("<br>------<br>", length(geneIdUniverse(hgOver, cond = conditional(hgOver))), testName(hgOver)[1], testName(hgOver)[2], "ids tested(", dim(summary(hgOver, pvalue = 0.05))[1], " p<0.05)<br>"), file = filename, append = app)
  write(paste("Gene universe size:", universeMappedCount(hgOver)), file = filename, append = TRUE)
  write(paste("<br>Selected gene set size:", geneMappedCount(hgOver)), file = filename, append = TRUE)
  write(paste("<br>Conditional:", conditional(hgOver)), file = filename, append = TRUE)
  write(paste("<br>Annotation:", annotation(hgOver), "<br><br>"), file = filename, append = TRUE)
  htmlReport(hgOver, summary.args = list(htmlLinks = TRUE, categorySize = categorySize, pvalue = pvalue), file = filename, append = TRUE)
}

htmlheader <- function(towrite, filename) {
  write(paste("<b", towrite, "</b>"), file = filename)
}



keggReport <-function(hgOver, fileout = "report.txt", pack.annot, pvalue = 0.05) {
  pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep="")))
  pack.annot.SYMBOL <- eval(as.name(paste(gsub(".db", "", pack.annot), "SYMBOL", sep="")))
  
  write(file = fileout, paste("KEGG HyperGTest report\n"), append = FALSE, sep = ",")
  entrez <- unlist(as.list(pack.annot.EID))

  a <- geneIdsByCategory(hgOver)
  b <- geneIdUniverse(hgOver, cond=conditional(hgOver))
  
  a <- a[sigCategories(hgOver, pvalue)]
  b <- b[sigCategories(hgOver, pvalue)]
  
  for (i in as.vector(unlist(attributes(a)))) {
    if (length(b[[i]])>10) {
      write(file = fileout, paste("\n",i,":", length(a[[i]]), " geneIds | ", length(b[[i]]), " universIds"), append = TRUE, sep = ",", ncolumns = length(a[[i]]))
      
      pbset <- unique(names(entrez[which(is.element(entrez, a[[i]]))]))
      pbset <- intersect(pbset, unique(names(entrez[is.element(entrez,geneIds(hgOver))])))
      write(file = fileout, paste(length(pbset),"ProbeSets :"), append = TRUE, sep = ",")
      write(file = fileout, pbset, append = TRUE, sep=",", ncolumns = length(pbset))
      
      gs <- unique(unlist(mget(pbset, pack.annot.SYMBOL)))
      write(file = fileout, paste(length(gs), "GeneSymbols :"), append=TRUE, sep=",")
      write(file = fileout, gs, append = TRUE, sep=",", ncolumns = length(gs))     
    }
  }
}

goReport <- function(hgOver, fileout = "report.txt", type = c("CC", "MF", "BP"), pack.annot, pvalue = 0.05,  categorySize = 1) {

  pack.annot.ALLPROBES <- eval(as.name(paste(gsub(".db", "", pack.annot), "GO2ALLPROBES", sep = "")))
  pack.annot.EID <- eval(as.name(paste(gsub(".db", "", pack.annot), "ENTREZID", sep = "")))
  pack.annot.SYMBOL <- eval(as.name(paste(gsub(".db", "", pack.annot), "SYMBOL", sep = "")))

  if (type == "CC") {
    write(file = fileout, paste("Component Cellular - GO HyperGTest report\n"), append = FALSE, sep = ",")
    probes <- mget('GO:0005575', pack.annot.ALLPROBES)
  } else if (type == "BP") {
    write(file = fileout, paste("Biological Process - GO HyperGTest report\n"), append = FALSE, sep = ",")
    probes <- mget('GO:0008150', pack.annot.ALLPROBES)
  } else if (type == "MF") {
    write(file = fileout, paste("Molecular Function - GO HyperGTest report\n"), append = FALSE, sep = ",")
    probes <- mget('GO:0003674', pack.annot.ALLPROBES)
  }
  entrez <- mget(unique(unlist(probes)), pack.annot.EID)

  a <- geneIdsByCategory(hgOver)
  b <- geneIdUniverse(hgOver, cond=conditional(hgOver))
  
  a <- a[sigCategories(hgOver, pvalue)]
  b <- b[sigCategories(hgOver, pvalue)]
  
  for (i in as.vector(unlist(attributes(a)))) {
    if (length(b[[i]]) > categorySize) {
      write(file = fileout, paste("\n",i,":", length(a[[i]]), " geneIds | ", length(b[[i]]), " universIds"), append = TRUE, sep=",", ncolumns = length(a[[i]]))
      
      pbset <- unique(names(unlist(entrez[which(entrez%in%a[[i]])])))
      pbset <- intersect(pbset, unique(names(entrez[is.element(entrez,geneIds(hgOver))])))
      write(file = fileout, paste(length(pbset),"ProbeSets :"), append = TRUE, sep = ",")
      write(file = fileout, pbset, append = TRUE, sep = ",",ncolumns = length(pbset))
      
      gs <- unique(unlist(mget(pbset, pack.annot.SYMBOL)))
      write(file = fileout, paste(length(gs), "GeneSymbols :"), append = TRUE, sep = ",")
      write(file = fileout, gs, append = TRUE, sep=",", ncolumns = length(gs))     
    }
  }
}


####
##
## Annotations
##
####

bioMartAnnot <-  function (data, inputTypeId, outputTypeId = c("entrezgene", "hgnc_symbol", "ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position", "band", "strand"), dataset= c("hsapiens_gene_ensembl"), database = "ensembl",  sort.by = NULL, outfile = NA) {

    if (missing(data) || missing(inputTypeId)){
      stop("Error : Data or inputTypeId not found")  
    }
    if (missing(database) || missing(dataset)){
      stop("Error : You have to define the database and the dataset you want to use")  
    }
    else{
      mart <- useMart(database, dataset = dataset)
    }
    
    ##Input - vector or matrix
    if (is.vector(data)) {	
      namesGenes <- data
      data <- data.frame(row.names = namesGenes)
      data[[inputTypeId]] <- namesGenes
    }
    else {
      namesGenes <- rownames(data)
      data[[inputTypeId]] <- rownames(data)
    }

    d<-getBM(attributes = unique(c(inputTypeId, outputTypeId, if (!("hgnc_symbol" %in% c(outputTypeId, inputTypeId))) "hgnc_symbol", if (!( "ensembl_gene_id" %in% c(outputTypeId, inputTypeId)))  "ensembl_gene_id")), filters = inputTypeId, values = namesGenes, mart = mart) 
  
    ## Unannotated Gene
    noAnnot <- namesGenes[!(namesGenes %in% d[[1]])]
    
    if (length(noAnnot) > 0) {
      nr <- nrow(d)
      d[(nr+1):(nr+length(noAnnot)),inputTypeId] = noAnnot
    }

    d <- merge(data, d, by = inputTypeId)

    ## Sort data.frame
    if (!is.null(sort.by)) {
      d <- d[order(d[[sort.by]], decreasing = TRUE),]
    }
    else{
      d <- d[unlist(sapply(namesGenes, function(x,col){which(col==x)},col=d[,inputTypeId])),]
    }

    if (!is.na(outfile)) {
      dhtml = d
      ##Add 2 columns :  link_geneCards = gene cards web address for the gene
      ##link_proteinAtlas = protein atlas web address for the gene
      dhtml$link_geneCards = paste("<p><a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=", dhtml$hgnc_symbol, "\">", dhtml$hgnc_symbol,"</a></p>",sep = "")
      dhtml$link_proteinAtlas = paste("<p><a href='http://www.proteinatlas.org/gene_info.php?ensembl_gene_id=", dhtml$ensembl_gene_id, "'>", dhtml$ensembl_gene_id,"</a>/</p/>",sep = "") 
      dhtml = dhtml[c(inputTypeId, setdiff(names(data),c(inputTypeId,outputTypeId)),"link_geneCards", "link_proteinAtlas", if (!is.null(outputTypeId)) outputTypeId[!(outputTypeId %in% inputTypeId)])]
      d = d[c(inputTypeId, setdiff(names(data),c(inputTypeId,outputTypeId)), if (!is.null(outputTypeId)) outputTypeId[!(outputTypeId %in% inputTypeId)])]
      
      write.table(d, file = paste(outfile,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE) 
      print(xtable(dhtml), type = "html", file = paste(outfile,sep=""), sanitize.text.function = force,include.rownames=FALSE)
    }
    return(d)
}

