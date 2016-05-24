matchFeatures <-
    function(CN=NULL, Exp=NULL, CN.by.feat=c("gene", "probe"),
                           Exp.by.feat=c("gene", "probe"), ref="hg19") {
    if (length(CN.by.feat)==2) {
      CN.by.feat <- "probe"
      }
    if (length(Exp.by.feat)==2) {
      Exp.by.feat <- "gene"
      }
    IntClustMemb <- NULL
    data(Map.CN, envir=environment())
    data(Map.Exp, envir=environment())
    data(train.CN, envir=environment())
    data(train.Exp, envir=environment())
    data(IntClustMemb, envir=environment())
    data(Map.All, envir=environment())
    if (!is.null(CN) & is.null(Exp)) {
        if (CN.by.feat=="probe") {
            Probes <- Map.All$Probe_ID
            Synonyms <- NULL
            train.CN <- as.matrix(train.CN)
        }
        else {
            Map.All <- Map.All[which(Map.All$Gene.Chosen=="YES"),]
            Probes <- as.character(Map.All$Gene_symbol)
            Synonyms <- as.character(Map.All$Synonyms_0)
            tmp <- rep(Probes, sapply(strsplit(Synonyms, ";"), length))
            Synonyms <- do.call("c", strsplit(Synonyms, ";"))
            Synonyms <- gsub(" ", "", Synonyms)
            names(Synonyms) <- tmp
            Synonyms <- Synonyms[which(!is.na(Synonyms))]
            train.CN <- as.matrix(train.CN[rownames(Map.All),])
            rownames(train.CN) <- Map.All$Gene_symbol
        }
        CN <- getCNfeatures(CN, Probes, Map.All, CN.by.feat, ref, Synonyms)
        all.na <- which(apply(CN, 1, function(x) mean(is.na(x)))<1)
        CN <- CN[all.na,]
        train.CN <- train.CN[all.na,]
        Map.CN <- Map.All[all.na,]
	Map.Exp <- NULL
        train.Exp <- NULL

        }
    if (!is.null(CN) & !is.null(Exp)) {

        if (Exp.by.feat=="probe") {
            Probes <- Map.Exp$Probe_ID
            Synonyms <- as.character(Map.Exp$Synonyms_0)
            tmp <- rep(Probes, sapply(strsplit(Synonyms, ";"), length))
            Synonyms <- do.call("c", strsplit(Synonyms, ";"))
            Synonyms <- gsub(" ", "", Synonyms)
            names(Synonyms) <- tmp
            Synonyms <- Synonyms[which(!is.na(Synonyms))]
            train.Exp <- train.Exp[Map.Exp$Probe_ID,]
        } else {
            Map.Exp <- Map.Exp[which(Map.Exp$Gene.Chosen=="YES"), ]
            Probes <- as.character(Map.Exp$Gene_symbol)
            Synonyms <- as.character(Map.Exp$Synonyms_0)
            tmp <- rep(Probes, sapply(strsplit(Synonyms, ";"), length))
            Synonyms <- do.call("c", strsplit(Synonyms, ";"))
            Synonyms <- gsub(" ", "", Synonyms)
            names(Synonyms) <- tmp
            Synonyms <- Synonyms[which(!is.na(Synonyms))]
            train.Exp <- train.Exp[Map.Exp$Probe_ID,]
            rownames(train.Exp) <- Map.Exp$Gene_symbol
        }
        if (CN.by.feat=="probe") {
            Probes.CN <- Map.CN$Probe_ID
            Synonyms.CN <- NULL
            train.CN <- as.matrix(train.CN[Map.CN$Probe_ID,])
        }
        else {
            tmp.map <- Map.All[,c('Probe_ID', 'Synonyms_0', 'Gene.Chosen')]
            Map.CN <- merge(Map.CN, tmp.map, sort=FALSE)
            rownames(Map.CN) <- Map.CN$Probe_ID
            Map.CN <- Map.CN[which(Map.CN$Gene.Chosen=="YES"), ]
            Probes.CN <- as.character(Map.CN$Gene_symbol)
            Synonyms.CN <- as.character(Map.CN$Synonyms_0)
            tmp <- rep(Probes.CN, sapply(strsplit(Synonyms.CN, ";"), length))
            Synonyms.CN <- do.call("c", strsplit(Synonyms.CN, ";"))
            Synonyms.CN <- gsub(" ", "", Synonyms.CN)
            names(Synonyms.CN) <- tmp
            Synonyms.CN <- Synonyms.CN[which(!is.na(Synonyms.CN))]
            train.CN <- as.matrix(train.CN[rownames(Map.CN),])
            rownames(train.CN) <- Map.CN$Gene_symbol
        }
        CN <- getCNfeatures(CN, Probes.CN, Map.CN, CN.by.feat, ref, Synonyms.CN)
        Exp <- getExpfeatures(Exp, Probes, Synonyms, Exp.by.feat)
	common.cols <- intersect(colnames(Exp), colnames(CN))
	Exp <- Exp[,which(colnames(Exp) %in% common.cols)]
	CN <- CN[,which(colnames(CN) %in% common.cols)]
        CN <- CN[,match(colnames(Exp), colnames(CN))]
        all.na <- which(apply(Exp, 1, function(x) mean(is.na(x)))<1)
        Exp <- Exp[all.na,]
        train.Exp <- train.Exp[all.na,]
        Map.Exp <- Map.Exp[all.na,]
        all.na <- which(apply(CN, 1, function(x) mean(is.na(x)))<1)
        CN <- CN[all.na,]
        train.CN <- train.CN[all.na,]
        Map.CN <- Map.CN[all.na,]
        train.Exp <- as.matrix(train.Exp)
        }
    if(!is.null(Exp) & is.null(CN)) {
        if (Exp.by.feat=="probe") {
            Probes <- Map.All$Probe_ID
            Synonyms <- as.character(Map.All$Synonyms_0)
            tmp <- rep(Probes, sapply(strsplit(Synonyms, ";"), length))
            Synonyms <- do.call("c", strsplit(Synonyms, ";"))
            Synonyms <- gsub(" ", "", Synonyms)
            names(Synonyms) <- tmp
            Synonyms <- Synonyms[which(!is.na(Synonyms))]
        } else {
            Map.All <- Map.All[which(Map.All$Gene.Chosen=="YES"), ]
            Probes <- as.character(Map.All$Gene_symbol)
            Synonyms <- as.character(Map.All$Synonyms_0)
            tmp <- rep(Probes, sapply(strsplit(Synonyms, ";"), length))
            Synonyms <- do.call("c", strsplit(Synonyms, ";"))
            Synonyms <- gsub(" ", "", Synonyms)
            names(Synonyms) <- tmp
            Synonyms <- Synonyms[which(!is.na(Synonyms))]
            train.Exp <- train.Exp[Map.All$Probe_ID,]
            rownames(train.Exp) <- Map.All$Gene_symbol
        }
        Exp <- getExpfeatures(Exp, Probes, Synonyms, Exp.by.feat)
        all.na <- which(apply(Exp, 1, function(x) mean(is.na(x)))<1)
        Exp <- Exp[all.na,]
        train.Exp <- train.Exp[all.na,]
        train.CN <- NULL
        train.Exp <- as.matrix(train.Exp)
        Map.CN <- NULL
        Map.Exp <- Map.All[all.na,]

    }
    res <- list(CN=CN, Exp=Exp, train.CN=train.CN, train.Exp=train.Exp, train.iC10=IntClustMemb,
                map.cn=Map.CN, map.exp=Map.Exp)
    attr(res, "CN.by.feat") <- CN.by.feat
    attr(res, "Exp.by.feat") <- Exp.by.feat
    attr(res, "ref") <- ref
    res
}
