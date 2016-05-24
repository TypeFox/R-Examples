fisher.multcomp <-
function(tab.cont,p.method="fdr") {
  if (is.matrix(tab.cont)) {tab.cont <- as.table(tab.cont)}
  call <- match.call()
  dname <- if(length(call$tab.cont)==1) {
    call$tab.cont
  } else {
    paste(call$tab.cont[1],"(",paste(call$tab.cont[-1],collapse=","),")",sep="")
  }
  if (!is.table(tab.cont)) {stop("'",dname,"' is not a \"table\" object",sep="")}
  if (is.null(colnames(tab.cont)) | any(0%in%nchar(colnames(tab.cont)))) {
    colnames(tab.cont) <- paste0("col",1:ncol(tab.cont))
  }
  if (is.null(rownames(tab.cont)) | any(0%in%nchar(rownames(tab.cont)))) {
    rownames(tab.cont) <- paste0("row",1:nrow(tab.cont))
  }
  colonnes <- combn(colnames(tab.cont),2)
  lignes <- combn(rownames(tab.cont),2)
  colonnes2 <- apply(colonnes,2,function(x) paste(x,collapse=":"))
  lignes2 <- apply(lignes,2,function(x) paste(x,collapse=":"))
  if (ncol(tab.cont)>2) {
    p.no.adjust <- matrix(0,nrow=ncol(lignes),ncol=ncol(colonnes),dimnames=list(lignes2,colonnes2))
    for (i in 1:ncol(colonnes)) {
	for (j in 1:ncol(lignes)) {
	  tab <- tab.cont[lignes[,j],colonnes[,i]]
	  p.no.adjust[j,i] <- fisher.test(tab)$p.value
	}
    }
    p.adjust <- matrix(p.adjust(p.no.adjust,method=p.method),nrow=nrow(p.no.adjust),ncol=ncol(p.no.adjust),
	dimnames=dimnames(p.no.adjust))
  } else {
    p.no.adjust <- matrix(NA,nrow=nrow(tab.cont),ncol=nrow(tab.cont),dimnames=list(rownames(tab.cont),rownames(tab.cont)))
    for (i in 1:ncol(lignes)) {
	tab.temp <- tab.cont[lignes[,i],]
	p.no.adjust[lignes[2,i],lignes[1,i]] <- fisher.test(tab.temp)$p.value
    }
    p.adjust <- matrix(p.adjust(p.no.adjust,method=p.method),nrow=nrow(p.no.adjust),ncol=ncol(p.no.adjust),
	dimnames = dimnames(p.no.adjust))
    p.adjust <- p.adjust[-1,-ncol(p.adjust)]
  }
  call <- match.call()
  dname <- if(length(call$tab.cont)==1) {call$tab.cont} else {paste(call$tab.cont[1],"(",paste(call$tab.cont[-1],collapse=","),")",sep="")}
  result<-list(method="Fisher's exact test for count data",data.name=dname,p.adjust.method=p.method,p.value=p.adjust)
  class(result) <- "RV.multcomp"
  return(result)
}
