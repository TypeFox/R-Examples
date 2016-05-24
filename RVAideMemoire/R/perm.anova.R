perm.anova <-
function(formula,nest.f2=c("fixed","random"),data,nperm=999) {
  if (missing(formula)||(length(formula)!=3)) {stop("missing or incorrect formula")}
  allnames <- all.names(formula)
  if (length(formula[[3]])==3) {
    if (formula[[3]][[1]]=="|") {formula[[3]][[1]] <- as.name("+")}
  }
  m <- match.call()
  m$formula <- formula
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$nest.f2 <- m$nperm <- NULL
  mf <- eval(m,parent.frame())
  resp <- mf[,1]
  fact1 <- mf[,2]
  variables <- names(mf)
  if (length(nest.f2)!=1) {nest.f2 <- "fixed"}
  data.name <- variables[1]
  tab <- NULL
  if (length(variables)==2) {
    tab <- .perm.anova.1way(resp,fact1,variables,nperm)$tab
  } else if (length(variables)==3) {
    fact2 <- mf[,3]
    if (allnames[3]=="+") {
	tab <- .perm.anova.2wayA(resp,fact1,fact2,variables,nperm)$tab
    } else if (allnames[3]=="*") {
	tab <- .perm.anova.2wayB(resp,fact1,fact2,variables,nperm)$tab
    } else if (allnames[3]=="/") {
	tab <- .perm.anova.2wayC(resp,fact1,fact2,nest.f2,variables,nperm)$tab
    } else if (allnames[3]=="|") {
	data.name <- paste(data.name,"\nBlock: ",variables[3],sep="")
	tab <- .perm.anova.2wayD(resp,fact1,fact2,variables,nperm)$tab
    }
  } else if (length(variables)==4) {
    if (allnames[3]!="|") {stop("incorrect 'formula'")}
    fact2 <- mf[,3]
    fact3 <- mf[,4]
    data.name <- paste(data.name,"\nBlock: ",variables[4],sep="")
    if (allnames[4]=="+") {
	tab <- .perm.anova.3wayA(resp,fact1,fact2,fact3,variables,nperm)$tab
    } else if (allnames[4]=="*") {
	tab <- .perm.anova.3wayB(resp,fact1,fact2,fact3,variables,nperm)$tab
    } else {
	stop("only additive and multiplicative models are permitted")
    }
  }
  h <- c("Permutational Analysis of Variance Table",paste("\nResponse: ",data.name,"\n",nperm," permutations",sep=""))
  structure(tab,heading=h,class=c("anova","data.frame"))
}

.perm.anova.1way <-
function(resp,fact1,variables,nperm) {
  anova.ref <- anova(lm(resp~fact1))
  F.ref <- anova.ref[1,"F value"]
  tab <- data.frame("Sum Sq"=anova.ref[,"Sum Sq"],"Df"=anova.ref[,"Df"],"Mean Sq"=anova.ref[,"Mean Sq"],
    "F value"=c(anova.ref[1,"F value"]," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
  rownames(tab) <- c(variables[2],"Residuals")
  F.perm <- numeric(nperm+1)
  F.perm[1] <- anova.ref[1,"F value"]
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    anova.perm <- anova(lm(sample(resp)~fact1))
    F.perm[i+1] <- anova.perm[1,"F value"]
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue <- length(which((F.perm+.Machine$double.eps/2) >= F.ref))/(nperm+1)
  tab[1,"Pr(>F)"] <- format(pvalue,digits=5,nsmall=5)
  return(list(tab=tab))
}

.perm.anova.2wayA <-
function(resp,fact1,fact2,variables,nperm) {
  if (any(diff(tapply(resp,list(fact2,fact1),length))!=0)) {stop("this function is not made for unbalanced design")}
  anova.ref <- anova(lm(resp~fact1+fact2))
  F1.ref <- anova.ref[1,"F value"]
  F2.ref <- anova.ref[2,"F value"]
  tab <- data.frame("Sum Sq"=anova.ref[,"Sum Sq"],"Df"=anova.ref[,"Df"],"Mean Sq"=anova.ref[,"Mean Sq"],
    "F value"=c(anova.ref[1:2,"F value"]," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
  rownames(tab) <- c(variables[2],variables[3],"Residuals")
  F1.perm <- numeric(nperm+1)
  F2.perm <- numeric(nperm+1)
  F1.perm[1] <- F1.ref
  F2.perm[1] <- F2.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    anova.perm <- anova(lm(sample(resp)~fact1+fact2))
    F1.perm[i+1] <- anova.perm[1,"F value"]
    F2.perm[i+1] <- anova.perm[2,"F value"]
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue1 <- length(which((F1.perm+.Machine$double.eps/2) >= F1.ref))/(nperm+1)
  pvalue2 <- length(which((F2.perm+.Machine$double.eps/2) >= F2.ref))/(nperm+1)
  tab[1:2,"Pr(>F)"] <- format(c(pvalue1,pvalue2),digits=5,nsmall=5)
  return(list(tab=tab))
}

.perm.anova.2wayB <-
function(resp,fact1,fact2,variables,nperm) {
  if (any(diff(tapply(resp,list(fact2,fact1),length))!=0)) {stop("this function is not made for unbalanced design")}
  if (tapply(resp,list(fact2,fact1),length)[1,1]==1) {
    warning("no repetition of ",variables[2],":",variables[3],", ANOVA without interaction")
    tab <- .perm.anova.2wayA(resp,fact1,fact2,variables,nperm)$tab
  } else {
    anova.ref <- anova(lm(resp~fact1*fact2))
    F1.ref <- anova.ref[1,"F value"]
    F2.ref <- anova.ref[2,"F value"]
    F3.ref <- anova.ref[3,"F value"]
    tab <- data.frame("Sum Sq"=anova.ref[,"Sum Sq"],"Df"=anova.ref[,"Df"],"Mean Sq"=anova.ref[,"Mean Sq"],
	"F value"=c(anova.ref[1:3,"F value"]," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
    rownames(tab) <- c(variables[2],variables[3],paste(variables[2],variables[3],sep=":"),"Residuals")
    F1.perm <- numeric(nperm+1)
    F2.perm <- numeric(nperm+1)
    F3.perm <- numeric(nperm+1)
    F1.perm[1] <- F1.ref
    F2.perm[1] <- F2.ref
    F3.perm[1] <- F3.ref
    pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    for (i in 1:nperm) {
	anova.perm <- anova(lm(sample(resp)~fact1*fact2))
	F1.perm[i+1] <- anova.perm[1,"F value"]
	F2.perm[i+1] <- anova.perm[2,"F value"]
	F3.perm[i+1] <- anova.perm[3,"F value"]
	setTxtProgressBar(pb,round(i*100/nperm,0))
    }
    cat("\n")
    pvalue1 <- length(which((F1.perm+.Machine$double.eps/2) >= F1.ref))/(nperm+1)
    pvalue2 <- length(which((F2.perm+.Machine$double.eps/2) >= F2.ref))/(nperm+1)
    pvalue3 <- length(which((F3.perm+.Machine$double.eps/2) >= F3.ref))/(nperm+1)
    tab[1:3,"Pr(>F)"] <- format(c(pvalue1,pvalue2,pvalue3),digits=5,nsmall=5)
  }
  return(list(tab=tab))
}

.perm.anova.2wayC <-
function(resp,fact1,fact2,nest.f2,variables,nperm) {
  if (any(diff(tapply(resp,fact2,length))!=0)) {stop("this function is not made for unbalanced design")}
  anova.ref <- anova(lm(resp~fact1/fact2))
  if (nest.f2=="random") {anova.ref[1,"F value"] <- anova.ref[1,"Mean Sq"]/anova.ref[2,"Mean Sq"]}
  F1.ref <- anova.ref[1,"F value"]
  F2.ref <- anova.ref[2,"F value"]
  tab <- data.frame("Sum Sq"=anova.ref[,"Sum Sq"],"Df"=anova.ref[,"Df"],"Mean Sq"=anova.ref[,"Mean Sq"],
    "F value"=c(anova.ref[1:2,"F value"]," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
  rownames(tab) <- c(variables[2],paste(variables[2],variables[3],sep=":"),"Residuals")
  F1.perm <- numeric(nperm+1)
  F2.perm <- numeric(nperm+1)
  F1.perm[1] <- F1.ref
  F2.perm[1] <- F2.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    anova.perm1 <- anova(lm(sample(resp)~fact1/fact2))
    if (nest.f2=="random") {
	F1.perm[i+1] <- anova.perm1[1,"Mean Sq"]/anova.perm1[2,"Mean Sq"]
    } else {
	F1.perm[i+1] <- anova.perm1[1,"F value"]
    }
    ordre <- order(fact1)
    ordre.new <- unlist(tapply(ordre,factor(rep(1:nlevels(fact1),each=length(fact1)/nlevels(fact1))),sample))
    anova.perm2 <- anova(lm(resp[ordre.new]~fact1/fact2))
    F2.perm[i+1] <- anova.perm2[2,"F value"]
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue1 <- length(which((F1.perm+.Machine$double.eps/2) >= F1.ref))/(nperm+1)
  pvalue2 <- length(which((F2.perm+.Machine$double.eps/2) >= F2.ref))/(nperm+1)
  tab[1:2,"Pr(>F)"] <- format(c(pvalue1,pvalue2),digits=5,nsmall=5)
  if (tapply(resp,fact2,length)[1]==1) {
    warning("only 1 observation per level of '",variables[3],"', permutation of '",variables[2],"' only")
    tab[2,"Pr(>F)"] <- "NA"
  }
  if (nest.f2=="random") {tab <- tab[-2,]}
  return(list(tab=tab))
}

.perm.anova.2wayD <-
function(resp,fact1,fact2,variables,nperm) {
  if (any(diff(tapply(resp,fact2,length))!=0)) {stop("this function is not made for unbalanced design")}
  if (tapply(resp,list(fact2,fact1),length)[1,1]==1) {stop("no repetition of ",variables[2],":",variables[3]," -> no interaction")}
  anova.ref <- anova(lm(resp~fact1*fact2))
  anova.ref[1,"F value"] <- anova.ref[1,"Mean Sq"]/anova.ref[3,"Mean Sq"]
  F1.ref <- anova.ref[1,"F value"]
  tab <- data.frame("Sum Sq"=anova.ref[c(1,4),"Sum Sq"],"Df"=anova.ref[c(1,4),"Df"],"Mean Sq"=anova.ref[c(1,4),"Mean Sq"],
    "F value"=c(anova.ref[1,"F value"]," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
  rownames(tab) <- c(variables[2],"Residuals")
  F1.perm <- numeric(nperm+1)
  F1.perm[1] <- F1.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    anova.perm <- anova(lm(sample(resp)~fact1*fact2))
    F1.perm[i+1] <- anova.perm[1,"Mean Sq"]/anova.perm[3,"Mean Sq"]
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue <- length(which((F1.perm+.Machine$double.eps/2) >= F1.ref))/(nperm+1)
  tab[1,"Pr(>F)"] <- format(pvalue,digits=5,nsmall=5)
  return(list(tab=tab))
}

.perm.anova.3wayA <-
function(resp,fact1,fact2,fact3,variables,nperm) {
  if (any(diff(table(fact1,fact2,fact3))!=0)) {stop("this function is not made for unbalanced design")}
  anova.ref <- anova(lm(resp~fact1*fact3+fact2*fact3))
  MSres <- sum(anova.ref[4:5,"Sum Sq"])/sum(anova.ref[4:5,"Df"]) 
  F1.ref <- anova.ref[1,"Mean Sq"]/MSres
  F2.ref <- anova.ref[3,"Mean Sq"]/MSres
  tab <- data.frame("Sum Sq"=anova.ref[c(1,3,6),"Sum Sq"],"Df"=anova.ref[c(1,3,6),"Df"],"Mean Sq"=anova.ref[c(1,3,6),"Mean Sq"],
    "F value"=c(c(F1.ref,F2.ref)," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
  rownames(tab) <- c(variables[2],variables[3],"Residuals")
  F1.perm <- numeric(nperm+1)
  F2.perm <- numeric(nperm+1)
  F1.perm[1] <- F1.ref
  F2.perm[1] <- F2.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    anova.perm <- anova(lm(sample(resp)~fact1*fact3+fact2*fact3))
    MSres.perm <- sum(anova.perm[4:5,"Sum Sq"])/sum(anova.perm[4:5,"Df"]) 
    F1.perm[i+1] <- anova.perm[1,"Mean Sq"]/MSres.perm
    F2.perm[i+1] <- anova.perm[3,"Mean Sq"]/MSres.perm
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue1 <- length(which((F1.perm+.Machine$double.eps/2) >= F1.ref))/(nperm+1)
  pvalue2 <- length(which((F2.perm+.Machine$double.eps/2) >= F2.ref))/(nperm+1)
  tab[1:2,"Pr(>F)"] <- format(c(pvalue1,pvalue2),digits=5,nsmall=5)
  return(list(tab=tab))
}

.perm.anova.3wayB <-
function(resp,fact1,fact2,fact3,variables,nperm) {
  if (any(diff(table(fact1,fact2,fact3))!=0)) {stop("this function is not made for unbalanced design")}
  if (table(fact1,fact2,fact3)[1,1,1]==1) {stop("no repetition of ",variables[2],":",variables[3],":",variables[4]," -> no interaction")}
  anova.ref <- anova(lm(resp~fact1*fact2*fact3))
  MSres <- sum(anova.ref[5:7,"Sum Sq"])/sum(anova.ref[5:7,"Df"]) 
  F1.ref <- anova.ref[1,"Mean Sq"]/MSres
  F2.ref <- anova.ref[2,"Mean Sq"]/MSres
  F1F2.ref <- anova.ref[4,"Mean Sq"]/MSres
  tab <- data.frame("Sum Sq"=anova.ref[c(1,2,4,8),"Sum Sq"],"Df"=anova.ref[c(1,2,4,8),"Df"],"Mean Sq"=anova.ref[c(1,2,4,8),"Mean Sq"],
    "F value"=c(c(F1.ref,F2.ref,F1F2.ref)," "),"Pr(>F)"=NA,stringsAsFactors=FALSE,check.names=FALSE)
  rownames(tab) <- c(variables[2],variables[3],paste(variables[2],":",variables[3],sep=""),"Residuals")
  F1.perm <- numeric(nperm+1)
  F2.perm <- numeric(nperm+1)
  F1F2.perm <- numeric(nperm+1)
  F1.perm[1] <- F1.ref
  F2.perm[1] <- F2.ref
  F1F2.perm[1] <- F1F2.ref
  pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  for (i in 1:nperm) {
    anova.perm <- anova(lm(sample(resp)~fact1*fact2*fact3))
    MSres.perm <- sum(anova.perm[5:7,"Sum Sq"])/sum(anova.perm[5:7,"Df"]) 
    F1.perm[i+1] <- anova.perm[1,"Mean Sq"]/MSres.perm
    F2.perm[i+1] <- anova.perm[2,"Mean Sq"]/MSres.perm
    F1F2.perm[i+1] <- anova.perm[4,"Mean Sq"]/MSres.perm
    setTxtProgressBar(pb,round(i*100/nperm,0))
  }
  cat("\n")
  pvalue1 <- length(which((F1.perm+.Machine$double.eps/2) >= F1.ref))/(nperm+1)
  pvalue2 <- length(which((F2.perm+.Machine$double.eps/2) >= F2.ref))/(nperm+1)
  pvalue3 <- length(which((F1F2.perm+.Machine$double.eps/2) >= F1F2.ref))/(nperm+1)
  tab[1:3,"Pr(>F)"] <- format(c(pvalue1,pvalue2,pvalue3),digits=5,nsmall=5)
  return(list(tab=tab))
}
