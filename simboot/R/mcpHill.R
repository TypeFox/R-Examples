mcpHill <- function(dataf, fact, align=FALSE, block, boots=5000, udmat=FALSE, usermat, mattype="Dunnett", dunbase=1, qval=seq(-1,3), opt="two.sided")
{
    # require(boot, quietly=TRUE)
    # require(multcomp, quietly=TRUE)
    if (!is.data.frame(dataf)) {stop("dataf must be an object of class 'data.frame'.")}
    if (length(fact) != dim(dataf)[1]) {stop("The length of fact must equal the number of rows in dataf.")}
    if (length(levels(fact)) <= 1) {stop("The factor variable fact should have at least 2 levels to be compared.")}
    ni <- as.vector(summary(fact))
                                        # Kontrastmatrix
    if(udmat==FALSE){
        if(mattype=="Dunnett")
        {cmat <- contrMat(ni, type="Dunnett", base=dunbase)}
        else
        {cmat <- contrMat(ni, type=mattype)}
    }
    else
    {cmat <- usermat}
                                        # Hillfunktion
    hillq <- function(dat,q)
    {apply(dat, 1, function(da)
       {p_i <- (da/sum(da))[(da/sum(da))>0]
        if(q==1)
        {return(exp(-1*sum(p_i*log(p_i))))}
        else
        {return((sum(p_i^q))^(1/(1-q)))}})}
                                        # Tabelle mit transformierten Indices und Gruppenzugehoerigkeit
    hilltab <- function(bui){sapply(bui, function(bae){hillq(dat=dataf, q=bae)})}
    tabtab <- hilltab(qval)
    group <- factor(fact)
    tabelle <- cbind(tabtab, group)
                                        # Blockalignment
    if(align==TRUE){
        alignfunc <- function(aliblock)
        {tfit <- lm(tabelle[, c(1:length(qval))] ~ aliblock)
         tabelleneu <- tabelle[, c(1:length(qval))]- predict(tfit)}
        tabelle2 <- alignfunc(block)
    }
    else
    {tabelle2 <- tabelle}
                                        # Tabelle mit Epsilons
    funcA <- function(f)
    {fit <- lm(tabelle2[, c(1:length(qval))] ~ f-1)
     epsilon <- residuals(fit)}
    epstabelle <- funcA(f=group)
                                        # Teststatistiken
    tstatshort <- function(mytab, i, f, cmat, ni)
    {mytab <- mytab[i]
     ni <- ni
     FIT <- lm(mytab ~ f-1)
     mi <- coefficients(FIT)
     res <- residuals(FIT)
     varpool <- sum(res^2)/(sum(ni)-length(mi))
     estC <- (cmat %*% mi)
     varC <- (cmat^2) %*% (varpool/ni)
     ti <- estC/sqrt(varC)
     return(ti)}

    if(length(qval)==1){
        funcfunc <- function(abc, i, f, cmat, ni){tstatshort(mytab=abc, f=f, i=i, cmat=cmat, ni=ni)}
    }
    else{
        funcfunc <- function(abc, i, f, cmat, ni){apply(abc, 2, FUN=function(xx){tstatshort(mytab=xx, f=f, i=i, cmat=cmat, ni=ni)})}
    }

                                        # Westfall-Young / Bootstrap
    wyboot <- boot(epstabelle, funcfunc, R=boots, stype="i", f=group, cmat=cmat, ni=ni)
                                        # p-Werte
    laenge <- 1:(length(qval) * dim(cmat)[1])
    tact <- as.vector(funcfunc(abc=tabelle2[,c(1:length(qval))], i=1:sum(ni), f=group, cmat=cmat, ni=ni))
                                        # two-sided
    bothfunc <- function(bobo){apply(bobo, 1, function(huhu)max(abs(huhu)))}
    tboth <- bothfunc(wyboot$t)
    bothpval <- function(k){
        tact <- tact[k]
        sapply(k, function(ka)sum(tboth > abs(tact[ka]))/boots)}
                                        # greater
    maxfunc <- function(mama){apply(mama, 1, max)}
    tmax <- maxfunc(wyboot$t)
    grpval <- function(k){
        tact <- tact[k]
        sapply(k, function(ka)sum(tmax > tact[ka])/boots)}
                                        # less
    minfunc <- function(mimi){apply(mimi, 1, min)}
    tmin <- minfunc(wyboot$t)
    lepval <- function(k){
        tact <- tact[k]
        sapply(k, function(ka)sum(tmin < tact[ka])/boots)}
                                        # switch options
    allqs <- rep(qval, each=dim(cmat)[1])
    names(allqs) <- rep(rownames(cmat), times=length(qval))
    switch(opt,
           two.sided = cbind("q"=allqs, "p-value"=bothpval(laenge)),
           greater = cbind("q"=allqs, "p-value"=grpval(laenge)),
           less = cbind("q"=allqs, "p-value"=lepval(laenge)))
}
