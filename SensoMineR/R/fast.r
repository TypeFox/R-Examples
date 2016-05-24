fast <- function (don, alpha = 0.05, sep.words = " ", word.min = 5, graph = TRUE,
    axes = c(1, 2), ncp = 5, B = 200, label.miss = NULL,ncp.boot=NULL)
{
    don = as.data.frame(don)
    I = nrow(don)
    J = ncol(don)
    for (i in 1:J) {
        don[, i] = as.factor(don[, i])
    }
    acm = MCA(don, graph = FALSE, ncp = ncp)
    if (graph) {
        plot.MCA(acm, choix = "ind", invisible = "var", axes = axes,new.plot=TRUE)
        plot.MCA(acm, choix = "ind", invisible = "ind", axes = axes,new.plot=TRUE)
        plot.MCA(acm, choix = "ind", axes = axes,new.plot=TRUE)
        plot.MCA(acm, choix = "var", axes = axes,new.plot=TRUE)
    }
    X = 0.05
    tab = as.matrix(acm$var$v.test)
    vtestmean1 = mean(tab[, axes[1]])
    vtestmean2 = mean(tab[, axes[2]])
    vtestsd1 = sd(tab[, axes[1]])
    vtestsd2 = sd(tab[, axes[2]])
    seuil1 = vtestmean1 + vtestsd1
    seuil2 = vtestmean2 + vtestsd2
    modext = which(abs(tab[, axes[1]]) >= seuil1 | abs(tab[,
        axes[2]]) >= seuil2)
    tab2 = which(abs(tab[, axes[1]]) < seuil1 & abs(tab[, axes[2]]) <
        seuil2)
    modmoy = tab2[sample(1:length(tab2), X * length(tab2))]
    mod_kept = cbind(t(modext), t(modmoy))
    res = acm
    res$var$coord = acm$var$coord[mod_kept, ]
    res$var$cos2 = acm$var$cos2[mod_kept, ]
    res$var$contrib = acm$var$contrib[mod_kept, ]
    res$var$v.test = acm$var$v.test[mod_kept, ]
    if (graph) {
        plot.MCA(res, choix = "ind", invisible = "ind", axes = axes,new.plot=TRUE)
    }
    lev = rep(NA, J)
    for (i in 1:J) {
        lev[i] = nlevels(don[, i])
    }
    lev2 = as.factor(lev)
    if (graph) {
        dev.new()
        plot(lev2, main = "Number of groups formed from sorting tasks")
    }
    nb_prod_grp = rep(NA, sum(lev))
    grp = 0
    for (i in 1:J) {
        nb_prod_grp[(grp + 1):(grp + lev[i])] = table(don[, i])
        grp = grp + lev[i]
    }
    nb_prod_grp2 = as.factor(nb_prod_grp)
    if (graph) {
        dev.new()
        plot(nb_prod_grp2, main = "Number of products per group")
    }
    tdc = tab.disjonctif(don)
    compte = tdc %*% t(tdc)
    ordre_prod = order(acm$ind$coord[, 1])
    compte2 = compte[ordre_prod, ]
    compte2 = compte2[, ordre_prod]
    chi2_t = function(x, y) {
        obs = table(x, y)
        chi = matrix(NA, length(levels(x)), length(levels(y)))
        for (i in 1:length(levels(x))) {
            for (j in 1:length(levels(y))) {
                chi[i, j] = (obs[i, j] - (sum(obs[i, ]) * sum(obs[,
                  j])/length(x)))^2/(sum(obs[i, ]) * sum(obs[,
                  j])/length(x))
            }
        }
        chi2 = sum(rowSums(chi))
        return(chi2)
    }
    cramer = function(x, y) {
        chi2 = chi2_t(x, y)
        n = length(x)
        p = length(levels(x))
        q = length(levels(y))
        m = min(p - 1, q - 1)
        V = sqrt(chi2/(n * m))
        return(V)
    }
    res = matrix(NA, J, J)
    for (i in 1:J) {
        for (j in i:J) {
            res[i, j] = res[j, i] = cramer(don[, i], don[, j])
        }
    }
    colnames(res) = rownames(res) = colnames(don)
    afc = CA(res, graph = FALSE)
    ord = order(afc$row$coord[, 1])
    res2 = res[ord, ]
    res2 = res2[, ord]
    out = matrix(NA, I, J)
    tdc = tab.disjonctif(don)
    gp = 0
    for (i in 1:J) {
        conc = cbind(c(1:lev[i]), acm$var$coord[(1 + gp):(gp +
            lev[i]), 1])
        o = order(conc[, 2])
        conc2 = cbind(conc[o, ], c(1:lev[i]))
        o2 = order(conc2[, 1])
        conc3 = conc2[o2, ]
        out[, i] = tdc[, (1 + gp):(gp + lev[i])] %*% conc3[,
            3]
        gp = gp + lev[i]
    }
    out = data.frame(out)
    for (i in 1:J) {
        out[, i] = as.factor(out[, i])
    }
    catego_num2 = out[ordre_prod, ]
    catego_num2 = catego_num2[, ord]
    rownames(catego_num2) = rownames(don)[ordre_prod]
    colnames(catego_num2) = colnames(don)[ord]
    
    if (graph) {
    #ellipses
    res.boot <- boot(don,method="sorting",ncp=ncp.boot,nbsim=B,level.conf=1-alpha,
    axes=axes)
    }

    texte = matrix(NA, (I * J), 3)
    texte = data.frame(texte)
    texte[, 1] = rep(rownames(don), J)
    texte[, 2] = rep(colnames(don), each = I)
    for (i in 1:J) {
        texte[((I * (i - 1)) + 1):(I * i), 3] = paste(don[, i])
    }
    restext = textual(texte, 3, 1, sep.word = sep.words, maj.in.min=TRUE)
    mod.suppr = paste("g", 1:99, sep = "")
    mod.suppr = intersect(colnames(restext$cont.table), mod.suppr)
    if (length(mod.suppr) != 0) {
        num.mod.suppr = which(colnames(restext$cont.table) %in%
            mod.suppr)
        restext$cont.table = restext$cont.table[, -num.mod.suppr]
        num.mod.suppr2 = which(rownames(restext$nb.words) %in%
            mod.suppr)
        restext$nb.words = restext$nb.words[-num.mod.suppr2,
            ]
    }
    nb_mot_diff = nrow(restext$nb.words)
    cat("Number of different words : ", nb_mot_diff, "\n")
    mots = rep(NA, sum(lev))
    grp = 0
    for (i in 1:J) {
        mots[(grp + 1):(grp + lev[i])] = levels(don[, i])
        grp = grp + lev[i]
    }
    mots_split = strsplit(mots, split = sep.words)
    nb_mots = rep(NA, length(mots_split))
    for (i in 1:length(mots_split)) {
        if (mots_split[[i]][1] %in% paste("G", 1:99, sep = "")) {
            nb_mots[i] = 0
        }
        else {
            nb_mots[i] = length(mots_split[[i]])
        }
    }
    nb_mots2 = as.factor(nb_mots)
    if (graph) {
        dev.new()
        plot(nb_mots2, main = "Number of words per group")
    }
    freq_min = which(apply(restext$cont.table, 2, sum) <= word.min)
    if (length(freq_min) != 0) {
        restext$cont.table = restext$cont.table[, -freq_min]
    }
    if (graph) {
        dev.new()
        compt.words=as.table(restext$nb.words[(min(20,nrow(restext$nb.words))):1,1])
        coordonnees = barplot(compt.words,las = 2, horiz = TRUE, names.arg=FALSE, main ="Occurrences of the most used words",
        sub="Number of products characterized by the words")
        text(x = 1, y = coordonnees, labels = rownames(restext$nb.words)[(min(20,nrow(restext$nb.words))):1], adj = 0)

    }

    caract_prod = descfreq(restext$cont.table)

    acm_call = list(X = acm$call$X, marge.col = acm$call$marge.col,
        marge.row = acm$call$marge.row, ncp = acm$call$ncp, quali = acm$call$quali,
        mca = acm,sep.words=sep.words)
    group_afm = list(coord = acm$var$eta2)

    res = list(eig = acm$eig, var = acm$var, ind = acm$ind, group = group_afm,
        acm = acm, call = acm_call, cooccur = compte2, reord = catego_num2,
        cramer = res2, textual = caract_prod, descriptor=data.frame(Occurrence=restext$nb.words[,1],row.names=rownames(restext$nb.words)))
    class(res) <- c("fast", "list ")
    return(res)
}
