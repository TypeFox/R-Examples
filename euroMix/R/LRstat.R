LRstat=function (ped_claim,ped_true, ids, alleles, afreq = NULL, known_genotypes = list(), 
loop_breakers = NULL, Xchrom = F, plot = T) {
    if (inherits(ped_claim, "linkdat")) 
        ped_claim = list(ped_claim)
    if (inherits(ped_true, "linkdat")) 
        ped_true = list(ped_true)
    ids_claim = lapply(ped_claim, function(x) ids[ids %in% x$orig.ids])
    ids_true = lapply(ped_true, function(x) ids[ids %in% x$orig.ids])
    loops_claim = lapply(ped_claim, function(x) {
        lb = x$orig.ids[x$orig.ids %in% loop_breakers]
        if (length(lb) == 0) 
            lb = NULL
        lb
    })
    loops_true = lapply(ped_true, function(x) {
        lb = x$orig.ids[x$orig.ids %in% loop_breakers]
        if (length(lb) == 0) 
            lb = NULL
        lb
    })
    N_claim = length(ped_claim)
    N_true = length(ped_true)
    N = N_claim + N_true
    if (length(alleles) == 1) 
        alleles = 1:alleles
    if (Xchrom) 
        chrom = 23
    else chrom = NA
    partial_claim = lapply(1:N_claim, function(i) {
        x = ped_claim[[i]]
        m = marker(x, alleles = alleles, afreq = afreq, chrom = chrom)
        for (tup in known_genotypes) if (tup[1] %in% x$orig.ids) 
            m = modifyMarker(x, m, ids = tup[1], genotype = tup[2:3])
        m
    })
    partial_true = lapply(1:N_true, function(i) {
        x = ped_true[[i]]
        m = marker(x, alleles = alleles, afreq = afreq, chrom = chrom)
        for (tup in known_genotypes) if (tup[1] %in% x$orig.ids) 
            m = modifyMarker(x, m, ids = tup[1], genotype = tup[2:3])
        m
    })
    if (isTRUE(plot) || plot == "plot_only") {
        op = par(oma = c(0, 0, 3, 0), xpd = NA)
        widths = ifelse(sapply(c(ped_claim, ped_true), inherits, 
            what = "singleton"), 1, 2)
        claim_ratio = sum(widths[1:N_claim])/sum(widths)
        layout(rbind(1:N), widths = widths)
        has_genotypes = length(known_genotypes) > 0
        for (i in 1:N) {
            if (i <= N_claim) {
                x = ped_claim[[i]]
                avail = ids_claim[[i]]
                mm = if (has_genotypes) 
                  partial_claim[[i]]
                else NULL
            }
            else {
                x = ped_true[[i - N_claim]]
                avail = ids_true[[i - N_claim]]
                mm = if (has_genotypes) 
                  partial_true[[i - N_claim]]
                else NULL
            }
            cols = ifelse(x$orig.ids %in% avail, 2, 1)
            plot(x, marker = mm, col = cols, margin = c(2, 4, 
                2, 4), title = "")
        }
        mtext("HP", outer = TRUE, at = claim_ratio/2)
        mtext("HD", outer = TRUE, at = 0.5 + claim_ratio/2)
        rect(grconvertX(0.02, from = "ndc"), grconvertY(0.02, 
            from = "ndc"), grconvertX(claim_ratio - 0.02, from = "ndc"), 
            grconvertY(0.98, from = "ndc"))
        rect(grconvertX(claim_ratio + 0.02, from = "ndc"), grconvertY(0.02, 
            from = "ndc"), grconvertX(0.98, from = "ndc"), grconvertY(0.98, 
            from = "ndc"))
        par(op)
        if (plot == "plot_only") 
            return()
    }
    p.g = Reduce('%o%', lapply(1:N_true, function(i) 
     oneMarkerDistribution(ped_true[[i]], ids = ids_true[[i]], partialmarker = partial_true[[i]], loop_breakers = loops_true[[i]], verbose = F)))
    p.c = Reduce('%o%', lapply(1:N_claim, function(i) 
      oneMarkerDistribution(ped_claim[[i]], ids = ids_claim[[i]], partialmarker = partial_claim[[i]], loop_breakers = loops_claim[[i]], verbose = F)))  
    I.g = Reduce('%o%', lapply(1:N_claim, function(i) 
      oneMarkerDistribution(ped_claim[[i]], ids = ids_claim[[i]], partialmarker = partial_claim[[i]], loop_breakers = loops_claim[[i]], verbose = F) == 0))
    PE = sum(p.g*I.g)
  
   LR= p.c/p.g #Pr(data|Hp)/PR(data|Hd), claim=HP;true=HD
   LRv = as.vector(LR)
   ord1=!is.na(LRv)
   LRv=LRv[ord1]
   ord = order(LRv)
   LRv = LRv[ord]
   p1 = as.vector(p.c)
   p1=p1[ord1]
   p1=p1[ord]
   Pr.LR.claim = sapply(split(p1, LRv), sum)
   p2 = as.vector(p.g)
   p2=p2[ord1]
   p2=p2[ord]
   Pr.LR.true = sapply(split(p2, LRv), sum)
   d1=as.double(names(Pr.LR.claim))
   ind=(1:length(d1))[d1>0]
   d1=d1[ind]
   z1 = sum(d1*Pr.LR.claim[ind])
   z2 = sum((1/d1)*Pr.LR.claim[ind])
   dum1=as.double(names(Pr.LR.claim[ind]))
   v1=sum(dum1^2*Pr.LR.claim[ind])-z1^2
   dum1=as.double(names(Pr.LR.true))
   z3 = sum(dum1*Pr.LR.true)
   v2=sum(dum1^2*Pr.LR.true)-z3^3
   z4 = sum((1/as.double(names(Pr.LR.true)))*Pr.LR.true)
   mainOutput=c(E.LR.HP=z1,E.LR.HP_1=z2,v1=v1,E.LR.HD=z3,E.LR.HD_1=z4,v2=v2,RMNE=1-PE)
   names(mainOutput)=c("E(LR(HP))","E(LR(HP)^(-1))","var(LR(HP))","E(LR(HD))","E(LR(HD)^(-1))","var(LR(HD))","RMNE")
   extraOutput=list(Pr.HD=p.g, PR.HP = p.c)
   LRdist=data.frame(LR=unique(round(LRv,10)),Pr.LR.HP = Pr.LR.claim, F.HP=cumsum(Pr.LR.claim),
   Pr.LR.HD = Pr.LR.true,F.HD=cumsum(Pr.LR.true))
   rownames(LRdist)=NULL
   list(main=mainOutput,extra=extraOutput,LRdist=LRdist)
}

