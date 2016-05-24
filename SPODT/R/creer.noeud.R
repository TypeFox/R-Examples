creer.noeud <-
function(data, data.sp, perm, vql, vqt, ordre.vqt, bord,
                        eff.init, nv, id,
                        ponderer, greffer,
                        nv.max, min.pere, min.fils, var.exp.min)
{
    n <- nrow(data)    
    m <- nrow(bord)

    v <- var(data$z)
        
    
    # Verification des 3 premieres conditions d'arret :
    #--------------------------------------------------
    
    if (nrow(data) <= min.pere  |  nv >= nv.max  |  v == 0)
    {
        feuille <- new("f.spodt")
        feuille@id <- id
        feuille@n <- n
        feuille@m <- mean(data$z)
        centre <- apply(data[,c("x","y")], MARGIN=2, sum)
        feuille@G <- centre / n
        if (n == 1)
        {
            feuille@v=0
        }
        else
        {
            feuille@v <- var(data$z)
        }

        partition <- rep(0, eff.init)
        partition[as.integer(rownames(data))] <- id
        
        bord <- cbind.data.frame(bord, id=rep(id, m))
        
        return(list(noeud=feuille, part=partition, bord=bord))
    }


    # Recherche du meilleur decoupage :
    #----------------------------------
    
    dcp <- list(vic=0)
    
    if (nrow(data.sp) > 1)
    {
        dcp <- part.obl(data, data.sp, perm, min.fils, ponderer)
        type.dcp <- "sp"
    }

    if (length(vqt) >= 1)
    {
        dcp.vqt <- part.vqt(data, vqt, ordre.vqt, min.fils, ponderer)
        
        if (dcp$vic < dcp.vqt$vic)
        {
            dcp <- dcp.vqt
            type.dcp="vqt"
        }
    }
   
    if (length(vql) >= 1)
    {
   
        dcp.vql <- part.vql(data, vql, min.fils, ponderer)
        
        if (length(dcp.vql$vql.s) >= 1)
        {
            data <- data[,which(!(colnames(data) %in% dcp.vql$vql.s))]
            vql <- vql[which(!(vql %in% dcp.vql$vql.s))]
            
        }
        if (dcp$vic < dcp.vql$vic)
        {
            dcp <- dcp.vql
            type.dcp="vql"
        }        
    }

    # Verification des deux autres conditions d'arret :
    #--------------------------------------------------
    
    var.exp <- dcp$vic / ((n-1)*v)
    
    if (var.exp <= var.exp.min  |  length(which(dcp$part==-1)) <= min.fils  |  length(which(dcp$part==1)) <= min.fils)
    {
        feuille <- new("f.spodt")
        feuille@id <- id
        feuille@n <- n
        feuille@m <- mean(data$z)
        centre <- apply(data[,c("x","y")], MARGIN=2, sum)
        feuille@G <- centre / n
        if (n == 1)
        {
            feuille@v=0
        }
        else
        {
            feuille@v <- var(data$z)
        }
        
        partition <- rep(0, eff.init)
        partition[as.integer(rownames(data))] <- id

        bord <- cbind.data.frame(bord, id=rep(id, m))

        return(list(noeud=feuille, part=partition, bord=bord))
    }


    # Analyse des differents types de decoupages, actualisation des donnees et creation du noeud :
    #---------------------------------------------------------------------------------------------

    if (type.dcp == "sp")
    {
        coup <- inter.sgmts(bord, dcp$coeff)
        bord <- scinder(coup, bord, dcp$coeff)

        loc.g <- data.sp$loc[which(dcp$part == -1)]
        loc.d <- data.sp$loc[which(dcp$part == 1)]
        
        data.g <- data[data$loc %in% loc.g,]
        data.d <- data[data$loc %in% loc.d,]
        
        data.sp.g <- data.sp[which(data.sp$loc %in% loc.g == TRUE),]
        data.sp.d <- data.sp[which(data.sp$loc %in% loc.d == TRUE),]
        rownames(data.sp.g) <- rownames(data.sp[which(data.sp$loc %in% loc.g == TRUE),])
        rownames(data.sp.d) <- rownames(data.sp[which(data.sp$loc %in% loc.d == TRUE),])
        
        if (nrow(data.sp.g) > 1)
        {
            perm.g <- perm[which((perm$loc1 %in% data.sp.g$loc  &  perm$loc2 %in% data.sp.g$loc) == TRUE),]
        }
        else
        {
            perm.g <- NULL
        }
        if (nrow(data.sp.d) > 1)
        {
            perm.d <- perm[which((perm$loc1 %in% data.sp.d$loc  &  perm$loc2 %in% data.sp.d$loc) == TRUE),]
        }
        else
        {
            perm.d <- NULL
        }
        
        if (length(vqt) >= 1)
        {
            ordre.vqt.g <- part.ordre.vqt(ordre.vqt, rownames(data.g))
            ordre.vqt.d <- part.ordre.vqt(ordre.vqt, rownames(data.d))
        }
        
        noeud <- new("sp.spodt")
        noeud@coeff <- dcp$coeff
        noeud@int <- coup$X
    }
    else
    {
        bord <- list(g=bord, d=bord)
        
        data.g <- data[which(dcp$part == -1),]
        data.d <- data[which(dcp$part == 1),]
        
        data.sp.g <- part.data.sp(data, data.sp, dcp$part, -1)
        data.sp.d <- part.data.sp(data, data.sp, dcp$part, 1)
        rownames(data.sp.g) <- rownames(data.sp[which(data.sp$loc %in% data.g$loc == TRUE),])
        rownames(data.sp.d) <- rownames(data.sp[which(data.sp$loc %in% data.d$loc == TRUE),])
        
        if (nrow(data.sp.g) > 1)
        {
            perm.g <- perm[which((perm$loc1 %in% data.sp.g$loc  &  perm$loc2 %in% data.sp.g$loc) == TRUE),]
        }
        else
        {
            perm.g <- NULL
        }
        if (nrow(data.sp.d) > 1)
        {
            perm.d <- perm[which((perm$loc1 %in% data.sp.d$loc  &  perm$loc2 %in% data.sp.d$loc) == TRUE),]
        }
        else
        {
            perm.d <- NULL
        }
        if (length(vqt) >=1)
        {
            ordre.vqt.g <- part.ordre.vqt(ordre.vqt, rownames(data.g))
            ordre.vqt.d <- part.ordre.vqt(ordre.vqt, rownames(data.d))
        }
        if (type.dcp == "vql")
        {
            noeud <- new("vql.spodt")
            noeud@vrbl <- dcp$vrbl
            noeud@mod <- dcp$mod
        }
        if (type.dcp == "vqt")
        {
            noeud <- new("vqt.spodt")
            noeud@vrbl <- dcp$vrbl
            noeud@seuil <- dcp$seuil
        }
    }
    
    noeud@R2 <- var.exp
    noeud@n <- n
    noeud@m <- mean(data$z)
    noeud@v <- var(data$z)
    
    
    # Recherche de decoupages dans chacune des deux classes de la partition :
    #------------------------------------------------------------------------
    
    partg <- creer.noeud(data.g, data.sp.g, perm.g, vql, vqt, ordre.vqt.g, bord$g,
                         eff.init, nv+1, 2*id,
                         ponderer, greffer,
                         nv.max, min.pere, min.fils, var.exp.min)        
    partd <- creer.noeud(data.d, data.sp.d, perm.d, vql, vqt, ordre.vqt.d, bord$d,
                         eff.init, nv+1, 2*id+1,
                         ponderer, greffer,
                         nv.max, min.pere, min.fils, var.exp.min)
                         
    # Collecte des resultats des decoupages inferieurs :
    #---------------------------------------------------

    noeud@fg <- partg[[1]]
    noeud@fd <- partd[[1]]
                         
    partition <- partg[[2]] + partd[[2]]

    bord <- rbind.data.frame(partg[[3]], partd[[3]])
    
    
    # Renvoie les resultats à la fonction "spodt" :
    #----------------------------------------------
    
    return(list(noeud=noeud, part=partition, bord=bord))
}
