sp1 <-
function(data,ql.fact=NULL, qt.fact=NULL, weight=FALSE, graft=FALSE,
                  level.max=5, min.parent=10, min.child=5, rtwo.min=0.001)
{   
    ponderer<-weight
    greffer<-graft
    nv.max<-level.max
    min.pere<-min.parent
    min.fils<-min.child
    var.exp.min<-rtwo.min
    vql<-ql.fact
    vqt<-qt.fact
    colnames(data)[1:4] <- c("loc", "x", "y", "z")
    n <- nrow(data)
    rownames(data) <- 1:n

    data.sp <- prep.data.sp(data)

    perm <- calculer.pentes(data, data.sp)
    bord <- prep.bord(data)
    
    if (length(vqt) > 1)
    {
        ordre.vqt <- apply(data[,vqt], MARGIN=2, order)
    }
    else
    {
        if (length(vqt) == 1)
        {
            ordre.vqt <- order(data[,vqt])
        }
        else
        {
            ordre.vqt <- NULL
        }
    }


    res <- creer.noeud(data, data.sp, perm, vql, vqt, ordre.vqt, bord,
                       n, 1, 1,
                       ponderer, greffer,
                       nv.max, min.pere, min.fils, var.exp.min)

    if (class(res$noeud) == "f.spodt")
    {
        warning("Root is a leaf")
    }

    arbre <- new("spodt")
    arbre@racine <- res$noeud
    arbre@partition <- res$part

    res <- classe.adj(res$bord)

    arbre@adj <- res$adj    
    arbre@brd <- as.matrix(res$bord)

    if (greffer)
    {
        if (is.numeric(greffer))
        {
            var.exp.min <- greffer
        }
        
        bord <- transformer.bord(res$bord)

        
        res <- realiser.greffe(data, bord, arbre@partition, ponderer, var.exp.min, arbre@adj)
        if (! is.null(res$bord)){
            arbre@brd <- as.matrix(res$bord)
        }
        if (!res$grf)
        {
            arbre@R2 <- R2.global(data$z, arbre@partition)
            return(arbre)
        }

        arbre@adj <- unlist(res$adj)
        arbre@R2 <- R2.global(data$z, arbre@partition)

        cl.grf <- matrix(as.integer(res$cl.grf ), ncol=3)
        dimnames(cl.grf) <- list(1:nrow(cl.grf), c("c1", "c2", "union(c1,c2)"))
        arbre@cl.grf <- cl.grf

        sgmts.grf <- matrix(res$sgmts.grf, ncol=4)
        dimnames(sgmts.grf) <- list(1:nrow(sgmts.grf), c("X1", "Y1", "X2", "Y2"))
        arbre@sgmts.grf <- sgmts.grf

        arbre@partition <- unlist(res$part)
        
        arbre@R2 <- R2.global(data$z, arbre@partition)
    }

    return(arbre)
}
