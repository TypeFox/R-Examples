actu.matrice.adj <-
function(mat.adj, cple)
{
    cple <- as.character(cple)

    m.cple <- mat.adj[which(!rownames(mat.adj)%in%cple), cple]

    if (is.matrix(m.cple))
    {
        v.cple <- apply(m.cple, MARGIN=1, max)                
    }
    else
    {    
        v.cple <- max(m.cple)
    }

    nm <- rownames(mat.adj)
    nm1 <- rownames(mat.adj)[which(!rownames(mat.adj)%in%cple)]    
    nm2 <- as.character(as.integer(nm[length(nm)])+ 2)
    
    nm <- sapply(c(nm1, nm2), as.character)

    mat.adj <- mat.adj[which(!rownames(mat.adj)%in%cple), which(!colnames(mat.adj)%in%cple)]

    mat.adj <- cbind(mat.adj, v.cple)
    mat.adj <- rbind(mat.adj, t(c(v.cple,1)))
    dimnames(mat.adj) <- list(nm, nm)      
    
    return(mat.adj)
}
