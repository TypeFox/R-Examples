# to get the segment extremities of the object sp

cpt <- 1
data.coord <- data.frame(matrix(, ncol = 4))
colnames(data.coord) <- c("X1", "Y1", "X2", "Y2")

recursif.carte.sp <- function(noeud){
    if (class(noeud) != "f.spodt")
    {
        if (class(noeud) == "sp.spodt")
        {
            X1 <- noeud@int[1]
            X2 <- noeud@int[2]
            Y1 <- noeud@coeff[1] * X1 + noeud@coeff[2]
            Y2 <- noeud@coeff[1] * X2 + noeud@coeff[2]

            data.coord[cpt, 1:2 ]   <- c(X1,Y1)
            data.coord[cpt, 3:4] <- c(X2,Y2)

            cpt <- cpt + 1
        }
        data.coord <- rbind(data.coord, recursif.carte.sp(noeud@fg))
        data.coord <- rbind(data.coord, recursif.carte.sp(noeud@fd))
    }
    return(data.coord)
}


