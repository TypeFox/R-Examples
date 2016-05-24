calculer.var.exp <-
function(cple, data, part, ponderer)
{
    c1 <- cple[1]
    c2 <- cple[2]

    n1 <- length(which(part == c1))
    n2 <- length(which(part == c2))

    v <- var(data$z[which(part == c1  |  part == c2)])*(n1 + n2 - 1)
    
    if (v == 0)
    {
        return(0)
    } 
    
    if (ponderer == TRUE)
    {
        pond1 <- ponderer.classe(data[which(part == c1), c("x", "y")])
        pond2 <- ponderer.classe(data[which(part == c2), c("x", "y")])
    }
    else
    {
        pond1 <- 1
        pond2 <- 1
    }

    
    m1 <- mean(data$z[which(part == c1)])
    m2 <- mean(data$z[which(part == c2)])
    
    m <- mean(data$z[which(part == c1  |  part == c2)])

    return((pond1 * n1 * (m1 - m)**2 + pond2 * n2 * (m2 - m)**2)  /  v)
}
