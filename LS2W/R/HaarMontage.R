HaarMontage <-
function (n=128,direction = "diagonal",sd=1) 

{
    m <- log2(n)-1

    if (direction == "diagonal") {
        x1 <- Haar2MA.diag(n = 2^m, order = 1,sd=sd)
        x2 <- Haar2MA.diag(n = 2^m, order = 2,sd=sd)
        x3 <- Haar2MA.diag(n = 2^m, order = 3,sd=sd)
        x4 <- Haar2MA.diag(n = 2^m, order = 4,sd=sd)
        temp1 <- cbind(x1, x2)
        temp2 <- cbind(x3, x4)
        Monty <- rbind(temp1, temp2)
        return(Monty)
    }
    if (direction == "vertical") {
        x1 <- Haar2MA.vert(n = 2^m, order = 1,sd=sd)
        x2 <- Haar2MA.vert(n = 2^m, order = 2,sd=sd)
        x3 <- Haar2MA.vert(n = 2^m, order = 3,sd=sd)
        x4 <- Haar2MA.vert(n = 2^m, order = 4,sd=sd)
        temp1 <- cbind(x1, x2)
        temp2 <- cbind(x3, x4)
        Monty <- rbind(temp1, temp2)
        return(Monty)
    }
    if (direction == "horizontal") {
        x1 <- Haar2MA.horiz(n = 2^m, order = 1,sd=sd)
        x2 <- Haar2MA.horiz(n = 2^m, order = 2,sd=sd)
        x3 <- Haar2MA.horiz(n = 2^m, order = 3,sd=sd)
        x4 <- Haar2MA.horiz(n = 2^m, order = 4,sd=sd)
        temp1 <- cbind(x1, x2)
        temp2 <- cbind(x3, x4)
        Monty <- rbind(temp1, temp2)
        return(Monty)
    }
    else stop("\nDirection can only take the values horizontal, vertical or diagonal!\n")
}

