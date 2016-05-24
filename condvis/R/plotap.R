plotap <-
function (k, pathindex = 1)
{
    rsk <- rowSums(k)
    plot(rsk, type = "l", xlab = "Path index", ylab = "sum of k")  
    abline(v = pathindex)
    structure(list(k = k, rsk = rsk, pathindex = pathindex, device = dev.cur(), 
        screen = screen(), mar = par()$mar, usr = par()$usr), class = "ap")
}

update.ap <- 
function (object, pathindex = NULL, ...)
{
    if (dev.cur() != object$device)
        dev.set(object$device)
    screen(n = object$screen, new = FALSE)
    par(mar = object$mar)
    par(usr = object$usr)    
    if (!is.null(pathindex)){
        abline(v = object$pathindex, col = "white")
        refreshindex <- max((object$pathindex - 5), 1):min((object$pathindex + 
            5), nrow(object$k))
        points(refreshindex, object$rsk[refreshindex], type = "l")
        abline(v = pathindex)
        box()
        object$pathindex <- pathindex
    }
    object
}