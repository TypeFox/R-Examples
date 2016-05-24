near1 <-
structure(function (x, data) 
{
    nd <- length(data[, 1])
    distall <- rep(0, nd)
    for (i in 1:nd) {
        distall[i] <- distancia(x, data[i, ])
    }
    ind1 <- order(distall)[1]
    near1 <- data[ind1, ]
    near1
}, source = c("function(x, data)", "{", "#****************************************", 
"# Esta funcion encuentra la observacion en el", "#conjunto de datos data que esta mas cerca a", 
"# la observacion x requiere la funcion distancia", "#*********************************************", 
"nd <- length(data[, 1])", "distall <- rep(0, nd)", "for(i in 1:nd) {", 
"distall[i] <- distancia(x, data[i,  ])", "}", "#print(sort(distall))", 
"ind1 <- order(distall)[1]", "near1 <- data[ind1,  ]", "near1", 
"}"))
