distancia1 <-
function (x, y) 
{#Manhattan distance
    if (class(y) == "matrix") {
#        print(class(x))
#        print(dim(t(y)))
        distancia = drop(colSums(abs(x - t(y))))
        distancia = t(distancia)
    }
    else distancia = sum(abs(x - y))
    distancia
}
