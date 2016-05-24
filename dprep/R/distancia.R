distancia <-
function (x, y) 
{#Find out the row in y that is closer to x using eucldean distance 
    if (class(y) == "matrix") {
        distancia = drop(sqrt(colSums((x - t(y))^2)))
        distancia = t(distancia)
    }
    else distancia = sqrt(sum((x - y)^2))
    distancia
}
