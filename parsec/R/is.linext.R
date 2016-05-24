is.linext <-
function(order, z)
    all(linzeta(order) - z >= 0)
