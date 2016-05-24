alphablend <-
function (col, alpha) 
{
    colrgb <- col2rgb(col)/255
    rgb(colrgb[1, ], colrgb[2, ], colrgb[3, ], alpha)
}
