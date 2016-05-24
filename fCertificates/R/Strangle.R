Strangle <- function (S, X1, X2, Time, r, r_d, sigma, ratio = 1) 
{
    if (Time == 0) {
        (pmax(S - X2, 0) + pmax(X1 - S, 0)) * ratio
    }
    else {
        call = GBSOption(TypeFlag = "c", S, X2, Time, r, b = r - 
            r_d, sigma)
        put = GBSOption(TypeFlag = "p", S, X1, Time, r, b = r - 
            r_d, sigma)
        price1 = pmax(attr(call, "price"), 0)
        price2 = pmax(attr(put, "price"), 0)
        (price1 + price2) * ratio
    }
}


