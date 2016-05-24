`gethats` <-
function (nobs, nvars, x.t, XWXinv, ww) 
{
    .C("hatsc", n = as.integer(nobs), p = as.integer(nvars), 
        x = as.double(x.t), invfisherinf = as.double(XWXinv), 
        w = as.double(ww), hat = double(nobs))$hat
}
