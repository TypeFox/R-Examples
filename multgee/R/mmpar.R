mmpar <-
function (LORem, LORstr, timepairs, homogeneous) 
{
    if (timepairs == 1) {
        LORem <- "2way"
        LORstr <- switch(LORstr, category.exch = "uniform", 
            RC = "time.exch",            
            uniform = "uniform", time.exch = "time.exch")
        if (LORstr == "uniform") 
            fmla <- counts ~ factor(x) + factor(y) + x:y
        if (LORstr == "time.exch") {
            fmla <- if (homogeneous) 
                counts ~ x + y + MultHomog(x,y)
            else counts ~ x + y + Mult(x,y)
        }
        
    }
    else if (LORem == "2way") {
        if (LORstr == "uniform") 
            fmla <- counts ~ (factor(x) + factor(y)) * factor(tp) + 
                factor(tp):x:y
        if (LORstr == "time.exch" | LORstr == "RC") {
            fmla <- if (homogeneous) 
                 counts ~ x + y + MultHomog(x,y)
            else counts ~ x + y + Mult(x,y)
        }
        
    }
    else {
        if (LORstr == "category.exch") 
            fmla <- counts ~ (factor(x) + factor(y)) * factor(tp) + 
                factor(tp):x:y
        if (LORstr == "uniform") 
            fmla <- counts ~ (factor(x) + factor(y)) * factor(tp) + 
                x:y
        if (LORstr == "time.exch") {
            fmla <- if (homogeneous) 
                counts ~ (x + y) * factor(tp) + 
                  MultHomog(x, y)
            else counts ~ (x + y) * factor(tp) + 
                Mult(x, y)
        }
        
    }
    list(LORem = LORem, LORstr = LORstr, fmla = fmla)
}

