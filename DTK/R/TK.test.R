TK.test <-
function (x = "data vector", f = "factor vector", a = "alpha level") 
{
    if (a == "alpha level") {
        a = 0.05
    }
    TukeyHSD(aov(x ~ f), conf.level = 1 - a)
}
