Print.Summary <-
function (XX) 
{
    XL <- XX
    fit1 <- survfitr(Survr(id, time, event) ~ as.factor(group), 
        data = XL, type = "pe")
    fit2 <- survfitr(Survr(id, time, event) ~ 1, data = XL, type = "pe")
    Plot.Surv.Rec(XL)
    Dif.Surv.Rec(XL, "all")
}
