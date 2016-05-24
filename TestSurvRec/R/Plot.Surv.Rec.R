Plot.Surv.Rec <-
function(XX)
{
    XL <- XX
    x <- factor(XL$group)
    Factores <- x
    x <- c(levels(x))
    Nivelesdefactores <- matrix(x)
    fit1 <- survfitr(Survr(id, time, event) ~ as.factor(group), 
        data = XL, type = "pe")
    fit2 <- survfitr(Survr(id, time, event) ~ 1, data = XL, type = "pe")
    print(summary(fit1))
    MedianaGPLECombinada <- q.search(fit2, q = 0.5)
    MedianaGPLGrupo01 <- q.search(fit1[[Nivelesdefactores[1, 
        1]]], q = 0.5)
    MedianaGPLGrupo02 <- q.search(fit1[[Nivelesdefactores[2, 
        1]]], q = 0.5)
    Nomb.grupos <- matrix(c("Pooled Group ", "1er Group ", "2do Group "))
    medianas <- matrix(c(MedianaGPLECombinada, MedianaGPLGrupo01, 
        MedianaGPLGrupo02))
    tabla <- data.frame(Group = Nomb.grupos, Median = medianas)
    print(tabla)
    dev.new()
    plot(fit2$time, fit2$survfunc, main = "Survival Curves of the groups", 
        xlab = "Time", ylab = "Probability of survival", xlim = c(0, 
            0.95 * max(fit2$time)), ylim = c(0, 1.05), type = "s", 
        col = "blue", lwd = 1, sub = R.version.string)
    lines(fit1[[Nivelesdefactores[1, 1]]]$time, fit1[[Nivelesdefactores[1, 
        1]]]$survfunc, type = "s", lty = 2, col = "red")
    lines(fit1[[Nivelesdefactores[2, 1]]]$time, fit1[[Nivelesdefactores[2, 
        1]]]$survfunc, type = "s", lty = 3, col = "black")
    legend("topright", c("Pooled Group", "First Group", "Second Group"), 
        col = c("blue", "red", "black"), lty = c(1, 2, 3))
}
