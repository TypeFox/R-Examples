boyce <-
function(prediction, prediction.map, categories = NULL, cv.sets = 10, type = "manual", 
    outcat = "cbi.results") {
    spearman.coef <- function(X, Y) {
        x <- rank(X)
        y <- rank(Y)
        x.mean <- mean(x)
        y.mean <- mean(y)
        dif.x <- x - x.mean
        dif.y <- y - y.mean
        dif.x.sq <- dif.x^2
        dif.y.sq <- dif.y^2
        prod.dif <- dif.x * dif.y
        return(sum(prod.dif)/sqrt(sum(dif.x.sq) * sum(dif.y.sq)))
    }
    plot.classification <- function(predictionc, prediction.mapc, categoriesc) {
        pi <- hist(predictionc, breaks = categoriesc, plot = FALSE)$counts
        Pi <- pi/sum(pi)
        ai.function <- function(prediction.map.subset) hist(prediction.map.subset, 
            breaks = categoriesc, plot = F)$counts
        ai <- matrix(unlist(tapply(prediction.mapc, cv.vector, ai.function)), cv.sets, 
            byrow = T)
        Ei <- ai/apply(ai, 1, sum)
        Fi <- t(apply(Ei, 1, function(x) Pi/x))
        Fi <- replace(Fi, is.nan(Fi), 0)
        spearman.function <- function(x) spearman.coef(x, categoriesc[-1])
        spea.coef <- apply(Fi, 1, spearman.function)
        median.Fi <- apply(Fi, 2, median)
        mean.Fi <- apply(Fi, 2, mean)
        sd.Fi <- apply(Fi, 2, sd)
        results[["coefficients"]] <- c(mean(spea.coef), summary(lm(median.Fi ~ categoriesc[-1]))$adj.r.squared)
        names(results[["coefficients"]]) <- c("spear.coef", "adj.r2")
        results[["intervals"]] <- categoriesc
        names(results[["intervals"]]) <- c("non", "unsuitable", "marginal", "suitable", 
            "optimal")
        plot((categoriesc[-5] + categoriesc[-1])/2, median.Fi, type = "l", xlab = "Habitat suitability", 
            ylab = "Predicted-expected ratio", xlim = c(categoriesc[1], categoriesc[5]), 
            ylim = c(0, max(mean.Fi + sd.Fi) + (max(mean.Fi + sd.Fi) * 0.25)))
        points((categoriesc[-5] + categoriesc[-1])/2, median.Fi, pch = 20)
        segments((categoriesc[-5] + categoriesc[-1])/2, mean.Fi - sd.Fi, (categoriesc[-5] + 
            categoriesc[-1])/2, mean.Fi + sd.Fi)
        abline(v = categoriesc[-1], lty = 2)
        text(categoriesc[-5] + ((categoriesc[-1] - categoriesc[-5])/2), rep(max(mean.Fi + 
            sd.Fi) + (max(mean.Fi + sd.Fi) * 0.25), 4), labels = c("unsuitable", 
            "marginal", "suitable", "optimal"), font = 3)
        text(x = 0, y = (max(mean.Fi + sd.Fi) + (max(mean.Fi + sd.Fi) * 0.25)) * 
            3/4, labels = paste("Spearman coefficient = ", round(mean(spea.coef), 
            2), sep = ""), font = 2, pos = 4)
        class(results) <- "CBI"
        exp.eval <- paste(outcat, " <<- results", sep="")
        eval(parse(text = exp.eval))
    }
    subsetting <- function(group, subgroups) c(sample(c(rep(1:subgroups, group%/%subgroups), 
        1:(group%%subgroups + 1)), group))
    cv.vector <- subsetting(length(prediction.map), cv.sets)
    results <- list()
    if (type == "none") {
        plot.classification(predictionc = prediction, prediction.mapc = prediction.map, 
            categoriesc = categories)
    }
    if (type == "manual") {
        manual.boyce.classification <- function(unsuitable = miniGUIscale(from = 0, 
            to = 1, by = 0.01), marginal = miniGUIscale(from = 0, to = 1, by = 0.01), 
            suitable = miniGUIscale(from = 0, to = 1, by = 0.01), origin = 0, optimal = 1) {
            plot.classification(predictionc = prediction, prediction.mapc = prediction.map, 
                categoriesc = c(origin, unsuitable, marginal, suitable, optimal))
        }
        hs.manual.classification <- makeWidgetCmd("Boyce classification", manual.boyce.classification)
        hs.manual.classification()
    }
}
