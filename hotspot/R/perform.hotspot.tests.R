perform.hotspot.tests <-
function (hotspot.data) 
{
    result <- c()
    for (project.name in names(hotspot.data)) {
        project.data <- hotspot.data[[project.name]]
        for (actual.version in names(project.data)) {
            version.data <- project.data[[actual.version]]
            cor.test.result <- cor.test(version.data$hotspot, 
                version.data$rmi, method = "spearman", alternative = "less")
            result <- rbind(result, c(project.name, actual.version, 
                "hotspot-rmi", cor.test.result$p.value, cor.test.result$estimate[["rho"]]))
            cor.test.result <- cor.test(version.data$hotspot, 
                version.data$bug, method = "spearman", alternative = "greater")
            result <- rbind(result, c(project.name, actual.version, 
                "hotspot-bug", cor.test.result$p.value, cor.test.result$estimate[["rho"]]))
        }
    }
    return(result)
}
