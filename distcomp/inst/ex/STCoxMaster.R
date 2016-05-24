library(distcomp)
defn <-
structure(list(id = "46c17593ba8a205f", compType = "StratifiedCoxModel", 
    projectName = "STCoxTest", projectDesc = "Stratified Cox Test", 
    formula = "Surv(time, censor) ~ age + becktota + ndrugfp1 + ndrugfp2 + ivhx3 + race + treat"), .Names = c("id", 
"compType", "projectName", "projectDesc", "formula"), row.names = c(NA, 
-1L), class = "data.frame")
sites <-
list(structure(list(name = "Site1", url = "http://localhost:8109/ocpu"), .Names = c("name", 
"url")), structure(list(name = "Site2", url = "http://localhost:8109/ocpu"), .Names = c("name", 
"url")))
master <- makeMaster(defn)
for (site in sites) {
   master$addSite(name = site$name, url = site$url)
}
result <- master$run()
print(master$summary())
