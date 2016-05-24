library(distcomp)
defn <-
structure(list(id = "a4ed44c4cbd39b0e", compType = "RankKSVD", 
    projectName = "SVDTest", projectDesc = "SVD Test Example", 
    rank = 2L, ncol = 5L), .Names = c("id", "compType", "projectName", 
"projectDesc", "rank", "ncol"), row.names = c(NA, -1L), class = "data.frame")
sites <-
list(structure(list(name = "Site1", url = "http://localhost:8109/ocpu"), .Names = c("name", 
"url")), structure(list(name = "Site2", url = "http://localhost:8109/ocpu"), .Names = c("name", 
"url")), structure(list(name = "Site3", url = "http://localhost:8109/ocpu"), .Names = c("name", 
"url")))
master <- makeMaster(defn)
for (site in sites) {
   master$addSite(name = site$name, url = site$url)
}
result <- master$run()
print(master$summary())
