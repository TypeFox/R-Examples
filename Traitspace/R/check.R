check <-
function(level_1,level_2,level_3, site, N){
# test data structure

level_1 <- as.factor(level_1)
level_2 <- data.frame(level_2)
level_3 <- data.frame(level_3)
site.name <- as.factor(site)
site <- data.frame(site)


if (!(length(level_1) == nrow(level_2) & nrow(level_2) == nrow(level_3) & nrow(level_3) == length(site.name))) stop("Please check your input data.")
if (!all(level_2 > 0)) stop("Negative value in level_2.")

result <- list(c, level_1 = as.factor(level_1), level_2 = data.frame(level_2), level_3 = data.frame(level_3), site.name = site.name, site = site)
return(result)
}
