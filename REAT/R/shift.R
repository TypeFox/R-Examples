shift <-
function (region_t, region_t1, nation_t, nation_t1) {
  
industries <- length(region_t)
sum_region_t <- sum(region_t)
sum_region_t1 <- sum(region_t1)
sum_nation_t <- sum(nation_t)
sum_nation_t1 <- sum(nation_t1)
 
nts <- sum_region_t1-(sum_region_t*(sum_nation_t1/sum_nation_t))

i <- 0
m_i <- vector()
for (i in 1:industries) {
m_i[i] <- (region_t[i]*(nation_t1[i]/nation_t[i]))
}
nps <- sum(m_i)-(sum_region_t*(sum_nation_t1/sum_nation_t))

nds <- nts-nps

results <- c(nps, nds, nts)
names(results) <- c("nps", "nds", "nts")
return(results)

}
