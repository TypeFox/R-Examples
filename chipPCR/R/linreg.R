#library(qpcR)
#linreg(fluo = reps[, 2])

#early implementation of linreg, still not done

# linreg <- function(cyc = 1L:length(fluo), fluo, max.it = 100) {
#   
#   if(abs(max(fluo)/min(fluo)) < 7)
#     stop("No amplification")
#   
#   baseline <- min(fluo)
#   
#   fluo_c <- fluo - baseline
#   #put here takeoff from qpcR as the alternative
#   #part I of figure 3A
#   sl <- get_slopes(cyc, fluo_c)
#   it1 <- 0
#   while(sl["up"] < sl["low"] && it1 < max.it) {
#     baseline <- 0.99*baseline
#     fluo_c <- fluo - baseline
#     sl <- get_slopes(cyc, fluo_c)
#     it1 <- it1 + 1
#   }
#   
#   #part II of figure 3A
#   step <- 0.5*baseline
#   baseline <- baseline + step
#   fluo_c <- fluo - baseline
#   sl <- get_slopes(cyc, fluo_c)
#   
#   it2 <- 0
#   while(abs(sl["up"] - sl["low"]) > 1e-5  && it2 < max.it) {
#     if(sl["up"] < sl["low"]) {
#       step <- step/2
#       baseline <- baseline - step
#       fluo_c <- fluo - baseline
#       sl <- get_slopes(cyc, fluo_c)
#     } else {
#       baseline <- baseline + step
#       fluo_c <- fluo - baseline
#       sl <- get_slopes(cyc, fluo_c)
#     } 
#     it2 <- it2 + 1
#   }
#   list(fluo = fluo_c, it1 = it1, it2 = it2)
# }
# 
# 
# get_slopes <- function(cyc, fluo_c) {
#   amp_range <- round(summary(bg.max(cyc, fluo_c), print = FALSE)[c("bg.stop", "amp.stop")], 0)
#   #start, midpoint, end
#   amp_points <- c(amp_range[1], round(mean(amp_range), 0), amp_range[2])
#   dat <- data.frame(cyc = cyc, fluo_c = fluo_c)
#   s_lower <- coef(lm(fluo_c ~ cyc, data = dat[amp_points[1]:amp_points[2], ]))["cyc"]
#   s_upper <- coef(lm(fluo_c ~ cyc, data = dat[amp_points[2]:amp_points[3], ]))["cyc"]
#   c(low = unname(s_lower), up = unname(s_upper))
# }
