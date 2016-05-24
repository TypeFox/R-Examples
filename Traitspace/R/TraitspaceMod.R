TraitspaceMod <-
function(level_1,level_2,level_3,site, N = 100){

#library(mclust)
mod = TRUE
input <- check(level_1,level_2,level_3,site, N) 

level_1 <- input$level_1
level_2 <- input$level_2
level_3 <- input$level_3
site.name <- input$site.name
unique.site <-unique(input$site)
pred.site <- unique.site[rep(seq_len(nrow(unique.site)), each=N),]


# step_1A 
step_1A <- Mstep(site.name, level_2)
par_2_3 <- step_1A$par

# step_1B 
step_1B <- Mstep(level_1, level_2)
par_2_1 <- step_1B$par

# step_2A 
step_2A <- step_2A(par_2_3, N, mod, level_2, level_3, step_1A$model, linear = FALSE)
level_2_sample <- step_2A$sample
P_level_2_level_3 <- step_2A$P_level_2_level_3

# step_2B 
step_2B <- step_2B(par_2_1, level_2_sample,level_1)
P_level_2_level_1 <- step_2B$P_level_2_level_1

# step_2C 
step_2C <- step_2C(P_level_2_level_1,level_1)
P_level_1_level_2_level_3 <- step_2C$P_level_1_level_2_level_3

# step_2D 
step_2D <- step_2D(P_level_2_level_3, P_level_1_level_2_level_3, N,level_1,level_3,pred.site)

#trueP
true.p <- trueP(level_1, site.name)

result <- list(check = input, step_1A = step_1A, step_1B = step_1B, step_2A = step_2A, step_2B = step_2B, 
step_2C = step_2C, step_2D = step_2D, BIC_TE = step_1A$summary.pdf, BIC_TS = step_1B$summary.pdf, predicted.p = step_2D, true.p = true.p)
class(result) <- "traitspace"
return(result)
}
