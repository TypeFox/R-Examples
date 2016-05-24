
data(x100)
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## out$blocks contains 3 data frames

## To illustrate two-level blocking, with multiple level two units per
##  level one unit:
x100.tmp <- x100
for(i in (1:nrow(x100.tmp))){if((i %% 2) == 0){x100.tmp$id[i] <- x100.tmp$id[i-1]}}
rm(i)

out2 <- block(x100.tmp, groups = "g", n.tr = 2, id.vars = c("id", "id2"),
              block.vars = c("b1", "b2"), algorithm="optGreedy",
              distance = "mahalanobis", level.two = TRUE, valid.var =
              "b1", valid.range = c(0,500), namesCol = 
              c("State 1", "City 1", "State 2", "City 2"), verbose = TRUE) 
