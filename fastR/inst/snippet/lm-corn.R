# the corn data frame has an inconvenient "shape" 
# (each type of corn is in its own column)
head(corn,3)                                   
# this puts all the yields in one column and type of seed in another
corn2 <- stack(corn)                         
corn2[c(1,2,12,13),]
# the default variable names aren't great, so we rename them
names(corn2) <- c('yield','treatment')       
corn2[c(1,2,12,13),]
summary(yield~treatment,data=corn2,fun=favstats)
corn.model <- lm(yield~treatment,data=corn2)
###hop:3-9
summary(corn.model)
