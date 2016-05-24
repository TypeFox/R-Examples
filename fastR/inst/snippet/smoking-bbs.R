smTab2 <- rbind(NonExperimenter=c(34,4,296), 
                   Experimenter=c(15,3,213))
colnames(smTab2) <- c('Never','Hardly ever','Sometimes or a lot')
smTab2
chisq.test(smTab2)
