## Create a sample dataset with 3 SNP markers
set.seed(125141)

g1 <- sample( x=c('C/C', 'C/T', 'T/T'), 
              prob=c(.6,.2,.2), 100, replace=T)
g2 <- sample( x=c('A/A', 'A/G', 'G/G'), 
              prob=c(.6,.1,.5), 100, replace=T)
g3 <- sample( x=c('C/C', 'C/T', 'T/T'), 
              prob=c(.2,.4, 4), 100, replace=T)

y <- rnorm(100) + (g1=='C/C')  +
     0.25 * (g2=='A/A' | g2=='A/G')

pid <- formatC( abs(rnorm(100))*1e6, format="d", flag=0, width=8)

## Form into a data frame
data <- data.frame( PID=pid, DELTA.BMI=y, c104t=g1, a1691g=g2, c2249t=g3)

## Save as a file
write.table(data, file="example_data.csv", sep=",", quote=F, row.names=F)
