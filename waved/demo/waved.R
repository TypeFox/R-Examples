

cat('Would you like to use the default parameter setting  \n')
cat('to reproduce the figures of: The WaveD Transform in R \n')
cat(' by Marc Raimondo and Michael Stewart, Journal of Statistical Software (2007).\n')
cat('Yes or No (Y/N)?\n')


answer <- readline()
pr=switch(answer,y=,Y=TRUE,FALSE)


data.demo=waved.example(pr)

