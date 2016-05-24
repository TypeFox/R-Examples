reproduction <- function(population,xfreq,mfreq,xoverpoints,pinvert,elitism){

        newpop <- array(NA,dim=dim(population))
        popsize <- dim(population)[1]

        if (elitism){
                xindex <- c(2:popsize)[runif((popsize-1))>xfreq]
                # ensure fittest solution from previous generation (population[1,]) is unchanged in next gen
                xindex <- c(1,xindex)
        }
        else {
                xindex <- c(1:popsize)[runif(popsize)>xfreq]
        }
        newpop[xindex,] <- population[xindex,]

        if(length(xindex)<popsize){
                xindex <- sample(setdiff(c(1:popsize),xindex))
                newpop[xindex,] <- crossover(population[xindex,],xoverpoints,pinvert)
        }
        if (elitism){
                newpop[c(2:popsize),] <- mutation(newpop[c(2:popsize),],mfreq)
        }
        else {
                newpop <- mutation(newpop,mfreq)
        }
        newpop

}
