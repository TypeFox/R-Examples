rm(list=ls())

abc_optim <- function(
  par,               # Vector de parametros a opti 
  fn,                # Funcion objetivo
  D=length(par),     # Numero de parametros
  ...,               # Argumentos de la funcion (M, x0, X, etc.)
  NP = 40,           # Numero de abejas
  FoodNumber = NP/2, # Fuentes de alimento 
  lb = -Inf,         # Limite inferior de recorrido
  ub = +Inf,         # Limite superior de recorrido
  limit = 100,        # Limite con que se agota una fuente de alimento
  maxCycle = 1000,   # Numero maximo de iteraciones 
  optiinteger=FALSE, # TRUE si es que queremos optimizar en [0,1] (binario)
  criter=50
)
{
  # Checking limits
  if (length(lb)>0) lb <- rep(lb, D)
  if (length(ub)>0) ub <- rep(ub, D)

  lb[is.infinite(lb)] <- -(.Machine$double.xmax*1e-10)
  ub[is.infinite(ub)] <- +(.Machine$double.xmax*1e-10)
  
  # Initial params
  Foods       <- matrix(double(FoodNumber*D), nrow=FoodNumber)
  f           <- double(FoodNumber)
  fitness     <- double(FoodNumber)
  trial       <- double(FoodNumber)
  prob        <- double(FoodNumber)
  solution    <- double(D)
  ObjValSol   <- double(1)
  FitnessSol  <- double(1)
  neighbour   <- integer(1)
  param2change<- integer(1)
  GlobalMin   <- double(1)
  GlobalParams<- double(D)
  #GlobalMins  <- double(runtime)
  r           <- integer(1)

  # Fun
  fun <- function(par) fn(par, ...)
  
  # Fitness function
  CalculateFitness <- function(fun)
  {
    if (fun >= 0) return(1/(fun + 1))
    else return(1 + abs(fun))
  }
  # CalculateFitness(f[1])
  
  # The best food source is memorized
  MemorizeBestSource <- function() 
  {
    change <- 0
    for(i in seq(1,FoodNumber)) {
      if (f[i] < GlobalMin) {
        change <- change + 1
        GlobalMin <<- f[i]
        
        # Replacing new group of parameters
        GlobalParams <<- Foods[i,]
      }
    }
    # Increasing persistance
    if (!change) persistance <<- persistance + 1
  }
  
  # Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub 
  # Counters of food sources are also initialized in this function
  
  init <- function(index, firstinit=FALSE, ...) {
    if (optiinteger) Foods[index,] <<- runif(D) > .5
    else {
      if (!firstinit) {
        Foods[index,] <<- sapply(1:D, function(k) runif(1,lb[k],ub[k]) )
      }
      else {
        # For the first initialization we set the bees at
        # specific places equaly distributed through the
        # bounds.
        Foods[index,] <<- 
          sapply(1:D, function(k) {
            seq(lb[k],ub[k],length.out=FoodNumber)[index]
          }
          )
      }
    }
    
    solution <<- Foods[index,]
    
    f[index] <<- fun(solution)

    fitness[index] <<- CalculateFitness(f[index])
    trial[index] <<- 0
    
  }
  # init(2)
  
  # All food sources are initialized
  initial <- function(firstinit=FALSE) {
    sapply(1:FoodNumber, init, firstinit=firstinit)
    
    GlobalMin <<- f[1]
    GlobalParams <<- Foods[1,]
  }
  
  # initial()
  
  
  SendEmployedBees <- function() {
    for (i in 1:FoodNumber) {
      r <- runif(1)
      # The parameter to be changed is determined randomly
      param2change <- floor(r*D) + 1 
      
      # A randomly chosen solution is used in producing a mutant solution of the solution i
      neighbour <- floor(r*FoodNumber) + 1 
 
      # Randomly selected solution must be different from the solution i
      while(neighbour==i)
        neighbour <- floor(runif(1)*FoodNumber) + 1
      
      solution <<- Foods[i,]
      
      # v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) 
      r <- runif(1)
      
      if (optiinteger) solution[param2change] <<- r > 0.5
      else {
        solution[param2change] <<- 
          Foods[i,param2change]+
          (Foods[i,param2change]-Foods[neighbour,param2change])*(r-0.5)*2

        # if generated parameter value is out of boundaries, it is shifted onto the boundaries
        if (solution[param2change]<lb[param2change])
          solution[param2change]<<-lb[param2change]
        
        if (solution[param2change]>ub[param2change])
          solution[param2change]<<-ub[param2change]
      }
      
      ObjValSol <<- fun(solution)
      FitnessSol <<- CalculateFitness(ObjValSol)
      
      # a greedy selection is applied between the current solution i and its mutant*/
      if (FitnessSol>fitness[i]) {
        # If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
        trial[i] <<- 0;
        #for(j in 1:D) Foods[i,j] <<- solution[j]
        Foods[i,] <<- solution
        f[i]<<- ObjValSol
        fitness[i]<<-FitnessSol
      }
      else {
        # the solution i can not be improved, increase its trial counter*/
        trial[i] <<- trial[i]+1
      }
    }
  }
  
  
  # A food source is chosen with the probability which is proportioal to its quality*/
  # Different schemes can be used to calculate the probability values*/
  # For example prob(i)=fitness(i)/sum(fitness)*/
  # or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
  # probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/
  CalculateProbabilities <- function() {
    maxfit <- fitness[1]
    for (i in 1:FoodNumber) 
      if (fitness[i] > maxfit) maxfit <- fitness[i]
    
    prob <<- .9*(fitness/maxfit) + .1
#     prob[is.nan(prob)]  <<- .1
  }
  
  SendOnlookerBees <- function()
  {
    # Onlooker Bee phase
    i <- 1
    t <- 0
    while (t < FoodNumber)
    {
      r <- runif(1)
      # choose a food source depending on its probability to be chosen
      if (r < prob[i]) {
        t <- t + 1
        r <- runif(1)
        
        # The parameter to be changed is determined randomly
        param2change <- floor(r*D) + 1
        
        # A randomly chosen solution is used in producing a mutant solution of the solution i
        neighbour <- floor(r*FoodNumber) + 1
        
        #Randomly selected solution must be different from the solution i*/        
        while(neighbour==i)
          neighbour <- floor(runif(1)*FoodNumber) + 1

        solution <<- Foods[i,]
        
        # v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
        r <- runif(1)
        
        if (optiinteger) solution[param2change] <<- r > .5
        else 
        {
          solution[param2change] <<- 
            Foods[i,param2change]+
            (Foods[i,param2change]-Foods[neighbour,param2change])*(r-0.5)*2
          
          # if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
          if (solution[param2change]<lb[param2change]) 
            solution[param2change] <<- lb[param2change]
          
          if (solution[param2change]>ub[param2change]) 
            solution[param2change] <<- ub[param2change]
          
        }
        
        ObjValSol <<- fun(solution)
        FitnessSol <<- CalculateFitness(ObjValSol)
        
        # a greedy selection is applied between the current solution i and its mutant*/
        if (FitnessSol>fitness[i])
        {
          # If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
          trial[i] <<- 0
          Foods[i,] <<- solution
          
          f[i]<<-ObjValSol
          fitness[i]<<-FitnessSol
        } #if the solution i can not be improved, increase its trial counter*/
        else trial[i] <<- trial[i]+1
      }
      i <- i + 1
      if (i==FoodNumber) i <- 1
      # end of onlooker bee phase
    }
  }
  
  # determine the food sources whose trial counter exceeds the "limit" value. In Basic ABC, only one scout is allowed to occur in each cycle*/
  
  SendScoutBees <- function() {
    maxtrialindex <- 1
    for (i in 1:FoodNumber) {
      if (trial[i] > trial[maxtrialindex]) maxtrialindex <- i
    }
    
    if (trial[maxtrialindex] >= limit) init(maxtrialindex)
  }
  
  persistance <- 0
  
  # Inicializa funcion
  initial(firstinit=T)
  
  # Memoriza la primera mejor solucion
  MemorizeBestSource() 
  
  iter <- 0
  # Comienza a iterar
  while ((iter <- iter + 1) < maxCycle)
  {
    SendEmployedBees()
    CalculateProbabilities()
    SendOnlookerBees() 
    MemorizeBestSource()
    if (persistance > criter) break
    SendScoutBees()
  }

  return(
    list(
      par=GlobalParams,
      value=fun(GlobalParams),
      counts=c("function"=iter)
      )
    )
  
}

################################################################################
# Ejemplos
################################################################################

X <- c(3,2,3,1)

# Funcion de matching
fun <- function(lambda, x0, X, M)
{
  norm((x0 - X)*lambda, type="2") + exp(abs(sum(lambda > 0) - M))
}

# Mejor vecino para
#  x0 = 2
#  X  = c(3,2,3,1)
#  M  = 1
# El mejor resultado debe ser [0,1,0,0]
x1 <- abc_optim(rep(0,4), fun, x0=2, X=X, M=1, lb=0, ub=1, optiinteger=T)
x1

# Mejores dos vecinos para
#  x0 = 3
#  X  = c(3,2,3,1)
#  M  = 2
# El mejor resultado debe ser [1,0,1,0]
x2 <- abc_optim(rep(0,4), fun, x0=3, X=X, M=2, lb=0, ub=1, optiinteger=T)
x2

################################################################################
# Definicion de la funcion
fun <- function(x) {
  -cos(x[1])*cos(x[2])*exp(-((x[1] - pi)^2 + (x[2] - pi)^2))
}

abc_optim(rep(0,2), fun, lb=-5, ub=5, criter=50)

optim(rep(0,2), fn=fun) #lower=-5,upper=5)

################################################################################
# Definicion de la funcion

fun <- function(x) {
  -4+(x[1]^2 + x[2]^2)
}

abc_optim(c(1,1), fn=fun, lb=-100000, ub=100000,criter=100)

################################################################################
# Definicion de la funcion

fun <- function(x) {
  -(x^4 - 2*x^2 - 8)
}

abc_optim(0, fn=fun, lb=-2, ub=2,criter=100)
# 
