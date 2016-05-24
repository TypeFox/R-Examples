## This example shows how to use frbs without 
 ## learning process.
 ## Note: some variables might be shared for other examples.

 ## Define shape and parameters of membership functions of input variables.
 ## Please see fuzzifier function to contruct the matrix.
 varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                       2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
                       2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                       nrow = 5, byrow = FALSE)

 ## Define number of fuzzy terms of input variables.
 ## Suppose, we have 3, 2, 3, and 3 numbers of fuzzy terms 
 ## for first, second, third and fourth variables, respectively.
 num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)
 
 ## Give the names of the fuzzy terms of each input variable.
 ## It should be noted that the names of the fuzzy terms must be unique,
 ## so we put a number for making it unique.
 varinput.1 <- c("a1", "a2", "a3")
 varinput.2 <- c("b1", "b2")
 varinput.3 <- c("c1", "c2", "c3")
 varinput.4 <- c("d1", "d2", "d3")
 names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)
 
 ## Define number of fuzzy terms of output variable.
 ## In this case, we set the number of fuzzy terms to 3.
 num.fvaloutput <- matrix(c(3), nrow=1)

 ## Give the names of the fuzzy terms of the output variable.
 ## Note: the names of the fuzzy terms must be unique.
 varoutput.1 <- c("e1", "e2", "e3")
 names.varoutput <- c(varoutput.1)

 ## Define the shapes and parameters of the membership functions of the output variables.
 varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                       nrow = 5, byrow = FALSE)

 ## Set interval of data.
 range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)

 ## Set weighted average method to be used as defuzzification method.
 type.defuz <- "WAM"
 ## We are using standard t-norm and s-norm.
 type.tnorm <- "MIN"
 type.snorm <- "MAX"
 type.implication.func <- "ZADEH"

 ## Give the name of simulation.
 name <- "Sim-0"

 ## Provide new data for testing. 
 newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow= 2, byrow = TRUE)
 ## the names of variables
 colnames.var <- c("input1", "input2", "input3", "input4", "output1")
 
 ######################
 ## 2. Using the same data as in the previous example, this example performs 
 ## step by step of the generation of a fuzzy rule-based system
 ######################
 ## Using Mamdani model
 type.model <- "MAMDANI"

 ## rules
 rule <- matrix(c("very a1","and","b1","and","slightly c1","and","very d1","->","e1",
                  "a2","and","extremely b2","and","c2","and","somewhat d2", "->", "e2", 
                  "slightly a3","and","b2","and","c2","and","dont_care", "->", "e3"), 
                  nrow=3, byrow=TRUE) 
 ## Check input data given by user.
 rule <- rulebase(type.model, rule, func.tsk = NULL)
 
 ## Fuzzification Module:
 ## In this function, we convert crisp values into fuzzy values 
 ## based on the data and the parameters of the membership function.
 ## The output: a matrix representing the degree of the membership of the data
 num.varinput <- ncol(num.fvalinput)
 MF <- fuzzifier(newdata, num.varinput, num.fvalinput, varinp.mf)
 
 ## Inference Module:
 ## In this function, we will calculate the confidence factor on the antecedent for each rule
 ## considering t-norm and s-norm.
 miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)

 ## Defuzzification Module
 ## In this function, we calculate and convert the fuzzy values back into crisp values. 
 range.output <- range.data[, ncol(range.data), drop = FALSE]
 result <- defuzzifier(newdata, rule, range.output, names.varoutput,
                   varout.mf, miu.rule, type.defuz, type.model, func.tsk = NULL)