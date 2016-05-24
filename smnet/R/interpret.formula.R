
interpret.formula <- function (formula.object) 
{
  formula.env      <- environment(formula.object)
  tf               <- terms.formula(formula.object, specials = c("m", "network"))
  response.term    <- attr(tf, "response") 
  variables.term   <- attr(tf, "variables")
  terms            <- attr(tf, "term.labels")
  specials         <- attr(tf, "specials")
  vtab             <- attr(tf, "factors")
  n.terms          <- length(terms)
  
  # what to do if a response variable is present
  if(response.term < 1) stop("response not provided in user formula")
  response <- as.character(variables.term[2])
  pf <- paste(response, "~", sep = "")

  #   identify smooth terms and network terms
  sp <- specials$m
  net <- specials$network

  # find the locations in the explanatory variables of the smooth terms
  if (length(sp) > 0) for(i in 1:length(sp)) sp[i] <- (1:n.terms)[as.logical(vtab[sp[i], ])]
  
  # find the locations in the explanatory variables of the network terms
  if (length(net) > 0) for (i in 1:length(net)) net[i] <- (1:n.terms)[as.logical(vtab[net[i], ])]

  # ns is the number of smooth terms; n.terms is total number of terms
  k <- l <- kp <- 1
  len.sp <- length(sp)
  len.net <- length(net)
  penalised.terms <- list()
  
  if (n.terms > 0){
    for (i in 1:n.terms){
      if (k <= len.sp  && sp[k] == i){
        penalised.terms[[k + l - 1]] <- eval(parse(text = terms[i]), envir = formula.env)
        k <- k + 1
      } else if (l <= len.net  && net[l] == i){
        penalised.terms[[l + k - 1]] <- eval(parse(text = terms[i]), envir = formula.env)
        l <- l + 1
      } else if (kp > 1) {
        # if its the first term then no + sign needed
        pf <- paste(pf, "+", terms[i], sep = "")
      } else {
        pf <- paste(pf, terms[i], sep = "")
        kp <- kp + 1
      }
    }
  } 
    
  # is the intercept omitted?  If so add a 1 or a -1 to pf
  if (attr(tf, "intercept") == 0) {
    pf   <- paste(pf, "-1", sep = "")
    pfok <- ifelse(kp > 1, 1, 0)
  } else {
    pfok <- 1
    if (kp == 1) pf <- paste(pf, "1")
  }
  
  fake.formula <- pf
  if (length(penalised.terms) > 0){
    for (i in 1:length(penalised.terms)){
      n.terms   <- length(penalised.terms[[i]]$term)
      ff1       <- paste(penalised.terms[[i]]$term[1:n.terms], collapse = "+")
      if(!ff1 == "") fake.formula <- paste(fake.formula, "+", ff1)
    }
  }
  fake.formula <- as.formula(fake.formula, formula.env)
  list(smooth.spec = penalised.terms, fake.formula = fake.formula, response = response)
}