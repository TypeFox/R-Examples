####################################################################
##########       Algorítmo EM para SNI            ###############
#                   alterado em 25/04/09                            #

smsn.nl <- function(y, x = NULL, z = NULL, betas = NULL, sigma2 = NULL, shape = NULL,  rho = NULL, 
           nu = NULL, nlf = NULL, rho.func = 1, reg.type = "Homoscedastic", criteria = FALSE, 
           family = "Skew.t", error = 0.00001, iter.max = 100){

  if(ncol(as.matrix(y)) > 1) stop("Only univariate non linear regression supported!")
  if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
  if( (length(x) == 0) | (length(betas) == 0) | (length(sigma2) == 0) | (length(shape) == 0) ) stop("All parameters must be provided.")
  if( !is.function(nlf) ) stop("nfl parameter must be a function!") 
  if( (family != "Skew.t") && (family != "t") && (family != "Skew.cn") && (family != "Skew.normal") && (family != "Normal") && (family != "Skew.slash"))
  stop("Distribution family not supported. Check documentation!") 
  if( (family == "Normal") | (family == "Skew.normal"))
  {
  if(length(nu) != 0) stop("For the Skew.Normal and the Normal family the Nu parameters must not be provided.")
  } 
  if( (family == "Skew.t") | (family == "t") | (family == "Skew.slash"))
  {
  if(length(nu) > 1) stop("Nu parameters must have only one parameter")
  if(length(nu) == 0) stop("Nu parameters must be provided.")
  if(nu <= 0) stop("Nu parameters must be positive.")
  } 
  if(family == "Skew.cn")
  {                               
  if(length(nu) !=2) stop("For the Skew.cn nu must have 2 parameters")
  if(nu[1] <=0 || nu[1] >= 1) stop("nu[1] must be in (0,1)")
  if(nu[2] <=0 || nu[2] >= 1) stop("nu[2] must be in (0,1)")
  }
  if( (reg.type != "Homoscedastic") && (reg.type != "Ho") && (reg.type != "Heteroscedastic") && (reg.type != "He") )
  stop("Regression type not supported")  
  if( (reg.type == "Homoscedastic") | (reg.type == "Ho") )

    out <- smsn.nl.intern(y, x, betas, sigma2, shape, nu, nlf, criteria, family, error, iter.max)
  if( (reg.type == "Heteroscedastic") | (reg.type == "He") ){
  if( (length(z) == 0) | (length(rho) == 0) ) stop("All parameters must be provided.")
  if((rho.func != 1) & (rho.func != 2)) stop("The non linear function for sigma can only have type 1 or 2. Check documentation!")

    out <- smsn.nlh.intern(y, x, z, betas, sigma2, shape, rho, nu, nlf, rho.func, criteria, family, error, iter.max)
  }
  
  out
}


