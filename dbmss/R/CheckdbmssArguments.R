CheckdbmssArguments <-
function() {

  # Get the list of arguments of the parent function
  ParentFunction <- sys.call(-1)[[1]]
  ErrorFunction <- paste("Error in ", ParentFunction, ":")
  Args <- formals(match.fun(ParentFunction))
  
  # Get the point pattern
  X <- eval(expression(X), parent.frame())
  
  # X 
  if (!is.na(names(Args["X"]))) {
    if (!inherits(X, "wmppp"))
      stop(paste(ErrorFunction, "X is not of class wmppp"))    
  }

  # r
  if (!is.na(names(Args["r"]))) {
    r <- eval(expression(r), parent.frame())
    if (!is.null(r)) {
      if (!is.numeric(r) && !is.vector(r)) 
        stop(paste(ErrorFunction, "r must be a numeric vector"))
      if (length(r) < 2) 
        stop(paste(ErrorFunction, "r has length", length(r), "- must be at least 2"))
      if (r[1] != 0) 
        stop(paste(ErrorFunction, "First r value must be 0"))
      if (any(diff(r) <= 0)) 
        stop(paste(ErrorFunction, "successive values of r must be increasing"))  
    }
  }
    
  # ReferenceType 
  if (!is.na(names(Args["ReferenceType"]))) {
    ReferenceType <- eval(expression(ReferenceType), parent.frame())
    if (ReferenceType!="" & !ReferenceType %in% X$marks$PointType)
      stop(paste(ErrorFunction, "ReferenceType must be a point type of the point pattern, it cannot be", sQuote(ReferenceType)))    
  }
  # NeighborType 
  if (!is.na(names(Args["NeighborType"]))) {
    NeighborType <- eval(expression(NeighborType), parent.frame())
    if (NeighborType!="" & !NeighborType %in% X$marks$PointType)
      stop(paste(ErrorFunction, "NeighborType must be a point type of the point pattern, it cannot be", sQuote(NeighborType)))    
  }
  # Cases 
  if (!is.na(names(Args["Cases"]))) {
    Cases <- eval(expression(Cases), parent.frame())
    if (!Cases %in% X$marks$PointType)
      stop(paste(ErrorFunction, "Cases must be a point type of the point pattern, it cannot be", sQuote(Cases)))    
  }
  # Controls 
  if (!is.na(names(Args["Controls"]))) {
    Controls <- eval(expression(Controls), parent.frame())
    if (!is.null(Controls)) {
      if (!(Controls %in% X$marks$PointType))
        stop(paste(ErrorFunction, "Controls must be a point type of the point pattern, it cannot be", sQuote(Controls)))
    }
  }
  
  # CaseControl 
  if (!is.na(names(Args["CaseControl"]))) {
    CaseControl <- eval(expression(CaseControl), parent.frame())
    if (!is.logical(CaseControl))
      stop(paste(ErrorFunction, "CaseControl must be TRUE or FALSE, it cannot be", sQuote(CaseControl)))    
  }
  # Intertype 
  if (!is.na(names(Args["Intertype"]))) {
    Intertype <- eval(expression(Intertype), parent.frame())
    if (!is.logical(Intertype))
      stop(paste(ErrorFunction, "Intertype must be TRUE or FALSE, it cannot be", sQuote(Intertype)))    
  }
  # Weighted 
  if (!is.na(names(Args["Weighted"]))) {
    Weighted <- eval(expression(Weighted), parent.frame())
    if (!is.logical(Weighted))
      stop(paste(ErrorFunction, "Weighted must be TRUE or FALSE, it cannot be", sQuote(Weighted)))    
  }
  # Original 
  if (!is.na(names(Args["Original"]))) {
    Original <- eval(expression(Original), parent.frame())
    if (!is.logical(Original))
      stop(paste(ErrorFunction, "Original must be TRUE or FALSE, it cannot be", sQuote(Original)))    
  }

  # lambda
  if (!is.na(names(Args["lambda"]))) {
    lambda <- eval(expression(lambda), parent.frame())
    if (!is.null(lambda)) {
      if (!inherits(lambda, "im") & !is.numeric(lambda))
        stop(paste(ErrorFunction, "lambda must be an image of class im or a numeric vector, it cannot be", sQuote(lambda)))
    }
  }
  
  # NumberOfSimulations 
  if (!is.na(names(Args["NumberOfSimulations"]))) {
    NumberOfSimulations <- eval(expression(NumberOfSimulations), parent.frame())
    if (!is.numeric(NumberOfSimulations))
      stop(paste(ErrorFunction, "NumberOfSimulations must be a number, it cannot be", sQuote(NumberOfSimulations)))    
    if (NumberOfSimulations <= 0)
      stop(paste(ErrorFunction, "NumberOfSimulations must be positive, it cannot be", sQuote(NumberOfSimulations)))    
  }
  # Alpha 
  if (!is.na(names(Args["Alpha"]))) {
    Alpha <- eval(expression(Alpha), parent.frame())
    if (!is.numeric(Alpha))
      stop(paste(ErrorFunction, "Alpha must be a number, it cannot be", sQuote(Alpha)))    
    if (Alpha<=0 | Alpha>=1)
      stop(paste(ErrorFunction, "Alpha must be strictly between 0 and 1, it cannot be", sQuote(Alpha)))    
  }
  
  # Adjust 
  if (!is.na(names(Args["Adjust"]))) {
    Adjust <- eval(expression(Adjust), parent.frame())
    if (!is.numeric(Adjust))
      stop(paste(ErrorFunction, "Adjust must be a number, it cannot be", sQuote(Adjust)))    
    if (Adjust<=0)
      stop(paste(ErrorFunction, "Adjust must be strictly positive, it cannot be", sQuote(Adjust)))    
  }
  
  # Approximate 
  if (!is.na(names(Args["Approximate"]))) {
    Approximate <- eval(expression(Approximate), parent.frame())
    if (!is.numeric(Approximate))
      stop(paste(ErrorFunction, "Approximate must be a number, it cannot be", sQuote(Approximate)))    
    if (Approximate < 0)
      stop(paste(ErrorFunction, "Approximate must be positive, it cannot be", sQuote(Approximate)))    
  }
  
  return (TRUE)
}
    