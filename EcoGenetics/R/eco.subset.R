#' Subsetting an ecogen object by group
#' 
#' @param eco Object of class "ecogen". 
#' @param fact The name of the S slot column with labels assigning individuals to groups.
#' @param grp Label for the subset of individuals, contained in fact. 
#' @param missing Missing data argument This can take three values ("0", "NA" or "MEAN"),
#' as described in  \code{\link{ecogen}}.
#' Missing elements are treated as zeros in the default option.
#' 
#' @examples
#' \dontrun{
#' data(eco3)
#' eco3
#' eco.sub <-eco.subset(eco3,"structure", 1) 
#' eco.sub
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.subset",
					 
					 function(eco, 
                    fact, 
                    grp, 
                    missing = c("0", "NA",  "MEAN"))  {
  
  grupo <- eco@S
  x <- match(fact, colnames(eco@S), nomatch = 0)
  x <- x[!x == 0]
  
  missing <- match.arg(missing)
  
  if(length(x) == 0) {
    stop("incorrect factor name")
  }
  
  if(grp > max(as.numeric(grupo[, x]))) {
    stop(sprintf("the number of groups (%d) exceeds the number of
                 groups in the data (%d)", grp,
                 max(as.numeric(grupo[, x]))))
  }
  
  grupo <- which(grupo[, x] == grp)
  z <- ecogen()
  z@P <- eco@P[grupo, ]
  z@G <- eco@G[grupo, ]
  z@A <- eco@A[grupo, ]
  z@E <- eco@E[grupo, ]
  z@XY <- eco@XY[grupo, ]
  
  z@S <- as.data.frame(eco@S[grupo, ])
  #all S columns of z as factors, removing unused levels
  if(dim(z@S)[1] != 0) {
    for(i in 1:(ncol(z@S))) {
      z@S[, i] <- factor(z@S[, i])
    }
  }
  
  temp <- int.df2genind(eco@G[grupo, ], 
                        missing = missing,
                        ncod = eco@INT@ncod,
                        ploidy = eco@INT@ploidy,
                        type = eco@INT@type)
  
z@A <- as.data.frame(temp@tab)
z@INT <- int.genind2gendata(temp)
z@G <- as.data.frame(int.genind2df(temp))
  
  colnames(z@S) <- colnames(eco@S)
  
  z@C <- eco@C[grupo, ]
  z@OUT <- list()
  
  z
  
})

