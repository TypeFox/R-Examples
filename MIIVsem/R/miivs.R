#' Model-implied instrumental variable (MIIV) search 
#'
#' A key step in the MIIV-2SLS approach is to transform the SEM by replacing the latent variables with their scaling indicators minus their errors.  Upon substitution the SEM is transformed from a model with latent variables to one containing observed variables with composite errors.  The miivs function automatically makes this transformation.
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param miivs.out A logical indicating whether or not to print the MIIVs as an object for later use. Default is \code{FALSE}.
#'
#' @return eqns
#' @return modeqns
#' 
#' @details 
#' \itemize{
#'  \item \code{LHS} The "dependent" variable.
#'  \item \code{RHS} The right hand side variables of the transformed equation.
#'  \item \code{Composite Disturbance}  Elements of the composite errors in the transformed equation.
#'  \item \code{MIIVs} The model implied instrumental variables for each equation.
#' }
#' 
#' @references 
#' 
#'  Bollen,	K. A. and	D. J.	Bauer.	2004.	Automating	the	Selection	of 
#' 	Model-Implied	Instrumental	Variables.	\emph{Sociological	Methods	and	
#' 	Research}, 32, 425-52.
#' 	
#' 	Bauldry, S.	2014.	miivfind: A command for identifying model-implied instrumental 
#' 	variables for structural equation models in Stata.	\emph{Stata Journal}, 14:4.
#' 	
#'
#' @examples
#'  bollen1989a_model <- '
#'
#'    Eta1 =~ y1 + y2  + y3  + y4  
#'    Eta2 =~ y5 + y6  + y7  + y8    
#'    Xi1  =~ x1 + x2 + x3 
#'
#'    Eta1 ~ Xi1  
#'    Eta2 ~ Xi1 
#'    Eta2 ~ Eta1 
#'
#'    y1   ~~ y5
#'    y2   ~~ y4
#'    y2   ~~ y6
#'    y3   ~~ y7
#'    y4   ~~ y8
#'    y6   ~~ y8 
#'  '
#'  
#'  miivs(bollen1989a_model)
#'  
#' @export
miivs <- function(model, miivs.out = FALSE) {

  fit <- lavaan(model,
                auto.fix.first = TRUE, 
                auto.var = TRUE, 
                auto.cov.lv.x = TRUE)

  pt  <- lavMatrixRepresentation(parTable(fit))

  if (any(!(is.na(pt$ustart) | pt$ustart ==1))) {
    stop(paste("Enter numerical constraints using labels."))
  }
  
  hof <- unique(pt[pt$op == "=~" & pt$rhs %in% pt[pt$op == "=~", "lhs"],"lhs"])
  
  if (length(hof) > 0) {higher_order <- TRUE} else {higher_order <- FALSE}
    
  if (higher_order == TRUE){
    fof <- unique(pt[pt$op == "=~" &!(pt$rhs %in% pt[pt$op == "=~", "lhs"]),"lhs"])
    sof <- unique(pt[pt$op == "=~" &  pt$rhs %in% fof, "lhs"])
    tof <- unique(pt[pt$op == "=~" &  pt$rhs %in% sof, "lhs"])
    
    if (length(tof) > 1){
      stop(paste("MIIVsem does not currently support factor orders 
                  greater than 2."))
    }
  
    latEnd <- setdiff(lavNames(fit, type = "lv.nox"), fof)
    latExo <- c(lavNames(fit, type = "lv.x"), fof) 
    latExo <- setdiff(latExo, sof) 
  }
  
  # determine maximum path length for
  
  if (higher_order == FALSE){
    latEnd <- lavNames(fit, type = "lv.nox")
    latExo <- lavNames(fit, type = "lv.x") 
  }
  
  y1 <- na.omit( pt[pt$mat    == "lambda" & 
                    pt$lhs   %in% latEnd & 
                    pt$ustart == 1, 
                    "rhs"] )
  
  y2 <- na.omit( pt[pt$op  == "=~" &
                    pt$lhs %in% latEnd  & 
                    is.na(pt$ustart), 
                    "rhs"] )
  
  y2 <- na.omit( unique( c(y2, 
                           intersect(
                             lavNames(fit, type = "ov.nox"), 
                             lavNames(fit, type = "eqs.y")))))
                    
  
  x1 <- na.omit( pt[pt$mat  == "lambda" &
                    pt$lhs %in% latExo &
                    pt$ustart == 1,
                    "rhs"] )
  
  x2 <- na.omit( pt[pt$mat == "lambda" & 
                    pt$lhs %in% latExo & 
                    (pt$ustart != 1 | is.na(pt$ustart)), 
                    "rhs"] )

  y1names <- na.omit( pt[pt$mat == "lambda" &
                      pt$lhs %in% latEnd  & 
                      pt$ustart == 1, 
                      c("lhs", "rhs")] )
  
  colnames(y1names) <- c("lat", "obs")

  x1names <- na.omit( pt[pt$mat == "lambda" & 
                         pt$lhs %in% latExo & 
                         pt$ustart == 1,  
                         c("lhs", "rhs")] )
  
  colnames(x1names) <- c("lat", "obs")
  
  n1 <- n2 <- c()
  
  if (higher_order == TRUE){ # Eta2
   
    n1names <- na.omit( pt[pt$mat == "beta" & pt$lhs %in% sof
                    & pt$ustart == 1, c("lhs", "rhs")] )
    
    for (i in 1:nrow(n1names)){
      n1names$rhs[i] <- x1names[x1names$lat == n1names$rhs[i], "obs"]
    }
    
    colnames(n1names) <- c("lat", "obs")
    
    n1 <- lavNames(fit, type = "lv.x")
    n1 <- n1names[n1names$lat %in% n1, "lat"]
    
    n2names <- na.omit( pt[pt$mat == "beta" & pt$lhs %in% sof
                    & is.na(pt$ustart), c("rhs", "lhs")] )
    
    for (i in 1:nrow(n2names)){
      n2names$lhs[i] <- x1names[x1names$lat == n2names$rhs[i], "obs"]
    }
    colnames(n2names) <- c("lat", "obs")
    n2 <- n2names$obs
  }
  
  if (!is.null(inspect(fit)$lambda)) {LA <- inspect(fit)$lambda}
  if (!is.null(inspect(fit)$beta))   {BA <- inspect(fit)$beta}
  if (!is.null(inspect(fit)$theta))  {TH <- inspect(fit)$theta}
  if (!is.null(inspect(fit)$psi))    {PS <- inspect(fit)$psi}
  
  if (!is.null(inspect(fit)$lambda)) {class(LA) <- "matrix"}
  if (!is.null(inspect(fit)$beta))   {class(BA) <- "matrix"}
  if (!is.null(inspect(fit)$theta))  {class(TH) <- "matrix"}
  if (!is.null(inspect(fit)$psi))    {class(PS) <- "matrix"}
  
  if (!is.null(inspect(fit)$lambda)) {
    nz <- na.omit( pt[pt$op  == "=~" &
                      pt$lhs %in% latEnd  & 
                      is.na(pt$ustart), 
                      c("lhs", "rhs")] )
    if (nrow(nz) > 0){
      for (i in 1:nrow(nz)){ # i =1
        LA[nz[i,"rhs"], nz[i,"lhs"]] <- nrow(pt)*2 + i
      }
    }
  }
  
  
  if (!is.null(inspect(fit)$lambda)) {
    LY2 <- LA[,latEnd, drop = FALSE]
    LY2 <- LY2[!apply(LY2, 1, function(row) all(row ==0 )),, drop = FALSE]
    LX2 <- LA[,latExo, drop = FALSE]
    LX2 <- LX2[!apply(LX2, 1, function(row) all(row ==0 )),, drop = FALSE]
  }
  
  if (!is.null(inspect(fit)$beta)) {
    BA1 <- BA[latEnd,latEnd, drop = FALSE]
    GA1 <- BA[latEnd,latExo, drop = FALSE]
  }
  
  if (!is.null(inspect(fit)$theta)) {
    TE1 <- as.matrix(forceSymmetric(Matrix(TH), "L"))
    TE1 <- TE1[c(y1, y2), c(y1, y2), drop = FALSE] 
    TD1 <- TH
    TD1 <- TD1[c(x1, x2), c(x1, x2), drop = FALSE] 
  }
  
  if (!is.null(inspect(fit)$psi)) {
    PS1 <- as.matrix(forceSymmetric(Matrix(PS), "L"))
    dimnames(PS1) <- dimnames(PS)
    PS1 <- PS1[latEnd, latEnd, drop = FALSE] 
  }
  
  # if there are no latent variables in the model
  # change the psi matrix to theta
  if (length(lavNames(fit, type = "lv")) == 0) {
    TE1 <- as.matrix(forceSymmetric(Matrix(PS), "L"))
    dimnames(TE1) <- dimnames(PS)
    TE1 <- TE1[y2, y2, drop = FALSE] 
  }
  
  
  dv <- c(y1, y2, x2)
  
  if (higher_order == TRUE){ 
    
    dv <- c(n2, y1, y2, x2)
    
    ET2 <- BA[!apply(BA, 1, function(row) all(row ==0 )),sof, drop = FALSE]
    BA2 <- BA[latExo,latExo]
    dv2 <- fof
  }


  eqns <- list(DVobs = "", DVlat = "", IVobs = "", IVlat = "", EQtype = "",
               CD = "", TE = "", PIV = "", IV = "", W = "", NOTE = "")
  eqns <- replicate(length(dv), eqns, simplify = FALSE)
  

  for (i in 1:length(dv)){ 
    
    eqns[[i]]$DVobs <- dv[i]
    
        if (eqns[[i]]$DVobs %in% n2) {
            eqns[[i]]$EQtype <- "n2"
            eqns[[i]]$DVlat <- n2names[n2names$obs == dv[i], "lat"]
            
            tmp <- ET2[rownames(ET2) %in% eqns[[i]]$DVlat,, drop=FALSE]
            eqns[[i]]$IVlat <- colnames(tmp)[which(tmp[1,] != 0)]
            eqns[[i]]$IVobs <- n1names[n1names$lat == eqns[[i]]$IVlat, "obs"]
            
            for ( m in 1:length(eqns[[i]]$IVobs)){
              c2 <- eqns[[i]]$IVobs[m]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
            
            for ( m in 1:length(eqns[[i]]$DVobs)){
              c2 <- eqns[[i]]$DVobs[m]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
            

            for ( m in 1:length(eqns[[i]]$DVlat)){
              c2 <- eqns[[i]]$DVlat[m]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
      
            for ( m in 1:length(eqns[[i]]$DVlat)){
              c2 <- setdiff(latExo,n2names$lat)
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
           
        } 
    
    
        if (eqns[[i]]$DVobs %in% y1) {
            eqns[[i]]$EQtype <- "y1"
            eqns[[i]]$DVlat <- y1names[y1names$obs == dv[i], "lat"]
            
            tmp <- BA1[rownames(BA1) %in% eqns[[i]]$DVlat,, drop=FALSE]
            eqns[[i]]$IVlat <- colnames(tmp)[which(tmp[1,] != 0)]
            
            tmp <- GA1[which(rownames(GA1) %in% eqns[[i]]$DVlat),,drop=FALSE]
            tmp <- colnames(tmp)[which(tmp[1,] != 0)]
            eqns[[i]]$IVlat <- c(eqns[[i]]$IVlat,tmp)
            eqns[[i]]$IVlat <- eqns[[i]]$IVlat[eqns[[i]]$IVlat != ""]
            
            for ( m in 1:length(eqns[[i]]$IVlat)){
              z <- eqns[[i]]$IVlat[m]
                if (z %in% y1names$lat){
                  c2 <- y1names[y1names$lat == z, "obs"]
                }
                if (z %in% x1names$lat){
                  c2 <- x1names[x1names$lat == z, "obs"]
                }
              eqns[[i]]$IVobs <- c(eqns[[i]]$IVobs, c2)
              eqns[[i]]$IVobs <- eqns[[i]]$IVobs[eqns[[i]]$IVobs!=""]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
            
            for ( m in 1:length(eqns[[i]]$DVobs)){
              c2 <- eqns[[i]]$DVobs[m]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
            
            for ( m in 1:length(eqns[[i]]$DVlat)){
              c2 <- eqns[[i]]$DVlat[m]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
          } 
    
        if (eqns[[i]]$DVobs %in% y2) {
            eqns[[i]]$EQtype <- "y2"
            eqns[[i]]$DVlat <- NA
            
            tmp <- na.omit( pt[pt$op  == "=~" &
                               pt$rhs == eqns[[i]]$DVobs  &  
                               is.na(pt$ustart),  "lhs"] )
            
            tmp <- na.omit( unique( c(tmp, 
                                      pt[pt$op  == "~" &
                                         pt$lhs == eqns[[i]]$DVobs  &  
                                         is.na(pt$ustart),  "rhs"])))
            eqns[[i]]$IVlat <- tmp
            
            for ( m in 1:length(eqns[[i]]$IVlat)){ 
              z <- eqns[[i]]$IVlat[m]
              
               c2 <- ""
                
                if (z %in% y1names$lat){
                  c2 <- y1names[y1names$lat == z, "obs"]
                }
                if (z %in% x1names$lat){
                  c2 <- x1names[x1names$lat == z, "obs"]
                }
                
                if (!(z %in% y1names$lat) & !(z %in% x1names$lat)){
                  eqns[[i]]$IVobs <- c(eqns[[i]]$IVobs, z)
                }
                
              eqns[[i]]$IVobs <- c(eqns[[i]]$IVobs, c2)
              eqns[[i]]$IVobs <- eqns[[i]]$IVobs[eqns[[i]]$IVobs!=""]
              
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
            
            for ( m in 1:length(eqns[[i]]$DVobs)){
              c2 <- eqns[[i]]$DVobs[m]
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
            
        } 
    
        if (eqns[[i]]$DVobs %in% x2) {
            eqns[[i]]$EQtype <- "x2"
            eqns[[i]]$DVlat = NA
            
            tmp <- LX2[which(rownames(LX2) %in% eqns[[i]]$DVobs),,drop=FALSE]
            eqns[[i]]$IVlat <- colnames(tmp)[which(tmp[1,] != 0)]
            eqns[[i]]$IVobs <- x1names[x1names$lat == eqns[[i]]$IVlat, "obs"]
           
            for ( m in 1:length(eqns[[i]]$IVlat)){
              z <- eqns[[i]]$IVlat[m]
                if (z %in% y1names$lat){
                  z2 <- y1names[y1names$lat == z, "obs"]
                  c2 <- c(eqns[[i]]$DVobs, z2)
                }
                if (z %in% x1names$lat){
                  z2 <- x1names[x1names$lat == z, "obs"]
                  c2 <- c(eqns[[i]]$DVobs, z2)
                }
              eqns[[i]]$CD <- c(eqns[[i]]$CD, c2)
              eqns[[i]]$CD <- eqns[[i]]$CD[eqns[[i]]$CD!=""]
            }
        } 
    
  }


  if (higher_order == FALSE){
    
    if (nrow(LY2) < 1) {TE_on_y2 <- NULL}
    
    if (nrow(LY2) > 1) {
      
      
      if (is.null(inspect(fit)$beta)) {
        TE_on_y2 <- LY2 
      }
      if (!is.null(inspect(fit)$beta)) {
        TE_on_y1 <- solve(diag(nrow(BA1)) - BA1) 
        TE_on_y2 <- LY2 %*% TE_on_y1 
      }
    }
    
  }
  
  
  if (higher_order == TRUE){
    TE_on_n2 <- solve(diag(nrow(BA2)) - BA2)
  }

  for (i in 1:length(eqns)){ 
    if (eqns[[i]]$EQtype == "n2") {
      eqns[[i]]$TE <- eqns[[i]]$DVobs
      tmp <- TE_on_n2[which(rownames(TE_on_n2) %in% eqns[[i]]$DVlat),,drop=FALSE]
      t2  <- colnames(tmp)[which(tmp[1,] != 0)]
      eqns[[i]]$TE <- c(eqns[[i]]$TE, t2)
      eqns[[i]]$TE <- eqns[[i]]$TE[eqns[[i]]$TE!=""]
    } # end  n2
    
    if (eqns[[i]]$EQtype == "y1") { 
      eqns[[i]]$TE <- eqns[[i]]$DVobs
      tmp <- TE_on_y1[which(rownames(TE_on_y1) %in% eqns[[i]]$DVlat),,drop=FALSE]
      t2  <- colnames(tmp)[which(tmp[1,] != 0)]
      eqns[[i]]$TE <- c(eqns[[i]]$TE, t2)
      eqns[[i]]$TE <- eqns[[i]]$TE[eqns[[i]]$TE!=""]
    } # end y1
    
    if (eqns[[i]]$EQtype == "y2") { 
      eqns[[i]]$TE <- eqns[[i]]$DVobs
      tmp <- TE_on_y2[which(rownames(TE_on_y2) %in% eqns[[i]]$DVobs),,drop=FALSE]
      t2  <- colnames(tmp)[which(tmp[1,] != 0)]
      eqns[[i]]$TE <- c(eqns[[i]]$TE, t2)
      eqns[[i]]$TE <- eqns[[i]]$TE[eqns[[i]]$TE!=""]
    } # end y2
    
    # total effects on x2
    if (eqns[[i]]$EQtype == "x2") {
      eqns[[i]]$TE <- eqns[[i]]$DVobs

      if (higher_order == TRUE){
        eqns[[i]]$TE <- c(eqns[[i]]$TE, eqns[[i]]$IVlat)
      }
      
      eqns[[i]]$TE <- eqns[[i]]$TE[eqns[[i]]$TE!=""]
      
    } # end x2
  }
  
  if (higher_order == TRUE)  { add_length <- ncol(ET2)}
  if (higher_order == FALSE & !is.null(inspect(fit)$beta)) {
    add_length <- ncol(GA1)
  }  
  if (is.null(inspect(fit)$beta))   {add_length <- 0}
  if (nrow(LY2) < 1) {add_length <- 0}
  
  if (higher_order == FALSE & is.null(inspect(fit)$beta)) {
    add_length <- length(x1)
  }  
    
  effects <- list(DVobs = "", TE = "")
  effects <- replicate(add_length, effects, simplify = FALSE)
  if (length(effects) > 0 ) {
    for (i in 1:(length(effects))) {
      if (higher_order == TRUE)  { 
        effects[[i]]$DVobs <- n1names$obs[i]
        effects[[i]]$TE <- n1names$obs[i]
      }
      if (higher_order == FALSE)  { 
        effects[[i]]$DVobs <- x1[i]
        effects[[i]]$TE <- x1[i]
      }
    }
  }
  effects <- append(lapply(eqns, "[", c("DVobs", "TE")), effects)
   
  for (i in 1:length(eqns)) { 
    eqns[[i]]$PIV <- lavNames(fit, type = "ov.x")
    eqns[[i]]$PIV <- eqns[[i]]$PIV[eqns[[i]]$PIV != ""]
  }
  
  # add higher order observed to MIIVs
  # needs to be expanded for multiple higher order MIIVs
  if (higher_order == TRUE)  { 
    obs_sof <- n1names$obs
    for (k in 1:length(obs_sof)){
      for (i in 1:length(eqns)) { 
        if (!(obs_sof[k] %in%  eqns[[i]]$IVobs) & !(obs_sof[k] %in%  eqns[[i]]$CD)) {
          eqns[[i]]$PIV <- c(eqns[[i]]$PIV,  obs_sof[k])
          eqns[[i]]$PIV <- eqns[[i]]$PIV[eqns[[i]]$PIV != ""]
        }
      }
    }
  }

  for (i in 1:length(eqns)) { 
    C <- unlist(eqns[[i]]$CD)
      for (j in 1:length(effects)) { 
        E <- unlist(effects[[j]]$TE)
        if (!any(C %in% E)) {
            eqns[[i]]$PIV <- c(eqns[[i]]$PIV, effects[[j]]$DVobs)
        }
      }
    eqns[[i]]$PIV <- eqns[[i]]$PIV[eqns[[i]]$PIV != ""]
  }
  
  ## remove and PIVs which are predicted by the DV
  for (i in 1:length(eqns)) { 
    if (eqns[[i]]$DVobs %in% lavNames(fit, type = "eqs.y")){
      tdv <- eqns[[i]]$DVobs
      for (j in 1:length(eqns)) {
        if (tdv %in% eqns[[j]]$IVobs) {
          eqns[[i]]$PIV <- setdiff(eqns[[i]]$PIV, eqns[[j]]$DVobs)
        }
      }
    }
  }
  
  for (i in 1:length(eqns)) { 
    C <- eqns[[i]]$CD
       for (j in 1:length(C)) { 
          W <- C[j]
          t2 <- c()
 
          if (W %in% c(y1, y2)){
            tmp <- TE1[which(rownames(TE1) %in% W),,drop=FALSE]
            t2  <- colnames(tmp)[which(tmp[1,] != 0)]
          }
          
          if (W %in% c(x1, x2)){
            tmp <- TD1[which(rownames(TD1) %in% W),,drop=FALSE]
            t2  <- colnames(tmp)[which(tmp[1,] != 0)]
          }
          
          eqns[[i]]$W <- c(eqns[[i]]$W, t2)
          eqns[[i]]$W <- eqns[[i]]$W[eqns[[i]]$W != ""]
        }
      eqns[[i]]$IV <- setdiff(eqns[[i]]$PIV, eqns[[i]]$W)
  }
  
 
  for (i in 1:length(eqns)) {
    C <- eqns[[i]]$CD
       for (j in 1:length(C)) { 
          W  <- C[j]
          t2 <- c()
          if (W %in% c(latEnd)){
            
            tmp <- PS1[which(rownames(PS1) %in% W),,drop=FALSE]
            t2  <- colnames(tmp)[which(tmp[1,] != 0)]
            
            for (p in 1:length(eqns)) { 
              t3 <- c()
              tmp <- intersect(t2,eqns[[p]]$CD)
              if (length(tmp) > 0){ 
                t3 <- eqns[[p]]$DVobs
              }
            eqns[[i]]$W <- c(eqns[[i]]$W, t3) 
            }
          }
        
          eqns[[i]]$W <- eqns[[i]]$W[eqns[[i]]$W != ""]
        }
      eqns[[i]]$IV <- setdiff(eqns[[i]]$IV, eqns[[i]]$W)
      eqns[[i]]$IV <- eqns[[i]]$IV[eqns[[i]]$IV != ""]
  }

    for (i in 1:length(eqns)) { 
    C <- eqns[[i]]$CD
       for (j in 1:length(C)) { 
          W  <- C[j]
          t2 <- c()
          if (W %in% c(latEnd)){
            
            tmp <- PS1[which(rownames(PS1) %in% W),,drop=FALSE]
            t2  <- colnames(tmp)[which(tmp[1,] != 0)]
            
            for (p in 1:length(effects)) { 
              t3 <- c()
              tmp <- intersect(t2,effects[[p]]$TE)
              if (length(tmp) > 0){ 
                t3 <- effects[[p]]$DVobs
              }
            eqns[[i]]$W <- c(eqns[[i]]$W, t3) 
            }
          }
          
          eqns[[i]]$W <- eqns[[i]]$W[eqns[[i]]$W != ""]
        }
      eqns[[i]]$IV <- setdiff(eqns[[i]]$IV, eqns[[i]]$W)
      eqns[[i]]$IV <- eqns[[i]]$IV[eqns[[i]]$IV != ""]
    }
  
    ## identified all descendants of parent and removes
    ## them from IVs
    for (i in 1:length(eqns)) { # i = 7;
      parent  <- eqns[[i]]$DVobs
      parents <- parent
      for (j in 1:length(eqns)){ 
        if ( length(intersect(parents,eqns[[j]]$IVobs)) > 0 ){ 
          parents <- c(parents, eqns[[j]]$DVobs)
          parents <- parents[parents != ""]
        }
      }
      for (k in length(eqns):1){
        if ( length(intersect(parents,eqns[[k]]$IVobs)) > 0 ){ 
          parents <- c(parents, eqns[[k]]$DVobs)
          parents <- parents[parents != ""]
        }
      }
      eqns[[i]]$IV <- setdiff(eqns[[i]]$IV, parents)
      eqns[[i]]$IV <- eqns[[i]]$IV[eqns[[i]]$IV != ""]
    }
        
    ## remove any remaining parents
    for (i in 1:length(eqns)) {
          eqns[[i]]$IV <- setdiff(eqns[[i]]$IV, eqns[[i]]$DVobs)
          eqns[[i]]$IV <- eqns[[i]]$IV[eqns[[i]]$IV != ""]
    }
  
    eq_plabels <- pt[is.na(pt$ustart) & 
                     !is.na(pt$plabel) &
                     pt$mat %in% c("lambda","beta"), 
                     c("plabel")]
    fix_labels <- pt$label[pt$label != ""]
    
    rs <- pt[pt$op == "==", ]
    

    rs <- rs[!(rs$lhs %in% fix_labels & rs$rhs %in% fix_labels), ]

    constr <- list(DV = "", SET = "", FIX = "", NAME = "")
    constr <- replicate(nrow(rs), constr, simplify = FALSE)
    
    if (nrow(rs) > 0){
      for (i in 1:nrow(rs)){
        
        if (!rs$lhs[i] %in% fix_labels & !rs$rhs[i] %in% fix_labels){
          
          if (pt[pt$plabel == rs$lhs[i], "op"] == "=~" &
              pt[pt$plabel == rs$lhs[i], "op"] == "=~" )
            {
              left  <- pt[pt$plabel == rs$lhs[i], "lhs"]
              right <- pt[pt$plabel == rs$rhs[i], "lhs"]
          }
        
          if (pt[pt$plabel == rs$lhs[i], "op"] == "~" &
              pt[pt$plabel == rs$lhs[i], "op"] == "~" )
          {
            left  <- pt[pt$plabel == rs$lhs[i], "rhs"]
            right <- pt[pt$plabel == rs$rhs[i], "rhs"]
          }
            
          if (left %in% y1names$lat){
            left <- y1names[y1names$lat == left, "obs"]
          }
          
          if (left %in% x1names$lat){
            left <- x1names[x1names$lat == left, "obs"]
          }
          
          if (right %in% y1names$lat){
            right <- y1names[y1names$lat == right, "obs"]
          }
          
          if (right %in% x1names$lat){
            right <- x1names[x1names$lat == right, "obs"]
          }
          
          
          constr[[i]]$DV   <- c(left, right)
          
          if (pt[pt$plabel == rs$lhs[i], "op"] == "=~" &
              pt[pt$plabel == rs$lhs[i], "op"] == "=~" )
          {
            constr[[i]]$SET  <- c(pt[pt$plabel == rs$lhs[i], "rhs"],
                                  pt[pt$plabel == rs$rhs[i], "rhs"])
            constr[[i]]$NAME <- paste(constr[[i]]$SET[1],
                                      " = ",
                                      constr[[i]]$SET[2],
                                      sep = "")
          }
          
          if (pt[pt$plabel == rs$lhs[i], "op"] == "~" &
              pt[pt$plabel == rs$lhs[i], "op"] == "~" )
          {
            constr[[i]]$SET  <- c(pt[pt$plabel == rs$lhs[i], "lhs"],
                                  pt[pt$plabel == rs$rhs[i], "lhs"])
            constr[[i]]$NAME <- paste(constr[[i]]$DV[1],
                                      " = ",
                                      constr[[i]]$DV[2],
                                      sep = "")
          }
          constr[[i]]$FIX  <- 0

        }
        
        if (rs$lhs[i] %in% fix_labels | rs$rhs[i] %in% fix_labels){
          
         constr[[i]]$DV   <- pt[pt$label == rs$lhs[i], "lhs"]
         
         if (constr[[i]]$DV  %in% y1names$lat){
            constr[[i]]$DV  <- y1names[y1names$lat == constr[[i]]$DV , "obs"]
          }
          
          if (constr[[i]]$DV  %in% x1names$lat){
            constr[[i]]$DV  <- x1names[x1names$lat == constr[[i]]$DV , "obs"]
          }
         
         constr[[i]]$SET  <- pt[pt$label == rs$lhs[i], "rhs"]
         constr[[i]]$FIX  <- as.numeric(rs[i, "rhs"])
         constr[[i]]$NAME <- paste(constr[[i]]$SET,
                                    " = ",
                                    constr[[i]]$FIX,
                                    sep = "")
        }
      }
    }
    
    for (i in 1:length(eqns)){
      LHS <- paste(eqns[[i]]$DVobs, collapse = ", ")
      RHS <- paste(eqns[[i]]$IVobs, collapse = ", ")
      Instruments <- paste(eqns[[i]]$IV, collapse = ", ")
      Disturbance <- paste("e.",eqns[[i]]$CD, collapse = ", ", sep="")
      modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
      colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
      if (i == 1) {modeqns <- modtemp }
      if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
    } 
  
  search <- list(eqns = eqns, df = modeqns, constr = constr, miivs.out = miivs.out)
  class(search) <- "miivs"
  search
}
