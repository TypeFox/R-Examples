isotab <-
function (ip, level = 1, phi.min = "auto", p.max = .05)
{ 
  if (class (ip) != 'isopam') stop ('Object does not seem to result from Isopam')  
  IO <- ip$dat
  IO [IO > 0] <- 1 
  N <- nrow (IO)
  SP <- ncol (IO)
  frq <- t (as.matrix (colSums (IO)))
  
  if (is.null (ip$hier)) 
  {
    if (level > 1) print ("No hierarchy levels available", quote = FALSE)
    tab <- t (aggregate (IO, by = list (ip$flat), FUN = sum))
  }
  else
  { 
    depth <- ncol (ip$hier)
    if (level > depth) 
    {
      level <- depth 
      print (paste ("Switching to lowest level", depth), 
        quote = FALSE)
    }
    tab <- t (aggregate (IO, by = list (ip$hier [,level]), FUN = sum))
  }
  
  ## Fisher's exact test for 2x2 tables
  fshtest <- function (x)
  {
    PVAL <- NULL
    m <- sum(x[, 1])
    n <- sum(x[, 2])
    k <- sum(x[1, ])
    x <- x[1, 1]
    lo <- max(0, k - n)
    hi <- min(k, m)
    support <- lo:hi
    logdc <- dhyper (support, m, n, k, log = TRUE)

    dnhyper <- function(ncp)
    {
      d <- logdc + log(ncp) * support
      d <- exp(d - max(d))
      d / sum(d)
    }

    pnhyper <- function(q, ncp = 1, upper.tail = FALSE)
    {
      if (ncp == 1)
      {
          if (upper.tail)
            return(phyper(x - 1, m, n, k, lower.tail = FALSE))
          else return(phyper(x, m, n, k))
      }
      if (ncp == 0)
      {
          if (upper.tail)
            return(as.numeric(q <= lo))
          else return(as.numeric(q >= lo))
      }
      if (ncp == Inf)
      {
          if (upper.tail)
            return(as.numeric(q <= hi))
          else return(as.numeric(q >= hi))
      }
      d <- dnhyper(ncp)
      if (upper.tail)
          sum(d[support >= q])
      else sum(d[support <= q])
    }
    PVAL <- switch ("two.sided", less = pnhyper(x, 1), greater = pnhyper (x, 1,
        upper.tail = TRUE), two.sided = {
          relErr <- 1 + 10^(-7)
          d <- dnhyper(1)
          sum(d [d <= d [x - lo + 1] * relErr])})
    return (PVAL)
  }
 
  ## 1) Contingency table  
  cnam <- tab [1,]
  tab <- tab [-1,]
  rnam <- rownames (tab)
  tab <- matrix (as.numeric (tab), nrow = nrow (tab))
  colnames (tab) <- cnam
  rownames (tab) <- rnam
  nc <- ncol(tab)
  
  ## 2) Calculate phi.min threshold
  if (phi.min == "auto") 
    phi.min <- round (0.483709 + nc * -0.003272 + N * -0.000489 + 
      SP * 0.000384 + sqrt (nc) * -0.01475, 2)   
  
  ## 3) Frequency table
  if (is.null (ip$hier)) siz <- table (ip$flat) ## Cluster sizes
  else siz <- table (ip$hier [,level]) ## Cluster sizes
  spc <- t (as.matrix (siz))[rep (1, nrow (tab)),]
  frq.2 <- tab / spc ## Frequency
  frq.2 <- round (frq.2 * 100, 0) ## as percentage
  
  ## 4) Fisher's significance table  
  ft <- tab
  for (fsp in 1:SP)         ## fsp-loop through species
  {                    
    spec_frq <- frq [fsp]
    spec_io <- IO [,fsp]
    
    for (fcl in 1:nc)       ## fcl-loop through clusters
    {                  
      Nj <- siz [fcl]
      insd <- tab [fsp, fcl]
      absci <- Nj - insd
      outs <- spec_frq - insd              ## occ. outside
      absco <- N - Nj - outs               ## abs. outside
      
      fshm <- matrix (c(insd, absci, outs, absco), 2, 2)    
      ft [fsp,fcl] <- fshtest (fshm)        
    }
  }      

  ## 5) Significance symbols
  ft.symb <- ft
  ft.symb [ft > 0.05] <- ""
  ft.symb [ft <= 0.05] <- "*"
  ft.symb [ft <= 0.01] <- "**"
  ft.symb [ft <= 0.001] <- "***"

  ## 6) Combined frequency table with significance symbols 
  frq.ft <- matrix (paste (frq.2, ft.symb, sep = ""), 
      nrow(frq.2), ncol(frq.2))
  frq.ft <- data.frame (frq.ft)
  colnames (frq.ft) <- colnames (ft.symb)
  rownames (frq.ft) <- rownames (ft.symb)

  ## 7) Standardized phi table 
  S <- 1 / nc                                   ## Constant s (Tichy '06)
  cs <- S * N                                   ## new cluster sizes
  phi <- tab

  for (i in 1:SP) 
  {
    for (j in 1:nc) 
    {
      insd <- tab [i, j]                        ## original n in cluster j
      outs <- sum (tab [i,-j])                  ## original n outside cluster j
      oc <- cs * (insd / siz [j])               ## new n in cluster j
      on <- (N - cs) * (outs / (N - siz [j]))   ## new n outside cluster j
      total <- oc + on                          ## new total value
      phi.1 <- nv <- (N * oc - total * cs) 
      phi.2 <- sqrt (total * cs * (N - total) * (N - cs))            
      nv <- phi.1 / phi.2
      phi [i,j] <- nv
    }
  }
  phi [is.na (phi)] <- 0  ## Replace NaN-values by 0

  ## Table sorting     
  ## 8) .... by phi and frequency
  phi.idx <- apply (phi, 1, which.max) ## Group affiliation by phi 
  frq.ord <- phi.idx  
  for (i in 1:length(frq.ord)) frq.ord [i] <- frq.2 [i, phi.idx [i]] 
  frq.top <- t (frq) [order (phi.idx, -frq.ord),] ## Sorting
  ord.top <- names (frq.top)
  frq.ft.top <- frq.ft [ord.top,]
  ft <- ft [ord.top,]
  phi <- phi [ord.top,]
  
  ## 9) Filter diagnostic species
  filter1 <- apply (ft, 1, min) <= p.max                                          
  filter2 <- apply (phi, 1, max) >= phi.min
  dia <- which (filter1 [filter2 == TRUE] == TRUE) ## diagnostic species 
  n.dia <- length (dia) ## how many diagnostic species
  if (n.dia == 0) diag <- "No diagnostic species with given thresholds." 
  if (n.dia > 0) diag <- frq.ft.top [names (dia),]
  
  ## 10) For later use in the bottom part of the tables
  ord.bot <- names (t (frq) [order (-frq),])
  frq.ft.b <- frq.ft [ord.bot,]

  ## 11) Move diagnostic species to top
  if (n.dia > 0) 
  {
    FRQ <- rbind (diag, frq.ft.b [rownames (frq.ft.b) %in% 
      rownames (diag) == FALSE,])
  }
  else FRQ <- frq.ft.b
  
  ## 12) Report cluster sizes
  siz <- t (as.matrix (siz))
  rownames (siz) <- "n"
  
  ## 13) Parameters
  param <- c(phi.min, p.max)
  names (param) <- c("phi.min", "p.max")

  ## 14) Print info about diagnostic species
  dig1 <- phi.idx [names (phi.idx) %in% names (dia)]
  dig2 <- dig1 [rownames (diag)]
  typ <- list ()
  for (i in 1:nc) 
  {
    if (length (names (dig2) [dig2 == i]) > 0)
      typ [i] <- paste (names (dig2) [dig2 == i], collapse = ', ')
    else typ [i] <- 'Nothing particularly typical'
  }
  names (typ) <- cnam
      
  ## 15) Output
  isotab.out <- list (
   tab = FRQ,  
   n = siz,
   thresholds = param,   
   typical = typ)
  
  ## 16) Output to screen
  
  cat ('$tab', fill = TRUE)
  print (FRQ)
  cat ('', fill = TRUE)
  cat ('$n', fill = TRUE)
  print (siz, quote = FALSE)
  cat ('', fill = TRUE)
  cat ('$thresholds', fill = TRUE)
  print (param, quote = FALSE)
  cat ('', fill = TRUE)
  cat ('$typical', fill = TRUE)
  
  for (i in 1:nc) 
  {
    cat (paste ('Typically found in ', cnam [i], ": ", sep =''), fill = TRUE)
    if (length (names (dig2) [dig2 == i]) > 0)
      cat (paste (names (dig2) [dig2 == i], collapse = ', '), fill = TRUE)
    else cat ('Nothing particularly typical', fill = TRUE)
    cat ('', fill = TRUE)                            
  }  

  invisible (isotab.out)
}

