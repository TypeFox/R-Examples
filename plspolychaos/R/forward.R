###############################################
# Functions for the option "forward'
###############################################

#######################################################
regsPC <- function(xin, dmax, yobs, struc) {
  #############################################
  # Compute the coefficient R2 of the simple regression
  # after X Legendre coding
  #############################################
  # INPUT
  # xin: matrix (nl x nvx) (number of rows x number of inputs)
  #  Calibrated inputs
  # dmax : scalar. Polynomial degree
  # yobs: vector (nl). The response
  # struc: matrix nmono+1 x nvx. Structure of the polynomial
  # RETURN
  # r2: vector (nmono-nvx). R2 of the simple regression
  # pour les monomes non unitaires

  nl <- nrow(xin)
  nvx <- ncol(xin)
  nmono <- nrow(struc) # nombre total de monomes+1

  # pas besoin de calculer le R2 des nvx+1 premiers monomes
  # on les garde toujours
  nr2 <- nmono - (1+nvx) # nombre de monomes a selectionner
  r2 <- rep(NA, nr2)
  ir2 <- 1
  ret <- rep(-99 , nl)
  
  for (imono in (2+nvx):nmono) {
    ## X Legendre coding and column multiplication
    xjt <- .C("Cpolleg2", as.double(xin), as.integer( struc[imono,]), as.integer(nl), as.integer(nvx),ret=as.double(ret))$ret

     # simple regression
    r <- lm(yobs~xjt)
   r2[ir2] <- summary(r)$r.squared
    ir2 <- ir2 + 1
  } # fin imono

  return(r2)
} # fin regsPC

#######################################################
## option forward
#######################################################
selexPC <- function(xin, dmax, Y, struc, forward) {
  #############################################
  # option forward (main function)
  #############################################
  # INPUT
  # xin: matrix (nl x nvx) (number of rows x number of inputs)
  #  Calibrated inputs
  # dmax : scalar. Polynomial degree
  # Y: vector (nl). The response
  # struc: matrix nmono+1 x nvx. Structure of the polynomial
  # forward: required number of monomials, not included
  # the constant term but with the unit monomials (forward>=0)
  # RETURN
  # the required number of monomials or NULL
  # and an object PCEpoly or NULL

  nvx <- ncol(xin)
  nl <- nrow(xin)
  nmono <- nrow(struc)
  nfo <- nmono - 1 # valeur max de forward
  if ( (forward < nvx) || (forward >= nfo)) {
    # on ne met pas le message en warning, car ceux-ci, par defaut,
    # n'apparaissent que quand le programme est termine
    # ce qui peut paraitre anormalement long a l'utilisateur
    # alors meme qu'il a utilise l'option forward
        cat(paste("The value of the option 'forward' should be greater or equal to the number of inputs,",
                      nvx,
                      "\n  and less than the total number of monomials,",
                   nfo,
                      "\n  Option 'forward' ignored.\n"))
        return(list(forward=NULL, object=NULL))
       } # fin warning
 

  
  # calcul des R2 de chaque monome (avec regression simple)
  r2 <- regsPC(xin, dmax, Y, struc)
  # r2 contient les R2 sans inclure  les nvx monomes unitaires
  r2tri <- sort(r2, decreasing = TRUE, index.return = TRUE)
  
  # on garde toujours le terme constant et les nvx monomes unitaires
  # sont comptes dans forward
  nmonosel <-1+forward 
  struc2 <- matrix(NA, nrow= nmonosel, ncol=nvx)
  struc2[1:(1+nvx),] <- struc[1:(1+nvx),]
  labelmono <- rownames(struc)[1:(1+nvx)]

  # on rajoute les monomes selectionnes a la matrice struc
  forward <- forward - nvx # nombre de monomes a selectionner
  
  if (forward >0) {
    nselect <- r2tri$ix[1:forward] # no des monomes selectionnes
    struc2[(2+nvx):nmonosel, ] <-  struc[(1+nvx+nselect),]
    labelmono <-c(labelmono, rownames(struc)[1+nvx+nselect])
  } # fin forward

                
  # calculer la matrice de Legendre des monomes selectionnes

  trav <- rep(99, nvx)
  XM <- rep(-99,  nl* nmonosel)
  XM <- .C("TCpolleg2", as.double(xin), as.integer(struc2),
           as.integer(nmonosel),
           as.integer(nl), as.integer(nvx),
           as.integer(trav),
           ret=as.double(XM))$ret
  XM <- matrix(XM,  nl, nmonosel)
  
  XMY <- cbind(XM, Y)
  dimnames(struc2) <- list(labelmono, colnames(struc))
  struc2 <- new("PCEdesign", .Data=struc2, degree=dmax, total.nmono=nfo)
  retour <- new("PCEpoly", .Data = XMY, STRUC = struc2,
                 nvx = nvx, call = match.call())
    return(list(forward=forward, object=retour))
} # fin selexPC


    
