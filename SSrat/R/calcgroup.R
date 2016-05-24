calcgroup <- function(schoolid = 1, groupid=1, dataframe, scalelength = c(5, 3, 7, 9), 
                      alpha = c(0.10, 0.05, 0.01), NBcriteria = F, printresult = F) {
  
  getsapi <- function(schoolid = 1, groupid=1, dataframe, scalelength , NBcriteria = FALSE) {
    #sc = as.numeric(match.arg(as.character(scalelength)))  #sc=7;schoolid=1 ;groupid=16
    sc=as.numeric(scalelength)
    sel = dataframe$schoolid == schoolid & dataframe$groupid == groupid  #dataframe=DF; school=3; group=2
    lastcol = ncol(dataframe)
    # count columns with name equal rn
    maxrated = sum(grepl("^r\\d", names(dataframe)))
    # ratings may be padded with NA
    maxnrrated = maxrated - 1
    # get ratings
    grr = dataframe[sel, (lastcol - maxnrrated):lastcol]
    lastrow = nrow(grr)
    if (lastrow == 0) {
      stop("Empty data selection. Please check school and group identification.")
    }

    grr=grr[ , ! apply( grr , 2 , function(x) all(is.na(x)) ) ] #remove cols with only NA
    rownames(grr) <- dataframe[sel, ]$respid
    
    if ((max(grr, na.rm = T) > sc) || (min(grr, na.rm = T) < 1) || (ncol(grr) == 
                                                                      0)) {
      stop("Invalid ratings found in dataframe")
    }
    nrresp = nrow(grr)
    if (nrresp <= 1) 
      stop(cat("Number of assessors (", nrresp, ") in this group <= 1.\n"))
    
    R = (sc + 1)/2
    
    P = grr
    nrAssessed <- apply(P, 1, function(x) {
      sum(!is.na(x))
    })  # nr ratings given
    nrAssessors <- apply(P, 2, function(x) {
      sum(!is.na(x))
    })  # nr ratings received
    
    if (NBcriteria) {
      P[(P > 1) & (P < sc)] <- 2
      P[P == sc] <- 3
      sc = 3
      R = 2
    }
    
    S = P - R
    S[S < 0] <- 0
    
    A = P - R
    A[A > 0] <- 0
    A = abs(A)
    
    I = A + S
    
    list(names = dataframe[sel, ]$names, S = S, A = A, P = P, I = I, 
         nrAssessed = nrAssessed, nrAssessors = nrAssessors)
  }  #end getsapi
  
  assessors <- function(subject) {
    if (subject <= ncol(SAPI$P)) {
      which(!is.na(SAPI$P[, subject]))
    } else NULL
  }
  # subject=1
  probtotrating <- function(subject, probmatrix) {
    a = assessors(subject)
    pr = probmatrix[a[1], ]
    for (i in a[2:length(a)]) {
      pr = convolve(pr, rev(probmatrix[i, ]), type = "o")
    }
    return(pr)
  }
  
 
#start
  sslabel = c("Popular", "Rejected", "Neglected", "Controversial", 
              "Average")
  sc = match.arg(as.character(scalelength), choices=c(3,5,7,9)) #sc=7 
  R = (sc + 1)/2
  if (!all(alpha < 1) || !all(alpha > 0)) 
    warning("Invalid number in alpha vector")
  sel = dataframe$schoolid == schoolid & dataframe$groupid == groupid
  SAPI = getsapi(schoolid, groupid, dataframe[sel, ], sc, NBcriteria=NBcriteria)  #dataframe=dataframe[sel,]
  
  # get probabilities Sp, Ap, Pp, Ip
  nrAors = length(SAPI$nrAssessed)
  Sp = matrix(0, nrAors, R)  #i=0
  for (i in 0:(R - 1)) {
    Sp[, i + 1] = apply(SAPI$S, 1, function(x) {
      mean(x == i, na.rm = T)
    })
  }
  Ap = matrix(0, nrAors, R)  #i=0
  for (i in 0:(R - 1)) {
    Ap[, i + 1] = apply(SAPI$A, 1, function(x) {
      mean(x == i, na.rm = T)
    })
  }
  
  Pp = matrix(0, nrAors, sc + 1)
  for (i in 0:sc) {
    Pp[, i + 1] = apply(SAPI$P, 1, function(x) {
      mean(x == i, na.rm = T)
    })
  }

Ip = matrix(0, nrAors, R)  #i=0
  for (i in 0:(R - 1)) {
    Ip[, i + 1] = apply(SAPI$I, 1, function(x) {
      mean(x == i, na.rm = T)
    })
  }
  
  # get total scores
  tr.S = colSums(SAPI$S, na.rm = T)
  tr.A = colSums(SAPI$A, na.rm = T)
  tr.P = colSums(SAPI$P, na.rm = T)
  tr.I = colSums(SAPI$I, na.rm = T)
  hscores = matrix(c(0:(R - 1)), ncol = 1, nrow = R)
  scores = matrix(c(0:(sc)), ncol = 1, nrow = sc + 1)
  
  # imr = structure for intermediate results: probabilty leftsided
  # testing (pl), probability rightsided testing (pr) and expected
  # value (e) for S, A, P and I
  nrAed = length(SAPI$nrAssessors)
  imr = matrix(0, nrow = nrAed, ncol = 12)
  colnames(imr) = c("pls", "prs", "es", "pla", "pra", "ea", "plp", 
                    "prp", "ep", "pli", "pri", "ei")
  
  status = matrix(0, nrow = nrAed, ncol = length(alpha))
  colnames(status) <- 1:dim(status)[2]
  
  # i=7; i=8; SAPI$assessed
  for (i in (1:nrAed)) {
    ass = assessors(i)
    if (is.null(ass)) {
      break
    }
    # sympathy
    cptr.Sp = cumsum(probtotrating(i, Sp))
    imr[i, "pls"] = cptr.Sp[tr.S[i] + 1]
    imr[i, "prs"] = ifelse(tr.S[i] == 0, 1, 1 - cptr.Sp[tr.S[i]])
    imr[i, "es"] = sum(Sp[ass, ] %*% hscores)
    
    # antipathy
    cptr.Ap = cumsum(probtotrating(i, Ap))
    imr[i, "pla"] = cptr.Ap[tr.A[i] + 1]
    imr[i, "pra"] = ifelse(tr.A[i] == 0, 1, 1 - cptr.Ap[tr.A[i]])
    imr[i, "ea"] = sum(Ap[ass, ] %*% hscores)
    
    # social preference
    cptr.Pp = cumsum(probtotrating(i, Pp))
    imr[i, "plp"] = cptr.Pp[tr.P[i] + 1]
    imr[i, "prp"] = ifelse(tr.P[i] == 0, 1, 1 - cptr.Pp[tr.P[i]])
    imr[i, "ep"] = sum(Pp[ass, ] %*% scores)
    
    # social impact
    cptr.Ip = cumsum(probtotrating(i, Ip))
    imr[i, "pli"] = cptr.Ip[tr.I[i] + 1]
    imr[i, "pri"] = ifelse(tr.I[i] == 0, 1, 1 - cptr.Ip[tr.I[i]])
    imr[i, "ei"] = sum(Ip[ass, ] %*% hscores)
    
    # j=2
    colnames(status) <- paste("SS", substr(sprintf("%.2f", alpha), 
                                           2, 4), sep = "")
    for (j in 1:length(alpha)) {
      if (NBcriteria) {
        status[i, j] = sslabel[5]
        if ((imr[i, "prs"] <= alpha[j]) & (tr.A[i] < imr[i, "ea"])) {
          status[i, j] = sslabel[1]
        } else if ((imr[i, "pra"] <= alpha[j]) & (tr.S[i] < imr[i, 
                                                                "es"])) {
          status[i, j] = sslabel[2]
        } else if (((imr[i, "prs"] <= alpha[j]) & (tr.A[i] > imr[i, 
                                                                 "ea"])) | ((imr[i, "pra"] <= alpha[j]) & (tr.S[i] > 
                                                                                                             imr[i, "es"]))) {
          status[i, j] = sslabel[4]
        } else if (imr[i, "pli"] <= alpha[j]) {
          status[i, j] = sslabel[3]
        }
      } else {
        status[i, j] = sslabel[5]
        if ((imr[i, "prp"] <= alpha[j]) & (tr.S[i] > imr[i, "es"]) & 
              (tr.A[i] < imr[i, "ea"])) {
          status[i, j] = sslabel[1]
        } else if ((imr[i, "plp"] <= alpha[j]) & (tr.S[i] < imr[i, 
                                                                "es"]) & (tr.A[i] > imr[i, "ea"])) {
          status[i, j] = sslabel[2]
        } else if ((imr[i, "pri"] <= alpha[j]) & (tr.S[i] > imr[i, 
                                                                "es"]) & (tr.A[i] > imr[i, "ea"])) {
          status[i, j] = sslabel[4]
        } else if ((imr[i, "pli"] <= alpha[j]) & (tr.S[i] < imr[i, 
                                                                "es"]) & (tr.A[i] < imr[i, "ea"])) {
          status[i, j] = sslabel[3]
        }
      }
    }
  }
  schoolid = rep(schoolid,nrAed) #create vectors
  groupid = rep(groupid, nrAed)
  respid = as.numeric(gsub("r",'',colnames(SAPI$P)))
  nrAss = SAPI$nrAssessors
  out = data.frame(schoolid, groupid, respid, nrAss, tr.S, 
                            tr.A, tr.P, tr.I, status, stringsAsFactors=FALSE)
  if ("resplabel" %in% colnames(dataframe[sel, ])) {
    resplabel = dataframe[sel, ]$resplabel
    out = cbind(resplabel, out)
  }
  row.names(out) <- NULL
  outlist = list(dataframe = out, S = SAPI$S, A = SAPI$A, P = SAPI$P, 
                 I = SAPI$I, intermediate = imr)
  if (printresult) {
    print(out)
  }
  invisible(outlist)
  
}
