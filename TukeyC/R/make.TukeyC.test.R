##
## Function to perform Tukey test
##

make.TukeyC.test <- function(r=r,
                             MSE=MSE,
                             m.inf=m.inf,
                             ord=ord,
                             sig.level=sig.level,
                             dfr=dfr,
                             bal=bal,
                             mt=mt,
                             round=round)
{
  if (length(r) < nrow(m.inf))  # expand n to the correct length if necessary
    r <- rep.int(r,
                 nrow(m.inf))

  r <- r[ord] 

  m.tmp <- m.inf[, 1]

  names(m.tmp) <- r

  vece <- outer(X=m.tmp,
                Y=m.tmp,  # (v)ariance (e)stimation of (c)ontrast (e)stimation
                function(X, Y) MSE * (1/as.numeric(names(X)) + (1/as.numeric(names(Y)))))

  qTukey <- qtukey(p=sig.level,
                   nmeans=nrow(m.inf),
                   df=dfr,
                   lower.tail=FALSE)

  if (!bal) {
    msd <- qTukey * sqrt(1/2 * vece)  # minimum significative difference

    diag(msd) <- 0

    dimnames(msd) <- list(rownames(m.inf),
                          rownames(m.inf))
  }
  else
    msd <- qTukey * sqrt(1/2 * vece)[1,1]

  m    <- m.inf[,1]

  difm <- abs(outer(m,
                    m,
                    "-"))  # means difference
  dif  <- difm - msd

  dif  <- ifelse(dif <= 0,
                 FALSE,
                 TRUE)

  res  <- make.TukeyC.groups(dif)

  res  <- cbind(format(round(m,
                             round),
                       nsmall=2),
                res)

  colnames(res) <- c('Means',
                     paste('G',
                           1:(ncol(res) - 1),
                           sep=''))

  if (bal) r <- r[1]

  # The below estimates the probability of observed difference betweeen means be significative
  # Matrix of the difference of means above diagonal and respective p-values of the Tukey test
  # below diagonal 
  difm[lower.tri(difm)] <- ptukey(q=difm[lower.tri(difm)] / sqrt(1/2 * vece[lower.tri(vece)]),
                                  nmeans=nrow(m.inf),
                                  df=dfr,
                                  lower.tail=FALSE)
  diag(difm) <- 0  # To be sure!

  out <- list(Table      = mt,
              Means      = m.inf,
              Result     = as.data.frame(res),
              Sig.level  = sig.level,
              Diff_Prob  = round(difm, 3),
              MSD        = round(msd, 3),
              Replicates = r)

  return(out)
}
