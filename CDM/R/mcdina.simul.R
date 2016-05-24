

#*******************************************************************
# simulation mcdina model
# input:
# alpha vectors
# pars_lc: Q-matrix, probabilities, ... , see below
# pars_lr: see below
#*******************************************************************

# sirt package is needed for loading this function

#################################################################
simul.mcdina <- function( alpha ,  pars_lc , pars_lr , skillcl ){
    # skills ... alpha vectors
    #   skillcl <- scan.vec( "P000 P100 P010 P110 P001 P101 P011 P111" )		
	base::requireNamespace("sirt")
    skills <- alpha 
    N <- length(alpha)
    I <- max( pars_lc$item )
    CC <- max( pars_lc$cats )
    dat <- matrix( NA , nrow=N , ncol=I )
    colnames(dat) <- paste0("I",1:I)
    # calculate probabilities and simulate
    for (ii in 1:I){
        # ii <- 4
        lc.ii <- pars_lc[ pars_lc$item == ii , ]
        lc.ii <- lc.ii[ lc.ii$sum == 1 , ]
        lr.ii <- pars_lr[ pars_lr$item == ii , ]
        lr.unique <- paste( unique( lr.ii$lr ) )
        # compute latent response pattern for item ii
        lr.ii <- paste(lr.ii[ match( skillcl[ skills ] , lr.ii$skillclass ) , "lr" ])
        probs <- lc.ii[ match( lr.ii , paste(lc.ii$lr) ) , grep( "Cat" , colnames(pars_lc ) ) ]
        Nc <- ncol(probs)
        rn <- stats::runif(N)	
        # probs1 <- rowCumsums.sirt(matr=as.matrix(probs)  )
		# eval( parse ( text = "probs1 <- sirt::rowCumsums.sirt(matr=as.matrix(probs)  )" ) )
		probs1 <- sirt::rowCumsums.sirt(matr=as.matrix(probs)  )
        # dat[,ii] <- rowIntervalIndex.sirt(matr= probs1 ,rn) 
		# eval( parse ( text = "dat[,ii] <- sirt::rowIntervalIndex.sirt(matr= probs1 ,rn) " ) )
		dat[,ii] <- sirt::rowIntervalIndex.sirt(matr= probs1 ,rn)
        print(paste0( "Item " ,ii )) ; utils::flush.console()
                }
      return(dat)
      }
#################################################################  

# Examples:

##   > pars_lc
##      item cats    lr max.cat lr_index    Q Cat0 Cat1 Cat2 Cat3 sum cat
##   1     1    1   LR0       0        1 Q000 0.80 0.10 0.05 0.05   1   1
##   2     1    2   LR1       0        2 Q100 0.10 0.60 0.10 0.20   1   2
##   3     1    3   LR2       0        3 Q010 0.25 0.10 0.60 0.05   1   3
##   4     1    4   LR3       1        4 Q110 0.02 0.02 0.16 0.80   1   4
##   5     2    1   LR0       0        1 Q000 0.70 0.10 0.15 0.05   1   1
##   6     2    2  LR12       0        2 Q100 0.10 0.35 0.50 0.05   1   2
##   7     2    3  LR12       0        2 Q100 0.00 0.00 0.00 0.00   0   3
##   8     2    4   LR3       1        3 Q110 0.05 0.10 0.20 0.65   1   4
##   9     3    1  LR01       0        1 Q000 0.25 0.55 0.10 0.10   1   1
##   10    3    2  LR01       0        1 Q000 0.00 0.00 0.00 0.00   0   2
##   11    3    3   LR2       1        2 Q110 0.04 0.01 0.85 0.10   1   3
##   12    3    4   LR3       0        3 Q100 0.10 0.10 0.20 0.60   1   4
##  [...]

##   > pars_lr
##       item skillclass skillclass_index    lr lr_index
##   1      1       P000                1   LR0        1
##   2      1       P100                2   LR1        2
##   3      1       P010                3   LR2        3
##   4      1       P110                4   LR3        4
##   5      1       P001                5   LR0        1
##   6      1       P101                6   LR1        2
##   7      1       P011                7   LR2        3
##   8      1       P111                8   LR3        4
##   9      2       P000                1   LR0        1
##   10     2       P100                2  LR12        2
##   11     2       P010                3   LR0        1
##   12     2       P110                4   LR3        3
##   13     2       P001                5   LR0        1
##   14     2       P101                6  LR12        2
##   15     2       P011                7   LR0        1
##   16     2       P111                8   LR3        3
##   17     3       P000                1  LR01        1
##   18     3       P100                2   LR3        3
##   19     3       P010                3  LR01        1
##   20     3       P110                4   LR2        2
##   21     3       P001                5  LR01        1
##   22     3       P101                6   LR3        3
##   23     3       P011                7  LR01        1
##   [...]

## pars_lr is created by .mcdina.prep.test.latent.response.
## The core of pars_lc is also created by this function, only
## probabilities must be specified by the user.