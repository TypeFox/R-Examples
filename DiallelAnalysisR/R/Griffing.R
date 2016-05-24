#' @name   Griffing
#' @aliases Griffing
#' @title Diallel Analysis using Griffing Approach
#' @description \code{Griffing} is used for performing Diallel Analysis using Griffing's Approach.
#' @param y Numeric Response Vector
#' @param Rep  Replicate as factor
#' @param Cross1 Cross 1 as factor
#' @param Cross2 Cross 2 as factor
#' @param data A \code{data.frame}
#' @param Method Method for Diallel Analysis using Griffing's approach.
#'        It can take \strong{1}, \strong{2}, \strong{3}, or \strong{4}
#'         as argument depending on the method being used.
#'        \enumerate{
#'        \item Method-I (Parents + \eqn{F_{1}}'s + reciprocals);
#'        \item Method-II (Parents and one set of \eqn{F_{1}}'s);
#'        \item Method-III (One set of \eqn{F_{1}}'s and reciprocals);
#'        \item Method-IV (One set of \eqn{F_{1}}'s only).
#'        }
#' @param Model Model for Diallel Analysis using Griffing's approach.
#'        It can take \strong{1} or \strong{2} as arguments depending on the model being used.
#'        \enumerate{
#'        \item Fixed  Effects Model;
#'        \item Random Effects Model.
#'        }
#' @details Diallel Analysis using Griffing's approach.
#' @return Means Means
#' @return ANOVA Analysis of Variance (ANOVA) table
#' @return Genetic.Components Genetic Components
#' @return Effects Effects of Crosses
#' @return StdErr Standard Errors of Crosses
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @references \enumerate{
#' \item Griffing, B. (1956) Concept of General and Specific Combining Ability
#'             in relation to Diallel Crossing Systems. \emph{Australian Journal of Biological Sciences},
#'              \strong{9(4)}, 463--493.
#' \item Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'  }
#' @seealso
#'    \code{\link{Hayman}}
#'  , \code{\link{GriffingData1}}
#'  , \code{\link{GriffingData2}}
#'  , \code{\link{GriffingData3}}
#'  , \code{\link{GriffingData4}}
#' @import ggplot2 stats
#' @export
#' @examples
#' #-------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 1 & Model 1
#' #-------------------------------------------------------------
#' Griffing1Data1 <-
#'  Griffing(
#'      y      = Yield
#'    , Rep    = Rep
#'    , Cross1 = Cross1
#'    , Cross2 = Cross2
#'    , data   = GriffingData1
#'    , Method = 1
#'    , Model  = 1
#'  )

#' names(Griffing1Data1)
#' Griffing1Data1

#' Griffing1Data1Means <- Griffing1Data1$Means
#' Griffing1Data1ANOVA <- Griffing1Data1$ANOVA
#' Griffing1Data1Genetic.Components <- Griffing1Data1$Genetic.Components
#' Griffing1Data1Effects <- Griffing1Data1$Effects
#' Griffing1Data1StdErr <- as.matrix(Griffing1Data1$StdErr)
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 1 & Model 2
#' #--------------------------------------------------------------
#' Griffing2Data1 <-
#'  Griffing(
#'      y      = Yield
#'    , Rep    = Rep
#'    , Cross1 = Cross1
#'    , Cross2 = Cross2
#'    , data   = GriffingData1
#'    , Method = 1
#'    , Model  = 2
#'  )
#' names(Griffing2Data1)
#' Griffing2Data1
#' Griffing2Data1Means <- Griffing2Data1$Means
#' Griffing2Data1ANOVA <- Griffing2Data1$ANOVA
#' Griffing2Data1Genetic.Components <- Griffing2Data1$Genetic.Components
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 2 & Model 1
#' #--------------------------------------------------------------
#' Griffing1Data2 <-
#'  Griffing(
#'      y      = Yield
#'    , Rep    = Rep
#'    , Cross1 = Cross1
#'    , Cross2 = Cross2
#'    , data   = GriffingData2
#'    , Method = 2
#'    , Model  = 1
#'  )
#' names(Griffing1Data2)
#' Griffing1Data2
#' Griffing1Data2Means <- Griffing1Data2$Means
#' Griffing1Data2ANOVA <- Griffing1Data2$ANOVA
#' Griffing1Data2Genetic.Components <- Griffing1Data2$Genetic.Components
#' Griffing1Data2Effects <- Griffing1Data2$Effects
#' Griffing1Data2StdErr <- as.matrix(Griffing1Data2$StdErr)
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 2 & Model 2
#' #--------------------------------------------------------------
#' Griffing2Data2 <-
#'  Griffing(
#'      y      = Yield
#'    , Rep    = Rep
#'    , Cross1 = Cross1
#'    , Cross2 = Cross2
#'    , data   = GriffingData2
#'    , Method = 2
#'    , Model  = 2
#'  )
#' names(Griffing2Data2)
#' Griffing2Data2
#' Griffing2Data2Means <- Griffing2Data2$Means
#' Griffing2Data2ANOVA <- Griffing2Data2$ANOVA
#' Griffing2Data2Genetic.Components <- Griffing2Data2$Genetic.Components
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 3 & Model 1
#' #--------------------------------------------------------------
#' Griffing1Data3 <-
#'  Griffing(
#'      y      = Yield
#'    , Rep    = Rep
#'    , Cross1 = Cross1
#'    , Cross2 = Cross2
#'    , data   = GriffingData3
#'    , Method = 3
#'    , Model  = 1
#'  )
#' names(Griffing1Data3)
#' Griffing1Data3
#' Griffing1Data3Means <- Griffing1Data3$Means
#' Griffing1Data3ANOVA <- Griffing1Data3$ANOVA
#' Griffing1Data3Genetic.Components <- Griffing1Data3$Genetic.Components
#' Griffing1Data3Effects <- Griffing1Data3$Effects
#' Griffing1Data3StdErr <- as.matrix(Griffing1Data3$StdErr)
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 3 & Model 2
#' #--------------------------------------------------------------
#' Griffing2Data3 <-
#'  Griffing(
#'      y       = Yield
#'    , Rep     = Rep
#'    , Cross1  = Cross1
#'    , Cross2  = Cross2
#'    , data    = GriffingData3
#'    , Method  = 3
#'    , Model   = 2
#'  )
#' names(Griffing2Data3)
#' Griffing2Data3
#' Griffing2Data3Means <- Griffing2Data3$Means
#' Griffing2Data3ANOVA <- Griffing2Data3$ANOVA
#' Griffing2Data3Genetic.Components <- Griffing2Data3$Genetic.Components
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 4 & Model 1
#' #--------------------------------------------------------------
#' Griffing1Data4 <-
#'  Griffing(
#'      y       = Yield
#'    , Rep     = Rep
#'    , Cross1  = Cross1
#'    , Cross2  = Cross2
#'    , data    = GriffingData4
#'    , Method  = 4
#'    , Model   = 1
#'  )
#' names(Griffing1Data4)
#' Griffing1Data4
#' Griffing1Data4Means <- Griffing1Data4$Means
#' Griffing1Data4ANOVA <- Griffing1Data4$ANOVA
#' Griffing1Data4Genetic.Components <- Griffing1Data4$Genetic.Components
#' Griffing1Data4Effects <- Griffing1Data4$Effects
#' Griffing1Data4StdErr <- as.matrix(Griffing1Data4$StdErr)
#'
#'
#' #--------------------------------------------------------------
#' ## Diallel Analysis with Griffing's Aproach Method 4 & Model 2
#' #--------------------------------------------------------------
#' Griffing2Data4 <-
#'  Griffing(
#'      y       = Yield
#'    , Rep     = Rep
#'    , Cross1  = Cross1
#'    , Cross2  = Cross2
#'    , data    = GriffingData4
#'    , Method  = 4
#'    , Model   = 2
#'  )
#' names(Griffing2Data4)
#' Griffing2Data4
#' Griffing2Data4Means <- Griffing2Data4$Means
#' Griffing2Data4ANOVA <- Griffing2Data4$ANOVA
#' Griffing2Data4Genetic.Components <- Griffing2Data4$Genetic.Components

Griffing <-
  function(y, Rep, Cross1, Cross2, data, Method, Model) {
  if(Method > 4 | Model > 2) {
    stop ("Either Method No. or Model No. is incorrect: \n Method=1, 2, 3, or 4 and Model=1, or 2:")
  }
  else if(Method==1 & Model==1)
  {
    #-----------------------------------------------------------------------------
    # Method-I(Parents + F1's + reciprocals)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-I (Fixed Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    DataTotals      <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = sum)
    n               <- nrow(DataTotals)
    r               <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    SS.gca          <- ((sum((colSums(Means)+rowSums(Means))^2))/(2*n)-(2*(sum(colSums(Means)))^2)/n^2)
    SS.sca          <- (sum(Means*(Means+t(Means)))/2-(sum((colSums(Means)+rowSums(Means))^2))/(2*n)+(sum(colSums(Means)))^2/n^2)
    SS.reciprocals  <- sum((Means-t(Means))^2)/4
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
          c(
              n-1
            , n*(n-1)/2
            , n*(n-1)/2
            , Simple.ANOVA["Residuals","Df"]
            )
    SS <-
          c(
              SS.gca
            , SS.sca
            , SS.reciprocals
            , Simple.ANOVA["Residuals", "Sum Sq"]/r
            )
    MS      <- SS/df
    F.Test  <- c(MS[-4]/MS[4], NA)
    P.Value <- c(pf(F.Test[-4], df[-4], df[4], lower.tail = FALSE, log.p = FALSE), NA)
    ANOVA <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
         "gca"
        , "sca"
        , "reciprocals"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v <- (MS[1]-MS[4])/(2*n)
    sca.v <- (MS[2]-MS[4])
    r.v   <- (MS[3]-MS[4])/2
    gca.v.to.sca.v <- gca.v/sca.v
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, r.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "Reciprocal", "gca/sca"))

    #-----------------------------------------------------------------------------
    # General Combining Ability (gca) Effects
    #-----------------------------------------------------------------------------
    G <- diag((colSums(Means)+rowSums(Means))/(2*n)-(sum(colSums(Means)))/n^2)

    #-----------------------------------------------------------------------------
    # Specific Combining Ability (sca) Effects
    #-----------------------------------------------------------------------------
    S <- matrix(NA, nrow=n, ncol=n, byrow=TRUE)

    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if (i >= j)
        {
          S[i, j] <- 0
        }
        else
        {
          S[i, j] <- (Means[i, j] + Means[j, i])/2-(rowSums(Means)[i] + colSums(Means)[i] + rowSums(Means)[j] + colSums(Means)[j])/(2*n) + (sum(colSums(Means)))/n^2
        }
      }
    }

    #-----------------------------------------------------------------------------
    # Reciprocal Effects
    #-----------------------------------------------------------------------------
    R               <- t(Means-t(Means))/2
    R[upper.tri(R)] <- 0
    Effects         <- G + S + R

    #-----------------------------------------------------------------------------
    # Standard Errors
    #-----------------------------------------------------------------------------
    MSE        <- Simple.ANOVA["Rep", "Mean Sq"]
    SE.gi      <- sqrt((n-1)/(2*n^2)*MSE/r)
    SE.sii     <- sqrt((n-1)^2/n^2*MSE/r)
    SE.sij     <- sqrt((n^2-2*n+2)/(2*n^2)*MSE/r)
    SE.rij     <- sqrt(1/2*MSE/r)
    SE.gi.gj   <- sqrt(1/n*MSE/r)
    SE.sii.sji <- sqrt(2*(n-2)/n*MSE/r)
    SE.sii.sij <- sqrt((3*n-2)/(2*n)*MSE/r)
    SE.sii.sjk <- sqrt(3*(n-2)/(2*n)*MSE/r)
    SE.sij.sik <- sqrt((n-1)/n*MSE/r)
    SE.sij.skl <- sqrt((n-2)/n*MSE/r)
    SE.rij.rkl <- sqrt(MSE/r)
    StdErr     <-
                  c(
                      SE.gi
                    , SE.sii
                    , SE.sij
                    , SE.rij
                    , SE.gi.gj
                    , SE.sii.sji
                    , SE.sii.sij
                    , SE.sii.sjk
                    , SE.sij.sik
                    , SE.sij.skl
                    , SE.rij.rkl
                    )

    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        , Effects            = Effects
        , StdErr             = StdErr
        )
    #-----------------------------------------------------------------------------
    }
  else
  if(Method==1 & Model==2)
  {
    #-----------------------------------------------------------------------------
    # Method-I(Parents + F1's + reciprocals)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-II (Random Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    DataTotals      <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = sum)
    n               <- nrow(DataTotals)
    r               <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    SS.gca          <- ((sum((colSums(Means)+rowSums(Means))^2))/(2*n)-(2*(sum(colSums(Means)))^2)/n^2)
    SS.sca          <- (sum(Means*(Means+t(Means)))/2-(sum((colSums(Means)+rowSums(Means))^2))/(2*n)+(sum(colSums(Means)))^2/n^2)
    SS.reciprocals  <- sum((Means-t(Means))^2)/4
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
      c(
        n-1
        , n*(n-1)/2
        , n*(n-1)/2
        , Simple.ANOVA["Residuals","Df"]
      )
    SS <-
      c(
        SS.gca
        , SS.sca
        , SS.reciprocals
        , Simple.ANOVA["Residuals", "Sum Sq"]/r
      )
    MS      <- SS/df
    F.Test  <- c(MS[1]/MS[2], MS[2:3]/MS[4], NA)
    P.Value <- c(pf(F.Test[-4], c(df[1], df[2:3]), c(df[2], df[4]), lower.tail = FALSE, log.p = FALSE), NA)

    ANOVA <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "reciprocals"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v <- (MS[1]-(MS[4]+n*(n-1)*MS[2])/(n^2-n+1))/(2*n)
    sca.v <- (MS[2]-MS[4])*(n^2)/(2*(n^2-n+1))
    r.v   <- (MS[3]-MS[4])/2
    e.v   <- MS[4]
    A.v   <- 2*gca.v
    D.v   <- sca.v
    gca.v.to.sca.v <- gca.v/sca.v
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, r.v, e.v,
                                                   A.v, D.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "Reciprocal", "Error",
                                                 "Additive", "Dominant", "gca/sca"))
    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        )
    #-----------------------------------------------------------------------------
  }
  else
  if(Method==2 & Model==1)
  {
    #-----------------------------------------------------------------------------
    # Method-II (Parents and one set of F1's)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-I (Fixed Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    Means1 <- Means
    Means1[lower.tri(Means1)] <- 0
    Means2 <- Means1 + t(Means1) - diag(diag(Means1))
    n <- nrow(Means2)
    r <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    SS.gca <- ((sum((rowSums(Means2) + diag(Means2))^2)-(4*(sum(rowSums(Means1)))^2)/n)/(n+2))
    SS.sca <- sum((Means1)^2)-sum((rowSums(Means2)+diag(Means2))^2)/(n+2)+(2*(sum(rowSums(Means1)))^2)/((n+1)*(n+2))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
      c(
          n-1
        , n*(n-1)/2
        , Simple.ANOVA["Residuals","Df"]
      )
    SS <-
      c(
          SS.gca
        , SS.sca
        , Simple.ANOVA["Residuals", "Sum Sq"]/r
      )
    MS      <- SS/df
    F.Test  <- c(MS[-3]/MS[3], NA)
    P.Value <- c(pf(F.Test[-3], df[-3], df[3], lower.tail = FALSE, log.p = FALSE), NA)
    ANOVA <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v <- (MS[1]-MS[3])/(n+2)
    sca.v <- (MS[2]-MS[3])
    gca.v.to.sca.v <- ((MS[1]-MS[3])/(n+2))/(MS[2]-MS[3])
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "gca/sca"))
    #-----------------------------------------------------------------------------
    # General Combining Ability (gca) Effects
    #-----------------------------------------------------------------------------
    G <- diag(((rowSums(Means2)+diag(Means2))-(2*(sum(rowSums(Means1))))/n)/(n+2))

    #-----------------------------------------------------------------------------
    # Specific Combining Ability (sca) Effects
    #-----------------------------------------------------------------------------
    S <- matrix(NA, nrow=n, ncol=n, byrow=TRUE)

    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if (i>=j)
        {
          S[i, j] <- 0
        }
        else
        {
          S[i, j] <- Means[i, j]-(rowSums(Means2)[i] + diag(Means2)[i] + rowSums(Means2)[j] + diag(Means2)[j])/(n+2) + (2*sum(rowSums(Means1)))/((n+1)*(n+2))
        }
      }
    }

    Effects                     <- G+S
    Effects[lower.tri(Effects)] <- NA
    dimnames(Effects)           <- dimnames(Means)
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Standard Errors
    #-----------------------------------------------------------------------------
    MSE        <- Simple.ANOVA["Rep", "Mean Sq"]
    SE.gi      <- sqrt((n-1)/(n*(n+2))*MSE/r)
    SE.sii     <- sqrt((n^2+n+2)/((n+1)*(n+2))*MSE/r)
    SE.sij     <- sqrt(n*(n-1)/((n+1)*(n+2))*MSE/r)
    SE.gi.gj   <- sqrt(2/(n+2)*MSE/r)
    SE.sii.sjj <- sqrt(2*(n-2)/(n+2)*MSE/r)
    SE.sij.sik <- sqrt(2*(n+1)/(n+2)*MSE/r)
    SE.sij.skl <- sqrt(2*n/(n+2)*MSE/r)
    StdErr     <-
                  c(
                      SE.gi
                    , SE.sii
                    , SE.sij
                    , SE.gi.gj
                    , SE.sii.sjj
                    , SE.sij.sik
                    , SE.sij.skl
                    )
    #-----------------------------------------------------------------------------

    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        , Effects            = Effects
        , StdErr             = StdErr
        )
    #-----------------------------------------------------------------------------
  }
  else
  if(Method==2 & Model==2)
  {
    #-----------------------------------------------------------------------------
    # Method-II (Parents and one set of F1's)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-II (Random Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    Means1 <- Means
    Means1[lower.tri(Means1)] <- 0
    Means2 <- Means1 + t(Means1) - diag(diag(Means1))
    n <- nrow(Means2)
    r <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))

    SS.gca <- ((sum((rowSums(Means2)+diag(Means2))^2)-(4*(sum(rowSums(Means1)))^2)/n)/(n+2))
    SS.sca <- sum((Means1)^2)-sum((rowSums(Means2)+diag(Means2))^2)/(n+2)+(2*(sum(rowSums(Means1)))^2)/((n+1)*(n+2))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
      c(
          n-1
        , n*(n-1)/2
        , Simple.ANOVA["Residuals","Df"]
      )
    SS <-
      c(
          SS.gca
        , SS.sca
        , Simple.ANOVA["Residuals", "Sum Sq"]/r
      )
    MS      <- SS/df
    F.Test  <- c(MS[1]/MS[2], MS[2]/MS[3], NA)
    P.Value <- c(pf(F.Test[-3], c(df[1], df[2]), c(df[2], df[3]), lower.tail = FALSE, log.p = FALSE),  NA)
    ANOVA   <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v <- (MS[1]-MS[2])/(n+2)
    sca.v <- (MS[2]-MS[3])
    e.v   <- MS[3]
    A.v   <- 2*gca.v
    D.v   <- sca.v
    gca.v.to.sca.v <- gca.v/sca.v
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, e.v, A.v, D.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "Error", "Additive", "Dominance", "gca/sca"))
    #-----------------------------------------------------------------------------
    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        )
    #-----------------------------------------------------------------------------
  }
  else
  if(Method==3 & Model==1)
  {
    #-----------------------------------------------------------------------------
    # Method-III (One set of F1's and reciprocals)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-I (Fixed Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    Means1 <- Means
    diag(Means1) <- 0
    n <- nrow(Means1)
    r <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    dimnames(Means1) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))

    SS.gca         <- ((sum((rowSums(Means1) + colSums(Means1))^2))/(2*(n-2))-(2*(sum(colSums(Means1)))^2)/(n*(n-2)))
    SS.sca         <- (sum(((Means1 + t(Means1))^2)/2)/2-(sum((rowSums(Means1) + colSums(Means1))^2))/(2*(n-2)) + (sum(colSums(Means1)))^2/((n-1)*(n-2)))
    SS.reciprocals <- sum((Means1 - t(Means1))^2)/4
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
         c(
            n-1
          , n*(n-3)/2
          , n*(n-1)/2
          , Simple.ANOVA["Residuals","Df"]
          )
    SS <-
          c(
              SS.gca
            , SS.sca
            , SS.reciprocals
            , Simple.ANOVA["Residuals", "Sum Sq"]/r
            )
    MS      <- SS/df
    F.Test  <- c(MS[-4]/MS[4], NA)
    P.Value <- c(pf(F.Test[-4], df[-4], df[4], lower.tail = FALSE, log.p = FALSE), NA)
    ANOVA   <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                          `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "reciprocals"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v              <- (MS[1]-MS[4])/(2*(n-2))
    sca.v              <- (MS[2]-MS[4])/2
    r.v                <- (MS[3]-MS[4])/2
    gca.v.to.sca.v     <- ((MS[1]-MS[4])/(2*(n-2)))/((MS[2]-MS[4])/2)
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, r.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "Reciprocal", "gca/sca"))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # General Combining Ability (gca) Effects
    #-----------------------------------------------------------------------------
    G <- diag((n*(colSums(Means1) + rowSums(Means1)) - (2*sum(colSums(Means1))))/(2*n*(n-2)))

    #-----------------------------------------------------------------------------
    # Specific Combining Ability (sca) Effects
    #-----------------------------------------------------------------------------
    S <- matrix(NA, nrow=n, ncol=n, byrow=TRUE)

    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if (i>=j)
        {
          S[i, j] <- 0
        }
        else
        {
          S[i, j] <- (Means1[i, j] + Means1[j, i])/2-(rowSums(Means1)[i] + colSums(Means1)[i] + rowSums(Means1)[j] + colSums(Means1)[j])/(2*(n-2)) + (sum(colSums(Means1)))/((n-1)*(n-2))
        }
      }
    }
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Reciprocal Effects
    #-----------------------------------------------------------------------------
    R               <- t(Means1-t(Means1))/2
    R[upper.tri(R)] <- 0
    Effects         <- G + S + R
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Standard Errors
    #-----------------------------------------------------------------------------
    MSE        <- Simple.ANOVA["Rep", "Mean Sq"]
    SE.gi      <- sqrt((n-1)/(2*n*(n-2))*MSE/r)
    SE.sii     <- sqrt((n-3)/(2*(n-1))*MSE/r)
    SE.rij     <- sqrt(1/2*MSE/r)
    SE.gi.gj   <- sqrt(1/(n-2)*MSE/r)
    SE.sij.sik <- sqrt((n-3)/(n-2)*MSE/r)
    SE.sij.skl <- sqrt((n-4)/(n-2)*MSE/r)
    StdErr     <-
                  c(
                      SE.gi
                    , SE.sii
                    , SE.rij
                    , SE.gi.gj
                    , SE.sij.sik
                    , SE.sij.skl
                    )
    #-----------------------------------------------------------------------------

    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        , Effects            = Effects
        , StdErr             = StdErr
        )
    #-----------------------------------------------------------------------------
  }
  else
  if(Method==3 & Model==2)
  {
    #-----------------------------------------------------------------------------
    # Method-III (One set of F1's and reciprocals)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-II (Random Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    Means1 <- Means
    diag(Means1) <- 0
    n <- nrow(Means1)
    r <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    dimnames(Means1) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))

    SS.gca         <- ((sum((rowSums(Means1) + colSums(Means1))^2))/(2*(n-2))-(2*(sum(colSums(Means1)))^2)/(n*(n-2)))
    SS.sca         <- (sum(((Means1 + t(Means1))^2)/2)/2-(sum((rowSums(Means1)+colSums(Means1))^2))/(2*(n-2))+(sum(colSums(Means1)))^2/((n-1)*(n-2)))
    SS.reciprocals <- sum((Means1 - t(Means1))^2)/4
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
      c(
          n-1
        , n*(n-3)/2
        , n*(n-1)/2
        , Simple.ANOVA["Residuals","Df"]
      )
    SS <-
      c(
          SS.gca
        , SS.sca
        , SS.reciprocals
        , Simple.ANOVA["Residuals", "Sum Sq"]/r
      )
    MS      <- SS/df
    F.Test  <- c(MS[1]/MS[2], MS[2:3]/MS[4], NA)
    P.Value <- c(pf(F.Test[-4], c(df[1], df[2:3]), c(df[2], df[4]), lower.tail = FALSE, log.p = FALSE), NA)
    ANOVA   <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                          `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "reciprocals"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v              <- (MS[1]-MS[2])/(2*(n-2))
    sca.v              <- (MS[2]-MS[4])/2
    r.v                <- (MS[3]-MS[4])/2
    e.v                <- MS[4]
    A.v                <- 2*gca.v
    D.v                <- sca.v
    gca.v.to.sca.v     <- gca.v/sca.v
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, r.v, e.v,
                                                   A.v, D.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "Reciprocal", "Error",
                                                 "Additive", "Dominant", "gca/sca"))
    #-----------------------------------------------------------------------------
    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        )
    #-----------------------------------------------------------------------------
  }
  else
  if(Method==4 & Model==1)
  {
    #-----------------------------------------------------------------------------
    # Method-IV (One set of F1's only)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-I (Fixed Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    Means1 <- rbind(cbind(rep(NA, ncol(Means)), Means), rep(NA, ncol(Means) + 1))
    Means1[lower.tri(Means1, diag=TRUE)] <- 0
    Means2 <- Means1 + t(Means1) - diag(diag(Means1))
    n <- ncol(Means2)
    r <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 2:n), paste0("Cross", 2:n))
    dimnames(Means1) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    dimnames(Means2) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))

    SS.gca <- ((sum((rowSums(Means2))^2))/(n-2)-(4*(sum(rowSums(Means2))/2)^2)/(n*(n-2)))
    SS.sca <- (sum((Means2)^2))/2-(sum((rowSums(Means2))^2))/(n-2)+(2*(sum(rowSums(Means2))/2)^2)/((n-1)*(n-2))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
      c(
          n-1
        , n*(n-3)/2
        , Simple.ANOVA["Residuals","Df"]
      )
    SS <-
      c(
          SS.gca
        , SS.sca
        , Simple.ANOVA["Residuals", "Sum Sq"]/r
      )
    MS      <- SS/df
    F.Test  <- c(MS[-3]/MS[3], NA)
    P.Value <- c(pf(F.Test[-3], df[-3], df[3], lower.tail = FALSE, log.p = FALSE), NA)
    ANOVA   <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                          `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v              <- (MS[1]-MS[3])/(n-2)
    sca.v              <- (MS[2]-MS[3])
    gca.v.to.sca.v     <- gca.v/sca.v
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "gca/sca"))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # General Combining Ability (gca) Effects
    #-----------------------------------------------------------------------------
    G <- diag((n*rowSums(Means2)-2*(sum(rowSums(Means2))/2))/(n*(n-2)))

    #-----------------------------------------------------------------------------
    # Specific Combining Ability (sca) Effects
    #-----------------------------------------------------------------------------
    S <- matrix(NA, nrow=n, ncol=n, byrow=TRUE)

    for(i in 1:n)
    {
      for(j in 1:n)
      {
        if (i>=j)
        {
          S[i, j] <- 0
        }
        else
        {
          S[i, j] <- Means2[i, j]-(rowSums(Means2)[i]+colSums(Means2)[j])/(n-2)+(2*sum(rowSums(Means2))/2)/((n-1)*(n-2))
        }
      }
    }
    #-----------------------------------------------------------------------------
    Effects                     <- G + S
    Effects[lower.tri(Effects)] <- NA
    dimnames(Effects)           <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Standard Errors
    #-----------------------------------------------------------------------------
    MSE        <- Simple.ANOVA["Rep", "Mean Sq"]
    SE.gi      <- sqrt((n-1)/(n*(n-2))*MSE/r)
    SE.sii     <- sqrt((n-3)/(n-1)*MSE/r)
    SE.sij     <- sqrt((n^2-2*n+2)/(2*n^2)*MSE/r)
    SE.gi.gj   <- sqrt(2/(n-2)*MSE/r)
    SE.sij.sik <- sqrt(2*(n-3)/(n-2)*MSE/r)
    SE.sij.skl <- sqrt(2*(n-4)/(n-2)*MSE/r)
    StdErr     <-
                  c(
                      SE.gi
                    , SE.sii
                    , SE.sij
                    , SE.gi.gj
                    , SE.sij.sik
                    , SE.sij.skl
                    )
    #-----------------------------------------------------------------------------
    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        , Effects            = Effects
        , StdErr             = StdErr
        )
    #-----------------------------------------------------------------------------
  }
  else
  {
    #-----------------------------------------------------------------------------
    # Method-IV (One set of F1's only)
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Model-II (Random Effects Model)
    #-----------------------------------------------------------------------------
    data$Trt <- paste(data[["Cross1"]], data[["Cross2"]], sep="-")

    y      <- deparse(substitute(y))
    Rep    <- deparse(substitute(Rep))
    Cross1 <- deparse(substitute(Cross1))
    Cross2 <- deparse(substitute(Cross2))
    Trt    <- deparse(substitute(Trt))

    Simple.ANOVA    <- summary(aov(data[[y]] ~ as.factor(data[[Rep]]) + data[[Trt]]), data=data)[[1]]
    rownames(Simple.ANOVA) <- c("Rep", "Trt", "Residuals")
    Means           <- tapply(X= data[[y]], INDEX = list(data[[Cross1]], data[[Cross2]]), FUN = mean)
    Means1 <- rbind(cbind(rep(NA, ncol(Means)), Means), rep(NA, ncol(Means) + 1))
    Means1[lower.tri(Means1, diag=TRUE)] <- 0
    Means2 <- Means1 + t(Means1) - diag(diag(Means1))
    n <- ncol(Means2)
    r <- length(levels(as.factor(data[[Rep]])))
    dimnames(Means) <- list(paste0("Cross", 2:n), paste0("Cross", 2:n))
    dimnames(Means1) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))
    dimnames(Means2) <- list(paste0("Cross", 1:n), paste0("Cross", 1:n))

    SS.gca <- ((sum((rowSums(Means2))^2))/(n-2)-(4*(sum(rowSums(Means2))/2)^2)/(n*(n-2)))
    SS.sca <- (sum((Means2)^2))/2-(sum((rowSums(Means2))^2))/(n-2)+(2*(sum(rowSums(Means2))/2)^2)/((n-1)*(n-2))
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
      c(
          n-1
        , n*(n-3)/2
        , Simple.ANOVA["Residuals","Df"]
      )
    SS <-
      c(
          SS.gca
        , SS.sca
        , Simple.ANOVA["Residuals", "Sum Sq"]/r
      )
    MS      <- SS/df
    F.Test  <- c(MS[1]/MS[2], MS[2]/MS[3], NA)
    P.Value <- c(pf(F.Test[-3], c(df[1], df[2]), c(df[2], df[3]), lower.tail = FALSE, log.p = FALSE), NA)
    ANOVA   <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                          `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
      c(
        "gca"
        , "sca"
        , "Error"
      )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Genetic Components
    #-----------------------------------------------------------------------------
    gca.v              <- (MS[1]-MS[2])/(n-2)
    sca.v              <- (MS[2]-MS[3])
    e.v                <- MS[3]
    A.v                <- 2*gca.v
    D.v                <- sca.v
    gca.v.to.sca.v     <- gca.v/sca.v
    Genetic.Components <- data.frame("Components"=c(gca.v, sca.v, e.v,
                                                   A.v, D.v, gca.v.to.sca.v),
                                  row.names=list("gca", "sca", "Error",
                                                 "Additive", "Dominant", "gca/sca"))
    #-----------------------------------------------------------------------------
    list(
          Means              = Means
        , ANOVA              = ANOVA
        , Genetic.Components = Genetic.Components
        )
    #-----------------------------------------------------------------------------
  }
}
