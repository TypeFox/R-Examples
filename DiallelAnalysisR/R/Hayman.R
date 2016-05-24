#' @name   Hayman
#' @aliases Hayman
#' @title Diallel Analysis using Hayman Approach
#' @description \code{Hayman} is used for performing Diallel Analysis using Hayman's Approach.
#' @param y Numeric Response Vector
#' @param Rep  Replicate as factor
#' @param Cross1 Cross 1 as factor
#' @param Cross2 Cross 2 as factor
#' @param data A \code{data.frame}
#' @details Diallel Analysis using Haymans's approach.
#' @return Means Means
#' @return ANOVA Analysis of Variance (ANOVA) table
#' @return Genetic.Components Genetic Components
#' @return Effects Effects of Crosses
#' @return StdErr Standard Errors of Crosses
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @references \enumerate{
#' \item Hayman, B. I. (1954 a) The Theory and Analysis of Diallel Crosses.
#'              \emph{Genetics},
#'              \strong{39}, 789--809.
#' \item Hayman, B. I. (1954 b) The Analysis of Variance of Diallel Tables.
#'              \emph{Biometrics},
#'              \strong{10}, 235--244.
#' \item Hayman, B. I. (1957) Interaction, Heterosis and Diallel Crosses.
#'              \emph{Genetics},
#'              \strong{42}, 336--355.
#' \item Singh, R. K. and Chaudhary, B. D. (2004) \emph{Biometrical Methods in Quantitative Genetic Analysis}.
#'              New Delhi: Kalyani.
#'  }
#' @seealso
#'    \code{\link{Griffing}}
#'  , \code{\link{HaymanData}}
#' @import ggplot2 stats
#' @export
#' @examples
#' #------------------------------------------
#' ## Diallel Analysis with Haymans's Aproach
#' #------------------------------------------
#'
#' Hayman1Data <-
#'  Hayman(
#'      y      = Yield
#'    , Rep    = Rep
#'    , Cross1 = Cross1
#'    , Cross2 = Cross2
#'    , data   = HaymanData
#'    )
#'
#' Hayman1Data
#' names(Hayman1Data)
#'
#' Hayman1DataMeans <- Hayman1Data$Means
#' Hayman1DataANOVA <- Hayman1Data$ANOVA
#' Hayman1DataWr.Vr.Table <- Hayman1Data$Wr.Vr.Table
#'
#' Hayman1DataComponents.of.Variation <- Hayman1Data$Components.of.Variation
#' Hayman1DataOther.Parameters <- Hayman1Data$Other.Parameters
#' Hayman1DataFr <- Hayman1Data$Fr
#'
#' #----------------
#' # Wr-Vr Graph
#' #----------------
#' VOLO     <- Hayman1Data$VOLO
#' In.Value <- Hayman1Data$In.Value
#' a        <- Hayman1Data$a
#' b        <- Hayman1Data$b
#' Wr.Vr    <- Hayman1Data$Wr.Vr.Table
#'
#'
#' library(ggplot2)
#' ggplot(data=data.frame(x=c(0, max(In.Value, Wr.Vr$Vr, Wr.Vr$Wr, Wr.Vr$Wrei))), aes(x)) +
#'   stat_function(fun=function(x) {sqrt(x*VOLO)}, color="blue") +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   geom_abline(intercept = a, slope = b) +
#'   geom_abline(intercept = mean(Wr.Vr$Wr)-mean(Wr.Vr$Vr), slope = 1) +
#'   geom_segment(aes(
#'       x     = mean(Wr.Vr$Vr)
#'     , y     = min(0, mean(Wr.Vr$Wr))
#'     , xend  = mean(Wr.Vr$Vr)
#'     , yend  = max(0, mean(Wr.Vr$Wr))
#'   )
#'   , color = "green"
#'   ) +
#'   geom_segment(aes(
#'       x     = min(0, mean(Wr.Vr$Vr))
#'     , y     = mean(Wr.Vr$Wr)
#'     , xend  = max(0, mean(Wr.Vr$Vr))
#'     , yend  = mean(Wr.Vr$Wr)
#'   )
#'   , color = "green"
#'   )  +
#'   lims(x=c(min(0, Wr.Vr$Vr, Wr.Vr$Wrei), max(Wr.Vr$Vr, Wr.Vr$Wrei)),
#'        y=c(min(0, Wr.Vr$Wr, Wr.Vr$Wrei), max(Wr.Vr$Wr, Wr.Vr$Wri))
#'   ) +
#'   labs(
#'          x = expression(V[r])
#'        , y = expression(W[r])
#'        , title = expression(paste(W[r]-V[r] , " Graph"))
#'        ) +
#'   theme_bw()

Hayman <-
  function(y, Rep, Cross1, Cross2, data) {
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

    SS.Total <- (sum(DataTotals^2)-(sum(colSums(DataTotals)))^2/(n^2))/r
    SS.a     <- ((sum((colSums(DataTotals) + rowSums(DataTotals))^2))/(2*n)-(2*(sum(colSums(DataTotals)))^2)/n^2)/r
    SS.b     <- ((sum((DataTotals+t(DataTotals))^2)/r)-sum((colSums(DataTotals)+rowSums(DataTotals))^2)/(2*n)+(sum(colSums(DataTotals)))^2/(n^2))/r
    SS.b1    <- ((sum(colSums(DataTotals))-n*sum(diag(DataTotals)))^2/(n^2*(n-1)))/r
    SS.b2    <- ((sum((colSums(DataTotals)+rowSums(DataTotals)-n*diag(DataTotals))^2))/(n*(n-2))-(2*sum(colSums(DataTotals))-n*sum(diag(DataTotals)))^2/(n^2*(n-2)))/r
    SS.b3    <- ((sum((DataTotals+t(DataTotals))^2)/r)-sum(diag(DataTotals)^2)-sum((colSums(DataTotals)+rowSums(DataTotals)-2*diag(DataTotals))^2)/(2*(n-2))+(sum(colSums(DataTotals))-sum(diag(DataTotals)))^2/((n-1)*(n-2)))/r
    SS.c     <- ((sum((colSums(DataTotals)-rowSums(DataTotals))^2))/(2*n))/r
    SS.d     <- ((sum((DataTotals-t(DataTotals))^2)/r)-(sum((colSums(DataTotals)-rowSums(DataTotals))^2))/(2*n))/r
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # ANOVA Table
    #-----------------------------------------------------------------------------
    df <-
          c(
              sum(Simple.ANOVA["Df"])
            , Simple.ANOVA["Rep","Df"]
            , Simple.ANOVA["Trt","Df"]
            , n-1
            , n*(n-1)/2
            , 1
            , n-1
            , n*(n-3)/2
            , n-1
            , (n-1)*(n-2)/2
            , Simple.ANOVA["Residuals","Df"]
            )

    SS <-
          c(
              sum(Simple.ANOVA["Sum Sq"])
            , Simple.ANOVA["Rep", "Sum Sq"]
            , Simple.ANOVA["Trt", "Sum Sq"]
            , SS.a
            , SS.b
            , SS.b1
            , SS.b2
            , SS.b3
            , SS.c
            , SS.d
            , Simple.ANOVA["Residuals", "Sum Sq"]
            )

    MS      <- c(NA, SS[-1]/df[-1])
    F.Test  <- c(NA, MS[2:10]/MS[11], NA)
    P.Value <- c(NA, pf(F.Test[2:10], df[2:10], df[11], lower.tail = FALSE, log.p = FALSE), NA)

    ANOVA <- data.frame(`Df` =  df, `Sum Sq` = SS, `Mean Sq` = MS,
                        `F value` = F.Test, `Pr(>F)` = P.Value, check.names = FALSE)
    rownames(ANOVA) <-
                      c(
                          "Total"
                        , "Rep"
                        , "Treatment"
                        , "  Additive"
                        , "  Non-Additive"
                        , "      b1"
                        , "      b2"
                        , "      b3"
                        , "  Maternal"
                        , "  Reciprocal"
                        , "Error"
                        )

    class(ANOVA) <- c("anova", "data.frame")
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    # Means, Varaince, and Covarainces
    #-----------------------------------------------------------------------------
    ParentalMean <- mean(diag(Means))
    VOLO         <- var(diag(Means))
    DataMeans1   <- Means
    DataMeans1[lower.tri(DataMeans1)] <- 0
    DataMeans2   <- DataMeans1 + t(DataMeans1) - diag(diag(DataMeans1))
    Vr           <- as.matrix(apply(DataMeans2, 2, var))
    V1L1         <- mean(Vr)
    V0L1         <- var(colMeans(DataMeans2))
    Wr           <- (Means %*% as.matrix(diag(Means)) - as.matrix(rowSums(Means)*sum(diag(Means))/n))/(n-1)
    WOLO1        <- mean(Wr)
    Wr.Vr        <- data.frame("Wr"=Wr, "Vr"=Vr, "Wr.Minus.Vr"=Wr-Vr,"Wr.Plus.Vr"=Wr+Vr, "Yr"=diag(Means))
    t.sq         <- ((n-2)/r)*(var(Wr.Vr$Vr)-var(Wr.Vr$Wr))^2/(var(Wr.Vr$Vr)*var(Wr.Vr$Wr)-(var(Wr.Vr$Vr, Wr.Vr$Wr))^2)
    b            <- var(Wr.Vr$Vr, Wr.Vr$Wr)/var(Wr.Vr$Vr)
    a            <- mean(Wr.Vr$Wr)-b*mean(Wr.Vr$Vr)
    In.Value     <- sqrt(V1L1*VOLO)
    Wr.Vr$Wri    <- sqrt(Wr.Vr$Vr*VOLO)
    Wr.Vr$Wrei   <- mean(Wr.Vr$Wr)-b*mean(Wr.Vr$Vr)+b*Wr.Vr$Vr
    Wr.Vr$Wreip  <- mean(Wr.Vr$Wr)-mean(Wr.Vr$Vr)+Wr.Vr$Vr
    Wr.Vr
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Estimation of Components of Variation
    #-----------------------------------------------------------------------------
    E        <- (Simple.ANOVA["Rep", "Sum Sq"] + Simple.ANOVA["Residuals", "Sum Sq"])/
                ((Simple.ANOVA["Rep", "Df"] + Simple.ANOVA["Residuals", "Df"])*r)
    D        <- VOLO - E
    F        <- 2*VOLO - 4*WOLO1 - 2*(n-2)*E/n
    H1       <- VOLO - 4*WOLO1 + 4*V1L1 - (3*n-2)*E/n
    H2       <- 4*V1L1 + 4*V0L1 - 2*E
    h2       <- 4*((sum(Means)/n-sum(diag(Means)))/n)^2-4*(n-1)*E/n^2
    s2       <- var(Wr.Vr$Wr.Minus.Vr)/2
    SE.E     <- sqrt(s2*(n^4/n^5))
    SE.D     <- sqrt(s2*(n^5+n^4)/n^5)
    SE.F     <- sqrt(s2*(4*n^5+20*n^4-16*n^3+16*n^2)/n^5)
    SE.H1    <- sqrt(s2*(n^5+41*n^4-12*n^3+4*n^2)/n^5)
    SE.H2    <- sqrt(s2*(36*n^4)/n^5)
    SE.h2    <- sqrt(s2*(16*n^4+16*n^3-32*n+16)/n^5)
    Estimate <- c(E, D, F, H1, H2, h2)
    StdErr   <- c(SE.E, SE.D, SE.F, SE.H1, SE.H2, SE.h2)
    tval     <- Estimate/StdErr
    Components.of.Variation <- data.frame("Estimate"=Estimate,
                                          "StdErr"=StdErr,
                                          "t.value"=tval,
                                          row.names=c("E", "D", "F", "H1", "H2", "h2"))

    Fr <- data.frame("Fr"=2*(VOLO-WOLO1+V1L1-(Wr.Vr$Wr+Wr.Vr$Vr))-2*(n-2)*E/n, row.names=paste("Fr", 1:n, sep=""))
    h2.ns <- (D+H1-H2-F)/(D+H1-0.5*H2-F+2*E)

    Other.Par <- data.frame(
                            "Other Parameters"=c(
                                                   sqrt(H1/D)
                                                 , H2/(4*H1)
                                                 , (sqrt(4*D*H1)+F)/(sqrt(4*D*H1)-F)
                                                 , cor(Wr.Vr$Wr.Plus.Vr, Wr.Vr$Yr)
                                                 , (cor(Wr.Vr$Wr.Plus.Vr,Wr.Vr$Yr))^2
                                                 , h2/H2
                                                 , h2.ns
                                                 )
                            )

    list(
           Means                   = Means
         , ANOVA                   = ANOVA
         , VOLO                    = VOLO
         , In.Value                = In.Value
         , a                       = a
         , b                       = b
         , Wr.Vr.Table             = Wr.Vr
         , Components.of.Variation = Components.of.Variation
         , Other.Parameters        = Other.Par
         , Fr                      = Fr
           )
    }
