mediation <- function (x, mediator, dv, S = NULL, N = NULL, x.location.S = NULL,  mediator.location.S = NULL, dv.location.S = NULL, mean.x = NULL, 
    mean.m = NULL, mean.dv = NULL, conf.level = 0.95, bootstrap = FALSE, B = 10000, which.boot="both", save.bs.replicates=FALSE, complete.set=FALSE) 
{
# Here is the internal mediation function.
    .mediation <- function(x = x, mediator = mediator, dv = dv, 
        S = NULL, N = N, x.location.S = x.location.S, mediator.location.S = mediator.location.S, 
        dv.location.S = dv.location.S, mean.x = mean.x, mean.m = mean.m, 
        mean.dv = mean.dv, conf.level = conf.level) {
        if (!is.null(S)) {
            These <- c(x.location.S, mediator.location.S, dv.location.S)
            Cov.Matrix <- as.matrix(S[These, These])
        }
        if (is.null(S)) {
            Data <- na.omit(cbind(x, mediator, dv))
            Cov.Matrix <- var(Data)
            N <- dim(Data)[1]
            mean.dv <- mean(dv)
            mean.x <- mean(x)
            mean.m <- mean(mediator)
            s.dv <- scale(dv)
            s.x <- scale(x)
            s.mediator <- scale(mediator)
            Y.on.X <- lm(dv ~ x)
            resid.Y.on.X <- resid(Y.on.X)
            standardized.Y.on.X <- lm(s.dv ~ s.x)
            standardized.resid.Y.on.X <- resid(standardized.Y.on.X)
            M.on.X <- lm(mediator ~ x)
            resid.M.on.X <- resid(M.on.X)
            standardized.M.on.X <- lm(s.mediator ~ s.x)
            standardized.resid.M.on.X <- resid(standardized.M.on.X)
            Y.on.X.and.M <- lm(dv ~ x + mediator)
            resid.Y.on.X.and.M <- resid(Y.on.X.and.M)
            standardized.Y.on.X.and.M <- lm(s.dv ~ s.x + s.mediator)
            standardized.resid.Y.on.X.and.M <- resid(standardized.Y.on.X.and.M)
            Y.on.M <- lm(dv ~ mediator)
            resid.Y.on.M <- resid(Y.on.M)
            standardized.Y.on.M <- lm(s.dv ~ s.mediator)
            standardized.resid.Y.on.M <- resid(standardized.Y.on.M)
            e.1M <- resid.M.on.X
            e.1Y <- resid.Y.on.X + resid.Y.on.M - resid.Y.on.X.and.M
            standardized.e.1M <- standardized.resid.M.on.X
            standardized.e.1Y <- standardized.resid.Y.on.X + 
                standardized.resid.Y.on.M - standardized.resid.Y.on.X.and.M
            e.0M <- mediator - mean.m
            e.0Y <- dv - mean.dv
            standardized.e.0M <- s.mediator - 0
            standardized.e.0Y <- s.dv - 0
      Residual.Based_Gamma <- as.numeric(1 -(sum(abs(e.1M) + 
                abs(e.1Y)))/(sum(abs(e.0M) + abs(e.0Y))))
                                
			Residual.Based.Standardized_gamma <- as.numeric(1 - 
              (sum(abs(standardized.e.1M) + abs(standardized.e.1Y)))/(sum(abs(standardized.e.0M) + 
              abs(standardized.e.0Y))))
        }
        Cor.Matrix <- cov2cor(Cov.Matrix)
        Dim.Cov.Matrix <- dim(Cov.Matrix)[1]
        s.XY <- Cov.Matrix[Dim.Cov.Matrix, -3]
        S.XX <- Cov.Matrix[1:(Dim.Cov.Matrix - 1), 1:(Dim.Cov.Matrix - 
            1)]
        B.Y_X <- solve(S.XX[1:1, 1:1]) %*% s.XY[1]
        B.Y_X <- cbind(mean.dv - mean.x * B.Y_X, B.Y_X)
        colnames(B.Y_X) <- c("Intercept.Y_X", "c (Regressor)")
        path.c <- B.Y_X[2]
        R2.Y_X <- (t(s.XY[1]) %*% solve(S.XX[1:1, 1:1]) %*% s.XY[1])/Cov.Matrix[Dim.Cov.Matrix, 
            Dim.Cov.Matrix]
        R2.Y_X.Adj <- 1 - ((1 - R2.Y_X) * ((N - 1)/(N - 1 - 1)))
        CI.R2.Y_X <- ci.R2(R2 = R2.Y_X, conf.level = conf.level, 
            Random.Predictors = TRUE, N = N, p = 1)
        Model.F.Y_X <- (R2.Y_X/1)/((1 - R2.Y_X)/(N - 1 - 1))
        
        MSE.Y_X <- (1 - R2.Y_X) * ((N - 1) * Cov.Matrix[3, 3]/(N - 2))
        RMSE.Y_X <- sqrt(round(MSE.Y_X, 10))
        
        SE.Y_X <- cbind(c(sqrt((1 - R2.Y_X) * ((N - 1) * Cov.Matrix[3, 
            3]/(N - 2))) * sqrt(1/N + mean.x^2/((N - 1) * S.XX[1, 
            1])), sqrt((1 - R2.Y_X)/(N - 1 - 1)) * sqrt(Cov.Matrix[3, 
            3]/S.XX[1, 1])))
        
        
        t.Y_X <- t(B.Y_X)/SE.Y_X
        p.Y_X <- 2 * (pt(-1 * abs(t.Y_X), df = N - 1 - 1))
        CL.Low.Y_X <- t(B.Y_X) - qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.Y_X
        CL.Up.Y_X <- t(B.Y_X) + qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.Y_X
        Values.Y_X <- cbind(t(B.Y_X), SE.Y_X, t.Y_X, p.Y_X, CL.Low.Y_X, 
            CL.Up.Y_X)
        colnames(Values.Y_X) <- c("Estimate", "Std. Error", "t value", 
            "p(>|t|)", "Low Conf Limit", "Up Conf Limit")
        Model.Fit.Y_X <- cbind(RMSE.Y_X, 1, N - 1 - 1, Model.F.Y_X, 
            1 - pf(Model.F.Y_X, 1, N - 1 - 1), R2.Y_X, R2.Y_X.Adj, 
            CI.R2.Y_X$Lower, CI.R2.Y_X$Upper)
        colnames(Model.Fit.Y_X) <- c("Residual standard error (RMSE)", 
            "numerator df", "denomenator df", "F-Statistic", 
            "p-value (F)", "R^2", "Adj R^2", "Low Conf Limit", 
            "Up Conf Limit")
        rownames(Model.Fit.Y_X) <- "Values"
        Regression.of.Y.on.X <- list(Regression.Table = Values.Y_X, 
            Model.Fit = Model.Fit.Y_X)
        B.M_X <- solve(S.XX[1:1, 1:1]) %*% S.XX[2, 1]
        B.M_X <- cbind(mean.m - mean.x * B.M_X, B.M_X)
        colnames(B.M_X) <- c("Intercept.M_X", "a (Regressor)")
        path.a <- B.M_X[2]
        R2.M_X <- (t(S.XX[2, 1]) %*% solve(S.XX[1:1, 1:1]) %*% 
            S.XX[2, 1])/S.XX[2, 2]
        R2.M_X.Adj <- 1 - ((1 - R2.M_X) * ((N - 1)/(N - 1 - 1)))
        CI.R2.M_X <- ci.R2(R2 = R2.M_X, conf.level = conf.level, 
            Random.Predictors = TRUE, N = N, p = 1)
        Model.F.M_X <- (R2.M_X/1)/((1 - R2.M_X)/(N - 1 - 1))
        RMSE.M_X <- sqrt((1 - R2.M_X) * ((N - 1) * S.XX[2, 2]/(N - 
            2)))
        SE.M_X <- cbind(c(sqrt((1 - R2.M_X) * ((N - 1) * S.XX[2, 
            2]/(N - 2))) * sqrt(1/N + mean.x^2/((N - 1) * S.XX[1, 
            1])), sqrt((1 - R2.M_X)/(N - 1 - 1)) * sqrt(S.XX[2, 
            2]/S.XX[1, 1])))
        t.M_X <- t(B.M_X)/SE.M_X
        p.M_X <- 2 * (pt(-1 * abs(t.M_X), df = N - 1 - 1))
        CL.Low.M_X <- t(B.M_X) - qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.M_X
        CL.Up.M_X <- t(B.M_X) + qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.M_X
        Values.M_X <- cbind(t(B.M_X), SE.M_X, t.M_X, p.M_X, CL.Low.M_X, 
            CL.Up.M_X)
        colnames(Values.M_X) <- c("Estimate", "Std. Error", "t value", 
            "p(>|t|)", "Low Conf Limit", "Up Conf Limit")
        Model.Fit.M_X <- cbind(RMSE.M_X, 1, N - 1 - 1, Model.F.M_X, 
            1 - pf(Model.F.M_X, 1, N - 1 - 1), R2.M_X, R2.M_X.Adj, 
            CI.R2.M_X$Lower, CI.R2.M_X$Upper)
        colnames(Model.Fit.M_X) <- c("Residual standard error (RMSE)", 
            "numerator df", "denomenator df", "F-Statistic", 
            "p-value (F)", "R^2", "Adj R^2", "Low Conf Limit", 
            "Up Conf Limit")
        rownames(Model.Fit.M_X) <- "Values"
        Regression.of.M.on.X <- list(Regression.Table = Values.M_X, 
            Model.Fit = Model.Fit.M_X)
        B.Y_XM <- solve(S.XX) %*% s.XY
        B.Y_XM <- t(cbind(c(mean.dv - (mean.x * B.Y_XM[1] + mean.m * 
            B.Y_XM[2]), cbind(B.Y_XM))))
        colnames(B.Y_XM) <- c("Intercept.Y_XM", "c.prime (Regressor)", 
            "b (Mediator)")
        path.c.prime <- B.Y_XM[2]
        path.b <- B.Y_XM[3]
        R2.Y_XM <- (t(s.XY) %*% solve(S.XX) %*% s.XY)/Cov.Matrix[Dim.Cov.Matrix, 
            Dim.Cov.Matrix]
        R2.Y_XM.Adj <- 1 - ((1 - R2.Y_XM) * ((N - 1)/(N - 2 - 
            1)))
        CI.R2.Y_XM <- ci.R2(R2 = R2.Y_XM, conf.level = conf.level, 
            Random.Predictors = TRUE, N = N, p = 2)
        Model.F.Y_XM <- (R2.Y_XM/2)/((1 - R2.Y_XM)/(N - 2 - 1))
        
        # Below is necessary so that rounding doesn't lead to a (very slightly) negative number (via rounding) that would not allow the square root to be taken.
        MSE.Y_XM <- (1 - R2.Y_XM) * ((N - 1) * Cov.Matrix[3, 3]/(N - 3))
        RMSE.Y_XM <- sqrt(round(MSE.Y_XM, 10))
        
        x.prime.x <- cbind(c(N, mean.x * N, mean.m * N), c(mean.x * 
            N, (S.XX[1, 1] * (N - 1) + N * mean.x^2), S.XX[1, 
            2] * (N - 1) + N * mean.x * mean.m), c(mean.m * N, 
            S.XX[1, 2] * (N - 1) + N * mean.x * mean.m, (S.XX[2, 
                2] * (N - 1) + N * mean.m^2)))
        SE.Y_XM <- cbind(sqrt(diag(solve(x.prime.x))) * RMSE.Y_XM)
        t.Y_XM <- t(B.Y_XM)/SE.Y_XM
        p.Y_XM <- 2 * (pt(-1 * abs(t.Y_XM), df = N - 2 - 1))
        CL.Low.Y_XM <- t(B.Y_XM) - qt(1 - (1 - conf.level)/2, 
            df = N - 2 - 1) * SE.Y_XM
        CL.Up.Y_XM <- t(B.Y_XM) + qt(1 - (1 - conf.level)/2, 
            df = N - 2 - 1) * SE.Y_XM
        Values.Y_XM <- cbind(t(B.Y_XM), SE.Y_XM, t.Y_XM, p.Y_XM, 
            CL.Low.Y_XM, CL.Up.Y_XM)
        colnames(Values.Y_XM) <- c("Estimate", "Std. Error", 
            "t value", "p(>|t|)", "Low Conf Limit", "Up Conf Limit")
        Model.Fit.Y_XM <- cbind(RMSE.Y_XM, 2, N - 2 - 1, Model.F.Y_XM, 
            1 - pf(Model.F.Y_XM, 2, N - 2 - 1), R2.Y_XM, R2.Y_XM.Adj, 
            CI.R2.Y_XM$Lower, CI.R2.Y_XM$Upper)
        colnames(Model.Fit.Y_XM) <- c("Residual standard error (RMSE)", 
            "numerator df", "denomenator df", "F-Statistic", 
            "p-value (F)", "R^2", "Adj R^2", "Low Conf Limit", 
            "Up Conf Limit")
        rownames(Model.Fit.Y_XM) <- "Values"
        Regression.of.Y.on.X.and.M <- list(Regression.Table = Values.Y_XM, 
            Model.Fit = Model.Fit.Y_XM)
        s2.X <- Cov.Matrix[1, 1]
        s2.M <- Cov.Matrix[2, 2]
        s2.Y <- Cov.Matrix[3, 3]
        s.YX <- Cov.Matrix[1, 3]
        s.XM <- Cov.Matrix[2, 1]
        s.YM <- Cov.Matrix[3, 2]
        
        ab <- path.a * path.b
        
a.contained <- c((s.YM * s.YX + sqrt(s2.M * s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y), 
                 (s.YM * s.YX - sqrt(s2.M * s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y))
        
b.contained <- c(sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * s2.M - s.XM^2), 
                -sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * s2.M - s.XM^2)) 

        
# Determins the M(.) from Preacher and Kelley (2011)

M.ab <- 0

if(round(ab, 10)!=0)
  {
  # To find M(a)
  if(path.a > 0) 
  {
to.get.M.a <- a.contained[a.contained > 0]
  M.a <- max(to.get.M.a)
  }
  if(path.a < 0) 
  {
to.get.M.a <- a.contained[a.contained < 0]
  M.a <- min(to.get.M.a)
  }
  if(path.a == 0) 
  {
  M.a <- 0
  }
  # Now to find M(b)  
  if(path.b > 0) 
  {
  to.get.M.b <- b.contained[b.contained > 0]
  M.b <- max(to.get.M.b)
  }
  if(path.b < 0) 
  {
  to.get.M.b <- b.contained[b.contained < 0]
  M.b <- min(to.get.M.b)
  }
  if(path.b == 0) 
  {
  M.b <- 0
  }  
  M.ab <- M.a*M.b
  } 
        
        Indirect.Effect <- c(Estimate = ab)
        
        Indirect.Effect.Partially.Standardized <- c(Estimate = ab/sqrt(Cov.Matrix[3, 3]))
        
        Index.of.Mediation <- c(Estimate = ab * (sqrt(Cov.Matrix[1, 1])/sqrt(Cov.Matrix[3, 3])))
        
        R2_4.5 <- c(Estimate = (Cor.Matrix[3, 2]^2) - (R2.Y_XM-R2.Y_X)) # Equation 13b from P&K
        
        R2_4.6 <- c(Estimate = ifelse(R2.Y_X==1, 0, (R2.M_X * (R2.Y_XM - R2.Y_X)/(1-R2.Y_X)))) # Equation 14b from P&K.
        
        R2_4.7 <- c(Estimate = ifelse(R2.Y_XM==0, 0, ((R2.M_X * (R2.Y_XM - R2.Y_X))/((1-R2.Y_X)*R2.Y_XM)))) # Equation 15b from P&K.

        Maximum.Possible.Mediation.Effect <- c(Estimate = M.ab)
        
        ab.to.Maximum.Possible.Mediation.Effect_kappa.squared <- c(Estimate = ifelse(M.ab==0, 0, ab/M.ab))
       
        Ratio.of.Indirect.to.Total.Effect <- c(Estimate = ifelse(path.c==0, 0, (1 - (path.c.prime/path.c))))
        
        Ratio.of.Indirect.to.Direct.Effect <- c(Estimate = ab/path.c.prime) # Under regular situations.
        if(path.a*path.b == 0) Ratio.of.Indirect.to.Direct.Effect <- c(Estimate = 0) # If the numerator is 0
        if(path.a*path.b != 0 & path.c.prime==0) Ratio.of.Indirect.to.Direct.Effect <- c(Estimate = sign(ab)*999999999999999999) # If the denominator is 0 but the numerator is not (+/- "infinity").

        Ratio.of.Indirect.to.Direct.Effect <- c(Estimate = ifelse(path.a*path.b == 0, 0, ab/path.c.prime))
        Ratio.of.Indirect.to.Direct.Effect <- c(Estimate = ifelse(path.a*path.b == 0, 0, ab/path.c.prime))
        
        
        Success.of.Surrogate.Endpoint <- c(Estimate = ifelse(path.a==0, 999999999999999999, path.c/path.a))
        
        SOS <- c(Estimate = ifelse(R2.Y_X==0, 0, c(Estimate = R2_4.5/R2.Y_X)))
        
        if (!is.null(S)) {
          
          Effect.Sizes <- rbind(
                Indirect.Effect=Indirect.Effect, 
                Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                Index.of.Mediation = Index.of.Mediation, 
                R2_4.5 = R2_4.5, 
                R2_4.6 = R2_4.6, 
                R2_4.7 = R2_4.7, 
                Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, 
                Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                SOS = SOS)
          
            Results.mediation <- list(
                Y.on.X = Regression.of.Y.on.X, 
                M.on.X = Regression.of.M.on.X, 
                Y.on.X.and.M = Regression.of.Y.on.X.and.M, 
                Effect.Sizes=Effect.Sizes)
        }
        if (is.null(S)) {
            if (sum(x == 0 | x == 1) != N) {
              
              
              Effect.Sizes <- rbind(
                  Indirect.Effect = Indirect.Effect, 
                  Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                  Index.of.Mediation = Index.of.Mediation, 
                  R2_4.5 = R2_4.5, 
                  R2_4.6 = R2_4.6, 
                  R2_4.7 = R2_4.7, 
                  Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                  ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                  Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, 
                  Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                  Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                  Residual.Based_Gamma = c(Estimate=Residual.Based_Gamma), 
                  Residual.Based.Standardized_gamma = c(Estimate=Residual.Based.Standardized_gamma), 
                  SOS = SOS)
              
                Results.mediation <- list(
                  Y.on.X = Regression.of.Y.on.X, 
                  M.on.X = Regression.of.M.on.X, 
                  Y.on.X.and.M = Regression.of.Y.on.X.and.M,
                  Effect.Sizes = Effect.Sizes)
                
                 
            }
            if (sum(x == 0 | x == 1) == N) {
              
                ES <- c(Estimate = ifelse((path.a==0 & path.b==0), 0, (((path.a * path.b)/(sqrt(SE.Y_XM[3]^2 * path.a^2 + SE.M_X[2]^2 * path.b^2))) * (sqrt(1/sum(x == 0) + 1/sum(x == 1))))))

                  Effect.Sizes <- rbind(
                  Indirect.Effect = Indirect.Effect, 
                  Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                  Index.of.Mediation = Index.of.Mediation, 
                  R2_4.5 = R2_4.5, 
                  R2_4.6 = R2_4.6, 
                  R2_4.7 = R2_4.7, 
                  Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                  ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                  Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, 
                  Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                  Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                  Residual.Based_Gamma = c(Estimate=Residual.Based_Gamma), 
                  Residual.Based.Standardized_gamma = c(Estimate=Residual.Based.Standardized_gamma), 
                  ES.for.two.groups = ES, 
                  SOS = SOS)
                
                Results.mediation <- list(
                  Y.on.X = Regression.of.Y.on.X, 
                  M.on.X = Regression.of.M.on.X, 
                  Y.on.X.and.M = Regression.of.Y.on.X.and.M, 
                  Effect.Sizes = Effect.Sizes)
                
                
                  
            }
        }
        
        return(Results.mediation)
    }
  
# Here is the internal bootstrap function.
.mediation.bs <- function(x, mediator, dv, S=S, conf.level = conf.level, B = B, which.boot="all", ...)
{
  
  if(!requireNamespace("boot", quietly = TRUE)) stop("The package 'boot' is needed; please install the package and try again.")
  
  if(!is.null(S)) stop("For the bootstrap procedures to be implemented, you must supply raw data (i.e., not a covariance matrix).")
       
Data <- na.omit(cbind(x, mediator, dv))
N <- dim(Data)[1]

Effect.Size.Names <- rownames(Result$Effect.Sizes)

       Boot.This <- function(Data, g) 
         {
          # Applies the internal mediation function to the data for the rows identified by g.
          Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3], S = NULL, conf.level = conf.level)
          as.numeric(Output.mediation$Effect.Sizes) # Returns only the values of effect size estimates
         }
        
print("Bootstrap resampling has begun. This process may take a considerable amount of time if the number of replications is large, which is optimal for the bootstrap procedure.")

boot.out <- boot::boot(data = Data, statistic=Boot.This, R = B, stype = "i")

if(save.bs.replicates==TRUE)
{
print("You requested saving the bootstrap replicates; they are available in the object \'Bootstrap.Replicates\' in the workspace.")
Bootstrap.Replicates <- as.data.frame(boot.out$t)
colnames(Bootstrap.Replicates) <- Effect.Size.Names
}
    
Number.of.Effect.Sizes <- length(Effect.Size.Names)
Get.CIs.for.These <- c(1:Number.of.Effect.Sizes) # All estimates of effect size (this is the index for each submitted to a loop)

if(which.boot=="both" | which.boot=="BOTH" | which.boot=="Both" | which.boot=="all" | which.boot=="ALL" | which.boot=="All") which.boot <- "Percentile.and.BCa"
if(which.boot=="perc" | which.boot=="Percentile" | which.boot=="Perc" | which.boot=="percentile" | which.boot=="PERCENTILE") which.boot <- "Percentile"     
if(which.boot=="bAc" | which.boot=="baC" | which.boot=="BCa" | which.boot=="Bca" | which.boot=="bca" | which.boot=="BCA") which.boot <- "BCa"     
 
Percentile.BS.Results <- matrix(NA, Number.of.Effect.Sizes, 3)
colnames(Percentile.BS.Results) <- c("Estimate", "CI.Lower_Percentile", "CI.Upper_Percentile") 
rownames(Percentile.BS.Results) <- Effect.Size.Names

BCa.BS.Results <- matrix(NA, Number.of.Effect.Sizes, 3)
colnames(BCa.BS.Results) <- c("Estimate", "CI.Lower_BCa", "CI.Upper_BCa") 
rownames(BCa.BS.Results) <- Effect.Size.Names

for(k in 1:Number.of.Effect.Sizes)
{
  
if(which.boot=="Percentile.and.BCa" | which.boot=="Percentile")
{
BS.Results <- try(boot::boot.ci(boot.out = boot.out, index=Get.CIs.for.These[k], conf = conf.level, type = c("perc"))$percent[4:5], silent=TRUE)
is.numeric(BS.Results)
if(is.numeric(BS.Results)==FALSE) BS.Results <- c(NA, NA)
Percentile.BS.Results[k, 1:3] <- c(Result$Effect.Sizes[k], BS.Results)
}


if(which.boot=="Percentile.and.BCa" | which.boot=="BCa") 
{
BS.Results <- c(NA, NA)
BS.Results <- try(boot::boot.ci(boot.out = boot.out, index=Get.CIs.for.These[k], conf = conf.level, type = c("bca"))$bca[4:5], silent=TRUE)
if(is.numeric(BS.Results)==FALSE) BS.Results <- c(NA, NA)
BCa.BS.Results[k, 1:3] <- c(Result$Effect.Sizes[k], BS.Results)
}
}

if(which.boot=="Percentile.and.BCa") BS.Output <- cbind(Percentile.BS.Results, BCa.BS.Results[,2:3])
if(which.boot=="Percentile") BS.Output <- Percentile.BS.Results
if(which.boot=="BCa") BS.Output <- BCa.BS.Results

return(BS.Output) # Bootstrap statistics are output here.
}
    
# Result if no bootstrapping is done.
Result <-   .mediation(x = x, mediator = mediator, dv = dv, S = S, N = N, x.location.S = x.location.S, mediator.location.S = mediator.location.S, 
            dv.location.S = dv.location.S, mean.x = mean.x, mean.m = mean.m, mean.dv = mean.dv, conf.level = conf.level)

# Result if bootstrapping is done.
if (bootstrap == TRUE) 
{
BS.Result <- .mediation.bs(x=x, mediator=mediator, dv=dv, S=S, conf.level=conf.level, B = B, which.boot=which.boot)
Result <- list(Y.on.X=Result$Y.on.X, M.on.X=Result$M.on.X, Y.on.X.and.M=Result$Y.on.X.and.M, Bootstrap.Results=BS.Result)
}


# Added to NOT report kappa squared by default (user has to request it specifically using the argument: kappa.squared=TRUE)
if(complete.set==FALSE)
{
if(bootstrap==FALSE)
{
Result <- cbind(Result$Effect.Sizes[!rownames(Result$Effect.Sizes)%in%c("Maximum.Possible.Mediation.Effect", "ab.to.Maximum.Possible.Mediation.Effect_kappa.squared"), 1])
colnames(Result) <- "Estimate"
}
if(bootstrap==TRUE)
{
col.names <- colnames(Result)
Result <- cbind(Result$Bootstrap.Results[!rownames(Result$Bootstrap.Results)%in%c("Maximum.Possible.Mediation.Effect", "ab.to.Maximum.Possible.Mediation.Effect_kappa.squared"), ])
}
}


return(Result)
}
