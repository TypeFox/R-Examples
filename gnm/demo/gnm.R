message("1. Set seed as gnm returns random parameterization")
set.seed(1)

{
    if (interactive()) {
        cat("\n3. Type <Return> to fit (linear) uniform association model,  ",
            "\n   using Diag() to fit diagonal effects: ")
        readline()
    }
    else
        message("2. Fit (linear) uniform association model, using Diag() to fit",
                "   diagonal effects")
}
Rscore <- scale(as.numeric(row(occupationalStatus)), scale = FALSE)
Cscore <- scale(as.numeric(col(occupationalStatus)), scale = FALSE)
Uniform <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               Rscore:Cscore, family = poisson, data = occupationalStatus)
summary(Uniform)
{
    if (interactive()) {
        cat("\n3. Type <Return> to fit an association model using Mult() to fit",
            "\n   separate row and column effects:")
        readline()
    }
    else message("3. Fit an association model using Mult() to fit separate row and ",
                 "column effects")
}
RC <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
          Mult(origin, destination), family = poisson,
          data = occupationalStatus)
summary(RC)
{
    if (interactive()) {
        cat("\n4. Type <Return> to fit an association model using MultHomog()",
            "\n   to fit homogeneous row-column effects:")
        readline()
    }
    else message("4. Fit an association model using MultHomog()\n",
                 "to fit homogeneous row-column effects")
}
RChomog <- gnm(Freq ~ origin + destination + Diag(origin, destination) +
               MultHomog(origin, destination), family = poisson,
               data = occupationalStatus)
summary(RChomog)
{
    if (interactive()) {
        cat("\n5. Type <Return> to compare models using anova:")
        readline()
    }
    else message("5. Compare models using anova")
}
anova(Uniform, RChomog, RC)

message("6. Produce diagnostic plots for RChomog")
plot(RChomog)

message("7. Get simple constrasts of homogeneous row-column effects")
getContrasts(RChomog, grep("MultHomog", names(coef(RChomog))))

message("End of demo. \n",
 "See vignette(\"gnmOverview\", package = \"gnm\") for full manual.")
