

cat("Running Paik-fit.R\n")

OncotypeRScutoffs = c(18, 30)  ### intermediate risk boundaries
TailorXRScutoffs = c(11, 25)  ### intermediate risk boundaries

pointStrings = grep(value=T, pattern = "<point",
                    read.csv(header=F,
                             "inst/Paik2006-10yearDFS.TAM.xml",
                             stringsAsFactors = F)[[1]]
)
dxBegin = regexpr("dx='", pointStrings)+4
dyBegin = regexpr("dy='", pointStrings)+4
dxEnd = dyBegin - 7
dyEnd = nchar(pointStrings) - 4
allX_TAM = as.numeric(substr(pointStrings, dxBegin, dxEnd))
allY_TAM = as.numeric(substr(pointStrings, dyBegin, dyEnd))


### Create pointStrings ####
pointStrings = grep(value=T, pattern = "<point",
                    read.csv(header=F,
                             "inst/Paik2006-10yearDFS.TAM+CHEMO.xml",
                             stringsAsFactors = F)[[1]]
)
dxBegin = regexpr("dx='", pointStrings)+4
dyBegin = regexpr("dy='", pointStrings)+4
dxEnd = dyBegin - 7
dyEnd = nchar(pointStrings) - 4
allX_TAM_CHEMO = as.numeric(substr(pointStrings, dxBegin, dxEnd))
allY_TAM_CHEMO = as.numeric(substr(pointStrings, dyBegin, dyEnd))
plot(allX_TAM, allY_TAM)
points(allX_TAM_CHEMO, allY_TAM_CHEMO)

lm_TenYearDFS_TAM = lm(allY_TAM ~ poly(allX_TAM,2))
lm_TenYearDFS_TAM_CHEMO = lm(allY_TAM_CHEMO ~ poly(allX_TAM_CHEMO,2))
tenYearDFS = data.frame(RS=c(allX_TAM, allX_TAM_CHEMO),
                        Recur=c(allY_TAM, allY_TAM_CHEMO),
                        group=rep(c('TAM', 'TAM_CHEMO'),
                                  times=c(length(allX_TAM), length(allX_TAM_CHEMO)))
)
tenYearDFS$RS = pmax(0, tenYearDFS$RS)
tenYearDFS$one = 1
tenYearDFS$RS2 = tenYearDFS$RS^2
lm_TenYearDFS = lm(data=tenYearDFS,
                   Recur ~ group*RS + group*RS2 - group
) ### shared intercept

predict(lm_TenYearDFS, newdata=
          data.frame(RS=c(0,0), RS2=c(0,0), group=c("TAM", "TAM_CHEMO"), one=c(1,1)))
RSvector = seq(0,100)
nRS = length(RSvector)
tenYearDFS_long = data.frame(RS=rep(RSvector,2), RS2=rep((RSvector)^2,2),
                             group=rep(c('TAM','TAM_CHEMO'), each=nRS), one=1)
predicted_TenYearDFS = predict( lm_TenYearDFS, newdata=tenYearDFS_long)
tenYearDFS_long$predicted = predicted_TenYearDFS

############# RS distribution modeling, it to match the categories   ##############
ProportionByRiskCat = c(51, 22, 27)/100
sampleSize = 668
### Let's try a gamma fit.
# shapescale = nlm(function(shapescale)
#   (pgamma(OncotypeRScutoffs[1], shape=shapescale[1], scale=shapescale[2]) - 0.51)^2 +
#       (pgamma(OncotypeRScutoffs[2], shape=shapescale[1], scale=shapescale[2]) -
#          pgamma(OncotypeRScutoffs[2], shape=shapescale[1], scale=shapescale[2]) - 0.22)^2 ,
#   p=c(shape=5, scale=17/5) # starting values
# )$estimate
# pgamma(OncotypeRScutoffs[1], shape=shapescale[1], scale=shapescale[2])
# pgamma(OncotypeRScutoffs[2], shape=shapescale[1], scale=shapescale[2])
#Perfect.  But generates long tail values.
### Let's try a beta fit.
abParam = nlm(function(abParam)
  (pbeta(OncotypeRScutoffs[1]/50, shape1=abParam[1], shape2=abParam[2]) - 0.51)^2 +
    (pbeta(OncotypeRScutoffs[2]/50, shape1=abParam[1], shape2=abParam[2]) -
       pbeta(OncotypeRScutoffs[1]/50, shape1=abParam[1], shape2=abParam[2]) - 0.22)^2 ,
  p=c(shape1=2, shape2=2)
)$estimate
pbeta(OncotypeRScutoffs[1]/50, shape1=abParam[1], shape2=abParam[2])
pbeta(OncotypeRScutoffs[2]/50, shape1=abParam[1], shape2=abParam[2])
#Perfect.

benefit =
  (-1)*apply(matrix(predicted_TenYearDFS, ncol=2), 1, diff)
tenYearDFS_long$benefit = c(rep(0, length(benefit)),
                            benefit)
tenYearDFS_long_treated = tenYearDFS_long[tenYearDFS_long$group=="TAM_CHEMO", ]
tenYearDFS_long_treated$control = tenYearDFS_long$predicted[1:101]
tenYearDFS_long_treated$benefit2 =  ### OK; the same.
  tenYearDFS_long_treated$control -
  tenYearDFS_long_treated$predicted

PaikSampleSize <- 668


RSsample = 50 * rbeta(sampleSize, shape1=abParam[1], shape2=abParam[2])
RSsample = sort(RSsample)
RSsampleBenefit = tenYearDFS_long_treated$benefit[
  match(ceiling(RSsample), tenYearDFS_long_treated$RS)]
RSsampleBenefit = pmax(0, RSsampleBenefit)
#   lines(RSsample, RSsampleBenefit, col=Benefitcolor, lwd=3)


# RSsampleBenefit = tenYearDFS_long_treated$benefit[]
whichBenefitted = which(1==rbinom(sampleSize, 1, RSsampleBenefit))
theseBenefitted = RSsample[whichBenefitted]
### Who benefits from T+C?   ##

benefitTable = table( RSsample %in% theseBenefitted,
                      cut(RSsample, c(0, OncotypeRScutoffs, Inf )))


# RSforNNTupper = NNTbiomarker::argmin(nnt, OncotypeNNTrange[1])
# RSforNNTlower = NNTbiomarker::argmin(nnt, OncotypeNNTrange[2])

#### NUMBER NEEDED TO TREAT   ####

nnt = 1/benefit
nnt[nnt<=0] = Inf
names(nnt) = RSvector
Paik_nnt = nnt

######   Marginal benefit  ####

RSsampleBtoT = rbinom(length(RSsampleBenefit), 1, RSsampleBenefit)
prevalence = sum(RSsampleBtoT)/length(RSsampleBtoT)  ### 3% benefit


benefitTable = table( RSsample %in% theseBenefitted,
                      cut(RSsample, c(0, OncotypeRScutoffs, Inf )))


