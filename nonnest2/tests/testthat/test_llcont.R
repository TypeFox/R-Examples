context("sum of llcont")

test_that("lavaan object", {
  if (isTRUE(require("lavaan"))) {
    HS.model <- 'visual  =~ x1 + x2 + x3
                 textual =~ x4 + x5 + x6
                 speed   =~ x7 + x8 + x9 '
    fit1 <- cfa(HS.model, data=HolzingerSwineford1939)
    fit2 <- cfa(HS.model, data=HolzingerSwineford1939, group="school")
    fit3 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
                group.equal = c("loadings"))
    fit4 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
                group.equal = c("loadings"),
                group.partial = c("visual=~x2", "x7~1"))

    expect_equal(sum(llcont(fit1)), as.numeric(logLik(fit1)))
    expect_equal(sum(llcont(fit2)), as.numeric(logLik(fit2)))
    expect_equal(sum(llcont(fit3)), as.numeric(logLik(fit3)))
    expect_equal(sum(llcont(fit4)), as.numeric(logLik(fit4)))
    
    HS.model2 <- 'visual =~ x1 + 0.5*x2 + c(0.6, 0.8)*x3
                  textual =~ x4 + start(c(1.2, 0.6))*x5 + a*x6
                  speed   =~ x7 + x8 + x9'
    fit5 <- cfa(HS.model2, data=HolzingerSwineford1939, group="school")
    
    expect_equal(sum(llcont(fit5)), as.numeric(logLik(fit5)))
    
    HS.model3 <- 'visual =~ x1 + x2 + c(v3,v3)*x3
                  textual =~ x4 + x5 + x6
                  speed   =~ x7 + x8 + x9'
    fit6 <- cfa(HS.model3, data=HolzingerSwineford1939, group="school")
    expect_equal(sum(llcont(fit6)), as.numeric(logLik(fit6)))
  }
})


test_that("glm object", {
  if (isTRUE(require("faraway")) && isTRUE(require("MASS"))) {
    ## binomial
    bin1 <- glm(formula=am ~ hp + wt, data=mtcars, family=binomial)
    bin2 <- glm(cbind(Menarche, Total-Menarche) ~ Age,
                family=binomial(logit), data=menarche)

    expect_equal(sum(llcont(bin1)), as.numeric(logLik(bin1)))
    expect_equal(sum(llcont(bin2)), as.numeric(logLik(bin2)))

    ## quasibinomial
    qbin1 <- glm(formula=am ~ hp + wt, data=mtcars, family=quasibinomial)
    qbin2 <- glm(cbind(Menarche, Total-Menarche) ~ Age,
                family=quasibinomial, data=menarche)

    expect_equal(sum(llcont(qbin1)), as.numeric(logLik(qbin1)))
    expect_equal(sum(llcont(qbin2)), as.numeric(logLik(qbin2)))

    ## gaussian
    gau1 <- glm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent,
                data=gala, family=gaussian)
    gau2 <- glm(Species ~ Area + Elevation + Nearest, data=gala,
                family=gaussian)

    expect_equal(sum(llcont(gau1)), as.numeric(logLik(gau1)))
    expect_equal(sum(llcont(gau2)), as.numeric(logLik(gau2)))

    ## inverse.gaussian
    invGau1 <- glm(actual ~ projected-1,
                   family=inverse.gaussian(link="identity"), cpd)

    expect_equal(sum(llcont(invGau1)), as.numeric(logLik(invGau1)))

    ## Gamma
    clotting <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
                           lot1 = c(118,58,42,35,27,25,21,19,18),
                           lot2 = c(69,35,26,21,18,16,13,12,12))
    gam1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma)

    expect_equal(sum(llcont(gam1)), as.numeric(logLik(gam1)))

    ## poisson
    counts <- c(18,17,15,20,10,20,25,13,12)
    outcome <- gl(3,1,9)
    treatment <- gl(3,3)
    d.AD <- data.frame(treatment, outcome, counts)
    pois1 <- glm(counts ~ outcome + treatment, family = poisson)
    pois2 <- glm(counts ~ outcome, family = poisson)

    expect_equal(sum(llcont(pois1)), as.numeric(logLik(pois1)))
    expect_equal(sum(llcont(pois2)), as.numeric(logLik(pois2)))

    ## quasipoisson
    qpois1 <- glm(counts ~ outcome + treatment, family = quasipoisson)
    qpois2 <- glm(counts ~ outcome, family = quasipoisson)

    expect_equal(sum(llcont(qpois1)), as.numeric(logLik(qpois1)))
    expect_equal(sum(llcont(qpois2)), as.numeric(logLik(qpois2)))

    ## negative-binomial: MASS
    nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)

    expect_equal(sum(llcont(nb1)), as.numeric(logLik(nb1)))
  }
})


test_that("clm object", {
  if (isTRUE(require("ordinal")) && isTRUE(require("MASS"))) {
    clm1 <- clm(rating ~ temp * contact, data = wine)
    clm2 <- update(clm1, ~.-temp:contact)
    clm3 <- update(clm1, link = "logit")
    clm4 <- update(clm1, link = "probit")
    clm5 <- update(clm1, link = "loglog")
    clm6 <- update(clm1, link = "cloglog")
    clm7 <- update(clm1, link = "cauchit")
    clm8 <- update(clm1, threshold = "symmetric")
    clm9 <- update(clm1, threshold = "equidistant")
    clm10 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

    expect_equal(sum(llcont(clm1)), as.numeric(logLik(clm1)))
    expect_equal(sum(llcont(clm2)), as.numeric(logLik(clm2)))
    expect_equal(sum(llcont(clm3)), as.numeric(logLik(clm3)))
    expect_equal(sum(llcont(clm4)), as.numeric(logLik(clm4)))
    expect_equal(sum(llcont(clm5)), as.numeric(logLik(clm5)))
    expect_equal(sum(llcont(clm6)), as.numeric(logLik(clm6)))
    expect_equal(sum(llcont(clm7)), as.numeric(logLik(clm7)))
    expect_equal(sum(llcont(clm8)), as.numeric(logLik(clm8)))
    expect_equal(sum(llcont(clm9)), as.numeric(logLik(clm9)))
    expect_equal(sum(llcont(clm10)), as.numeric(logLik(clm10)))
  }
})


test_that("hurdle object", {
  if (isTRUE(require("pscl"))) {
    hurdle1 <- hurdle(formula = art ~ ., data = bioChemists)
    hurdle2 <- hurdle(formula = art ~ ., data = bioChemists, separate=FALSE)
    hurdle3 <- hurdle(art ~ ., data = bioChemists, zero = "geometric")
    hurdle4 <- hurdle(art ~ fem + ment, data = bioChemists,
                      dist = "negbin", zero = "negbin")
    hurdle5 <- hurdle(art ~ ., data = bioChemists, dist = "negbin")

    expect_equal(sum(llcont(hurdle1)), as.numeric(logLik(hurdle1)))
    expect_equal(sum(llcont(hurdle2)), as.numeric(logLik(hurdle2)))
    expect_equal(sum(llcont(hurdle3)), as.numeric(logLik(hurdle3)))
    expect_equal(sum(llcont(hurdle4)), as.numeric(logLik(hurdle4)))
    expect_equal(sum(llcont(hurdle5)), as.numeric(logLik(hurdle5)))
  }
})


test_that("zeroinfl object", {
  if (isTRUE(require("pscl"))) {
    zi1 <- zeroinfl(art ~ . | 1, data = bioChemists)
    zi2 <- zeroinfl(art ~ . | 1, data = bioChemists, dist = "negbin")
    zi3 <- zeroinfl(art ~ . | ., data = bioChemists)
    zi4 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")

    expect_equal(sum(llcont(zi1)), as.numeric(logLik(zi1)))
    expect_equal(sum(llcont(zi2)), as.numeric(logLik(zi2)))
    expect_equal(sum(llcont(zi3)), as.numeric(logLik(zi3)))
    expect_equal(sum(llcont(zi4)), as.numeric(logLik(zi4)))
  }
})


test_that("lm object", {
  ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
  trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
  group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
  weight <- c(ctl, trt)
  lm1 <- lm(weight ~ group)

  expect_equal(sum(llcont(lm1)), as.numeric(logLik(lm1)))
})


test_that("mlogit object", {
  if (isTRUE(require("mlogit")) & isTRUE(require("AER"))) {
    data("Fishing", package = "mlogit")
    Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide",
                        choice = "mode")
    mlog1 <- mlogit(mode ~ price + catch, data = Fish)
    mlog2 <- mlogit(mode ~ 0 | income, data = Fish)
    mlog3 <- mlogit(mode ~ price+ catch | income, data = Fish)
    mlog4 <- mlogit(mode ~ price+ catch | income, data = Fish,
                    reflevel = "charter")
    mlog5 <- mlogit(mode ~ price+ catch | income, data = Fish,
                    alt.subset = c("charter", "pier", "beach"))
    Fishing2 <- Fishing
    Fishing2[1, "price.pier"] <- Fishing2[3, "price.beach"] <- NA
    mlog6 <- mlogit(mode~price+catch|income, Fishing2, shape="wide",
                    choice="mode", varying = 2:9)
    data("TravelMode", package = "AER")
    Tr2 <- TravelMode[-c(2, 7, 9),]
    mlog7 <- mlogit(choice~wait+gcost|income+size, Tr2, shape = "long",
                    chid.var = "individual", alt.var="mode", choice = "choice")
    data("TravelMode", package = "AER")
    mlog8 <- mlogit(choice ~ wait + travel + vcost, TravelMode,
                    shape = "long", chid.var = "individual",
                    alt.var = "mode",
                    method = "bfgs", heterosc = TRUE, tol = 10)
    TravelMode$avincome <- with(TravelMode, income * (mode == "air"))
    TravelMode$time <- with(TravelMode, travel + wait)/60
    TravelMode$timeair <- with(TravelMode, time * I(mode == "air"))
    TravelMode$income <- with(TravelMode, income / 10)
    TravelMode$incomeother <- with(TravelMode,
                                   ifelse(mode %in% c('air', 'car'),
                                          income, 0))
    mlog9 <- mlogit(choice~gcost+wait+incomeother, TravelMode,
                    shape='long', alt.var='mode',
                    nests=list(public=c('train', 'bus'),
                        other=c('car','air')))
    data("Game", package = "mlogit")
    mlog10 <- mlogit(ch~own|hours, Game, choice='ch', varying = 1:12,
                     ranked=TRUE, shape="wide", reflevel="PC")

    expect_equal(sum(llcont(mlog1)), as.numeric(logLik(mlog1)))
    expect_equal(sum(llcont(mlog2)), as.numeric(logLik(mlog2)))
    expect_equal(sum(llcont(mlog3)), as.numeric(logLik(mlog3)))
    expect_equal(sum(llcont(mlog4)), as.numeric(logLik(mlog4)))
    expect_equal(sum(llcont(mlog5)), as.numeric(logLik(mlog5)))
    expect_equal(sum(llcont(mlog6)), as.numeric(logLik(mlog6)))
    expect_equal(sum(llcont(mlog7)), as.numeric(logLik(mlog7)))
    expect_equal(sum(llcont(mlog8)), as.numeric(logLik(mlog8)))
    expect_equal(sum(llcont(mlog9)), as.numeric(logLik(mlog9)))
    expect_equal(sum(llcont(mlog10)), as.numeric(logLik(mlog10)))
  }
})


test_that("nls object", {
  DNase1 <- subset(DNase, Run == 1)
  nls1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
  nls2 <- nls(density ~ 1/(1 + exp((xmid - log(conc))/scal)),
              data = DNase1,
              start = list(xmid = 0, scal = 1),
              algorithm = "plinear")
  nls3 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
              data = DNase1,
              start = list(Asym = 3, xmid = 0, scal = 1))
  nls4 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
              data = DNase1,
              start = list(Asym = 3, xmid = 0, scal = 1),
              algorithm = "port")
  Treated <- Puromycin[Puromycin$state == "treated", ]
  weighted.MM <- function(resp, conc, Vm, K) {
    pred <- (Vm * conc)/(K + conc)
    (resp - pred) / sqrt(pred)
  }
  nls5 <- nls( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
              start = list(Vm = 200, K = 0.1))
  lisTreat <- with(Treated,
                   list(conc1 = conc[1], conc.1 = conc[-1], rate = rate))
  weighted.MM1 <- function(resp, conc1, conc.1, Vm, K) {
    conc <- c(conc1, conc.1)
    pred <- (Vm * conc)/(K + conc)
    (resp - pred) / sqrt(pred)
  }
  nls6 <- nls( ~ weighted.MM1(rate, conc1, conc.1, Vm, K),
              data = lisTreat, start = list(Vm = 200, K = 0.1))
  weighted.MM.grad <- function(resp, conc1, conc.1, Vm, K) {
    conc <- c(conc1, conc.1)
    K.conc <- K+conc
    dy.dV <- conc/K.conc
    dy.dK <- -Vm*dy.dV/K.conc
    pred <- Vm*dy.dV
    pred.5 <- sqrt(pred)
    dev <- (resp - pred) / pred.5
    Ddev <- -0.5*(resp+pred)/(pred.5*pred)
    attr(dev, "gradient") <- Ddev * cbind(Vm = dy.dV, K = dy.dK)
    dev
  }
  nls7 <- nls( ~ weighted.MM.grad(rate, conc1, conc.1, Vm, K),
              data = lisTreat, start = list(Vm = 200, K = 0.1))
  if(isTRUE(require("MASS"))){
      utils::data(muscle, package = "MASS")
      nls9 <- nls(Length ~ cbind(1, exp(-Conc/th)), muscle,
                  start = list(th = 1), algorithm = "plinear")
      b <- coef(nls9)
      nls10 <- nls(Length ~ a[Strip] + b[Strip]*exp(-Conc/th), muscle,
                   start = list(a = rep(b[2], 21), b = rep(b[3], 21),
                       th = b[1]))

      expect_equal(sum(llcont(nls9)), as.numeric(logLik(nls9)))
      expect_equal(sum(llcont(nls10)), as.numeric(logLik(nls10)))
  }
      
  expect_equal(sum(llcont(nls1)), as.numeric(logLik(nls1)))
  expect_equal(sum(llcont(nls2)), as.numeric(logLik(nls2)))
  expect_equal(sum(llcont(nls3)), as.numeric(logLik(nls3)))
  expect_equal(sum(llcont(nls4)), as.numeric(logLik(nls4)))
  expect_equal(sum(llcont(nls5)), as.numeric(logLik(nls5)))
  expect_equal(sum(llcont(nls6)), as.numeric(logLik(nls6)))
  expect_equal(sum(llcont(nls7)), as.numeric(logLik(nls7)))
})


test_that("polr object", {
  if (isTRUE(require("MASS"))) {
    options(contrasts = c("contr.treatment", "contr.poly"))
    polr1 <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    polr2 <- update(polr1, method = "probit", Hess = TRUE)
    polr3 <- update(polr1, method = "loglog", Hess = TRUE)
    polr4 <- update(polr1, method = "cloglog", Hess = TRUE)

    expect_equal(sum(llcont(polr1)), as.numeric(logLik(polr1)))
    expect_equal(sum(llcont(polr2)), as.numeric(logLik(polr2)))
    expect_equal(sum(llcont(polr3)), as.numeric(logLik(polr3)))
    expect_equal(sum(llcont(polr4)), as.numeric(logLik(polr4)))
  }
})


test_that("rlm object", {
  if (isTRUE(require("MASS"))) {
    rlm1 <- rlm(stack.loss ~ ., stackloss)
    rlm2 <- rlm(stack.loss ~ ., stackloss, psi = psi.hampel, init = "lts")
    rlm3 <- rlm(stack.loss ~ ., stackloss, psi = psi.bisquare)

    expect_equal(sum(llcont(rlm1)), as.numeric(logLik(rlm1)))
    expect_equal(sum(llcont(rlm2)), as.numeric(logLik(rlm2)))
    expect_equal(sum(llcont(rlm3)), as.numeric(logLik(rlm3)))
  }
})


# test_that("SingleGroupClass", {
#   if (isTRUE(require("mirt"))) {
#     data <- expand.table(LSAT7)
#     mirts1 <- mirt(data, 1)
#     mirts2 <- mirt(data, 1, SE = TRUE)
#     mirts3 <- mirt(data, 1, SE = TRUE, SE.type = 'SEM')
#     mirts4 <- mirt(data, 1, SE = TRUE, SE.type = 'BL')
#     mirts5 <- mirt(data, 1, itemtype = c('2PL', '2PL', '2PL', '2PL', '3PL'))
#     model <- 'F = 1-5
#               PRIOR = (5, g, norm, -1.5, 3)'
#     mirts6 <- mirt(data, model, itemtype = c('2PL', '2PL', '2PL', '2PL', '3PL'))
#     mirts7 <- mirt(data, 2)
#     
#     expect_equal(sum(llcont(mirts1)), as.numeric(mirts1@logLik))
#     expect_equal(sum(llcont(mirts2)), as.numeric(mirts2@logLik))
#     expect_equal(sum(llcont(mirts3)), as.numeric(mirts3@logLik))
#     expect_equal(sum(llcont(mirts4)), as.numeric(mirts4@logLik))
#     expect_equal(sum(llcont(mirts5)), as.numeric(mirts5@logLik))
#     expect_equal(sum(llcont(mirts6)), as.numeric(mirts6@logLik))
#     expect_equal(sum(llcont(mirts7)), as.numeric(mirts7@logLik))
#   }
# })
