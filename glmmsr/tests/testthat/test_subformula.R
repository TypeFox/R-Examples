library(glmmsr)
context("Parsing subformula")


set.seed(1)
player <- 1:10
player1 <- factor(rep(2:10, 10), levels = player)
player2 <- factor(rep(1:9, 10), levels = player)

x <- rnorm(length(player))
u0 <- rnorm(length(player))

beta0 <- 1
ability0 <- beta0*x + u0
p0 <- pnorm(ability0[player1] - ability0[player2])
y <- rbinom(length(p0), 1, p0)

formula <- y ~ 0 + Sub(ability[player1] - ability[player2])
subform <- ability[player] ~ 0 + x[player] + (1 | player)
data <- list(x = x, player1 = player1, player2 = player2)

test_that("passes to glFormula if no Sub() terms", {
    form <- y ~ 0 + (1 | player1) + (1 | player2)
    data <- list(y = y, player1 = player1, player2 = player2)
    modfr1 <- find_modfr_glmm(form, data = data, family = binomial)
    modfr2 <- lme4::glFormula(form, data = data, family = binomial)
    expect_equal(modfr1$reTrms, modfr2$reTrms)
})

test_that("splits up formula correctly", {
  expect_equal(length(split_formula(formula)$subexprs), 1)
  form0 <- formula(y ~ 0 + x)
  form1 <- formula(y ~ 0 + x + Sub(ability[player1] - ability[player2]))
  form2 <- formula(y ~ 0 + x + Sub(ability[player1] - ability[player2])
                  + Sub(stuff[other1] - stuff[other2]))
  expect_equal(length(split_formula(form0)$subexprs), 0)
  expect_equal(length(split_formula(form1)$subexprs), 1)
  expect_equal(length(split_formula(form2)$subexprs), 2)
})

test_that("finds subvar from subform", {
  subform0 <- formula(ability[player] ~ x[player] + (1 | player))
  expect_equal(find_subvar(subform0), "ability")
  subform1 <- formula(ability[player] + ability[player2] ~ (1 | player))
  expect_error(find_subvar(subform1), "single variable on LHS")
})


test_that("matches subforms and subexprs correctly", {
  subforms <- list(formula(ability[player] ~ 0 + (1 | player)),
                   formula(stuff[i, j] ~ x[i] + (1 | j)))
  subexprs <- list(quote(ability[player1] - ability[player2]),
                   quote(stuff[team1, team2] - 3*stuff[team2, team1]))
  data = list(player1 = c(1, 2), player2 = c(2, 1),
              team1 = c(1, 2, 3), team2 = c(4, 2, 1), x = c(1, 2, 3, 4))
  out <- match_subform_subexpr(subforms, subexprs, data)
  expect_equal_to_reference(out, "subs_ability_stuff.rds")
  subforms1 <- list(formula(ability[player] ~ 0 + (1 | player)),
                  formula(stuff[i, j] ~ x[i] + (1 | j)),
                  other[k] ~ (1 | k))
  expect_warning(match_subform_subexpr(subforms1, subexprs, data),
                "No subexpressions involving")
}
)
