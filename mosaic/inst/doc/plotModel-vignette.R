## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(mosaic)

## ------------------------------------------------------------------------
mod <- lm( mpg ~ factor(cyl), data = mtcars)
plotModel(mod)

# SLR
mod <- lm( mpg ~ wt, data = mtcars)
plotModel(mod, pch = 19)

# parallel slopes
mod <- lm( mpg ~ wt + factor(cyl), data=mtcars)
plotModel(mod)

# multiple categorical vars
mod <- lm( mpg ~ wt + factor(cyl) + factor(vs) + factor(am), data = mtcars)
plotModel(mod)
plotModel(mod, mpg = ~am)

# interaction
mod <- lm( mpg ~ wt + factor(cyl) + wt:factor(cyl), data = mtcars)
plotModel(mod)

# polynomial terms
mod <- lm( mpg ~ wt + I(wt^2), data = mtcars)
plotModel(mod)

# polynomial terms and covariate
mod <- lm( mpg ~ wt + I(wt^2) * factor(cyl), data = mtcars)
plotModel(mod)
plotModel(mod, formula = mpg ~ wt | cyl)

# GLM
mod <- glm(vs ~ wt, data=mtcars, family = 'binomial')
plotModel(mod)

# GLM with interaction
mod <- glm(vs ~ wt + factor(cyl), data=mtcars, family = 'binomial')
plotModel(mod, jitter.y = TRUE)

# GLM adjusted for factor response

mod <- glm( sex ~ width, data = KidsFeet, family = binomial())
plotModel(mod, jitter.y = TRUE)

# facets
mod <- lm(mpg ~ poly(wt,2) + factor(cyl) * factor(am), data=mtcars)
plotModel(mod, mpg ~ wt | factor(cyl))

## ------------------------------------------------------------------------
knitr::opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
#  mod <- lm( mpg ~ factor(cyl), data = mtcars)
#  plotModel(mod, system="g")
#  
#  # SLR
#  mod <- lm( mpg ~ wt, data = mtcars)
#  plotModel(mod, system="g")
#  
#  # parallel slopes
#  mod <- lm( mpg ~ wt + factor(cyl), data=mtcars)
#  plotModel(mod, system = "g")
#  
#  # multiple categorical vars
#  mod <- lm( mpg ~ wt + factor(cyl) + factor(vs) + factor(am), data = mtcars)
#  plotModel(mod, system="g")
#  plotModel(mod, mpg ~ am, system="g")
#  
#  # interaction
#  mod <- lm( mpg ~ wt + factor(cyl) + wt:factor(cyl) + factor(am), data = mtcars)
#  plotModel(mod, system="g")
#  plotModel(mod, system="g") + facet_grid( . ~ cyl)
#  
#  # polynomial terms
#  mod <- lm( mpg ~ wt + I(wt^2), data = mtcars)
#  plotModel(mod, system="g")
#  
#  # GLM
#  mod <- glm(vs ~ wt, data=mtcars, family = 'binomial')
#  plotModel(mod, system="g")
#  
#  # GLM with interaction
#  mod <- glm(vs ~ wt + factor(cyl), data=mtcars, family = 'binomial')
#  plotModel(mod, system="g")
#  plotModel(mod, system="g") + facet_grid(. ~ cyl)

