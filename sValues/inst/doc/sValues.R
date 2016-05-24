## ----include=FALSE, echo=FALSE-------------------------------------------
library(knitr)
render_sweave()
# hook2 <- function(x){ gsub("```\n+```\n", "", x) }
# knit_hooks$set(document = hook2)
opts_chunk$set(size="small", comment = "")

## ----size='small'--------------------------------------------------------
library(sValues) # loads package
data("economic_growth") # loads data
eg <- sValues(economic_growth) # runs analysis
eg # prints basic results

## ------------------------------------------------------------------------
# gets a complete table like in Leamer[3] 
full_table <- coef(eg) 
full_table[1:5, 1:5] # showing only first five columns and rows

# gets just the s_values 
just_svalues <- coef(eg, type = "s_values")
just_svalues[1:5, ] # showing only first five rows

## ------------------------------------------------------------------------
extreme_bounds(eg)[c("GOVNOM1", "IPRICE1"),
                   c("R2_0.5_1.low", "R2_0.5_1.up")]

## ----t_s_plot, size='small', fig.align='center', fig.lp="fig:", fig.cap = 't-statistics vs s-values'----
plot(eg, type = "t_s_plot", R2_bounds = c(0.5, 1))

## ----beta_plot, size='small', fig.align='center', fig.lp="fig:", fig.cap = 'Bayesian estimates for GOVNOM1, with error bars and extreme bounds (shaded areas).'----
plot(eg, type = "beta_plot", variables = "GOVNOM1",
     error_bar = TRUE, ext_bounds_shades = TRUE)

## ----size = 'small'------------------------------------------------------
favorites <- c("GDPCH60L", "OTHFRAC", "ABSLATIT", 
               "LT100CR", "BRIT", "GOVNOM1", 
               "WARTIME", "SCOUT","P60", "PRIEXP70", 
               "OIL", "H60", "POP1560", "POP6560")
eg_fav <- sValues(economic_growth, R2_bounds = c(0.5, 1),
                  favorites = favorites, R2_favorites = c(0.4, 0.8))
eg_fav

