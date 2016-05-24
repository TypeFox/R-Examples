pal <- palettes$palettes; dim(pal) <- c(5,4); pal
palperc <- 100 * row.perc(pal); palperc
palettes$palperc <- as.vector(palperc)
anova(lm(palperc~employee,palettes))
confint(glht(lm(palperc~employee,palettes),mcp(employee='Tukey')))
anova(lm(palperc~employee+day,palettes))
confint(glht(lm(palperc~employee+day,palettes),mcp(employee='Tukey')))
