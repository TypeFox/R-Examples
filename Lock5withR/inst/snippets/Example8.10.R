head(TextbookCosts)
Books.Model <- lm(Cost ~ Field, data = TextbookCosts)
anova(Books.Model)
summary(Books.Model)

