
lymph3.glm <- glm(nodes ~ X.ray + stage + grade + acid.ph,
                  data=lymph, family=binomial)
anova(lymph3.glm, test="Chisq")
summary(lymph3.glm)
