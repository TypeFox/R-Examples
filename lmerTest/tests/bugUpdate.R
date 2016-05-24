require(lmerTest)
#### ASSERTS ERRORS
m <- lmer(Coloursaturation ~ TVset*Picture +
            (1|Assessor)+(1|Assessor:TVset), data=TVbo)
m2 <- update(m)
tools::assertError(stopifnot(class(m)==class(m2), TRUE), TRUE)




