require(lmerTest)

m <- lmer(Coloursaturation ~ TVset*Picture+
            (1|Assessor)+(1|Assessor:TVset), data=TVbo)
step(m)

#does not with lmerTest attached
m.lm <- lm(Coloursaturation ~ TVset*Picture, data=TVbo)
tools::assertError(step(m.lm))

#stats package needs to implicitely specified then
stats::step(m.lm)


tools::assertError(lsmeans(m.lm))

