require(lmerTest)

## for carrots data
modelCarrots <- lmer(Preference ~  sens2*sens1*Homesize*Age
                      + (1 + sens1 | Consumer),
                     data=carrots)

st1 <- step(modelCarrots, keep.effs = c("sens1","sens2"))
st1

st2 <- step(modelCarrots, keep.effs = c("sens1", "Homesize:Age"))
st2

step(modelCarrots, keep.effs = c("sens1:Consumer"))

step(modelCarrots, keep.effs = c("dfgsdfg"))



## for TV data

m <- lmer(Dimglasseffect ~ TVset*Picture+
            (1|Assessor)+(1|Assessor:TVset), data=TVbo)

step(m, keep.effs = "TVset:Picture")
