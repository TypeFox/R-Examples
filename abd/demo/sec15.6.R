data(WalkingStickFemurs)
str(WalkingStickFemurs)

WalkingStickFemurs$specimenF <- factor(WalkingStickFemurs$specimen)
fm <- aov(femur.length ~ specimenF, data = WalkingStickFemurs)
summary(fm)

if (require(ICC)){
  ICCest(specimenF,
         femur.length,
         data = WalkingStickFemurs,
         alpha = 0.05,
         CI.type = "Smith")
}
