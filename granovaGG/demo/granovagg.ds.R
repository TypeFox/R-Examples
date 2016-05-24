### Using granovagg.ds to examine trends or effects for repeated measures data.

# This example corresponds to case 1b in Pruzek and Helmreich (2009). In this
# graphic we're looking for the effect of Family Treatment on patients with anorexia.

data(anorexia.sub)

granovagg.ds(anorexia.sub,
             revc = TRUE,
             main = "Assessment Plot for weights to assess \
                    Family Therapy treatment for Anorexia Patients",
             xlab = "Weight after therapy (lbs.)",
             ylab = "Weight before therapy (lbs.)"
)

### Using granovagg.ds to compare two experimental treatments (with blocking)

# This example corresponds to case 2a in Pruzek and Helmreich (2009). For this
# data, we're comparing the effects of two different virus preparations on the
# number of lesions produced on a tobacco leaf.

data(tobacco)
granovagg.ds(tobacco[, c("prep1", "prep2")],
             main = "Local Lesions on Tobacco Leaves",
             xlab = "Virus Preparation 1",
             ylab = "Virus Preparation 2"
)

### Using granovagg.ds to compare two experimental treatments (with blocking)

# This example corresponds to case 2a in Pruzek and Helmreich (2009). For this
# data, we're comparing the wear resistance of two different shoe sole
# materials, each randomly assigned to the feet of 10 boys.

data(shoes)
granovagg.ds(shoes,
             revc = TRUE,
             main = "Shoe Wear",
             xlab = "Sole Material B",
             ylab = "Sole Material A",
)

### Using granovagg.ds to compare matched individuals for two treatments

# This example corresponds to case 2b in Pruzek and Helmreich (2009). For this
# data, we're examining the level of lead (in mg/dl) present in the blood of
# children. Children of parents who had worked in a factory where lead was used
# in making batteries were matched by age, exposure to traffic, and neighborhood
# with children whose parents did not work in lead-related industries.

data(blood_lead)
granovagg.ds(blood_lead,
             sw = .1,
             main = "Dependent Sample Assessment Plot
             Blood Lead Levels of Matched Pairs of Children",
             xlab = "Exposed (mg/dl)",
             ylab = "Control (mg/dl)"
)
