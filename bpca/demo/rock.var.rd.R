oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp <- bpca(rock,
           var.rb=TRUE,
           var.rd=TRUE)
summary(bp)

# The total variance explained is satisfactory (>= .80)!

plot(bp)

# A more accurate diagnostic
bp$var.rd

# It is possible to observe that the variable 'perm'
# has not a good representation (bpca.2d)!

# Observed correlations:
cor(rock)

# Projected correlations:
bp$var.rb

# Aditional diagnostic
plot(qbpca(rock,
           bp))

# This variable reamains as important in a dimension not contemplated
# by the biplot reduction (PC3):

bp$eigenvectors

bp1 <- bpca(rock,
            d=3:4)

summary(bp1)

plot(bp1)

# The recommendation, knowing that this variable has a poor
# representaion is:
# 1- Avoid to discut it;
# 2- Consider to incorporate the information with a bpca.3d

bp2 <- bpca(rock,
            d=1:3,
            var.rb=TRUE,
            var.rd=TRUE)

summary(bp2)

plot(bp2)           # Static

plot(bp2,
     rgl.use=TRUE)  # Dinamic

bp2$var.rd          # Nice!

# Aditional diagnostic
plot(qbpca(rock,
           bp2))

devAskNewPage(oask)

