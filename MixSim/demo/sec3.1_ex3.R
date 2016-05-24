### Example 3 of Section 3.1
set.seed(1234)
(ex.3 <- MixSim(BarOmega = 0.05, K = 2, p = 4, sph = TRUE, hom = TRUE,
                int = c(0, 10), eps = 1e-10))
summary(ex.3)
