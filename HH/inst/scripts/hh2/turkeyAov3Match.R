
                    summary.lm(turkey3.aov)
               coef(summary.lm(turkey3.aov))
      dimnames(coef(summary.lm(turkey3.aov)))
      dimnames(coef(summary.lm(turkey3.aov)))[[1]]

               turkey3.aov$x
      dimnames(turkey3.aov$x)
      dimnames(turkey3.aov$x)[[2]]

match(dimnames(coef(summary.lm(turkey3.aov)))[[1]],
      dimnames(turkey3.aov$x)[[2]])
