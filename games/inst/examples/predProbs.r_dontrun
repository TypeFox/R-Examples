data("war1800")
f1 <- esc + war ~ s_wt_re1 + revis1 | 0 | regime1 | balanc + regime2
m1 <- egame12(f1, data = war1800, boot = 10)

pp1 <- predProbs(m1, x = "s_wt_re1", n = 5)
print(pp1)  ## Hypothetical observations and their predicted probs
plot(pp1, which = 2)  ## See ?plot.predProbs for more plot examples

## Changing the profile used
pp2 <- predProbs(m1, x = "s_wt_re1", n = 5, revis1 = 1, balanc = 0.7)
pp3 <- predProbs(m1, x = "s_wt_re1", n = 5, regime1 = "dem")
pp4 <- predProbs(m1, x = "s_wt_re1", n = 5, balanc = median(balanc))

## Variable names (other than `x`) must match exactly!
\dontrun{
    pp5 <- predProbs(m1, x = "s_wt_re1", bal = 0.7)  ## Error will result
}

## `x` can be a factor too
pp6 <- predProbs(m1, x = "regime1")

## Action probabilities
pp7 <- predProbs(m1, x = "regime1", type = "action")
