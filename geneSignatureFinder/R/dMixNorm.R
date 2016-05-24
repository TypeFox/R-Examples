dMixNorm <-
function(x, mmu1, ssd1, mmu2, ssd2, ttau)
  ttau * dnorm(x, mmu1, ssd1) + (1 - ttau) * dnorm(x, mmu2, ssd2)
