# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

normalize.compare = function (normd){
  for (j in 1:length(normd)){
    if (length(normd) > 1) tmp = normd[[j]]
    else tmp = normd
    between = 0
    within = 0
    withintot = 0
    total = 0
    formants = as.matrix (cbind(tmp$f1, tmp$f2))
    vowels = tmp$vowel
    vowelsf = levels(vowels)
    n = length (vowelsf)
    
    d = 0
    for (i in 1:(n-1))
      for (k in (i+1):n)     
        d = d + hotelling.test(formants[vowels == vowelsf[i],],
                               formants[vowels == vowelsf[k],])$f.value
    
    d = (d / (gamma (n+1) / (gamma (n-1) * 2)))^.5

    template = createtemplate (formants, tmp$vowel)
    winners = ldclassify (formants, vowels)
    correct = mean (winners == tmp$vowel)
    
    cat('\n   Method ', j, '\n\n')
    cat('   Average between-category distance: ', d, '\n\n')
    cat('   Correct Classification: ', correct*100, '%\n\n')
  }
  cat ('\n\n')
}
