#------ Funtion 1: Generate random numbers for experiments ------

gen.numbers = function(RNG, m, len){
  if (RNG==1){
    RNGkind("Wichmann-Hill")
    set.seed(283158301, kind = "Wichmann-Hill")
  }else if (RNG==2){
    RNGkind("Marsaglia-Multicarry")
    set.seed(283158301, kind = "Marsaglia-Multicarry")
  }else if (RNG==3){
    RNGkind("Super-Duper")
    set.seed(283158301, kind = "Super-Duper")
  }else if (RNG==4){
    RNGkind("Mersenne-Twister")
    set.seed(283158301, kind = "Mersenne-Twister")
  }else if (RNG==5){
    RNGkind("Knuth-TAOCP-2002")
    set.seed(283158301, kind = "Knuth-TAOCP-2002")
  }else if (RNG==6){
    RNGkind("Knuth-TAOCP")
    set.seed(283158301, kind = "Knuth-TAOCP")
  }else if (RNG==7){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(283158301, kind = "L'Ecuyer-CMRG")
  }
  numbers = array(NA, dim = len)
  mult = seq(0, m-1, 1)
  for (i in 1:len){
    say=round(runif(m,0,1))
    while (sum(say)==0){
      say=round(runif(m,0,1))
    }
    numbers[i] = round(runif(m,0,1))%*%(2^mult)
  }
  dat = digitsBase(numbers, base = 2, m)
  dat = as.vector(dat)
  bits = matrix(data = dat[1:(floor(length(dat)/m)*m)], nrow = m,
                ncol = floor(length(dat)/m), byrow = FALSE)
  result = list(numbers = numbers, bits = bits)
  return(result)
}

#----------------------------------------------------------------

#------------- Function 2: Carry on experiments -----------------

experiments = function(RNG, m, len, alpha, cv.TBT, mu.GCD, sd.GCD){
  result = array(NA, dim = 21)
  sim.data = gen.numbers(RNG, m, len)$numbers
  sim.data.bit = gen.numbers(RNG, m, len)$bits

  if (m<=64){
    sim.data = matrix(data = sim.data, ncol = len, byrow = FALSE)
  }else{
    sim.data = mpfrArray(sim.data, prec = m)
  }

  result[1] = topological.binary(x = sim.data.bit, B = m, critical.value = cv.TBT)$result.TBT

  result[2] = adaptive.chi.square(A = sim.data, B = m, S = 4, alpha)$result.acsq
  if (m<=32){
    mm = 8
    n = 2^m
    BDS = birthday.spacings(x = sim.data, m = mm, n = n, alpha=alpha, lambda = ((mm^3)/(4*n)))
    result[3] = BDS$AD.result
    result[4] = BDS$KS.result
    result[5] = BDS$CS.result
  }

  exc = FALSE
  expn = FALSE
  hei = FALSE
  do.test.g = TRUE
  if (m>=16){
    exc = TRUE
  }
  if (m>=32){
    hei = TRUE
  }
  if (m>=64){
    expn = TRUE
    do.test.g = FALSE
  }
  RWT = random.walk.tests(x = sim.data.bit, B = m, Excursion = exc,
                          Expansion = expn, Height = hei, alpha = alpha)
  result[6] = RWT$AD.result.Excursion
  result[7] = RWT$AD.result.Expansion
  result[8] = RWT$AD.result.Height
  result[9] = RWT$KS.result.Excursion
  result[10] = RWT$KS.result.Expansion
  result[11] = RWT$KS.result.Height
  result[12] = RWT$CS.result.Excursion
  result[13] = RWT$CS.result.Expansion
  result[14] = RWT$CS.result.Height

  if (m<=16){
    n = m*(2^(m/2)) #optimal length of sample that is composed of B-bit words
                    #given by Ryabko and Monarev (2005)
    dat.BS = sim.data[1:round(n/m)]
    BS = book.stack(A = dat.BS, B = m, k = n/m, alpha = alpha, bit = FALSE)
    result[15] = BS$BS.result
  }

  if (len%%2==1){
    len = len-1
  }
  len2 = len/2
  if (m<=64){
    dat.new = array(NA, dim = c(len2,2))
    dat.new[1:len2,1] = sim.data[1:len2]
    dat.new[1:len2,2] = sim.data[(len2+1):len]
  }else {
    dat.new = mpfrArray(NA, prec = m, dim = c(len2,2))
    dat.new[1:len2,1] = sim.data[1:len2]
    dat.new[1:len2,2] = sim.data[(len2+1):len]
  }
  EBOB=GCD.test(x = dat.new, B = m, mu = mu.GCD, sd = sd.GCD, test.g = do.test.g)
  result[16] = EBOB$KS.result.k
  result[17] = EBOB$CSQ.result.k
  result[18] = EBOB$AD.result.k
  result[19] = EBOB$JB.result.k
  if (m<=64){
    result[20] = EBOB$KS.result.g
    result[21] = EBOB$CSQ.result.g
  }
  return(result)
}
#----------------------------------------------------------------