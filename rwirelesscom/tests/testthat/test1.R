


context("BPSK error rate, Eb/No = 4 dB")
test_that("Test BPSK Eb/No = 4 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=2
  Eb=1
  Es = log2(M)*Eb
  Nsymbols=10000
  Ns=1
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- fbpskmod(bits)

  EbNodB=4
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No)
  r <- s+n
  bitsr <- fbpskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.03, info="BPSK EbNodb=4, Pberr should be < 0.03")
  expect_true(Pberr > 0.008, info="PBSK EbNodb=4, Pberr should be > 0.008")

} )

context("BPSK error rate, Eb/No = 8 dB")
test_that("Test BPSK EbNo_dB= 4, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=2
  Eb=1
  Es = log2(M)*Eb
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- sqrt(Eb)*fbpskmod(bits)

  EbNodB=8
    No = Eb/(10^(EbNodB/10))
    n <- fNo(Nsymbols,No)
    r <- s+n
    bitsr <- fbpskdemod(r)
    biterrs<-bits[bitsr!=bits]
    b<-factor(bits)
    Pberr=length(biterrs)/length(bits)

    # str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
    #print("",quote=FALSE)
    #print(str,quote=FALSE)
    expect_true(Pberr < 0.00032, info="BPSK EbNodb=8, Pberr should be < 0.0032")
    expect_true(Pberr > 0.0001, info="PBSK EbNodb=8, Pberr should be > 0.001")

} )

context("QPSK error rate, Eb/No = 4 dB")
test_that("Test BPSK Eb/No = 4 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=4
  Es=1
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- sqrt(Es)*fqpskmod(bits)

  EbNodB=4
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- fqpskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.02, info="BPSK EbNodb=4, Pberr should be < 0.015")
  expect_true(Pberr > 0.007, info="PBSK EbNodb=4, Pberr should be > 0.012")

} )

context("QPSK error rate, Eb/No = 8 dB")
test_that("Test BPSK EbNo_dB= 4, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=4
  Es=1
  Eb = Es/log2(M)
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)

  s <- sqrt(Es)*fqpskmod(bits)

  EbNodB=8
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- fqpskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: BPSK EbNo_dB = %d, Bits = %g, bit errors = %g, Pberr=%f",EbNodB, length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)
  expect_true(Pberr < 0.00032, info="BPSK EbNodb=8, Pberr should be < 0.0032")
  expect_true(Pberr > 0.0001, info="PBSK EbNodb=8, Pberr should be > 0.001")

} )


context("8-PSK error rate, Eb/No = 7 dB")
test_that("Test 8-PSK Eb/No = 11 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=8
  Es=1
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- f8pskmod(bits)

  EbNodB=7
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f8pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  # str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  # print("",quote=FALSE)
  # print(str,quote=FALSE)

  expect_true(Pberr < 0.015, info="8-PSK EbNodb=7, Pberr should be < 0.015")
  expect_true(Pberr > 0.01, info="8-BSK EbNodb=7, Pberr should be > 0.01")

} )

context("8-PSK error rate, Eb/No = 10 dB")
test_that("Test 8-PSK Eb/No = 11 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=8
  Es=1
  Eb = Es/log2(M)
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- f8pskmod(bits)

  EbNodB=10
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f8pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.002, info="8-PSK EbNodb=7, Pberr should be < 0.015")
  expect_true(Pberr > 0.0005, info="8-BSK EbNodb=7, Pberr should be > 0.01")

} )

context("16-PSK error rate, Eb/No = 12 dB")
test_that("Test 16-PSK Eb/No = 12 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=16
  Es=1
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)
  s <- f16pskmod(bits)

  EbNodB=12
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f16pskdemod(r)
  biterrs<-bits[bitsr!=bits]
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.01, info="16-PSK EbNodb=12, Pberr should be < 0.01")
  expect_true(Pberr > 0.006, info="16-BSK EbNodb=12, Pberr should be > 0.006")

} )

context("16-QAM error rate, Eb/No = 8 dB")
test_that("Test 16-QAM Eb/No = 8 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=16
  Es=10
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)

  s <- f16qammod(bits)
  EbNodB=8
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f16qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

   #str<-sprintf("Test: %d-QAM EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
   #print("",quote=FALSE)
   #print(str,quote=FALSE)

  expect_true(Pberr < 0.012, info="16-QAM EbNodb=8, Pberr should be < 0.012")
  expect_true(Pberr > 0.004, info="16-QAM EbNodb=8, Pberr should be > 0.004")

} )

context("16-QAM error rate, Eb/No = 10 dB")
test_that("Test 16-QAM Eb/No = 10 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=16
  Es=10
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)

  s <- f16qammod(bits)
  EbNodB=10
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f16qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-QAM EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.003, info="16-QAM EbNodb=10, Pberr should be < 0.005")
  expect_true(Pberr > 0.001, info="16-QAM EbNodb=10, Pberr should be > 0.002")

} )

context("64-QAM error rate, Es/No = 12 dB")
test_that("Test 64-QAM Eb/No = 12 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  M=64
  Es=42
  Eb = Es/log2(M)
  Nsymbols=10000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)


  #  Nbits=24
  #  Nsymbols=24/log2(M)
  #  bits=c(0,1,0,1,0,0,0,0,1,0,1,1, 1,0,0,0,1,1, 0,1,0,1,0,1)

  s <- f64qammod(bits)

  EbNodB=12
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f64qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  # str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  # print("",quote=FALSE)
  # print(str,quote=FALSE)

  expect_true(Pberr < 0.012, info="8-PSK EbNodb=12, Pberr should be < 0.012")
  expect_true(Pberr > 0.008, info="8-BSK EbNodb=12, Pberr should be > 0.008")

} )

context("64-QAM error rate, Es/No = 14 dB")
test_that("Test 64-QAM Eb/No = 14 dB, Modulator and Demodulator in AWGN, Bit Error Rate", {
  skip_on_cran()
  M=64
  Es=42
  Eb = Es/log2(M)
  Nsymbols=100000
  Nbits=log2(M)*Nsymbols
  bits <- sample(0:1,Nbits, replace=TRUE)

  #  Nbits=24
  #  Nsymbols=24/log2(M)
  #  bits=c(0,1,0,1,0,0,0,0,1,0,1,1, 1,0,0,0,1,1, 0,1,0,1,0,1)

  s <- f64qammod(bits)

  EbNodB=14
  No = Eb/(10^(EbNodB/10))
  n <- fNo(Nsymbols,No,type="complex")
  r <- s+n
  bitsr <- f64qamdemod(r)
  biterrs<-bits[bitsr!=bits]
  b<-factor(bits)
  Pberr=length(biterrs)/length(bits)

  #str<-sprintf("Test: %d-PSK EbNo_dB = %d, EsNo_dB = %g, Bits = %g, bit errors = %g, Pberr=%f",M,EbNodB, 10*log10(Es/No), length(bits), length(biterrs),Pberr)
  #print("",quote=FALSE)
  #print(str,quote=FALSE)

  expect_true(Pberr < 0.003, info="8-PSK EbNodb=14, Pberr should be < 0.003")
  expect_true(Pberr > 0.001, info="8-BSK EbNodb=14, Pberr should be > 0.001")

} )

context("Modulation w/ pulse shape")
test_that("Modulation pulse shaping", {
  skip_on_cran()

  Nsymbols=1
  M=2
  Ns=8
  hx=seq(-4*Ns,4*Ns,by=1)
  h1=rcosine(hx,B=1,Ns=Ns)
  h0_25=rcosine(hx,B=0.25,Ns=Ns)
  h0_5=rcosine(hx,B=0.5,Ns=Ns)
  bpsk_bits <- c(0,1)
  s_bpsk <- fbpskmod(bpsk_bits,Ns,h0_5)
  x_bpsk= c(-1.776357e-16,-2.199782e-03,-4.245149e-03,-5.556007e-03,-5.716294e-03,-4.658960e-03,-2.771533e-03,-8.675462e-04,0.000000e+00,1.062116e-03,-5.276463e-04,-5.010069e-03,-1.143259e-02,-1.752667e-02,-1.997811e-02,-1.512443e-02,8.073768e-17,2.646658e-02,6.238553e-02,1.022829e-01,1.371911e-01,1.556013e-01,1.452514e-01,9.551705e-02,-1.665335e-16,-1.412083e-01,-3.201165e-01,-5.206996e-01,-7.202531e-01,-8.921642e-01,-1.009738e+00,-1.050467e+00,-1.000000e+00,-8.550623e-01,-6.247323e-01,-3.297658e-01,1.712515e-16,3.297658e-01,6.247323e-01,8.550623e-01,1.000000e+00,1.050467e+00,1.009738e+00,8.921642e-01,7.202531e-01,5.206996e-01,3.201165e-01,1.412083e-01,2.220446e-16,-9.551705e-02,-1.452514e-01,-1.556013e-01,-1.371911e-01,-1.022829e-01,-6.238553e-02,-2.646658e-02,2.638617e-17,1.512443e-02,1.997811e-02,1.752667e-02,1.143259e-02,5.010069e-03,5.276463e-04,-1.062116e-03,-1.887379e-16,8.675462e-04,2.771533e-03,4.658960e-03,5.716294e-03,5.556007e-03,4.245149e-03,2.199782e-03,0.000000e+00,0.000000e+00,0.000000e+00,-6.661338e-17,-1.192185e-17,-8.881784e-17,0.000000e+00,-1.776357e-16)
  bpsk_err=max(abs(s_bpsk-x_bpsk))

  M=4
  Ns=4
  hx=seq(-2*Ns,2*Ns,by=1)
  h1=rcosine(hx,B=1,Ns=Ns)
  s_qpsk <- fqpskmod(c(0,1),Ns,h1)
  x_qpsk=c(2.220446e-17+2.220446e-17i,-5.716294e-03+5.716294e-03i,-4.636454e-17+1.059389e-16i,1.714888e-02-1.714888e-02i,-1.056949e-16+6.345947e-17i,-1.200422e-01+1.200422e-01i,-3.535535e-01+3.535535e-01i,-6.002109e-01+6.002109e-01i,-7.071068e-01+7.071068e-01i,-6.002109e-01+6.002109e-01i,-3.535535e-01+3.535535e-01i,-1.200422e-01+1.200422e-01i,-1.312711e-16+7.169683e-17i,1.714888e-02-1.714888e-02i,-4.973636e-17+9.197175e-17i,-5.716294e-03+5.716294e-03i,0.000000e+00-4.440892e-17i,-1.332268e-16+8.881784e-17i,0.000000e+00+0.000000e+00i,0.000000e+00+8.881784e-17i)
  qpsk_err=max(abs(s_qpsk-x_qpsk))

  s_8psk <- f8pskmod(c(0,1,0),3)
  x_8psk <- c(complex(real=-1/sqrt(2),imaginary = 1/sqrt(2)),0,0)
  err_8psk = Re(sum( (s_8psk - x_8psk) * Conj(s_8psk - x_8psk)))

  s_16psk <- f16pskmod(c(0,0,1,1),4)
  x_16psk <- c(complex(real=1/sqrt(2),imaginary = 1/sqrt(2)),0,0,0)
  err_16psk <- Re(sum( (s_16psk - x_16psk) * Conj(s_16psk - x_16psk)))

  s_16qam <- f16qammod(c(0,0,1,1),5)
  x_16qam <- c(complex(real=1,imaginary = -3),0,0,0,0)
  err_16qam <- Re(sum( (s_16qam - x_16qam)*Conj(s_16qam - x_16qam) ))

  s_64qam <- f64qammod(c(0,0,0,1,1,1),5)
  x_64qam <- c(complex(real=3,imaginary = -7),0,0,0,0)
  err_64qam <- Re(sum( (s_64qam - x_64qam)*Conj(s_64qam - x_64qam)))

  expect_true(bpsk_err  < 1e-6, info="BPSK w/ rcosine unit test error should be < 1e-6")
  expect_true(qpsk_err < 1e-7, info="QPSK w/ rcosine unit test error should be < 1e-7")
  expect_true(err_8psk < 1e-10, info="8-PSK test error")
  expect_true(err_16psk < 1e-10, info="16-PSK test error")
  expect_true(err_16qam < 1e-10, info="16-QAM test error")
  expect_true(err_64qam < 1e-10, info="64-QAM test error")

} )

context("rcosine and sqrtrcosine test")
test_that("unit test basic rcosine and sqrtrcosine", {
  skip_on_cran()
  Ns=8
  B=0.5
  hx=seq(-4*Ns,4*Ns,by=1)
  h05=rcosine(hx,B,Ns)
  x05= c(0.0000000000,0.0021997815,0.0042451487,0.0055560072,0.0057162941,0.0046589601,0.0027715325,0.0008675462,0.0000000000,0.0011376658,0.0047727950,0.0105660767,0.0171488822,0.0221856300,0.0227496429,0.0159919804,0.0000000000,-0.0253289130,-0.0576127319,-0.0917168572,-0.1200421755,-0.1334156763,-0.1225017380,-0.0795250670,0.0000000000,0.1158793833,0.2625037243,0.4289827129,0.6002108774,0.7587485457,0.8872360718,0.9709416669,1.0000000000,0.9709416669,0.8872360718,0.7587485457,0.6002108774,0.4289827129,0.2625037243,0.1158793833,0.0000000000,-0.0795250670,-0.1225017380,-0.1334156763,-0.1200421755,-0.0917168572,-0.0576127319,-0.0253289130,0.0000000000,0.0159919804,0.0227496429,0.0221856300,0.0171488822,0.0105660767,0.0047727950,0.0011376658,0.0000000000,0.0008675462,0.0027715325,0.0046589601,0.0057162941,0.0055560072,0.0042451487,0.0021997815,0.0000000000)
  h05_err=max(abs(h05-x05))
  h05one= h05[c(seq(1,65,by=8))]
  h05one[5]

  Ns=4
  B=0.25
  hx=seq(-2*Ns,2*Ns,by=1)
  h025=sqrtrcosine(hx,B,Ns)
  x025=c(0.05305165,-0.05499532,-0.17029763,-0.19871737,-0.06423716,0.23786183,0.62179741,0.94316532,1.06830989,0.94316532,0.62179741,0.23786183,-0.06423716,-0.19871737,-0.17029763,-0.05499532,0.05305165)
  h025_err=max(abs(h025-x025))

  expect_true(sum(h05one) == 1, info="rcosine zeros and one not = 1")
  expect_true(h05one[5] ==1, info="rcosine at x[n]=0 not = 1")
  expect_true(h05_err < 1e-10, info="rcosine error")
  expect_true(h025_err < 1e-8, info="sqrt rcsoine error")



} )
