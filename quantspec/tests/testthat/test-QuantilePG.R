context("QuantilePG")

test_that("quantilePG works as expected for various levels",{
      
      source("load-ref.R")
      
      lev.ok.all <- c(0.25,0.5,0.75)
      lev.ok.1   <- c(0.25)
      lev.ok.2   <- c(0.5,0.75)
      lev.err.1  <- c("non numeric",0.5)
      lev.err.2  <- c(0.5,1.5)

      # Check whether it works for all levels:
      V.qr.ref.1 <- array(V.qr.ref[,,,1], dim=c(64,3,3,1))
      V.fft.ref.1 <- array(V.fft.ref[,,,1], dim=c(64,3,3,1))
      
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.all, type="qr")
      V.qr <- getValues(qPG.qr)
      expect_that(dim(V.qr),equals(c(64,3,3,1)))
      expect_that(V.qr,equals(V.qr.ref.1))
      
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.all, type="clipped")
      V.fft <- getValues(qPG.fft)
      expect_that(dim(V.fft),equals(c(64,3,3,1)))
      expect_that(V.fft,equals(V.fft.ref.1))
      
      # Now, check whether it works for only one level:
      V.qr.ref.1 <- array(V.qr.ref[,1,1,1], dim=c(64,1,1,1))
      V.fft.ref.1 <- array(V.fft.ref[,1,1,1], dim=c(64,1,1,1))
      
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.1, type="qr")
      V.qr.1 <- getValues(qPG.qr)      
      expect_that(dim(V.qr.1),equals(c(64,1,1,1)))
      expect_that(V.qr.1,equals(V.qr.ref.1))
      
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.1, type="clipped")
      V.fft.1 <- getValues(qPG.fft)      
      expect_that(dim(V.fft.1),equals(c(64,1,1,1)))
      expect_that(V.fft.1,equals(V.fft.ref.1))
      
      # Now, check whether it works for two different sized levels:
      V.qr.ref.1 <- array(V.qr.ref[,1,2:3,1], dim=c(64,1,2,1))
      V.fft.ref.1 <- array(V.fft.ref[,1,2:3,1], dim=c(64,1,2,1))
      
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.1, levels.2=lev.ok.2, type="qr")
      V.qr <- getValues(qPG.qr)      
      expect_that(dim(V.qr),equals(c(64,1,2,1)))
      expect_that(V.qr,equals(V.qr.ref.1))
      
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.1, levels.2=lev.ok.2, type="clipped")
      V.fft <- getValues(qPG.fft)      
      expect_that(dim(V.fft),equals(c(64,1,2,1)))
      expect_that(V.fft,equals(V.fft.ref.1))
      
      # Now, check whether it fails for incorrect input:
      expect_that(quantilePG(Y,levels.1=lev.err.1,type="qr"),throws_error())
      expect_that(quantilePG(Y,levels.1=lev.err.2,type="clipped"),throws_error())
      
    })

test_that("quantilePG works as expected for various frequencies",{
      
      source("load-ref.R")
      
      lev.ok.all <- c(0.25,0.5,0.75)
      
      # Test whether computation on all Fourier frequencies works
      freq.init.all  <- 2*pi*(0:63)/64
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.all,
          frequencies=freq.init.all, type="qr")
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.all,
          frequencies=freq.init.all, type="clipped")
      
      # Call some frequencies that are just a little bit off
      V.qr.ref.1 <- array(V.qr.ref[,,,1], dim=c(64,3,3,1))
      V.fft.ref.1 <- array(V.fft.ref[,,,1], dim=c(64,3,3,1))
      
      freq.call.part <- 2*pi*((0:63)/64+1/256)
      expect_that(V.qr <- getValues(qPG.qr, frequencies=freq.call.part),
          gives_warning())
      expect_that(dim(V.qr),equals(c(64,3,3,1)))
      expect_that(V.qr, equals(V.qr.ref.1))
      
      expect_that(V.fft <- getValues(qPG.fft, frequencies=freq.call.part),
          gives_warning())
      expect_that(dim(V.fft),equals(c(64,3,3,1)))
      expect_that(V.fft, equals(V.fft.ref.1))
      
      # Now getValues for every second of those frequencies  
      # - should not give a warning and yield the correct numbers!!
      freq.call.part <- 2*pi*(0:31)/32
      
      V.qr <- getValues(qPG.qr, frequencies=freq.call.part)
      expect_that(dim(V.qr),equals(c(32,3,3,1)))
      expect_that(V.qr[,,,1], equals(V.qr.ref[1+2*(0:31),,,1]))

      V.fft <- getValues(qPG.fft, frequencies=freq.call.part)
      expect_that(dim(V.fft),equals(c(32,3,3,1)))
      expect_that(V.fft[,,,1], equals(V.fft.ref[1+2*(0:31),,,1]))      

      
      # Now the other way around (init with every second, call all!)
      freq.init.part  <- 2*pi*(0:31)/32
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.all,
          frequencies=freq.init.part, type="qr")
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.all,
          frequencies=freq.init.part, type="clipped")
      
      # Now getValues for every second of those frequencies  
      # - gives warning, but yields the correct numbers!!
      freq.call.all <- 2*pi*(0:63)/64
      expect_that(V.qr <- getValues(qPG.qr, frequencies=freq.call.all),
          gives_warning())
      expect_that(dim(V.qr),equals(c(64,3,3,1)))
      expect_that(V.qr[1+2*(0:31),,,1], equals(V.qr.ref[1+2*(0:31),,,1]))

      expect_that(V.fft <- getValues(qPG.fft, frequencies=freq.call.all),
          gives_warning())
      expect_that(dim(V.fft),equals(c(64,3,3,1)))
      expect_that(V.fft[1+2*(0:31),,,1], equals(V.fft.ref[1+2*(0:31),,,1]))
      
      
      # Now check initializing only in the beginning 
      freq.init.beg <- 2*pi*(0:15)/64
      freq.call.all <- 2*pi*(0:63)/64       
      
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.all, frequencies=freq.init.beg, type="qr")
      expect_that(V.qr <- getValues(qPG.qr, frequencies = freq.call.all),
          gives_warning())
      expect_that(dim(V.qr),equals(c(64,3,3,1)))
      expect_that(V.qr[1:33,,,1], equals(V.qr.ref[c(1:16,rep(16,17)),,,1]))
      expect_that(V.qr[34:64,,,1], equals(Conj(V.qr.ref[c(rep(16,16),16:2),,,1])))
      
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.all, frequencies=freq.init.beg, type="clipped")
      expect_that(V.fft <- getValues(qPG.fft, frequencies = freq.call.all),
          gives_warning())
      expect_that(dim(V.fft),equals(c(64,3,3,1)))
      expect_that(V.fft[1:33,,,1], equals(V.fft.ref[c(1:16,rep(16,17)),,,1]))
      expect_that(V.fft[34:64,,,1], equals(Conj(V.fft.ref[c(rep(16,16),16:2),,,1])))
      
      # Now check frequencies not from [0,2pi) and in various orders!

      freq.init <- 2*pi*(0:63)/64  
      freq.call <- 2*pi*c(64,32,128)/64  

      qPG.qr <- quantilePG(Y, levels.1=lev.ok.all, frequencies=freq.init, type="qr")
      V.qr <- getValues(qPG.qr, frequencies = freq.call)
      expect_that(dim(V.qr),equals(c(3,3,3,1)))
      expect_that(V.qr[,,,1], equals(V.qr.ref[c(1,33,1),,,1]))
      
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.all, frequencies=freq.init, type="clipped")
      V.fft <- getValues(qPG.fft, frequencies = freq.call)
      expect_that(dim(V.fft),equals(c(3,3,3,1)))
      expect_that(V.fft[,,,1], equals(V.fft.ref[c(1,33,1),,,1]))
    })
    
    
test_that("quantilePG works as expected with bootstrapping",{
      
      source("load-ref.R")
      
      lev.ok.all <- c(0.25,0.5,0.75)
      set.seed(2581)
      
      # Check whether it works for all levels, with bootstrapping:
      qPG.qr <- quantilePG(Y, levels.1=lev.ok.all, type="qr", B=1, l=8, type.boot="mbb")
      V.qr <- getValues(qPG.qr)
      expect_that(dim(V.qr),equals(c(64,3,3,2)))
      expect_that(V.qr,equals(V.qr.ref))
      
      # Check whether it works for all levels, with bootstrapping:
      qPG.fft <- quantilePG(Y, levels.1=lev.ok.all, type="clipped", B=1, l=8, type.boot="mbb")
      V.qr <- getValues(qPG.qr)
      expect_that(dim(V.qr),equals(c(64,3,3,2)))
      expect_that(V.qr,equals(V.qr.ref))
      
    })