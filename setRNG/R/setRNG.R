  
###########################################################

#    Random number generation.

###########################################################

#  Prior to R 0.99 Wichmann-Hill was the default and the DSE version of
#  Box-Muller was used for rnorm.
  
setRNG <- function(kind=NULL, seed=NULL, normal.kind=NULL)
      {# with a null argument this also serves as getRNG 
       #The next line means that setRNG with null args does not truly
       #  return the state of the RNG in the case when it has not been 
       #  initialize. It first initializes it and then returns the state. The
       #  rational is that querying the state is usually for the purpose of
       #  reproducing it, so it must first be initialized to put it in a 
       #  reproducible state.
	if (!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE))
	     z <- runif(1)
	old <- list(kind=RNGkind()[1],normal.kind=RNGkind()[2], 
	     seed=get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)[-1])

        if (is.null(kind) & is.null(seed) & is.null(normal.kind)) return (old)
	if (is.list(kind)) 
          {seed        <- kind$seed
	   normal.kind <- kind$normal.kind
	   kind        <- kind$kind
	  }
       #Something like this might be nice, but there does not seem to be a way
       #  to associate rng versions other than using R versions, and 
       #  that is too strict. Perhaps RNGversion() could return something.
       #v <- some number associated with kind and normal.kind
       #if (is.null(v)) warning("version cannot be verified.")
       #if (!all(unlist(RNGversion()) == unlist(v)))
       # warning("rng used but version does not correspond to original. Differences may occur.")
	k <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)[1]
	remove(".Random.seed", envir=.GlobalEnv) # otherwise RNGkind complains
	RNGkind(kind=kind, normal.kind=normal.kind) # this sets .Random.seed
	seed <- as.integer(seed) #fix case where an integer is stored as double
	if ( 1==length(seed)) set.seed(seed) 
        else assign(".Random.seed", 
	    c(get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)[1], seed), 
	    envir=.GlobalEnv)
	invisible(old)
      }


getRNG <- function(e=NULL)UseMethod("getRNG")
getRNG.default <- function(e=NULL)
  {if (is.null(e))             k <- setRNG()
   else if (!is.null(e$rng))   k <- e$rng
   else if (!is.null(e$noise)) k <- e$noise$rng
   else stop("RNG information not found.")
   k
  }


#########################################################

#   test function (This test is run by other tests so it should
#   be defined here and not moved to tests directory.

#########################################################


random.number.test <- function()
 {cat("Random number generator tests ...")
  if (is.R())  
     {test.seed<- 979   #previous to R 1.0.0: c( 979, 1479, 1542) 
      # values from 0.49 beta
      #test.valueU <-c(5.693354055333957e-01,1.051357751852140e-01,
      #    5.846933178718317e-02, 7.537960906527452e-02, 7.043734921992200e-01)
      #test.valueN <-c(-5.559389931781886e-01,
      #                   -1.902431069568611e+00,  1.524595894866778e+00,
      #                   -7.863494805034426e-01, 1.328128164898773e-01)
      # values from 0.99.0a
      # test.valueU <-c(0.25603057077527752, 0.07879165329010961,
      # 		 0.60394682330171257, 0.20843868707503158, 0.97636939375111098)
      # test.valueN <-c( -1.39965726956837644, -0.24025807684466990,
      #          2.34362137305187446, -0.66321208109989371, -0.71183559894654214)
      # values from 1.1.0 (devel) should also work in 1.0.1
      test.valueU.R1.6.2 <-c(0.59132864479704950, 0.76406894316060192,
                  0.18576870606880833, 0.81087542344137897, 0.05835627439859235)
      test.valueN.R1.7.1 <-c( 0.959409416512901569, 0.046546246157109825,
            -0.775306427558083988, -0.777761120323404387, -1.363043207314097227)
      test.valueN.R1.6.2 <-c( 0.959409416509782953, 0.046546246156130192,
            -0.775306427558391964, -0.777761120325662803, -1.363043207314267313)
      test.valueU <-c(0.6159259434789419, 0.2200001638848335,
             0.7024894717615098, 0.4869760072324425, 0.2284618646372110)
      test.valueN <- c(0.2947980960879867, 0.5315740209738240,
            -0.7439218830522146, -0.9861002105572579, -0.3542773118879623)
     }
  if (!is.R()) 
     {# above should really be if (is.Splus()) 
      test.seed<- c(37, 39, 39, 4, 7, 2, 27, 58, 38, 15, 32, 2)
      test.valueU.S <- c(0.4299328043125570, 0.3092006836086512,
            0.5808096211403608, 0.3061958812177181, 0.8137333435006440)
      test.valueN.S <- c( -0.7613318231781665, -0.5724360196433543,
            0.8536399448225964, -0.2269096022522968, -0.8126790170570223)
     }

  old.seed <- setRNG(kind="default", seed=test.seed, normal.kind="default")
  on.exit(setRNG(old.seed))

  ok <- TRUE
if ( !is.R())
 {if (1e-14 < max(abs(runif(5)-test.valueU.S)))
    {warning("The default runif number generator has been changed.")
     ok <- FALSE
    }

   setRNG(kind="default", seed=test.seed, normal.kind="default")

   if (1e-14  < max(abs(rnorm(5)-test.valueN.S)))
    {warning("The default rnorm number generator has been changed.")
     ok <- FALSE
    }
 }
if ( is.R())
 {if (as.numeric(version$major)+0.1*as.numeric(version$minor) > 1.62 )
  {if (1e-14 < max(abs(runif(5)-test.valueU)))
    {warning("The default runif number generator has been changed.")
     ok <- FALSE
    }

   setRNG(kind="default", seed=test.seed, normal.kind="default")

   if (1e-14  < max(abs(rnorm(5)-test.valueN)))
    {warning("The default rnorm number generator has been changed.")
     ok <- FALSE
    }
  } 
# As of R 1.7.0 the default generator changed. For R 0.99 to 1.6.2 the
#   following corresponded to kind="default",normal.kind="default"
  setRNG(kind="Marsaglia-Multicarry", seed=test.seed, normal.kind="Kinderman-Ramage")

  if (1e-14 < max(abs(runif(5)-test.valueU.R1.6.2)))
    {warning("The Marsaglia-Multicarry runif number generator has been changed.")
     ok <- FALSE
    }
 
# As of R 1.7.1 a bug was recognized in Kinderman-Ramage and it was changed. The
#   old version was named "Buggy Kinderman-Ramage"
  setRNG(kind="Marsaglia-Multicarr", seed=test.seed, normal.kind="Kinderman-Ramage")

if (as.numeric(version$major)+0.1*as.numeric(version$minor) >= 1.71 )
  {if (1e-14  < max(abs(rnorm(5) -  test.valueN.R1.7.1  )))
    {warning("The Kinderman-Ramage rnorm number generator has been changed.")
     ok <- FALSE}
   setRNG(kind="Marsaglia-Multicarr", seed=test.seed, normal.kind="Buggy Kinderman-Ramage")
   if (1e-14  < max(abs(rnorm(5) - test.valueN.R1.6.2 )))
    {warning("The Buggy Kinderman-Ramage rnorm number generator has been changed.")
     ok <- FALSE
    }
  }
 else if (1e-14  < max(abs(rnorm(5) - test.valueN.R1.6.2 )))
    {warning("The Kinderman-Ramage rnorm number generator (pre R 1.7.1)has been changed.")
     ok <- FALSE
    }
 }

  setRNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
  if (setRNG()$kind        != "Wichmann-Hill" ||
      setRNG()$normal.kind != "Box-Muller"    ||
      all(setRNG()$seed    != c(979,1479,1542) )) 
     {warning("RNG is not being set properly")
      ok <- FALSE
     }
  if (1e-14 < max(abs(runif(5) -
      c(0.56933540553339546, 0.10513577518521355, 0.05846933178718317,
        0.07537960906527452, 0.70437349219921996))))
    {warning("The Wichmann-Hill runif number generator has been changed.")
     ok <- FALSE
    }

# for the next R 1.0.0 
# setRNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
# rnorm gives
#[1] -1.92425218107175877 -0.89568905204068128  2.12213361588187510
#[4]  0.81669202948845299 -0.13569189805256629
# as does
  setRNG(kind="Wichmann-Hill", seed=c(979,1479,1542), normal.kind="Box-Muller")
  if (1e-14 < max(abs(rnorm(5) -
      c(-1.92425218107175877, -0.89568905204068128,  2.12213361588187510,
         0.81669202948845299, -0.13569189805256629)
      # pre R 1,0,0 c(0.4605069059114530, 0.7685565310963474, -0.3737680932387061,
      # 0.5926372779538560, 1.4995245125275518)
	)))
    {warning("The Box-Muller rnorm number generator has been changed.")
     ok <- FALSE
    }
  # this needs to be done twice with odd and even n to chech completely
  if (1e-14 < max(abs(rnorm(5) -
      c(-0.4602838255495997,  0.2761541652763959,  1.3265434523731297,
        0.6856247181400722, -1.8336523890846541) )))
    {warning("The Box-Muller rnorm state is not properly preserved.")
     ok <- FALSE
    }
  if (1e-14 < max(abs(rnorm(6) -
      c(1.9850437759531543,  0.6107700961454843, -0.9419893721776470,
        1.1031328847642050,  0.4184702210057414,  0.9167797157851526) )))
    {warning("The Box-Muller rnorm state is not properly preserved.")
     ok <- FALSE
    }

  if (1e-14 < max(abs(rnorm(6) -
      c(-0.724539745179790251, -0.439138566092752758,  1.466237618877826110,
         0.289289597398559639,  0.003007778996985022,  1.008712871048744297) )))
    {warning("The Box-Muller rnorm state is not properly preserved.")
     ok <- FALSE
    }

  if (ok) {cat("ok\n"); invisible(TRUE) } else { cat("failed!\n"); stop("failed")}
 }
  
