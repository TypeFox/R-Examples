### ============================================================================
### Butcher tables for selected explicit ODE solvers of Runge-Kutta type
### Note that for fixed step methods A is a vector (the subdiagonal of matrix A)
###   For variable time step methods, A must be strictly lower triangular.
###   The underlying rk codes support explicit methods
###   and (still experimentally) some implicit methods.
### ============================================================================

rkMethod <- function(method = NULL, ...) {
  methods <- list(
    euler = list(ID = "euler",
        varstep = FALSE,
          A      = c(0),
          b1     = c(1),
          c      = c(0),
          stage  = 1,
          Qerr   = 1
    ),
    ## Heun's method
    rk2 = list(ID = "rk2",
        varstep = FALSE,
          A      = c(0, 1),
          b1     = c(0.5, 0.5),
          c      = c(0, 1),
          stage  = 2,
          Qerr   = 1
    ),
    ## classical Runge-Kutta 4th order method
    rk4 = list(ID = "rk4",
        varstep = FALSE,
          A      = c(0, .5, .5, 1),
          b1     = c(1/6, 1/3, 1/3, 1/6),
          c      = c(0, .5, .5, 1),
          stage  = 4,
          Qerr   = 4
    ),
    ## One of the numerous RK23 formulae
    rk23 = list(ID = "rk23",
      varstep = TRUE,
      FSAL    = FALSE,
      A  = matrix(c(0, 0, 0,
                  1/2, 0, 0,
                  -1, 2, 0), 3, 3, byrow = TRUE),
      b1 = c(0, 1, 0),
      b2 = c(1/6, 2/3, 1/6),
      c  = c(0, 1/2, 2),
      stage = 3,
      Qerr  = 2
    ),
    ## Bogacki & Shampine
    rk23bs = list(ID = "rk23bs",
      varstep = TRUE,
      FSAL    = TRUE,
      A  = matrix(c(0, 0, 0, 0,
                  1/2, 0, 0, 0,
                  0, 3/4, 0, 0,
                  2/9, 1/3, 4/9, 0), 4, 4, byrow = TRUE),
      b1 = c(7/24, 1/4, 1/3, 1/8),
      b2 = c(2/9, 1/3, 4/9, 0),
      c  = c(0, 1/2, 3/4, 1),
      stage = 4,
      Qerr  = 2
    ),
    ## RK-Fehlberg 34
    rk34f = list(ID = "rk34f",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0,
                      2/7, 0, 0, 0,
                      77/900, 343/900, 0, 0,
                      805/1444, -77175/54872, 97125/54872, 0,
                      79/490, 0, 2175/3626, 2166/9065),
                      5, 4, byrow = TRUE),
         b1 = c(79/490, 	0, 2175/3626, 2166/9065, 	0),
         b2 = c(229/1470, 0, 1125/1813, 13718/81585, 1/18),
         c  = c(0,	2/7, 	7/15, 35/38, 	1),
         stage = 5,
         Qerr  = 3
    ),
    ## RK-Fehlberg Method 45
    rk45f = list(ID = "rk45f",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/4, 0, 0, 0, 0,
                      3/32, 9/32, 0, 0, 0,
                      1932/2197, -7200/2197, 7296/2197, 0, 0,
                      439/216, -8, 3680/513, -845/4104, 0,
                      -8/27, 2, -3544/2565, 1859/4104, -11/40),
                      6, 5, byrow = TRUE),
         b1 = c(25/216, 	0, 	1408/2565, 	2197/4104, 	-1/5, 	0),
         b2 = c(16/135, 	0, 	6656/12825, 	28561/56430, 	-9/50, 	2/55),
         c  = c(0,	1/4, 	3/8, 	12/13, 	1, 	1/2),
         stage = 6,
         Qerr  = 4
    ),
    ## Cash-Karp method
    rk45ck = list(ID = "rk45ck",
         varstep = TRUE,
         FSAL = TRUE,
         A = matrix(c(0,    0,       0,         0,            0,
                      1/5,  0,       0,         0,            0,
                      3/40, 9/40,    0,         0,            0,
                      3/10, -9/10,   6/5,       0,            0,
                    -11/54, 5/2,    -70/27,     35/27,        0,
                   1631/55296, 175/512, 575/13824, 44275/110592, 253/4096),
                      6, 5, byrow = TRUE),
         b1 = c(2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4),
         b2 = c(37/378, 0, 250/621, 125/594, 0, 512/1771),
         c = c(0, 1/5, 3/10, 3/5,  1, 7/8),
         densetype = 2, # special dense output type 2
         stage = 6,
         Qerr = 4),
    ## England Method
    rk45e = list(ID = "rk45e",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/2, 0, 0, 0, 0,
                      1/4, 1/4, 0, 0, 0,
                      0, -1, 2, 0, 0,
                      7/27, 10/27, 0, 1/27, 0,
                      28/625, -125/625, 546/625, 54/625, -378/625),
                      6, 5, byrow = TRUE),
         b1 = c(1/6, 	0, 4/6, 1/6, 	0, 	0),
         b2 = c(14/336, 0, 0,	35/336, 162/336, 125/336),
         c  = c(0,	1/2, 	1/2, 	1, 	2/3, 	1/5),
         stage = 6,
         Qerr  = 4
    ),
    ## Prince-Dormand 5(4)6m
    rk45dp6 = list(ID = "rk45dp6",
         varstep = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0,
                      1/5, 0, 0, 0, 0,
                      3/40, 9/40, 0, 0, 0,
                      3/10, -9/10, 6/5, 0, 0,
                      226/729, -25/27, 880/729, 55/729, 0,
                      -181/270, 5/2, -266/297, -91/27, 189/55),
                      6, 5, byrow = TRUE),
         b1 = c(31/540, 	0, 190/297, -145/108, 351/220, 1/20),
         b2 = c(19/216, 0, 1000/2079,	-125/216, 81/88, 5/56),
         c  = c(0,	1/5, 	3/10, 3/5, 	2/3, 	1),
         stage = 6,
         Qerr  = 4
    ),
    ## Prince-Dormand 5(4)7m -- recommended by the Octave developers
    rk45dp7 = list(ID = "rk45dp7",
         varstep = TRUE,
         FSAL    = TRUE,
         A  = matrix(c(0, 0, 0, 0, 0, 0,
                      1/5, 0, 0, 0, 0, 0,
                      3/40, 9/40, 0, 0, 0, 0,
                      44/45, -56/15, 32/9, 0, 0, 0,
                      19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,
                      9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,
                      35/384, 0, 500/1113, 125/192, -2187/6784, 11/84),
                      7, 6, byrow = TRUE),
         b1 = c(5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40),
         b2 = c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0),
         c  = c(0, 1/5, 3/10, 4/5, 8/9, 1, 1),
         d  = c(-12715105075.0/11282082432.0, 0, 87487479700.0/32700410799.0,
                -10690763975.0/1880347072.0, 701980252875.0/199316789632.0,
                -1453857185.0/822651844.0, 69997945.0/29380423.0),
         densetype = 1, # default type of dense output formula, if available
         stage = 7,
         Qerr  = 4
    ),
    ## Prince-Dormand 78 method
    rk78dp = list(ID = "rk78dp",
      varstep = TRUE,
       FSAL = FALSE,
       A = matrix(c(
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0,
         3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0,
         29443841/614563906, 0, 0, 77736538/692538347,
         -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0,
         16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777,
         545815736/2771057229,  -180193667/1043307555, 0, 0, 0, 0, 0,
         39632708/573591083, 0, 0, -433636366/683701615,
         -421739975/2616292301, 100302831/723423059, 790204164/839813087,
         800635310/3783071287, 0, 0, 0, 0,
         246121993/1340847787, 0, 0, -37695042795/15268766246,
         -309121744/1061227803, -12992083/490766935, 6005943493/2108947869,
         393006217/1396673457, 123872331/1001029789, 0, 0, 0,
         -1028468189/846180014, 0, 0, 8478235783/508512852,
         1311729495/1432422823, -10304129995/1701304382,
         -48777925059/3047939560, 15336726248/1032824649,
         -45442868181/3398467696, 3065993473/597172653, 0, 0,
         185892177/718116043, 0, 0, -3185094517/667107341,
         -477755414/1098053517, -703635378/230739211,
         5731566787/1027545527, 5232866602/850066563,
         -4093664535/808688257, 3962137247/1805957418,
         65686358/487910083, 0,
         403863854/491063109, 0, 0, -5068492393/434740067,
         -411421997/543043805, 652783627/914296604, 11173962825/925320556,
         -13158990841/6184727034, 3936647629/1978049680,
         -160528059/685178525, 248638103/1413531060, 0),
         nrow = 13, ncol = 12 , byrow = TRUE),
      b1 = c(13451932/455176623, 0, 0, 0, 0, -808719846/976000145,
         1757004468/5645159321, 656045339/265891186,
         -3867574721/1518517206, 465885868/322736535,
         53011238/667516719, 2/45, 0),
      b2 = c(14005451/335480064, 0, 0, 0, 0, -59238493/1068277825,
            181606767/758867731, 561292985/797845732, -1041891430/1371343529,
            760417239/1151165299, 118820643/751138087, -528747749/2220607170,
            1/4),
      c = c(0, 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200,
        5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1),
      stage = 13,
      Qerr = 7
    ),

    ## Runge-Kutta-Fehlberg 78 method
    rk78f = list(ID = "rk78f",
        varstep = TRUE,
        FSAL    = FALSE,
        A  = matrix(
         c(rep(0,12),
         2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0,
         0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0,
         -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0,
         31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0,
         2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0,
         -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0,
         2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0,
         3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0,
         -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1
        ), nrow=13, ncol=12, byrow = TRUE),
        b1 = c(41/840, 0,0,0,0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840, 0, 0),
        b2 = c(0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840),
        c  = c(0, 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1),
        stage = 13,
        Qerr  = 7
    ),
    ## -------------------------------------------------------------------------
    ## Implicit methods; experimental!
    ## -------------------------------------------------------------------------

    ## Radau order 3
    irk3r = list(ID = "irk3r",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(
            c(5/12, -1/12,
              3/4,   1/4),
             nrow = 2, ncol = 2, byrow = TRUE),
      b1 = c(3/4, 1/4) ,
      c = c(1/3, 1/4),
      stage = 2,
      Qerr = 3
    ),

    ## Radau IIA order 5
    irk5r = list(ID = "irk5r",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(
            c((88-7*sqrt(6))/360, (296-169*sqrt(6))/1800, (-2+3*sqrt(6))/225,
              (296+169*sqrt(6))/1800, (88+7*sqrt(6))/360, (-2-3*sqrt(6))/225,
             (16-sqrt(6))/36,         (16+sqrt(6))/36,     1/9),
             nrow = 3, ncol = 3, byrow = TRUE),
      b1 = c((16-sqrt(6))/36, (16+sqrt(6))/36, 1/9),
      c  = c(0.4-sqrt(6)/10, 0.4+sqrt(6)/10, 1),
      stage = 3,
      Qerr = 5
    ),

    ## Hammer - Hollingsworth coefficients , order 4
    irk4hh = list(ID = "irk4hh",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(
            c(1/4,           1/4-sqrt(3)/6,
              1/4+sqrt(3)/6, 1/4),
             nrow = 2, ncol = 2, byrow = TRUE),
      b1 = c(1/2, 1/2),
      c  = c(0.5-sqrt(3)/6, 0.5+sqrt(3)/6),
      stage = 2,
      Qerr = 4
    ),

    ## Kuntzmann and Butcher order 6
    irk6kb = list(ID = "irk6kb",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(5/36,       2/9-sqrt(15)/15, 5/36 - sqrt(15)/30,
             5/36+sqrt(15)/24, 2/9,             5/36-sqrt(15)/24,
             5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36),
             nrow = 3, ncol = 3, byrow = TRUE),
      b1 = c(5/18, 4/9, 5/18),
      c  = c(1/2-sqrt(15)/10, 1/2, 1/2+sqrt(15)/10),
      stage = 3,
      Qerr = 6
    ),

    ## Lobatto order 4
    irk4l = list(ID = "irk4l",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(c(0,  0,  0,
                   1/4,1/4,0,
                   0,  1,  0),
                   nrow=3, ncol=3, byrow = TRUE),
      b1 = c(1/6, 2/3, 1/6) ,
      c  = c(0,   1/2, 1),
      stage = 3,
      Qerr = 4
    ),

    ## Lobatto order 6
    irk6l = list(ID = "irk6l",
      varstep = FALSE,
      implicit = TRUE,
      A = matrix(
            c(0,              0,                  0,                0,
              (5+sqrt(5))/60, 1/6,               (15-7*sqrt(5))/60, 0,
              (5-sqrt(5))/60, (15+7*sqrt(5))/60,  1/6,              0,
              1/6,            (5-sqrt(5))/12,    (5+sqrt(5))/12,    0),
              nrow = 4, ncol = 4, byrow = TRUE),
      b1 = c(1/12, 5/12, 5/12, 1/12) ,
      c = c(0,(5-sqrt(5))/10, (5+sqrt(5))/10, 1),
      stage = 4,
      Qerr = 6
    )
  )
  ## ---------------------------------------------------------------------------
  ## look if the method is known; ode23 and ode45 are used as synonyms
  ## ---------------------------------------------------------------------------
  knownMethods <- c(lapply(methods,"[[", "ID"), "ode23", "ode45")

  if (!is.null(method)) {
    method <- unlist(match.arg(method, knownMethods))
    if (method == "ode23")
      method <- "rk23bs"
    else if (method == "ode45")
      method <- "rk45dp7"

    out <- methods[[method]]
  } else {
    out <- vector("list", 0)
  }

  ## modify a known or add a completely new method)
  ldots <- list(...)
  out[names(ldots)] <- ldots

  ## return the IDs of the methods if called with an empty argument list
  if (is.null(method) & length(ldots) == 0) {
    out <- as.vector(unlist(knownMethods))
  } else {
    ## check size consistency of parameter sets
    sl    <- lapply(out, length)
    stage <- out$stage
    if (is.matrix(out$A)) {
      if (nrow(out$A) != stage | ncol(out$A)  < stage -1 | ncol(out$A) > stage)
        stop("Size of matrix A does not match stage")
    } else {
      if (length(out$A) != stage) stop("Size of A does not match stage")
    }
    if (stage != sl$b1 | stage != sl$c)
      stop("Wrong rkMethod, length of parameters do not match")
    if (out$varstep & is.null(out$b2))
      stop("Variable stepsize method needs non-empty b2")
    if (!is.null(out$b2))
      if (sl$b2 != stage)
        stop("Wrong rkMethod, length of b2 must be empty or equal to stage")
    if (!is.null(out[["d"]])) # exact argument matching!
      if (sl[["d"]] != stage)
        stop("Wrong rkMethod, length of d must be empty or equal to stage")
    
    ## check densetype
    if (!is.null(out$densetype)) {
      if (out$densetype == 1)
        if (!(out$ID %in% c("rk45dp7", "ode45")))
          stop("densetype = 1 not implemented for this method")
  
      if (out$densetype == 2)
        if (!(out$ID %in% c("rk45ck")))
          stop("densetype = 2 not implemented for this method")
    }    
    class(out) <- c("list", "rkMethod")
  }

  out
}
