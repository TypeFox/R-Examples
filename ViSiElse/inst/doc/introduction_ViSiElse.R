## ----fig.show='asis',fig.width=7-----------------------------------------
    rm(list = ls())
    library(ViSiElse)
    vars <- c( "Action1" , "Action2" , "Action3" , "Action4" , "Action5" , 
               "Action6" , "Action7","Action8")
    label <- c( "Action 1", "Action 2", "Action 3", "Action 4", "Action 5", 
                "Action 6", "Action 7","Action 8" )
    typeA <- c( "p" , "l" , "p" , "p" , "p", "l" , "p" , "p" )
    deb <- rep(NA,8)
    deb[2] <- vars[1]
    deb[6] <- vars[4]
    fin <- rep( NA, 8)
    fin[2] <- vars[3]
    fin[6] <- vars[7]
    showorder <- c(1, 2 ,3 ,NA ,4 , 5, NA, 6 )
    book_df <- data.frame(vars ,label ,typeA ,showorder, deb,fin)
    book <-  ConvertoViSibook(book_df)
    plot(book)
    print(book)

## ----fig.show='asis',fig.width=7-----------------------------------------
    
  GZDeb <- rep(NA,8)
  GZDeb[1] <- 10
  GZDeb[5] <- 35
  GZDeb[8] <- 60
  
  GZFin <- rep(NA,8)
  GZFin[1] <- 25
  GZFin[5] <- 40
  GZFin[8] <- 80
  
  Repetition <- rep(NA,8)
  Repetition[5] <- 20
  
  BZBeforeDeb <- rep(NA,8)
  BZBeforeDeb[ 8 ] <- 0
  
  BZBeforeFin <- rep(NA,8)
  BZBeforeFin[ 8 ] <- 55
  
  BZAfterDeb <- rep(NA,8)
  BZAfterDeb[ 8 ] <- 95
  BZAfterDeb[ 1 ] <- 25
  
  BZAfterFin <- rep(NA,8)
  BZAfterFin[ 8 ] <- Inf
  BZAfterFin[ 1 ] <- Inf
  
  BZLong <- rep(NA,8)
  BZLong[2] <- 5
  BZLong[6] <- 40
  
  BZLtype <- rep(NA,8)
  BZLtype[2] <- "span"
  BZLtype[6] <- "time"
  
  book_df <- data.frame(vars ,label ,typeA ,showorder, deb,fin,
                        GZDeb,GZFin, Repetition,BZBeforeDeb,BZBeforeFin,
                        BZAfterDeb,BZAfterFin, BZLong, BZLtype)
    book <-  ConvertoViSibook(book_df)
    plot(book)
    print(book)

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
Action1 <- rbinom(50, 25, 0.5)
Action3 <- Action1 + rbinom(50, 10, 0.9)
Action4 <- Action3 + rbinom(50, 20, 1/2)
Action5 <- Action4 + rbinom(50, 5, 1/3)
Action7 <- Action5 + rbinom(50, 100, 1/5)
Action8 <- Action7 + rbinom(50, 110, 1/3)
X <- data.frame(id = seq(1,50), Action1, Action3, Action4, Action5, Action7, Action8)
head(X)

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
Action5sup1 <- Action5 + rep( 20, 50)
Action5sup2 <- Action5sup1[ seq( 1,25)] +  rep( 20, 25)
Action5sup3 <- Action5sup2 +  rep( 5, 25)
Action5sup <- c(Action5sup1, Action5sup2, Action5sup3)
id <- c( seq( 1,50), seq( 1,25), seq( 1,25))
Action1sup <- rep( NA, 50 + 25  + 25  )
Action3sup <- rep( NA, 50 + 25  + 25  )
Action4sup <- rep( NA, 50 + 25  + 25  )
Action7sup <- rep( NA, 50 + 25  + 25  )
Action8sup <- rep( NA, 50 + 25  + 25  )
Xsup <- data.frame(id , Action1sup, Action3sup, 
                   Action4sup, Action5sup, Action7sup, Action8sup)
colnames(Xsup) <- colnames(X)
head(Xsup)

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup )
plot(x)


## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, pixel = 5, Xsup = Xsup )
plot(x, main = "PIXEL = 5")
x <-  buildViSiGrid( X = X, book = book, pixel = 20, Xsup = Xsup )
plot(x, main = "PIXEL = 20")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup )
plot(x, main = "PIXEL = 1")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup )
plot(x, main = "informer =  Median")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup,  informer = "mean" )
plot(x, main = "informer =  Mean")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup,  informer = NULL )
plot(x, main = "informer =  NULL")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup )
plot(x, main = "t_0 =  0")
x <-  buildViSiGrid( X = X, book = book, pixel = 1,  t_0 = 10 )
plot(x, main = "t_0 = 10")
x <-  buildViSiGrid( X = X, book = book, pixel = 1,  t_0 = "Action1" )
plot(x, main = "t_0 = Action1")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
book_change <- changeShoworder( book, c(2,3,4,5,6) )
book_change[5,7] <- NA
book_change[5,8] <- NA
book_change[5,9] <- NA
x <-  buildViSiGrid( X = X, book = book_change, pixel = 1,  t_0 = "Action1" )
plot(x, main = "t_0 = Action1")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup )
plot(x, main = " quantity = N")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup,  quantity = "dens" )
plot(x, main = " quantity = dens")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup )
plot(x, main = " sorted.line = TRUE")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, Xsup = Xsup, sorted.line = FALSE )
plot(x, main = "sorted.line = FALSE")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
group = c(rep( "G1", 15),rep("G2",35))
x <-  buildViSiGrid( X = X, book = book, pixel = 1, group = group, method = "cut" )
plot(x, main = " Method = cut ")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, group = group, method = "join" )
plot(x, main = " Method = join ")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, group = group, method = "within", 
                     grwithin = "G1"  )
plot(x, main = " Method = within ; group1 ")
x <-  buildViSiGrid( X = X, book = book, pixel = 1, group = group, method = "within", 
                     grwithin = "G2", quantity = "dens"  )
plot(x, main = " Method = within ; group2 / density of ind ")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, group = group, method = "cut" , pixel = 1, 
                     tests = TRUE )
plot(x, main = " test, threshold p-value 0.01")
x <-  buildViSiGrid( X = X, book = book, group = group, method = "cut" , pixel = 1, 
                     tests = TRUE, threshold.test = 0.5  )
plot(x, main = " test, threshold p-value 0.5")

## ----fig.show='asis',fig.width=7, fig.height=7---------------------------
x <-  buildViSiGrid( X = X, book = book, group = group, method = "cut" , pixel = 1, 
                     Xsup = Xsup )
plot(x, main = " decrgr2=FALSE")
x <-  buildViSiGrid( X = X, book = book, group = group, method = "cut" , pixel = 1, 
                     Xsup = Xsup, decrgr2 = TRUE )
plot(x, main = "decrgr2 = TRUE")

