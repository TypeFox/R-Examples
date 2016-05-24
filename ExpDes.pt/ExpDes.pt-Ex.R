pkgname <- "ExpDes.pt"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ExpDes.pt')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ccF")
### * ccF

flush(stderr()); flush(stdout())

### Name: ccF
### Title: Teste de comparacoes multiplas de Calinski & Corsten baseado na
###   distribuicao F
### Aliases: ccF

### ** Examples

data(ex2)
attach(ex2)
dbc(trat, provador, aparencia, quali = TRUE, mcomp='ccf', sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("ccboot")
### * ccboot

flush(stderr()); flush(stdout())

### Name: ccboot
### Title: Comparacao multipla: Bootstrap
### Aliases: ccboot

### ** Examples

data(ex1)
attach(ex1)
dic(trat, ig, quali = TRUE, mcomp='ccboot', sigF = 0.05)



cleanEx()
nameEx("dbc")
### * dbc

flush(stderr()); flush(stdout())

### Name: dbc
### Title: Delineamento em Blocos Casualizados
### Aliases: dbc

### ** Examples

data(ex2)
attach(ex2)
dbc(trat, provador, aparencia, quali = TRUE, mcomp="lsd", sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("dic")
### * dic

flush(stderr()); flush(stdout())

### Name: dic
### Title: Delineamento Inteiramente Casualizado Simples
### Aliases: dic

### ** Examples

data(ex1)
attach(ex1)
dic(trat, ig, quali = FALSE, sigF = 0.05)




cleanEx()
nameEx("dql")
### * dql

flush(stderr()); flush(stdout())

### Name: dql
### Title: Delineamento em Quadrado Latino
### Aliases: dql

### ** Examples

data(ex3)
attach(ex3)
dql(trat, linha, coluna, resp, quali = TRUE, mcomp = "snk", sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("est21Ad")
### * est21Ad

flush(stderr()); flush(stdout())

### Name: est21Ad
### Title: Percevejos no milho: tratamento adicional
### Aliases: est21Ad

### ** Examples

data(est21Ad)
## maybe str(est21Ad) ; plot(est21Ad) ...



cleanEx()
nameEx("ex")
### * ex

flush(stderr()); flush(stdout())

### Name: ex
### Title: Videiras: parcelas subdivididas em DBC
### Aliases: ex

### ** Examples

data(ex)



cleanEx()
nameEx("ex1")
### * ex1

flush(stderr()); flush(stdout())

### Name: ex1
### Title: Yacon: DIC
### Aliases: ex1

### ** Examples

data(ex1)



cleanEx()
nameEx("ex2")
### * ex2

flush(stderr()); flush(stdout())

### Name: ex2
### Title: Barras alimenticias: DBC
### Aliases: ex2

### ** Examples

data(ex2)



cleanEx()
nameEx("ex3")
### * ex3

flush(stderr()); flush(stdout())

### Name: ex3
### Title: Forrageiras: DQL
### Aliases: ex3

### ** Examples

data(ex3)



cleanEx()
nameEx("ex4")
### * ex4

flush(stderr()); flush(stdout())

### Name: ex4
### Title: Compostagem: fatorial duplo em DIC
### Aliases: ex4

### ** Examples

data(ex4)



cleanEx()
nameEx("ex5")
### * ex5

flush(stderr()); flush(stdout())

### Name: ex5
### Title: Barras alimenticias: fatorial duplo em DBC
### Aliases: ex5

### ** Examples

data(ex5)



cleanEx()
nameEx("ex6")
### * ex6

flush(stderr()); flush(stdout())

### Name: ex6
### Title: Dados ficticios 1
### Aliases: ex6

### ** Examples

data(ex6)



cleanEx()
nameEx("ex7")
### * ex7

flush(stderr()); flush(stdout())

### Name: ex7
### Title: Estatura de plantas de milho 21 dias apos a emergencia.
### Aliases: ex7

### ** Examples

data(ex7)



cleanEx()
nameEx("ex8")
### * ex8

flush(stderr()); flush(stdout())

### Name: ex8
### Title: Compostagem: fatorial duplo com um tratamento adicional em DIC
### Aliases: ex8

### ** Examples

data(ex8)



cleanEx()
nameEx("ex9")
### * ex9

flush(stderr()); flush(stdout())

### Name: ex9
### Title: Coberturas vegetais: parcelas subdivididas em DIC
### Aliases: ex9

### ** Examples

data(ex9)



cleanEx()
nameEx("fat2.ad.dbc")
### * fat2.ad.dbc

flush(stderr()); flush(stdout())

### Name: fat2.ad.dbc
### Title: Fatorial duplo com um tratamento adicional em DBC
### Aliases: fat2.ad.dbc

### ** Examples

data(ex7)
attach(ex7)
data(est21Ad)
fat2.ad.dbc(periodo, nivel, bloco, est21, est21Ad, quali = c(TRUE, FALSE), mcomp = "sk", fac.names = c("Periodo", "Nivel"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat2.ad.dic")
### * fat2.ad.dic

flush(stderr()); flush(stdout())

### Name: fat2.ad.dic
### Title: Fatorial duplo com um tratamento adicional em DIC
### Aliases: fat2.ad.dic

### ** Examples

data(ex8)
attach(ex8)
data(secaAd)
fat2.ad.dic(inoculante, biodiesel, vaso, seca, secaAd, quali = c(TRUE,FALSE), mcomp = "tukey", fac.names = c("Inoculante", "Biodiesel"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat2.dbc")
### * fat2.dbc

flush(stderr()); flush(stdout())

### Name: fat2.dbc
### Title: Fatorial duplo em DBC
### Aliases: fat2.dbc

### ** Examples

data(ex5)
attach(ex5)
fat2.dbc(trat, genero, bloco, sabor, quali=c(TRUE,TRUE), mcomp="lsd", fac.names=c("Amostras","Genero"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat2.dic")
### * fat2.dic

flush(stderr()); flush(stdout())

### Name: fat2.dic
### Title: Fatorial duplo em DIC
### Aliases: fat2.dic

### ** Examples

data(ex4)
attach(ex4)
fat2.dic(revol,esterco,zn,quali=c(FALSE,TRUE),mcomp="tukey",fac.names=c("Revolvimento","Esterco"),sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.ad.dbc")
### * fat3.ad.dbc

flush(stderr()); flush(stdout())

### Name: fat3.ad.dbc
### Title: Fatorial triplo com um tratamento adicional em DBC
### Aliases: fat3.ad.dbc

### ** Examples

data(ex6)
attach(ex6)
data(respAd)
fat3.ad.dbc(fatorA, fatorB, fatorC, rep, resp, respAd, quali = c(TRUE, TRUE, TRUE), mcomp = "snk", fac.names = c("Fator A", "Fator B", "Fator C"),
 sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.ad.dic")
### * fat3.ad.dic

flush(stderr()); flush(stdout())

### Name: fat3.ad.dic
### Title: Fatorial triplo com um tratamento adicional em DIC
### Aliases: fat3.ad.dic

### ** Examples

data(ex6)
attach(ex6)
data(respAd)
fat3.ad.dic(fatorA, fatorB, fatorC, rep, resp, respAd, quali = c(TRUE, TRUE, TRUE), mcomp = "duncan", fac.names = c("Fator A", "Fator B", "Fator C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.dbc")
### * fat3.dbc

flush(stderr()); flush(stdout())

### Name: fat3.dbc
### Title: Fatorial triplo em DBC
### Aliases: fat3.dbc

### ** Examples

data(ex6)
attach(ex6)
fat3.dbc(fatorA, fatorB, fatorC, rep, resp, quali = c(TRUE, TRUE, TRUE), mcomp = "tukey", fac.names = c("Fator A", "Fator B", "Fator C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.dic")
### * fat3.dic

flush(stderr()); flush(stdout())

### Name: fat3.dic
### Title: Fatorial triplo em DIC
### Aliases: fat3.dic

### ** Examples

data(ex6)
attach(ex6)
fat3.dic(fatorA, fatorB, fatorC, resp, quali = c(TRUE, TRUE, TRUE), mcomp = "lsdb", fac.names = c("Fator A", "Fator B", "Fator C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("ginv")
### * ginv

flush(stderr()); flush(stdout())

### Name: ginv
### Title: Inversa generalizada
### Aliases: ginv

### ** Examples

## Not run: 
# The function is currently defined as
function(X, tol = sqrt(.Machine$double.eps))
{
## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X, dimnames = dnx[2:1])
}

## End(Not run)



cleanEx()
nameEx("lastC")
### * lastC

flush(stderr()); flush(stdout())

### Name: lastC
### Title: Setting the last character of a chain
### Aliases: lastC

### ** Examples

x<-c("a","ab","b","c","cd")
lastC(x)
# "a" "b" "b" "c" "d"



cleanEx()
nameEx("psub2.dbc")
### * psub2.dbc

flush(stderr()); flush(stdout())

### Name: psub2.dbc
### Title: Parcelas subdivididas em DBC
### Aliases: psub2.dbc

### ** Examples

data(ex)
attach(ex)
psub2.dbc(trat, dose, rep, resp, quali = c(TRUE, FALSE), mcomp = "tukey", fac.names = c("Tratamento", "Dose"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("psub2.dic")
### * psub2.dic

flush(stderr()); flush(stdout())

### Name: psub2.dic
### Title: Parcelas subdivididas em DIC
### Aliases: psub2.dic

### ** Examples

data(ex9)
attach(ex9)
psub2.dic(cobertura, prof, rep, pH, quali = c(TRUE, TRUE), mcomp = "lsd", fac.names = c("Cobertura", "Profundidade"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("respAd")
### * respAd

flush(stderr()); flush(stdout())

### Name: respAd
### Title: Dados ficticios: tratamento adicional
### Aliases: respAd

### ** Examples

data(respAd)
## maybe str(respAd) ; plot(respAd) ...



cleanEx()
nameEx("secaAd")
### * secaAd

flush(stderr()); flush(stdout())

### Name: secaAd
### Title: Compostagem: tratamento adicional
### Aliases: secaAd

### ** Examples

data(secaAd)
## maybe str(secaAd) ; plot(secaAd) ...



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
