### R code from vignette source 'apprex.Rnw'

###################################################
### code chunk number 1: loadlibs
###################################################
library("approximator")


###################################################
### code chunk number 2: datagen_design_matrix
###################################################
n <- 16
set.seed(0)
D1.vig <- latin.hypercube(n,2)


###################################################
### code chunk number 3: datagen_subsets
###################################################
subsets.vig <- subsets.fun(n,levels=4,prob=0.6)
names(subsets.vig) <- paste("level",1:4,sep=".")


###################################################
### code chunk number 4: headD1vig
###################################################
head(D1.vig)
nrow(D1.vig)


###################################################
### code chunk number 5: showsubsets
###################################################
subsets.vig


###################################################
### code chunk number 6: basis.vig.function.def
###################################################
basis.vig <- 
function (x) 
{
    if (is.vector(x)) {
        stopifnot(length(x) == 2)
        out <- c(1, x , x[1]*x[2])
        names(out) <- c("const", LETTERS[1:2], "interaction")
        return(out)
    }
    else {
        return(t(apply(x, 1, match.fun(sys.call()[[1]]))))
    }
}


###################################################
### code chunk number 7: basis.vig.in.action
###################################################
head(basis.vig(D1.vig))


###################################################
### code chunk number 8: datagen_generate_vig_fun
###################################################
"generate.vig.observations" <-
function (D1, subsets, basis.fun, hpa, betas = NULL, export.truth=FALSE)
{

      if(is.null(betas)){
        betas <-
          rbind(c(1, 2, 3, 4),
                c(1, 1, 3, 4),
                c(1, 1, 1, 4),
                c(1, 1, 1, 1))
        colnames(betas) <- c("const", LETTERS[1:2], "interaction")
        rownames(betas) <- paste("level", 1:4, sep = "")
      }
      
      if(export.truth){
        return(list(
                    hpa=hpa,
                    betas=betas
                    )
               )
      }


    sigma_squareds <- hpa$sigma_squareds
    B <- hpa$B
    rhos <- hpa$rhos
    delta <- function(i) {
      out <- rmvnorm(n = 1,
                     mean = basis.fun(D1[subsets[[i]], , drop =
                       FALSE]) %*% betas[i, ],
                     sigma = sigma_squareds[i] * corr.matrix(xold = D1[subsets[[i]], , drop = FALSE], pos.def.matrix = B[[i]])
                     )
      out <- drop(out)
      names(out) <- rownames(D1[subsets[[i]], , drop = FALSE])
      return(out)
    }
    
    use.clever.but.untested.method <- FALSE
    
    if(use.clever.but.untested.method){
      z1 <- delta(1)
      z2 <- delta(2) + rhos[1] * z1[match(subsets[[2]], subsets[[1]])]
      z3 <- delta(3) + rhos[2] * z2[match(subsets[[3]], subsets[[2]])]
      z4 <- delta(4) + rhos[3] * z3[match(subsets[[4]], subsets[[3]])]
      return(list(z1 = z1, z2 = z2, z3 = z3, z4 = z4))
    } else {
      out <- NULL
      out[[1]] <- delta(1)
      for(i in 2:length(subsets)){
        out[[i]] <- delta(i) + rhos[i-1] *
    out[[i-1]][match(subsets[[i]], subsets[[i-1]])]
      }
      return(out)
    }
  }


###################################################
### code chunk number 9: datagen_hpa.fun
###################################################
"hpa.fun.vig" <-
function (x) 
{
    if (length(x) != 15) {
        stop("x must have 19 elements")
    }
    "pdm.maker" <- function(x) {
        jj <- diag(x[1:2],nrow=2)
        rownames(jj) <- LETTERS[1:2]
        colnames(jj) <- LETTERS[1:2]
        return(jj)
    }
    sigma_squareds <- x[1:4]
    names(sigma_squareds) <- paste("level", 1:4, sep = "")
    B <- list()
    B[[1]] <- pdm.maker(x[05:06])
    B[[2]] <- pdm.maker(x[07:08])
    B[[3]] <- pdm.maker(x[09:10])
    B[[4]] <- pdm.maker(x[11:12])
    rhos <- x[13:15]
    names(rhos) <- paste("level", 1:3, sep = "")
    return(list(sigma_squareds = sigma_squareds, B = B, rhos = rhos))
}


###################################################
### code chunk number 10: datagen_hpa.vig_creator
###################################################
hpa.vig <- hpa.fun.vig(c(rep(0.01,4),rep(20,8),rep(1,3)))


###################################################
### code chunk number 11: datagen_generate_vig_obs
###################################################
z.vig <- 
generate.vig.observations(D1=D1.vig,subsets=subsets.vig, basis.fun=basis.vig,hpa=hpa.vig)


###################################################
### code chunk number 12: lapplyzvig
###################################################
lapply(z.vig,head)
lapply(z.vig,length)


###################################################
### code chunk number 13: hpafig
###################################################
hpa.vig


###################################################
### code chunk number 14: morestuff
###################################################
     jj <- list(trace=100,maxit=10)
     hpa.vig.level1 <- opt.1(D=D1.vig, z=z.vig, basis=basis.vig, subsets=subsets.vig, hpa.start=hpa.vig,control=jj)
     hpa.vig.level1


###################################################
### code chunk number 15: fourlevels
###################################################
     jj <- list(trace=0,maxit=4)
     hpa.vig.level2 <- opt.gt.1(level=2, D=D1.vig, z=z.vig, basis=basis.vig, subsets=subsets.vig, hpa.start=hpa.vig.level1,control=jj)
     hpa.vig.level3 <- opt.gt.1(level=3, D=D1.vig, z=z.vig, basis=basis.vig, subsets=subsets.vig, hpa.start=hpa.vig.level2,control=jj)
     hpa.vig.level4 <- opt.gt.1(level=4, D=D1.vig, z=z.vig, basis=basis.vig, subsets=subsets.vig, hpa.start=hpa.vig.level3,control=jj)
     hpa.vig.level4


###################################################
### code chunk number 16: betahatapp
###################################################
betahat.app(D1=D1.vig,subsets=subsets.vig,basis=basis.vig, hpa=hpa.vig.level4, z=z.vig)


###################################################
### code chunk number 17: mdashfun
###################################################
mdash.fun(x=c(0.5,0.5),D1=D1.vig, subsets=subsets.vig,hpa=hpa.vig.level4,z=z.vig,basis=basis.vig)


###################################################
### code chunk number 18: error
###################################################
cdash.fun(x=c(0.5,0.5), D1=D1.vig, subsets=subsets.vig,
  basis=basis.vig, hpa=hpa.vig)


###################################################
### code chunk number 19: showbutdonotevaldatagen (eval = FALSE)
###################################################
## n <- 16
## set.seed(0)
## D1.vig <- latin.hypercube(n,2)


###################################################
### code chunk number 20: subsetshow (eval = FALSE)
###################################################
## subsets.vig <- subsets.fun(n,levels=4,prob=0.6)
## names(subsets.vig) <- paste("level",1:4,sep=".")


###################################################
### code chunk number 21: hpafunvig
###################################################
hpa.fun.vig


###################################################
### code chunk number 22: callthis (eval = FALSE)
###################################################
## hpa.vig <- hpa.fun.vig(c(rep(0.01,4),rep(20,8),rep(1,3)))


###################################################
### code chunk number 23: hpauser (eval = FALSE)
###################################################
## 


###################################################
### code chunk number 24: datagen_hpa.fun
###################################################



###################################################
### code chunk number 25: gendat2 (eval = FALSE)
###################################################
## "generate.vig.observations" <-
## function (D1, subsets, basis.fun, hpa, betas = NULL, export.truth=FALSE)
## {
## 
##       if(is.null(betas)){
##         betas <-
##           rbind(c(1, 2, 3, 4),
##                 c(1, 1, 3, 4),
##                 c(1, 1, 1, 4),
##                 c(1, 1, 1, 1))
##         colnames(betas) <- c("const", LETTERS[1:2], "interaction")
##         rownames(betas) <- paste("level", 1:4, sep = "")
##       }
##       
##       if(export.truth){
##         return(list(
##                     hpa=hpa,
##                     betas=betas
##                     )
##                )
##       }
## 
## 
##     sigma_squareds <- hpa$sigma_squareds
##     B <- hpa$B
##     rhos <- hpa$rhos
##     delta <- function(i) {
##       out <- rmvnorm(n = 1,
##                      mean = basis.fun(D1[subsets[[i]], , drop =
##                        FALSE]) %*% betas[i, ],
##                      sigma = sigma_squareds[i] * corr.matrix(xold = D1[subsets[[i]], , drop = FALSE], pos.def.matrix = B[[i]])
##                      )
##       out <- drop(out)
##       names(out) <- rownames(D1[subsets[[i]], , drop = FALSE])
##       return(out)
##     }
##     
##     use.clever.but.untested.method <- FALSE
##     
##     if(use.clever.but.untested.method){
##       z1 <- delta(1)
##       z2 <- delta(2) + rhos[1] * z1[match(subsets[[2]], subsets[[1]])]
##       z3 <- delta(3) + rhos[2] * z2[match(subsets[[3]], subsets[[2]])]
##       z4 <- delta(4) + rhos[3] * z3[match(subsets[[4]], subsets[[3]])]
##       return(list(z1 = z1, z2 = z2, z3 = z3, z4 = z4))
##     } else {
##       out <- NULL
##       out[[1]] <- delta(1)
##       for(i in 2:length(subsets)){
##         out[[i]] <- delta(i) + rhos[i-1] *
##     out[[i-1]][match(subsets[[i]], subsets[[i-1]])]
##       }
##       return(out)
##     }
##   }


###################################################
### code chunk number 26: calld (eval = FALSE)
###################################################
## z.vig <- 
## generate.vig.observations(D1=D1.vig,subsets=subsets.vig, basis.fun=basis.vig,hpa=hpa.vig)


###################################################
### code chunk number 27: zfigshow
###################################################
z.vig


###################################################
### code chunk number 28: SetTheBib
###################################################
bib <- system.file( "doc", "uncertainty.bib", package = "emulator" )
bib <- sub('.bib$','',bib)


###################################################
### code chunk number 29: usethebib
###################################################
cat( "\\bibliography{",bib,"}\n",sep='')


