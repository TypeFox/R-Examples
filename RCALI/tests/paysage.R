# source("moninit")
 library("RCALI")    
file <- "paysage.data"
dispf <-  4
param <-
list(input=2,output=3,warn.conv=FALSE, send.and.receive=TRUE,
     poly=list(c(1,2), c(3,4), c(11,13)))

#calcul avec RCALI cub
# -----------------------
param$method <- "cub"
resfile <- "paysage.cubres"

start.run <- Sys.time()
califlopp(file=file, dispf=dispf, param = param, resfile = resfile)
end.run <- Sys.time()
duration <- end.run - start.run
# print(duration)


#calcul avec RCALI grid sur un couple
# -----------------------
param$method <- "grid"
param$poly <- c(1,1)
param$send.and.receive=FALSE
resfile <- "paysage.grid1.1res"


start.run <- Sys.time()
califlopp(file=file, dispf=dispf, param = param, resfile = resfile)
end.run <- Sys.time()
duration <- end.run - start.run
# print(duration)


# Charger le fichier de résultats cub
# ---------------------------------
#avec la fonction getRes de Rcali
cat("avec la fonction getRes paysage.cubres de Rcali\n")
start.run <- Sys.time()
dispcub.res <- getRes("paysage.cubres")
end.run <- Sys.time()
duration <- end.run - start.run
# print(duration)

#avec un read.table
cat("avec un read.table\n")
start.run <- Sys.time()
dispcubtab.res <- read.table("paysage.cubres", skip=1)
end.run <- Sys.time()
duration <- end.run - start.run
# print(duration)
