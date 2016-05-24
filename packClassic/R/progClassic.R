######################################################
###                  ProgClassic.R                 ###
######################################################

publicA     <- function(x){plot(x+1)}
privateA    <- function(x){x+2}
.publicB    <- function(x){x+10}
.privateB   <- function(x){x+20}
publicC     <- function(x){publicA(privateA(x))}
privateC    <- function(x){.publicB(.privateB(x))}
