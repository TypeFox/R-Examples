setClass("virtual.spodt")

setClass("f.spodt",
         representation(id="numeric", n="numeric", m="numeric", v="numeric", G="numeric"),
         prototype(id=0, n=0, m=0, v=0, G=c(0,0)),
         contains="virtual.spodt")

setClass("n.spodt",
         representation(n="numeric", m="numeric", v="numeric", R2="numeric",
                        fg="virtual.spodt", fd="virtual.spodt"),
         prototype(n=0, m=0, v=0, R2=0,
                   fg=new("f.spodt"), fd=new("f.spodt")),
         contains="virtual.spodt")
                             
setClass("sp.spodt",
         representation(coeff="numeric", int="numeric"),
         prototype(coeff=c(0,0), int=c(0,0)),
         contains=c("virtual.spodt", "n.spodt"))

setClass("vql.spodt",
         representation(vrbl="character", mod="character"),
         prototype(vrbl="", mod=""),
         contains=c("virtual.spodt", "n.spodt"))
         
setClass("vqt.spodt",
         representation(vrbl="character", seuil="numeric"),
         prototype(vrbl="", seuil=0),
         contains=c("virtual.spodt", "n.spodt"))

setClass("spodt",
         representation(racine="virtual.spodt",
                        R2="numeric", partition="vector", adj="matrix", cl.grf="matrix", sgmts.grf="matrix", brd="matrix"),
         prototype(racine=new("f.spodt"),
                   R2=0, partition=0, adj=matrix(0), cl.grf=matrix(ncol=3), sgmts.grf=matrix(ncol=4), brd=matrix(ncol=5)),
         contains="virtual.spodt")
         