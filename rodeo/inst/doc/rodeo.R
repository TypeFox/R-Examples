### R code from vignette source 'rodeo.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rodeo.rnw:35-36
###################################################
options(continue=" ", prompt=" ")


###################################################
### code chunk number 2: rodeo.rnw:120-123
###################################################
library(rodeo, quietly=TRUE)
data(exampleIdentifiers)
print(format(exampleIdentifiers, justify="right"), row.names=FALSE)


###################################################
### code chunk number 3: rodeo.rnw:133-136
###################################################
library(rodeo, quietly=TRUE)
data(exampleProcesses)
print(format(exampleProcesses, justify="right"), row.names=FALSE)


###################################################
### code chunk number 4: rodeo.rnw:146-149
###################################################
library(rodeo, quietly=TRUE)
data(exampleStoichiometry)
print(format(exampleStoichiometry, justify="right"), row.names=FALSE)


###################################################
### code chunk number 5: rodeo.rnw:167-177
###################################################
library(rodeo, quietly=TRUE)

# Load sample data frames (contents shown above)
data(exampleIdentifiers, exampleProcesses, exampleStoichiometry)

# Instantiate new object
model= new("rodeo", vars=subset(exampleIdentifiers,type=="v"),
  pars=subset(exampleIdentifiers,type=="p"),
  funs=subset(exampleIdentifiers,type=="f"),
  pros=exampleProcesses, stoi=exampleStoichiometry)


###################################################
### code chunk number 6: rodeo.rnw:182-187 (eval = FALSE)
###################################################
## # Built-in method
## model$show()
## 
## # Show stoichiometry information as a matrix
## print(model$stoichiometry())


###################################################
### code chunk number 7: rodeo.rnw:194-204
###################################################
# 'normal' functions
O2sat= function(t) {
  14.652 - 0.41022*t + 0.007991*t^2 - 0.000077774*t^3
}
ka= function(u, d) {
  (0.728*sqrt(u) - 0.317*u + 0.0372*u^2) / d / 86400
}
monod= function(s,h) {
  s / (s + h)
}


###################################################
### code chunk number 8: rodeo.rnw:207-214
###################################################
# forcings are functions of special variable 'time'
c_z_in= function(seconds) {
  0.1 * seconds/(7*86400 + seconds)
}
c_do_in= function(seconds) {
  10.  # taken as constant
}


###################################################
### code chunk number 9: setDataSingleBox
###################################################
pars= list(kd=5.78e-7, h_do=0.5, s_do_z=2.76, wind=1, depth=2,
 temp=20, q_in=1, q_ex=1)
vars= list(c_z=1, c_do=9.022, v=1.e6)
p= model$arrangePars(pars)
v= model$arrangeVars(vars)


###################################################
### code chunk number 10: rodeo.rnw:231-233
###################################################
m= model$stoichiometry(c(v, p, time=0))
print(signif(m, 3))


###################################################
### code chunk number 11: rodeo.rnw:244-246
###################################################
code= model$generate(name="derivs",lang="r")
derivs= eval(parse(text=code))


###################################################
### code chunk number 12: rodeo.rnw:255-261
###################################################
library(deSolve)
t= seq(0, 30*86400, 3600)
out= ode(y=v, times=t, func=derivs, parms=p, NLVL=1)
layout(matrix(1:9, ncol=3, byrow=TRUE))
plot(out, mfrow=NULL)
layout(1)


###################################################
### code chunk number 13: setDataMultiBox (eval = FALSE)
###################################################
## nbox= 3
## pars= list(kd=rep(5.78e-7, nbox), h_do=0.5, s_do_z=2.76, wind=1,
##   depth=2, temp=20, q_in=1, q_ex=1)
## vars= list(c_z=seq(from=0, to=50, length.out=nbox), c_do=9.022,
##   v=1.e6)
## p= model$arrangePars(pars)
## v= model$arrangeVars(vars)


###################################################
### code chunk number 14: rodeo.rnw:291-295
###################################################
nbox= 3
pars= list(kd=rep(5.78e-7, nbox), h_do=0.5, s_do_z=2.76, wind=1,
  depth=2, temp=20, q_in=1, q_ex=1)
vars= list(c_z=seq(from=0, to=50, length.out=nbox), c_do=9.022,
  v=1.e6)
p= model$arrangePars(pars)
v= model$arrangeVars(vars)
out= ode(y=v, times=t, func=derivs, parms=p, NLVL=nbox)
layout(matrix(1:nbox, nrow=1))
plot(out, which=paste("c_do",1:nbox,sep="."), mfrow=NULL)


###################################################
### code chunk number 15: rodeo.rnw:307-310
###################################################
code= model$generate(name="derivs",lang="f95")
# Optionally display generated code
#cat(code)


###################################################
### code chunk number 16: rodeo.rnw:338-339
###################################################
lib= model$compile(fileFun="functionsCode.f95", NLVL=nbox)


###################################################
### code chunk number 17: rodeo.rnw:347-351
###################################################
file_ffuns= "functionsCode.f95"
text= readLines(file_ffuns, n=-1L, ok=TRUE, warn=TRUE, encoding="unknown", skipNul=FALSE)
text= paste(text,"\n")
cat(text)


###################################################
### code chunk number 18: rodeo.rnw:359-366
###################################################
nbox= 3
pars= list(kd=rep(5.78e-7, nbox), h_do=0.5, s_do_z=2.76, wind=1,
  depth=2, temp=20, q_in=1, q_ex=1)
vars= list(c_z=seq(from=0, to=50, length.out=nbox), c_do=9.022,
  v=1.e6)
p= model$arrangePars(pars)
v= model$arrangeVars(vars)
dyn.load(lib["libFile"])
out= ode(y=v, times=t, func=lib["libFunc"], parms=p,
  dllname=lib["libName"], initfunc="initmod", nout=model$lenPros()*nbox)
layout(matrix(1:nbox, nrow=1))
dyn.unload(lib["libFile"])
plot(out, which=paste("c_do",1:nbox,sep="."), mfrow=NULL)


###################################################
### code chunk number 19: rodeo.rnw:369-371
###################################################
# Clean up dll file
invisible(file.remove(lib["libFile"]))


###################################################
### code chunk number 20: rodeo.rnw:392-397
###################################################
dat= data.frame(time=1:10, temp=round(rnorm(n=10, mean=20, sd=3)),
  humid=round(runif(10)*100))
write.table(x=dat, file="meteo.txt", col.names=TRUE,
  row.names=FALSE, sep="\t", quote=FALSE)
print(dat)


###################################################
### code chunk number 21: rodeo.rnw:402-408
###################################################
dat= data.frame(name=c("temp","humid"),
  column=c("temp","humid"), file="meteo.txt", mode=-1, default=FALSE)
code= forcingFunctions(dat)
write(x=code, file="forc.f95")
# Optionally inspect generated code
# cat(code)


###################################################
### code chunk number 22: rodeo.rnw:428-431
###################################################
text= readLines("fortranForcingsTest.f95", n=-1L, ok=TRUE, warn=TRUE, encoding="unknown", skipNul=FALSE)
text= paste(text,"\n")
cat(text)


###################################################
### code chunk number 23: exportTex (eval = FALSE)
###################################################
## # Select columns to export
## df= model$getVars()[,c("tex","unit","description")]
## # Define formatting functions
## bold= function(x){paste0("\\textbf{",x,"}")}
## mathmode= function(x) {paste0("$",x,"$")}
## # Export
## tex= exportDF(x=df, tex=TRUE,
##   colnames=c(tex="symbol"),
##   funHead=setNames(replicate(ncol(df),bold),names(df)),
##   funCell=list(tex=mathmode)
## )
## cat(tex)


###################################################
### code chunk number 24: rodeo.rnw:462-463
###################################################
# Select columns to export
df= model$getVars()[,c("tex","unit","description")]
# Define formatting functions
bold= function(x){paste0("\\textbf{",x,"}")}
mathmode= function(x) {paste0("$",x,"$")}
# Export
tex= exportDF(x=df, tex=TRUE,
  colnames=c(tex="symbol"),
  funHead=setNames(replicate(ncol(df),bold),names(df)),
  funCell=list(tex=mathmode)
)
cat(tex)


###################################################
### code chunk number 25: exportMarkdown (eval = FALSE)
###################################################
## to_markdown= function(dat, which_cols){
##   cols= which(names(dat) %in% which_cols)
##   for(i in cols){
##     dat[, i]= ifelse(dat[, i] != "", paste0("$", dat[, i], "$"), "")
##     }
##   return(dat)
## } 
## 
## ids= model$getVars()[,c("tex", "unit", "description")]
## names(ids)= c("Symbol", "Unit", "Description")
## kable(to_markdown(ids, which_cols=c("Symbol")),
##   caption = "State variables")


###################################################
### code chunk number 26: rodeo.rnw:500-502
###################################################
pars= list(kd=5.78e-7, h_do=0.5, s_do_z=2.76, wind=1, depth=2,
 temp=20, q_in=1, q_ex=1)
vars= list(c_z=1, c_do=9.022, v=1.e6)
p= model$arrangePars(pars)
v= model$arrangeVars(vars)
model$plotStoichiometry(values=c(v, p, time=0), cex=0.3)


###################################################
### code chunk number 27: rodeo.rnw:508-528
###################################################
pars= list(kd=5.78e-7, h_do=0.5, s_do_z=2.76, wind=1, depth=2,
 temp=20, q_in=1, q_ex=1)
vars= list(c_z=1, c_do=9.022, v=1.e6)
p= model$arrangePars(pars)
v= model$arrangeVars(vars)
signsymbol= function(x) {
  if (as.numeric(x) > 0) return("\\textcolor{orange}{$\\blacktriangle$}")
  if (as.numeric(x) < 0) return("\\textcolor{cyan}{$\\blacktriangledown$}")
  return("")
}
rot90= function(x) { paste0("\\rotatebox{90}
  {$",gsub(pattern="*", replacement="\\cdot ", x=x, fixed=TRUE),"$}") }
m= model$stoichiometry(c(v, p, time=0))
tbl= cbind(data.frame(process=rownames(m), stringsAsFactors=FALSE),
  as.data.frame(m))
tex= exportDF(x=tbl, tex=TRUE,
  colnames= setNames(c("",model$getVars()$tex[match(colnames(m),
    model$getVars()$name)]), names(tbl)),
  funHead= setNames(replicate(ncol(m),rot90), colnames(m)),
  funCell= setNames(replicate(ncol(m),signsymbol), colnames(m)),
  lines=TRUE
)
tex= paste0("%\n% THIS IS A GENERATED FILE\n%\n", tex)
# write(tex, file="/home/dkneis/temp/stoichiometry.tex")


###################################################
### code chunk number 28: rodeo.rnw:533-534
###################################################
write(tex, file="./stoiGraphics.tex")


###################################################
### code chunk number 29: rodeo.rnw:544-559
###################################################
signsymbol= function(x) {
  if (as.numeric(x) > 0) return("&#9651;")
  if (as.numeric(x) < 0) return("&#9661;")
  return("")
}
m= model$stoichiometry(c(v, p, time=0))
tbl= cbind(data.frame(process=rownames(m), stringsAsFactors=FALSE),
  as.data.frame(m))
html= exportDF(x=tbl, tex=FALSE,
  colnames= setNames(c("Process",model$getVars()$html[match(colnames(m),
    model$getVars()$name)]), names(tbl)),
  funCell= setNames(replicate(ncol(m),signsymbol), colnames(m))
)
html= paste("<html>", html, "</html>", sep="\n")
# write(html, file="/home/dkneis/temp/stoichiometry.html")


###################################################
### code chunk number 30: stoiMarkdown (eval = FALSE)
###################################################
## signsymbol= function(x) {
##   if (as.numeric(x) > 0) return("$\\blacktriangle$")
##   if (as.numeric(x) < 0) return("$\\blacktriangledown$")
##   return("")
## }
## stoi_mat= model$stoichiometry(c(v, p, time=0))
## stoi_mat= data.frame(apply(stoi_mat, MARGIN = c(1, 2), signsymbol))
## stoi_mat= setNames(stoi_mat, paste0("$",     
##   model$getVars()$tex[match(colnames(stoi_mat),
##   model$getVars()$name)], "$"))
## stoi_mat= cbind(Process=rownames(stoi_mat), stoi_mat)
## 
## kable(stoi_mat, row.names= FALSE, caption= "Stoichiometric matrix")


