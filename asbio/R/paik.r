paik <- function(formula, counts, resp.lvl = 2, data, circle.mult = .4, xlab = NULL, ylab = NULL, leg.title = NULL, leg.loc = NULL, show.mname = FALSE,...){

vars <- as.character(attr(terms(formula), "variables")[-1])
cond.var = vars[3]
rv <- data[,names(data)==vars[1]]
cv <- data[,names(data)==cond.var]
ov <- vars[vars!=vars[1]&vars!=cond.var]
or <- data[,names(data)==ov]

new.formula <- formula(counts ~ rv + ov + cv)
cl <- levels(data[,names(data) == cond.var])
ol <- levels(data[,names(data) == ov])
rl <- levels(data[,names(data) == vars[1]])

x <- xtabs(count ~ rv + or + cv, data = data)
xm <- xtabs(count ~ rv + or, data = data)
m.prop <- apply(xm, 2, function(x)x[resp.lvl]/sum(x))


r.prop <- matrix(nrow=length(cl), ncol=length(ol), dimnames=list(paste(cond.var,cl,sep="."),paste(ov,ol,sep=".")))
r.sum <- matrix(nrow=length(cl), ncol=length(ol), dimnames=list(paste(cond.var,cl,sep="."),paste(ov,ol,sep=".")))

for(i in 1:length(cl)){
tab <- x[,,cl = cl[i]]
r.prop[i,] <- apply(tab, 2, function(x)x[resp.lvl]/sum(x))
r.sum[i,] <- apply(tab, 2, sum)
}

b<-barplot(1:(length(ol)+2), plot = FALSE)
pts <- b[-c(1,(length(ol)+2))]
y.loc <- stack(as.data.frame(r.prop))[,1]
x.loc <- rep(pts,each=length(cl))
temp.ylab <- bquote(paste("Proportion  ", .(vars[1]), "  =  ", .(rl[resp.lvl])))

p <- plot(x.loc, y.loc, xlim = c(1.5,(length(ol)+1.5)), type = "n", xaxt = "n", xlab = ifelse(is.null(xlab), ov, xlab), ylab = ifelse(is.null(ylab), eval(temp.ylab), ylab))
grid(p)
axis(1, at = pts , labels = ol)

tprop <- stack(as.data.frame(t(r.prop)))[,1]
tx <- rep(pts,length(cl))
tpropm <- matrix(ncol = length(pts), nrow = length(cl), data = tprop, byrow = TRUE)
txm <- matrix(ncol = length(pts), nrow = length(cl), data = tx, byrow = TRUE)

col <- gray(seq(1:length(cl))/length(cl))
circle.col <- rep(col, length(cl))
radii <- r.sum/sum(r.sum)
radii <- stack(as.data.frame(radii))[,1]*circle.mult
for(i in 1: length(radii))draw.circle(x.loc[i], y.loc[i], radii[i], col = circle.col[i])


if(length(pts)==2){
for(i in 1:length(cl)){
segments(txm[i,][1], tpropm[i,][1], txm[i,][2], tpropm[i,][2])
}
segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, lwd = 2)
}

if(length(pts)==3){
for(i in 1:length(cl)){
segments(txm[i,][1], tpropm[i,][1], txm[i,][2], tpropm[i,][2]);segments(txm[i,][2], tpropm[i,][2], txm[i,][3], tpropm[i,][3])
}
segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, lwd = 2); segments(pts[2], m.prop[2], pts[3], m.prop[3], lty = 2, lwd = 2)
}

if(length(pts)==4){
for(i in 1:length(cl)){
segments(txm[i,][1], tpropm[i,][1], txm[i,][2], tpropm[i,][2]);segments(txm[i,][2], tpropm[i,][2], txm[i,][3], tpropm[i,][3]);segments(txm[i,][3], tpropm[i,][3], txm[i,][4], tpropm[i,][4]) 
}
segments(pts[1], m.prop[1], pts[2], m.prop[2], lty = 2, lwd = 2); segments(pts[2], m.prop[2], pts[3], m.prop[3], lty = 2, lwd = 2); segments(pts[3], m.prop[3], pts[4], m.prop[4], lty = 2, lwd = 2)
}

if(length(pts)==5)stop("Number of rows in table must be less than 5")


points(x.loc, y.loc, pch = 19, cex=.6)
legend(ifelse(is.null(leg.loc),"topright",leg.loc),pch=rep(21,length(cl)),pt.bg = col, bg = "white", pt.cex = 1.5, title = ifelse(is.null(leg.title), cond.var, leg.title), legend = cl) 


degree <- diff(m.prop)/diff(pts)*360

if(show.mname==TRUE)text(mean(pts),mean(m.prop) + 0.03*max(y.loc),srt=degree, "Marginal prop.")
res <- invisible(list(marginal.prop = r.prop, group.prop = r.sum/sum(r.sum))) 
invisible(res)
}

