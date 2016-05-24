## ----libraries, message = FALSE------------------------------------------
library(ltbayes)
library(ggplot2)
library(gridExtra)
library(reshape2)

## ----fig1----------------------------------------------------------------
samp <- 5000
burn <- 1000

alph <- c(1.27,1.34,1.14,1,0.67)
beta <- c(1.19,0.59,0.15,-0.59,-2)
gamm <- c(0.1,0.15,0.15,0.2,0.1)

zeta <- postsamp(fmodel3pl, c(0,0,1,1,1), 
	apar = alph, bpar = beta, cpar = gamm,
	control = list(nbatch = samp + burn, scale = 3))

zeta <- data.frame(sample = 1:samp, 
	zeta = zeta$batch[(burn + 1):(samp + burn)])

p <- ggplot(zeta, aes(x = sample, y = zeta))
p <- p + geom_line() + ylab(expression(zeta)) + xlab("Sample")
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(axis.text.x = element_text(size = 10, color = "black"))
p <- p + theme(text = element_text(size = 15))
p <- p + theme(axis.text.y = element_text(color = "black"))
p <- p + ggtitle("(a)\n")

p1 <- p

p <- ggplot(zeta, aes(x = zeta)) + geom_density(adjust = 2) + coord_flip()
p <- p + xlab(expression(zeta)) + ylab("Density")
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(axis.text.x = element_text(size = 10, color = "black"))
p <- p + theme(text = element_text(size = 15))
p <- p + theme(axis.text.y = element_text(color = "black"))
p <- p + ggtitle("(b)\n")

p2 <- p

grid.arrange(p1, p2, nrow = 1)

## ----fig2----------------------------------------------------------------
y <- patterns(5, 2, 0:5)
zeta <- vector("list", length = 32)
for (i in 1:32) {
  zeta[[i]] <- postsamp(fmodel3pl, y[i,], apar = alph, bpar = beta, cpar = gamm,
    control = list(nbatch = samp + burn))$batch[(burn + 1):(samp + burn)]
}
zeta <- data.frame(zeta = unlist(zeta))
zeta$pattern <- rep(apply(y, 1, paste, collapse = ""), each = samp)
zeta$sum <- rep(apply(y, 1, sum), each = samp)

p <- ggplot(zeta, aes(x = pattern, y = zeta))
p <- p + geom_violin(trim = FALSE)
p <- p + facet_grid(. ~ sum, scales = "free_x", space = "free_x")
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 15)) + xlab("Response Pattern")
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, 
	vjust = 0.5, color = "black"))
p <- p + theme(axis.text.y = element_text(color = "black")) 
p <- p + ylab(expression(zeta)) + xlab("Reponse Pattern")
p <- p + stat_summary(fun.y = mean, 
   fun.ymin = function(x) quantile(x, probs = 0.025), 
   fun.ymax = function(x) quantile(x, probs = 0.975), size = 0.5)

plot(p)

## ----fig3, message = FALSE-----------------------------------------------
zeta <- vector("list", length = 6)
for (i in 1:6) {
	y <- patterns(5, 2, i-1)
	zeta[[i]] <- postsamp(fmodel3pl, y, apar = alph, bpar = beta, cpar = gamm,
		control = list(nbatch = samp + burn))
	zeta[[i]] <- zeta[[i]]$batch[(burn+1):(samp + burn)]
}
zeta <- data.frame(zeta = unlist(zeta))
zeta$sum <- factor(rep(0:5, each = samp))

p <- ggplot(zeta, aes(x = zeta, color = sum, linetype = sum)) + stat_ecdf()
p <- p + ylab("") + xlim(c(-3,3)) + guides(color = guide_legend(title = "Sum"))
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(axis.text.x = element_text(size = 10, color = "black"))
p <- p + theme(text = element_text(size = 15))
p <- p + theme(axis.text.y = element_text(color = "black")) 
p <- p + ylab("") + xlab(expression(zeta))
p <- p + guides(linetype = guide_legend(title = "Sum"))
p <- p + ggtitle("(b)\n")

p1 <- p

p <- ggplot(zeta, aes(x = sum, y = zeta))
p <- p + geom_violin(trim = FALSE)
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 15)) + xlab("Sum Score")
p <- p + theme(axis.text.x = element_text(size = 10, color = "black"))
p <- p + theme(axis.text.y = element_text(color = "black")) + ylab(expression(zeta))
p <- p + stat_summary(fun.y = mean, 
	fun.ymin = function(x) quantile(x, probs = 0.025), 
	fun.ymax = function(x) quantile(x, probs = 0.975), size = 0.5)
p <- p + ggtitle("(a)\n")

p2 <- p

grid.arrange(p2, p1, nrow = 1)

## ----fig4, warning = FALSE-----------------------------------------------
zetal <- vector("list", length = 5)
zetau <- vector("list", length = 5)
for (i in 1:5) {
	zetal[[i]] <- postsamp(fmodel3pl, patterns(5, 2, 0:(i-1)), 
		apar = alph, bpar = beta, cpar = gamm,
		control = list(nbatch = samp + burn))
	zetau[[i]] <- postsamp(fmodel3pl, patterns(5, 2, i:5), 
		apar = alph, bpar = beta, cpar = gamm,
		control = list(nbatch = samp + burn))
	zetal[[i]] <- zetal[[i]]$batch[(burn+1):(samp + burn)]
	zetau[[i]] <- zetau[[i]]$batch[(burn+1):(samp + burn)]
}
zetal <- data.frame(zeta = unlist(zetal))
zetau <- data.frame(zeta = unlist(zetau))
zetal$interval <- "Lower"
zetau$interval <- "Upper"
zetal$cut <- rep(c("y1","y2","y3","y4","y5"), each = samp)
zetau$cut <- rep(c("y1","y2","y3","y4","y5"), each = samp)
zeta <- rbind(zetal, zetau)

p <- ggplot(zeta, aes(x = zeta)) + geom_density(aes(linetype = interval), 
	adjust = 2, show_guide = FALSE)
p <- p + facet_wrap(~ cut, ncol = 5) + xlim(c(-4,4)) + ylab("") + xlab(expression(zeta))
p <- p + geom_hline(yintercept = 0) + stat_ecdf(aes(linetype = interval))
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 15))
p <- p + theme(axis.text.x = element_text(size = 10, color = "black"))
p <- p + theme(axis.text.y = element_text(color = "black"))
p <- p + scale_linetype_manual(values = c(6,1), name = "Interval")

plot(p)

## ----fig5----------------------------------------------------------------
zeta <- vector("list", length = 6)
for (i in 1:6) {
	y <- patterns(5, 2, i-1)
	zeta[[i]] <- postsamp(fmodel3pl, y, apar = alph, bpar = beta, cpar = gamm,
		control = list(nbatch = samp + burn))
	zeta[[i]] <- zeta[[i]]$batch[(burn + 1):(samp + burn)]
	zeta[[i]] <- table(cut(zeta[[i]], c(-Inf, qnorm(c(0.25,0.5,0.75)), Inf)))
}
zeta <- data.frame(zeta = unlist(zeta))
zeta$quartile = factor(rep(1:4, 6))
levels(zeta$quartile) <- c("First","Second","Third","Fourth")
zeta$zeta <- zeta$zeta/samp * 100
zeta$sum <- factor(rep(0:5, each = 4))
zeta$text <- paste(round(zeta$zeta, 1), "%", sep = "")
zeta$side <- factor(ifelse(zeta$zeta > 50, 1, 0))

p <- ggplot(zeta, aes(x = sum, y = quartile, fill = zeta)) + geom_tile()
p <- p + scale_fill_gradient(name = "Probability (%)", low = grey(0.9), 
	high = grey(0.1), breaks = c(0, 50, 100), limits = c(0, 100))
p <- p + geom_text(aes(label = text, color = side)) 
p <- p + scale_color_manual(values = c("black","white"), guide = FALSE)
p <- p + theme_bw() + coord_fixed()
p <- p + theme(text = element_text(size = 10))
p <- p + theme(axis.text = element_text(color = "black"))
p <- p + ylab("") + xlab("Sum")

plot(p)

## ----fig6----------------------------------------------------------------

prob <- matrix(NA, 6, 6)
for (j in 1:6) {
	for (i in j:6) {
		if (i == j) {
			prob[i,j] <- 0.5
		}      
		else {
			zeta1 <- postsamp(fmodel3pl, patterns(5, 2, i-1), 
				apar = alph, bpar = beta, cpar = gamm,
				control = list(nbatch = samp + burn))
			zeta2 <- postsamp(fmodel3pl, patterns(5, 2, j-1), 
				apar = alph, bpar = beta, cpar = gamm,
				control = list(nbatch = samp + burn))
      zeta1 <- zeta1$batch[(burn + 1):(samp + burn)]
      zeta2 <- zeta2$batch[(burn + 1):(samp + burn)]
      prob[i,j] <- mean(zeta1 > zeta2)
      prob[j,i] <- 1 - prob[i,j]
		}
	}
}
prob <- data.frame(p = as.vector(prob), sumr = factor(rep(0:5, 6)), 
	sumc = factor(rep(0:5, each = 6)))
prob$p <- prob$p * 100
prob$text <- paste(round(prob$p, 1), "%", sep = "")
prob$side <- factor(ifelse(prob$p > 50, 1, 0))

p <- ggplot(prob, aes(x = sumc, y = sumr, fill = p)) + geom_tile()
p <- p + scale_fill_gradient(name = "Probability (%)", low = grey(0.9), 
	high = grey(0.1), breaks = c(0, 50, 100), limits = c(0, 100))
p <- p + geom_text(aes(label = text, color = side)) 
p <- p + scale_color_manual(values = c("black","white"), guide = FALSE)
p <- p + theme_bw() + coord_fixed()
p <- p + theme(text = element_text(size = 10))
p <- p + theme(axis.text = element_text(color = "black"))
p <- p + ylab("Examinee A Sum") + xlab("Examinee B Sum")

plot(p)

## ----fig7----------------------------------------------------------------
post <- vector("list", length = 6)
for (i in 1:6) {
	y <- patterns(5, 2, i-1)
	post[[i]] <- posttrace(fmodel3pl, y, apar = alph, bpar = beta, cpar = gamm, 
		zmin = -5, zmax = 5, length = 100)
	post[[i]] <- cbind(post[[i]]$zeta, post[[i]]$post)
}
post <- data.frame(do.call("rbind", post))
names(post) <- c("zeta","post")
post$sum <- factor(rep(0:5, each = 100))

post.norm <- post
post.norm$prior <- "normal"
	
post <- vector("list", length = 6)
for (i in 1:6) {
	y <- patterns(5, 2, i-1)
	post[[i]] <- posttrace(fmodel3pl, y, apar = alph, bpar = beta, cpar = gamm,
		zmin = -5, zmax = 5, prior = function(z) 1, length = 100)
	post[[i]] <- cbind(post[[i]]$zeta, post[[i]]$post)
}
post <- data.frame(do.call("rbind", post))
names(post) <- c("zeta","post")
post$sum <- factor(rep(0:5, each = 100))
	
post.unif <- post
post.unif$prior = "uniform"
	
post <- rbind(post.norm, post.unif)

names(post)[4] <- "Prior"

p <- ggplot(post, aes(x = zeta, y = post, linetype = Prior))
p <- p + geom_line() + facet_wrap(~ sum)
p <- p + xlab(expression(zeta)) + ylab("Log-Posterior Density")
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 15))
p <- p + theme(axis.text.x = element_text(size = 10, color = "black"))
p <- p + theme(axis.text.y = element_text(color = "black"))

plot(p)

## ----fig8----------------------------------------------------------------
post <- post[post$Prior == "uniform",]
post <- post[post$sum %in% c(2,3,4),]

y2 <- profileci(fmodel3pl, patterns(5, 2, 2), apar = alph, 
	bpar = beta, cpar = gamm, lower = FALSE)
y3 <- profileci(fmodel3pl, patterns(5, 2, 3), apar = alph, 
	bpar = beta, cpar = gamm)
y4 <- profileci(fmodel3pl, patterns(5, 2, 4), apar = alph, 
	bpar = beta, cpar = gamm)

names(post)[3] <- "Sum"

p <- ggplot(post, aes(x = zeta, y = post, linetype = Sum))
p <- p + geom_line() + ylim(c(-14,0))
p <- p + xlab(expression(zeta)) + ylab("Log-Posterior Density")
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 20))
p <- p + theme(axis.text.x = element_text(size = 15, color = "black"))
p <- p + theme(axis.text.y = element_text(color = "black"))
p <- p + annotate("segment", x = y2$zeta, xend = y2$zeta, 
	y = -14, yend = y2$post, linetype = 1)
p <- p + annotate("segment", x = y2$upper, xend = y2$upper, 
	y = -14, yend = y2$f.upper, linetype = 1)
p <- p + annotate("segment", x = y3$zeta, xend = y3$zeta, 
	y = -14, yend = y3$post, linetype = 3)
p <- p + annotate("segment", x = y3$lower, xend = y3$lower, 
	y = -14, yend = y3$f.lower, linetype = 3)
p <- p + annotate("segment", x = y3$upper, xend = y3$upper, 
	y = -14, yend = y3$f.upper, linetype = 3)
p <- p + annotate("segment", x = y4$zeta, xend = y4$zeta, 
	y = -14, yend = y4$post, linetype = 2)
p <- p + annotate("segment", x = y4$lower, xend = y4$lower, 
	y = -14, yend = y4$f.lower, linetype = 2)
p <- p + annotate("segment", x = y4$upper, xend = y4$upper, 
	y = -14, yend = y4$f.upper, linetype = 2)

plot(p)

## ----fig9----------------------------------------------------------------
zeta <- seq(-3, 3, length = 100)
info <- matrix(NA, 100, 5)
for (j in 1:100) {
	info[j,] <- information(fmodel3pl, c(0,1,1,1,1), zeta = zeta[j], 
		apar = alph, bpar = beta, cpar = gamm)$item
}

matplot(zeta, info, type = "l", ylab = "Information", bty = "n", xlab = expression(zeta))
legend(-3, 0.3, paste("Item", 1:5), lty = 1:5, col = 1:5)

## ----fig10---------------------------------------------------------------
fmodelhcm <- function(zeta, y, alph, beta, prior = dnorm, ...) {
	if (is.vector(y)) y <- matrix(y, 1, length(y))
	m <- ncol(y)
	n <- nrow(y)
	prob <- matrix(NA, m, 2)
  yprb <- matrix(NA, n, m)
	prob[,1] <- 2*cosh(zeta - beta)
	prob[,2] <- exp(alph)
	prob <- sweep(prob, 1, apply(prob, 1, sum), "/") 
	for (i in 1:n) {
		yprb[i,] <- prob[row(prob) == 1:m & col(prob) == y[i,] + 1]
	}
	return(list(post = log(sum(apply(yprb, 1, prod))) 
		+ log(prior(zeta, ...)), prob = prob))
}

dmixnorm <- function(z, pi, m1, m2, s1, s2) {
   return(pi*dnorm(z, m1, s1) + (1 - pi)*dnorm(z, m2, s2))
}

apar <- c(1,1,1,1,1)
bpar <- c(-2,-1,0,1,2)
n <- 10000

tmp.1 <- postsamp(fmodelhcm, patterns(5,2,1), alph = apar, beta = bpar, prior = dmixnorm, 
   m1 = -2, m2 = 2, s1 = 2, s2 = 2, pi = 0.5, control = list(nbatch = n))

tmp.2 <- postsamp(fmodelhcm, c(1,0,0,0,0), alph = apar, beta = bpar, prior = dmixnorm, 
   m1 = -2, m2 = 2, s1 = 2, s2 = 2, pi = 0.5, control = list(nbatch = n))

tmp.3 <- postsamp(fmodelhcm, c(1,1,1,0,0), alph = apar, beta = bpar, prior = dmixnorm, 
   m1 = -2, m2 = 2, s1 = 2, s2 = 2, pi = 0.5, control = list(nbatch = n))

tmp.4 <- postsamp(fmodelhcm, c(0,1,1,1,0), alph = apar, beta = bpar, prior = dmixnorm, 
   m1 = -2, m2 = 2, s1 = 2, s2 = 2, pi = 0.5, control = list(nbatch = n))
	
tmp <- data.frame(zeta = c(tmp.1$batch, tmp.2$batch, tmp.3$batch, tmp.4$batch), 
	Pattern = rep(letters[1:4], times = c(n,n,n,n*5)))	

p <- ggplot(tmp, aes(x = zeta, color = Pattern, linetype = Pattern)) 
p <- p + geom_density(adjust = 2)
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 20))
p <- p + theme(axis.text.x = element_text(size = 15, color = "black"))
p <- p + theme(axis.text.y = element_text(color = "black"))
p <- p + xlab(expression(zeta)) + ylab("Posterior Density")
plot(p)

## ----fig11---------------------------------------------------------------
zeta <- seq(-3, 3, length = 100)
info <- matrix(NA, 100, 5+1)
for (j in 1:100) {
	tmp <- information(fmodelhcm, c(0,0,0,0,0), zeta = zeta[j], alph = apar, beta = bpar)
	info[j,] <- c(tmp$test, tmp$item)
}

info <- data.frame(cbind(zeta, info))
names(info) <- c("zeta","test",paste("item", 1:5))
info <- melt(info, id.vars = "zeta", variable.name = "Type", value.name = "information")

p <- ggplot(info, aes(x = zeta, y = information, color = Type, linetype = Type))
p <- p + geom_line() + ylab("Information") + xlab(expression(zeta))
p <- p + theme_bw() + theme(legend.key = element_blank())
p <- p + theme(text = element_text(size = 15))
p <- p + theme(axis.text = element_text(size = 15)) 

plot(p)

