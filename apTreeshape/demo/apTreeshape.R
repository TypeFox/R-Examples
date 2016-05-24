if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))


## universal tree - test

data(universal.treeshape)
tree <- universal.treeshape
plot(tree)
summary(tree)

likelihood.test(tree,  model = "yule", alternative = "two.sided")
likelihood.test(tree,  model = "pda", alternative = "two.sided")
colless.test(tree, model = "yule", alternative = "greater")

# aldous test

#aldous.test(tree)
aldous.test(rtreeshape(n=1, tip.number=2000, model="yule"))


#resolution de polytomies par simulation

phy=dbtrees(db="treebase", tree=741, class="phylo")
plot(phy)
tree=dbtrees(db="treebase", tree=741, class="treeshape")
tree=dbtrees(db="treebase", tree=741, class="treeshape", model="yule")
plot(tree)

 
### familles d'oiseaux

data(bird.families)
tree.b <- bird.families
tree.birds <- as.treeshape(tree.b, model ="yule")
plot(tree.birds)
print(tree.birds$names)
class(tree.birds) <- "treeshape"
likelihood.test(tree.birds, alternative="two.sided")
likelihood.test(tree.birds, model ="pda", alternative="two.sided")

### test du modele d'Aldous
## Exemple

plot(rtreeshape(n=1, tip.number=40, model="aldous"))

## simulation

trees.ab <- rtreeshape(n=500, tip.number=137, model="aldous")

shape.lst <- sapply(trees.ab, FUN = shape.statistic)
quantile(shape.lst, prob = c(0.025, 0.975))
shape.statistic(tree.birds)

## scan pandit

trees <- dbtrees(db="pandit", tree=200:400, class="treeshape")
#on ne garde que les arbre qui ont plus de 20 feuilles

tr <- trees[sapply(trees, FUN = function(x){length(x$merge)>38} )]

#on regarde les tailles
sapply(tr, FUN = function(x){length(x$merge)/2} )


# Pda N(0,1) ?

hist( sapply( tr, FUN = shape.statistic, norm="pda"))
qqnorm(sapply( tr, FUN = shape.statistic, norm="pda"))

# un arbre particulier et des sous arbres de taille 25
tree <- trees[[85]]
tr.25 <- NULL; for (i in 1:200) tr.25[[i]] <- tipsubtree(tree, sample(1:256, 25), numeric = TRUE )
hist( sapply( tr.25, FUN = shape.statistic, norm="pda"))

# simulation
summary(tree)
tr.s <- rtreeshape(500, 256, FUN = function(n,i){if ((i>0)&(i<n)) 1/i/(n-i) else 0})
quantile( sapply(tr.s, FUN = colless), prob = c(0.025, 0.5, 0.975))

hist( sapply( tr.25, FUN = shape.statistic, norm="pda"))
tr.25s <- rtreeshape(200, 25, FUN = function(n,i){if ((i>0)&(i<n)) 1/i/(n-i) else 0})
x11(); hist( sapply( tr.25s, FUN = shape.statistic, norm="pda"))


#on ne garde que les arbre qui ont 50-100 feuilles
tr <- trees[sapply(trees, FUN = function(x){(length(x$merge)>97)& (length(x$merge)<199)} )]
 hist( sapply( tr, FUN = shape.statistic, norm="pda"))
