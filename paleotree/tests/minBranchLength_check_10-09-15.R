library(paleotree)

# tests to make sure results of minimum branch length algorithm are
	# not stochastic when no stochastic element exists in the input

# a small portion of the tree from Benson et al. 2014 (PLOS Biology)
	# courtesy of Manabu Sakamoto

tree<-read.tree(text="(((Spinosaurus,Angaturama,Irritator),(Suchomimus,Baryonyx)),(Eustreptospondylus,((Afrovenator,Magnosaurus,Dubreuillosaurus,Leshansaurus,Piveaeausaurus),(Duriavenator,(Megalosaurus,Torvosaurus)))));")
 
ages<-cbind(FAD=c(98.5, 106.75, 106.75, 119.65, 129.675, 163.75, 168.8,
		169.85, 167.2, 160.4, 163.75, 168.85, 167.2, 154.7),
	LAD=c(98.5, 106.75, 106.75, 119.65, 129.675, 163.75, 154.25,
		169.85, 167.2, 156, 163.75, 168.85, 167.2, 148.55))
rownames(ages)<-c('Spinosaurus','Angaturama','Irritator',
	'Suchomimus','Baryonyx','Eustreptospondylus','Afrovenator',
	'Magnosaurus','Dubreuillosaurus','Leshansaurus','Piveaeausaurus',
	'Duriavenator','Megalosaurus','Torvosaurus')

# test to make sure that time-scaling is not stochastic under MBL
	# based on code from Manabu Sakamoto

timetreesMBL <- list()
for(i in 1:100){
	timetreesMBL[[i]] <- timePaleoPhy(tree, ages,
		type="mbl", vartime = 1)
	}
# compare first tree with all others
testEqual<-all(sapply(timetreesMBL,all.equal,timetreesMBL[[1]]))
if(!testEqual){stop("Not all MBL trees are identical")}
