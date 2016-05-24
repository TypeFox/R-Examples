pkgname <- "vwr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('vwr')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ald")
### * ald

flush(stderr()); flush(stdout())

### Name: ald
### Title: Compute average Levenshtein distances
### Aliases: ald old20

### ** Examples

data(basque.words)
ald(basque.words[1:10],basque.words,20)
old20(basque.words[1:10],basque.words)



cleanEx()
nameEx("basque.words")
### * basque.words

flush(stderr()); flush(stdout())

### Name: basque.words
### Title: A list of Basque Words
### Aliases: basque.words

### ** Examples

data(basque.words)



cleanEx()
nameEx("coltheart.N")
### * coltheart.N

flush(stderr()); flush(stdout())

### Name: coltheart.N
### Title: Compute Coltheart's N
### Aliases: coltheart.N

### ** Examples

data(spanish.words)
sample.words<-sample(spanish.words,20)
coltheart.N(sample.words,spanish.words)
coltheart.N(sample.words,spanish.words, method='Levenshtein')



cleanEx()
nameEx("dutch.words")
### * dutch.words

flush(stderr()); flush(stdout())

### Name: dutch.words
### Title: A list of Dutch Words
### Aliases: dutch.words

### ** Examples

data(dutch.words)



cleanEx()
nameEx("english.words")
### * english.words

flush(stderr()); flush(stdout())

### Name: english.words
### Title: A list of English Words
### Aliases: english.words

### ** Examples

data(english.words)



cleanEx()
nameEx("french.words")
### * french.words

flush(stderr()); flush(stdout())

### Name: french.words
### Title: A list of French Words
### Aliases: french.words

### ** Examples

data(french.words)



cleanEx()
nameEx("german.words")
### * german.words

flush(stderr()); flush(stdout())

### Name: german.words
### Title: A list of German Words
### Aliases: german.words

### ** Examples

data(german.words)



cleanEx()
nameEx("hamming.distance")
### * hamming.distance

flush(stderr()); flush(stdout())

### Name: hamming.distance
### Title: Compute Hamming distances
### Aliases: hamming.distance

### ** Examples

data(english.words)
hamming.distance('electroencephalogram',english.words)



cleanEx()
nameEx("hamming.neighbors")
### * hamming.neighbors

flush(stderr()); flush(stdout())

### Name: hamming.neighbors
### Title: Compute Hamming neighbors
### Aliases: hamming.neighbors

### ** Examples

data(english.words)
hamming.neighbors('electroencephalogram',english.words)
hamming.neighbors('hello',english.words)



cleanEx()
nameEx("ldknn")
### * ldknn

flush(stderr()); flush(stdout())

### Name: ldknn
### Title: Run the ldknn algorithm
### Aliases: ldknn ld1nn plot.ld1nn.run print.ld1nn.run

### ** Examples

data(english.words)
data(basque.words)
# set up a mock experiment: English stimuli are words, Basque stimuli are nonwords
experiment<-data.frame(stimulus=c(sample(english.words,500),
 sample(basque.words,500)),
 type=factor(rep(c('Word','Nonword'),each=500),levels=c('Word','Nonword')))
# randomize the trials
experiment<-experiment[sample(1:1000,1000),]
# run the ldknn algorithm
results<-ldknn(experiment$stimulus,experiment$type,'Word')
print(results)
plot(results)



cleanEx()
nameEx("levenshtein.distance")
### * levenshtein.distance

flush(stderr()); flush(stdout())

### Name: levenshtein.distance
### Title: Compute Levenshtein distances
### Aliases: levenshtein.distance

### ** Examples

data(french.words)
levenshtein.distance('pourquoi',sample(french.words,20))



cleanEx()
nameEx("levenshtein.neighbors")
### * levenshtein.neighbors

flush(stderr()); flush(stdout())

### Name: levenshtein.neighbors
### Title: Compute Levenshtein neighbors
### Aliases: levenshtein.neighbors

### ** Examples

data(serbian_latin.words)
levenshtein.neighbors('pola',serbian_latin.words)[1:2]



cleanEx()
nameEx("serbian_cyrillic.words")
### * serbian_cyrillic.words

flush(stderr()); flush(stdout())

### Name: serbian_cyrillic.words
### Title: A list of Serbian Words in Cyrillic alphabet
### Aliases: serbian_cyrillic.words

### ** Examples

data(serbian_cyrillic.words)



cleanEx()
nameEx("serbian_latin.words")
### * serbian_latin.words

flush(stderr()); flush(stdout())

### Name: serbian_latin.words
### Title: A list of Serbian Words in Latin alphabet
### Aliases: serbian_latin.words

### ** Examples

data(serbian_latin.words)



cleanEx()
nameEx("spanish.words")
### * spanish.words

flush(stderr()); flush(stdout())

### Name: spanish.words
### Title: A list of Spanish Words
### Aliases: spanish.words

### ** Examples

data(spanish.words)



cleanEx()
nameEx("vietnamese.words")
### * vietnamese.words

flush(stderr()); flush(stdout())

### Name: vietnamese.words
### Title: A list of Vietnamese Words
### Aliases: vietnamese.words

### ** Examples

data(vietnamese.words)



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
