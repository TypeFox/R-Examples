### 
### testing routines for LSA package
### 

Sys.setenv(NOAWT=TRUE)

library(lsa)

lsatest <- function(test, description) {
    if (!test) stop(description)
}

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# testing cosine

a = c(2,1,1,1,0)
b = c(0,0,0,1,0)
lsatest( (round(cosine(a,b),3) == 0.378), "[cosine] - vector comparison");

m = cbind(a,b,a,a);
simm = round(cosine(m),2);
sims = c(1, 0.38, 1, 1, 0.38, 1, 0.38, 0.38, 1, 0.38, 1, 1, 1, 0.38, 1, 1);
dim(sims) = c(4,4);
lsatest( all(simm==sims), "[cosine] - matrix comparison");

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# routines for get/set/delTriples

myTextMatrix = matrix(2,2,3);
colnames(myTextMatrix) = c("c1","c2","c3");
rownames(myTextMatrix) = c("dog","cat");
environment(myTextMatrix) = new.env();
class(myTextMatrix) = "textmatrix";

setTriple(myTextMatrix, "c1", "has_category", 15)
setTriple(myTextMatrix, "c1", "has_category", 11)
setTriple(myTextMatrix, "c3", "has_category", 20)
setTriple(myTextMatrix, "c2", "has_category", 20)
setTriple(myTextMatrix, "c1", "has_category", 20)

lsatest( all( (getTriple(myTextMatrix)[[1]] == as.vector(c("1", "1", "3", "2", "1"))) == TRUE), "[triples] - getTriple(all) subjects");
lsatest( all( (getTriple(myTextMatrix)[[2]] == as.vector(rep("has_category", 5))) == TRUE), "[triples] - getTriple(all) predicates");
lsatest( all( (getTriple(myTextMatrix)[[3]] == as.vector(c("15","11","20","20","20"))) == TRUE), "[triples] - getTriple(all) objects");

lsatest( all( (getTriple(myTextMatrix, "c1")[[1]] == as.vector(rep("has_category",3))) == TRUE), "[triples] - getTriple(c1) predicates" )
lsatest( all( (getTriple(myTextMatrix, "c1")[[2]] == as.vector(c("15","11","20"))) == TRUE), "[triples] - getTriple(c1) objects")

lsatest( getTriple(myTextMatrix, "c2")[[1]][1] == "has_category", "[triples] - getTriple(c2) predicates")
lsatest( getTriple(myTextMatrix, "c2")[[2]][1] == "20", "[triples] - getTriple(c2) objects")

lsatest( all( (getTriple(myTextMatrix, "c1", "has_category") == c("15","11","20")) ), "[triples] - getTriple(c1,has_category) objects")

delTriple(myTextMatrix, "c1", "has_category", 11)

lsatest( all( (getTriple(myTextMatrix, "c1")[[1]] == as.vector(rep("has_category",2))) == TRUE), "[triples] - deletion, predicates")
lsatest( all( (getTriple(myTextMatrix, "c1")[[2]] == as.vector(c("15","20"))) == TRUE), "[triples] - deletion, objects")

setTriple(myTextMatrix, "c1", "has_category", 17)

lsatest( all( (getTriple(myTextMatrix)[[1]] == as.vector(c("1", "3", "2", "1", "1"))) == TRUE), "[triples] - insertion after deletion, subjects")
lsatest( all( (getTriple(myTextMatrix)[[2]] == as.vector(rep("has_category", 5))) == TRUE), "[triples] - insertion after deletion, predicates" )
lsatest( all( (getTriple(myTextMatrix)[[3]] == as.vector(c("15","20","20","20","17"))) == TRUE), "[triples] - insertion after deletion, objects")

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# routines for textmatrix

# create landauer example with files
td = tempfile()
dir.create(td)
write( c("human", "interface", "computer"), file=paste(td,"/c1", sep=""))
write( c("survey", "user", "computer", "system", "response", "time"), file=paste(td,"/c2", sep=""))
write( c("EPS", "user", "interface", "system"), file=paste(td,"/c3", sep=""))
write( c("system", "human", "system", "EPS"), file=paste(td,"/c4", sep=""))
write( c("user", "response", "time"), file=paste(td,"/c5", sep=""))
write( c("trees"), file=paste(td,"/m1", sep=""))
write( c("graph", "trees"), file=paste(td,"/m2", sep=""))
write( c("graph", "minors", "trees"), file=paste(td,"/m3", sep=""))
write( c("graph", "minors", "survey"), file=paste(td,"/m4", sep=""))

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# test normal matrix

dtm = textmatrix(td)
lsatest( all( rownames(dtm) == c("computer", "human",  "interface", "response", "survey", "system", "time", "user", "eps", "trees", "graph", "minors")), "[textmatrix] - landauer, terms")
lsatest( all( colnames(dtm) == c("c1","c2","c3","c4","c5","m1","m2","m3","m4") ), "[textmatrix] - landauer, docs")

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# test with reduced vocabulary (replaces former function pseudo_docs)

dtm2 = textmatrix(td, vocabulary = rownames(dtm)[-(3:7)])
lsatest( all( rownames(dtm2) == c("computer", "human",  "user", "eps", "trees", "graph", "minors")), "[textmatrix] - controlled vocabulary (terms)")
lsatest( all( colnames(dtm2) == c("c1","c2","c3","c4","c5","m1","m2","m3","m4") ), "[textmatrix] - controlled vocabulary (docs)")

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# test with stemming

dtm = textmatrix(td, stemming=TRUE, language="english")
lsatest( all( rownames(dtm) == c("comput", "human", "interfac", "respons", "survey", "system", "time", "user", "ep", "tree", "graph", "minor")), "[textmatrix] - stemming")

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# test with stopping

write( c("the", "das", "minor", "die", "it"), file=paste(td,"/stopwords", sep=""))
data(stopwords_en)
dtm = textmatrix(td, stopwords=stopwords_en, minDocFreq=1, minWordLength=1)
lsatest( all( rownames(dtm) == c("computer", "human", "interface", "response", "survey", "system", "time", "user", "eps", "trees", "graph", "minors", "das", "die", "minor")), "[textmatrix] - stopping english")
data(stopwords_de)
dtm2 = textmatrix(td, stopwords=stopwords_de, minDocFreq=1, minWordLength=1)
lsatest( all( rownames(dtm2) == c("computer", "human", "interface", "response", "survey", "system", "time", "user", "eps", "trees", "graph", "minors", "it", "minor", "the" )), "[textmatrix] - stopping german")

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# clean up
unlink(td, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# test query

lsatest( all( query("response.interface human", rownames(dtm2)) == c(0,1,1,1,0,0,0,0,0,0,0,0,0,0,0)), "[textmatrix] - query")

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# test for word order sensitivity

td1 = tempfile()
dir.create(td1)
write( c("word4", "word3", "word2"), file=paste(td1,"/c1", sep=""))
write( c("word1", "word4", "word2"), file=paste(td1,"/c2", sep=""))
td2 = tempfile()
dir.create(td2)
write( c("word1", "word2", "word3"), file=paste(td2,"/c1", sep=""))
write( c("word1", "word2", "word3"), file=paste(td2,"/c2", sep=""))
dtm1 = textmatrix(td1)
dtm2 = textmatrix(td2, vocabulary=rownames(dtm1))
dtm3 = textmatrix(td2)
lsatest( length(rownames(dtm1)) == length(rownames(dtm2)), "[textmatrix] - word order (vocabulary by number)")
lsatest( all(rownames(dtm1) == rownames(dtm2)), "[textmatrix] - word order (vocabulary by content)")
lsatest( length(rownames(dtm2)) != length(rownames(dtm3)), "[textmatrix] - word order (vocabulary test 2)")
unlink(td1, recursive=TRUE)
unlink(td2, recursive=TRUE)


# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# word length

td1 = tempfile()
write( c("123", "1234", "12345", "1234567"), file=td1)
dtm = textmatrix(td1, minWordLength=4, maxWordLength=5)
lsatest( length(rownames(dtm)) == 2, "[textmatrix] - word length")
unlink(td1, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# frequency cut-offs

td1 = tempfile()
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td1)
dtm = textmatrix(td1, minWordLength=0, maxWordLength=F, minDocFreq=3, maxDocFreq=4)
lsatest( length(rownames(dtm)) == 2, "[textmatrix] - document frequencies")
unlink(td1, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# input files: file, files, dir, dirs, files and dirs

td1 = tempfile()
dir.create(td1)

td2 = paste(td1,"A1.txt",sep="/")
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td2)

td3 = paste(td1,"A2.txt",sep="/")
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td3)

td4 = paste(td1,"B",sep="/")
dir.create(td4)
td5 = paste(td4, "B1.txt", sep="/")
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td5)
td6 = paste(td4, "B2.txt", sep="/")
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td6)

td7 = paste(td1,"C",sep="/")
dir.create(td7)
td8 = paste(td7, "C1.txt", sep="/")
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td8)
td9 = paste(td7, "C2.txt", sep="/")
write( c("A", "A", "A", "A", "B", "B", "B", "B", "B", "C", "C", "C", "D"), file=td9)

dtm1 = textmatrix(td2, minWordLength=0)
dtm2 = textmatrix(c(td2, td3), minWordLength=0)
dtm3 = textmatrix(td4, minWordLength=0)
dtm4 = textmatrix(c(td4,td7), minWordLength=0)
dtm5 = textmatrix(td1, minWordLength=0)

lsatest( length(colnames(dtm1)) == 1, "[textmatrix] - input files: file ")
lsatest( length(colnames(dtm2)) == 2, "[textmatrix] - input files: files ")
lsatest( length(colnames(dtm3)) == 2, "[textmatrix] - input files: dir")
lsatest( length(colnames(dtm4)) == 4, "[textmatrix] - input files: dirs ")
lsatest( length(colnames(dtm5)) == 6, "[textmatrix] - input files: files + dirs ")

unlink(td1, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# XML/HTML parsing

td1 = tempfile()
write( c("<html>\n<head></head><body background=\"#FFFFFF\">&auml;bc<h1>test></h1></body></html>"), file=td1)
lsatest( length(rownames(textmatrix(td1, removeXML=TRUE, language="german"))) == 2, "[textmatrix] - XML removal")
unlink(td1)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# remove numbers

td1 = tempfile()
write( c("abc def 123 324 def12 ab1 defg abc"), file=td1)
lsatest( length(rownames(textmatrix(td1, removeNumbers=TRUE))) == 5, "[textmatrix] - remove numbers")
unlink(td1)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# arabic buckwalter support

td1 = tempfile()
dir.create(td1)
write( "HalAwaY_1    jan~ap_1    muriyd_1   taHoDiyr_1      |l_2", file=paste(td1,"A1",sep="/") )
write( "HalAwap_2    jaraH-a_1   muro$id_1  taHoDiyriy~_1   |laY_1", file=paste(td1,"A2",sep="/") )
write( "HalAyib_2    jaraY-i_1   muroDiy_1  taHoSiyl_1      |lam_1 ", file=paste(td1,"A3",sep="/") )
#lcc = Sys.getlocale("LC_ALL")
#Sys.setlocale("LC_ALL", "C")
lsatest("tahodiyriy~_1" %in% rownames(textmatrix(td1, language="arabic")),
        "[textmatrix] - arabic")
unlink(td1, recursive=TRUE)
#Sys.setlocale("LC_ALL", lcc)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# global frequency boundaries

td1 = tempfile()
dir.create(td1)
write( "hund hund katze maus", file=paste(td1,"A1",sep="/") )
write( "hund katze birne apfel banana birne", file=paste(td1,"A2",sep="/") )
write( "hund birne", file=paste(td1,"A3",sep="/") )
lsatest( length(rownames(textmatrix(td1, minGlobFreq=2)))==3, "[textmatrix] - global frequency (minimum)" )
lsatest( length(rownames(textmatrix(td1, maxGlobFreq=2, minGlobFreq=2)))==2, "[textmatrix] - global frequency (minimum)" )
unlink(td1, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# stemming in german

td1 = tempfile()
dir.create(td1)
write( "systeme", file=paste(td1,"A1",sep="/") )
write( "hunde", file=paste(td1,"A2",sep="/") )
td2 = tempfile()
dir.create(td2)
write( "von systemen von hunden", file=paste(td2,"A3",sep="/") )
tm1 = textmatrix(td1, stemming=T, language="german")
tm2 = textmatrix(td2, stemming=T, language="german", vocabulary=rownames(tm1) )
lsatest( all( (rownames(tm1) == c("system", "hund")) == TRUE), "[textmatrix] - stemming (german)" )
lsatest( all( (rownames(tm2) == c("system", "hund")) == TRUE), "[textmatrix] - fold-in with stemming (german)" )
unlink(td2, recursive=TRUE)
unlink(td1, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# one term matrices

td1 = tempfile()
dir.create(td1)
write( "systeme", file=paste(td1,"A1",sep="/") )
tm1 = textmatrix(td1, stemming=T, language="german")
lsatest( (rownames(tm1) == "system"), "[textmatrix] - one term matrix" )
unlink(td1, recursive=TRUE)

# -  -  -  -  -  -  -  -  -  -  -  -  -  -  
# polish language support

td1=tempfile()
load("polski.RData")
writeLines(polski, con=td1, useBytes=TRUE)
polski2 = textvector(td1, language="polish")
polski2['terms']
Encoding(as.character(polski2[,"terms"]))
polvec = "test\xc4\x85\xc4\x85\xc4\x99\xc4\x99\xc3\xb3\xc3\xb3\xc4\x87\xc4\x87\xc5\x82\xc5\x82\xc5\x84\xc5\x84\xc5\x9b\xc5\x9b\xc5\xba\xc5\xba\xc5\xbc\xc5\xbc"
# enc2utf
Encoding(polvec) = "UTF-8"
polvec
Encoding(polvec)
lsatest( polski2['terms'] == polvec, "[textvector] - polish language support" )
unlink(td1)
