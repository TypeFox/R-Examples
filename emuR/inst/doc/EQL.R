## ----results='hide', message=FALSE, warning=FALSE------------------------
# load the package
library(emuR)

# create demo data in folder provided by the tempdir() function
create_emuRdemoData(dir = tempdir())

# get the path to emuDB called 'ae' that is part of the demo data
path2folder = file.path(tempdir(), "emuR_demoData", "ae_emuDB")

# load emuDB into current R session
ae = load_emuDB(path2folder)

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phonetic == m]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phonetic == m | n]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phonetic != m | n]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable =~ .*]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ a.*]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text !~ a.*]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
# NOTE: all row entries in the resulting segment list have the start time of "@", the end time of "n" and their labels will be "@->n"
query(ae, "[Phonetic == @ -> Phonetic == n]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
# NOTE: all row entries in the resulting segment list have the start time of "@", the end time of "@" and their labels will also be "@"
query(ae, "[#Phonetic == @ -> Phonetic == n]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
# NOTE: all row entries in the resulting segment list have the start time of "n", the end time of "n" and their labels will also be "n"
query(ae, "[Phonetic == @ -> #Phonetic == n]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Phonetic == @ -> Phonetic == n ] -> Phonetic =~ s]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Text == to -> Text == offer ] -> Text == any]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[[Text =~ offer -> Text =~ .*] -> Text =~ .* ] -> Text == resistance]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text == always & Word == C]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ .* & Word == F]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ .* & Word == C & Accent == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == p ^ Syllable == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable =~ .* ^ Phoneme == p]")
# or
query(ae, "[Phoneme == p ^ #Syllable =~ .*]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable =~ .* ^ Phoneme != p | t | k]")
# or
query(ae, "[Phoneme != p | t | k ^ #Syllable =~ .*]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[#Phonetic =~ .* ^ Syllable == S] ^ Text == amongst | beautiful]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Phonetic =~ .* ^ Syllable == S] ^ #Text == amongst | beautiful]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Start(Word, Syllable) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Start(Word, Phoneme) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Start(Word, Syllable) == 0]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[End(Word, Syllable) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Medial(Word, Syllable) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == n & Start(Syllable, Phoneme) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == m & End(Word, Phoneme) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable == S & End(Word, Syllable) == 0]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == p ^ Start(Word, Syllable) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme =~ .* ^ End(Word, Syllable) == 0]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Num(Word, Syllable) == 4]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Num(Syllable, Phoneme) > 6]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ .* & Num(Text, Phoneme) > 5 ]")
# or
query(ae, "[Text =~ .* & Num(Word, Phoneme) > 5]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable == S & Num(Syllable, Phoneme) == 5]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == m ^ Num(Word, Syllable) == 3]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable = W ^ Num(Word, Syllable) < 3]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ .* ^ Num(Syllable, Phoneme) == 4]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Phoneme == m -> Phoneme =~ p] ^ Syllable == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == s -> [Phoneme == t ^ Syllable == W]]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[#Syllable == S ^ Phoneme == s] -> Syllable == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Syllable == S ^ #Phoneme == s] -> Syllable == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ .* ^ Phoneme == @ & Start(Text, Phoneme) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Phoneme == m & Start(Word, Phoneme) == 1 -> Phoneme == o:] ^ Syllable == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[[Phoneme == m & Start(Word, Phoneme) == 1 -> Phoneme == o:] ^ Syllable == S] ^ #Text =~ .*]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Text =~ .* & Num(Text, Syllable) == 3 ^ [Phoneme == @ ^ Start(Word, Syllable) == 1]] -> Text == his]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phoneme == m | n & Medial(Word, Phoneme) == 1]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[Phonetic == H -> Phonetic =~ .*] -> Phonetic == I | U ]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable =~ .* & Medial(Word, Syllable) == 0]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text =~ .* & Num(Text, Syllable) == 2]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Text == the -> #Text =~ .* & Accent == S]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable = S ^ Num(Word, Phoneme) == 5]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable == W ^ Phoneme == @]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Phonetic =~ .* ^ #Syllable == W]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[Syllable == W & End(Word, Syllable) == 1 ^ Num(Word, Syllable) == 3]")

## ----results='hide', message=FALSE, warning=FALSE------------------------
query(ae, "[[[Phoneme =~ .* ^ Phonetic == H] ^ Start(Word, Syllable) == 1] ^ Accent == S]")

## ----echo=FALSE, results='hide', message=FALSE, warning=FALSE------------
# disconnect to avoid file locking to sqliteDB that causes unlink
# to fail under windows
DBI::dbDisconnect(ae$connection)

## ----results='hide', message=FALSE, warning=FALSE------------------------
# remove emuR_demoData as we will not be needing it 
# throughout the rest of this vignette
unlink(file.path(tempdir(), "emuR_demoData"), recursive = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  # query all "p" ITEMs on the "Phoneme" level that are dominated by "S" (strong) syllables
#  query(emuDBhandle = andosl,
#        query = "[Phoneme == p ^ Syllable == S]")
#  
#  # same query as before but this time using
#  # the sessionPattern and bundlePattern arguments
#  # to only select specific sessions / bundles
#  # using regular expressions (RegEx)
#  query(emuDBhandle = andosl,
#        query = "[Phoneme == p ^ Syllable == S]",
#        sessionPattern = "000.", # RegEx that matches session 0000; 000a; 0001; ...
#        bundlePattern = "msajc0[1-2].") # RegEx that matches bundles msajc01a; msajc02a; msajc021; ...

## ----eval=FALSE----------------------------------------------------------
#  emu.query(template = "andosl",
#            pattern = "*",
#            query = "[Text=spring & #Accent=S]")

## ----eval=FALSE----------------------------------------------------------
#  emu.query(template = "andosl",
#            pattern = "*",
#            query = "[#Text=spring & #Accent=S]")

## ----eval=FALSE----------------------------------------------------------
#  > query(emuDBhandle = andosl,
#          query = "[Text == spring & #Accent == S]",
#          resultType == "emusegs")

## ----eval=FALSE----------------------------------------------------------
#  > query(dbName = "andosl",
#          query = "[#Text == spring & #Accent == S]")

## ----eval=FALSE----------------------------------------------------------
#  emu.query(template = "ae",
#            pattern = "*",
#            query = "[Text!=beautiful|futile ^ Phoneme=u:]")

## ----eval=FALSE----------------------------------------------------------
#  query(dbName = "ae",
#        query = "[Text != beautiful | futile ^ Phoneme == u:]",
#        resultType = "emusegs")

