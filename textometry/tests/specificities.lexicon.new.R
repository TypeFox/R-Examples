library(textometry);
#               A  D  G  J  B  C
sublexicon <- c(1, 2, 1, 3, 2, 1);
names(sublexicon) <- c("A", "D", "G", "J", "B", "C");
#            A  B  C  D  E  F  G  H  I  J  K  L  M
lexicon <- c(1, 2, 1, 3, 2, 1, 2, 2, 2, 4, 1, 3, 4);
names(lexicon) <- LETTERS[1:length(lexicon)];

specificities.lexicon.new(lexicon, sublexicon);
