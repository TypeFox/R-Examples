# match.hash is an infomal (S3) class representing the
# chain of hash tables stored in the .match.hash attribute
# of tables that have been hashed

# we provide a (sort of dummy) print method so
# the output is not as ugly
print.match.hash <- function(x, ...) { cat("<hash table>\n"); x }
