library("aroma.affymetrix")

## Empty CEL file
df <- AffymetrixCelFile()
print(df)

## Missing CEL file
df <- AffymetrixCelFile(NA_character_, mustExist=FALSE)
print(df)

