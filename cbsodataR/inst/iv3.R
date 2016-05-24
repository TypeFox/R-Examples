
# installeer cbsodataR package als dat nog niet gebeurd is
if (!requireNamespace("cbsodataR")){
  install.packages( "cbsodataR", repos=c("http://cran.rstudio.com", 
                    'http://edwindj.github.io/drat'))
}

library(cbsodataR)
IV3 <- "http://dataderden.cbs.nl"

tl <- get_table_list(base_url = IV3)
# tl bevat nu een lijst met tabellen uit de IV3

meta <- get_meta("45001NED", base_url = IV3)
# meta bevat metadata informatie van de tabellen

# Zoek in meta$gemeente de sleutel van je gemeente op bijv. Alkmaar. 
# NB: let op de spaties: sleutels moeten exact overgenomen worden!

# Voor kleine tabellen kun je get_data gebruiken.
data <- get_data( "45001NED"              # tabel nummer (uit 2012)
                , Gemeenten = "GM0361   "       # filter op Gemeente Alkmaar
                , Verslagsoort = "2012X005"     # pak het jaarverslag
                , base_url = IV3                # haal deze tabel uit IV3
                )

# gebruik anders download_table
download_table("45001NED"                     # tabel nummer (uit 2012)
              , Gemeenten = "GM0361   "       # filter op Gemeente Alkmaar
              , Verslagsoort = "2012X005"     # pak het jaarverslag
              , base_url = IV3                # haal deze tabel uit IV3
              )