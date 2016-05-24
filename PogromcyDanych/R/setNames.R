pushGlobal <- function(name, value) {
  tf <- tempfile()
  assign(name, value = value)
  save(list=name, file=tf)
  load(tf, envir = .GlobalEnv)
  unlink(tf)
}

setLang <- function(lang = "eng") {
  auta2012_eng <- PogromcyDanych::auta2012
  pearson_eng <- PogromcyDanych::pearson
  galton_eng <- PogromcyDanych::galton
  WIG_eng <- PogromcyDanych::WIG
  TCGA_BRCA_eng <- PogromcyDanych::TCGA_BRCA
  diagnosis <- PogromcyDanych::diagnoza
  diagnosisDict <- PogromcyDanych::diagnozaDict
  mandatySejmik2014_eng <- PogromcyDanych::mandatySejmik2014
  imiona_warszawa_eng <- PogromcyDanych::imiona_warszawa
  seriale_eng <- PogromcyDanych::serialeIMDB
  cats_birds <- PogromcyDanych::koty_ptaki
  if (lang == "eng") {
    colnames(auta2012_eng) <- c("Price", "Currency", "Price.in.PLN", "Gross.Net", "HP", "kW",
                                "Brand", "Model", "Version", "Nubmer.of.doors", "Engine.cubic.capacity",
                                "Mileage", "Type.of.fuel", "Year", "Color", "Country.of.current.registration",
                                "Country.of.origin", "Is.damaged", "Transmission", "Is.imported",
                                "Accessories")
  
    colnames(TCGA_BRCA_eng) <- c("TP53", "gender", "vital.status", "days.to.death", "new.tumor")
    colnames(seriale_eng) <- c("id", "series", "name", "season","part","note","votes","imdbId")
    colnames(pearson_eng) <- c("son", "father")
    colnames(galton_eng) <- c("son", "mid_parent")
    colnames(imiona_warszawa_eng) <- c("name", "sex", "year", "month", "count")
    colnames(mandatySejmik2014_eng) <- c("Voivodeship", "PSL", "PiS", "PO", "SLD", "Other", "Prc_valid_votes", 
                               "long", "lat")
    
    colnames(WIG_eng) <- c("Date", "Name", "Opening Price", "Max Price", "Min Price",
                           "Closing Price", "Change", "Turnover")
    
    colnames(cats_birds) <- c("species", "weight", "length", "speed", "habitat", "lifespan", "group" )
    cats_birds$species <- c("Tiger", "Lion", "Cheetah", "Jaguar", "Puma", "Leopard", "Irbis", "Swift", 
                            "Ostrich", "Golden Eagle", "Peregrine Falcon", "Falcon Norwegian", "Albatros")
    cats_birds$group <- c(rep("Cat", 7), rep("Bird", 6))
    cats_birds$habitat <- c("Asia", "Africa", "America", "America", "Asia", "Africa", "Asia", "Eurasia", "Africa", "North", "North", "North", "South")

    levels(diagnosis$plec) = c("MAN", "WOMAN")
    for (i in 1:ncol(diagnosis)) {
      if(class(diagnosis[,i]) == "factor")
        diagnosis[,i] <- droplevels(diagnosis[,i])
    }
    levels(diagnosis$eduk4_2013) = c("PRIMARY/NO EDUCATION", "VOCATIONAL/GRAMMAR", 
                   "SECONDARY", "HIGHER AND POST-SECONDARY")
    levels(diagnosis$status9_2013) = c("EMPLOYEES IN PUBLIC SECTOR", 
                     "EMPLOYEES IN PRIVATE SECTOR", "ENTREPRENEUR/SELF-EMPLOYED", 
                     "FARMERS", "PENSIONERS", "RETIREES", "PUPILS AND STUDENTS", 
                     "UNEMPLOYED", "OTHER PROFESSIONALLY INACTIVE")
    levels(diagnosis$gp3) = c("DELIGHTED", "PLEASED", 
            "MOSTLY SATISFIED", "MIXED", "MOSTLY DISSATISFIED", "UNHAPPY", 
            "TERRIBLE")
    levels(diagnosis$gp29) = c("FUN, WELL-BEING, LACK OF STRESS", 
             "SENSE OF PURPOSE, ACHIEVING IMPORTANT GOALS DESPITE DIFFICUL")
    levels(diagnosis$gp54_01) = c("DEFINITELY AGREE", "AGREE", "RATHER AGREE", 
                "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", "DISAGREE", 
                "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_02) = c("DEFINITELY YES", "YES", 
                "RATHER YES", "NEITHER YES OR NO", "PROBABLY NOT", "NO", 
                "DEFINITELY NOT")
    levels(diagnosis$gp54_03) = c("DEFINITELY YES", "YES", "RATHER YES", 
                "NEITHER YES OR NO", "PROBABLY NOT", "NO", "DEFINITELY NOT")
    levels(diagnosis$gp54_04) = c("DEFINITELY AGREE", "AGREE", "RATHER AGREE", 
                "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", "DISAGREE", 
                "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_05) = c("DEFINITELY AGREE", "AGREE", 
                "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_06) = c("DEFINITELY AGREE", 
                "AGREE", "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_07) = c("ZDECYDOWANIE TAK", 
                "TAK", "RACZEJ TAK", "ANI TAK, ANI NIE", "RACZEJ NIE", "NIE", 
                "ZDECYDOWANIE NIE") 
    levels(diagnosis$gp54_08) = c("DEFINITELY YES", "YES", 
                "RATHER YES", "NEITHER YES OR NO", "PROBABLY NOT", "NO", 
                "DEFINITELY NOT")
    levels(diagnosis$gp54_09) = c("DEFINITELY YES", "YES", "RATHER YES", 
                "NEITHER YES OR NO", "PROBABLY NOT", "NO", "DEFINITELY NOT")
    levels(diagnosis$gp54_10) = c("DEFINITELY YES", "YES", "RATHER YES", "NEITHER YES OR NO", 
                "PROBABLY NOT", "NO", "DEFINITELY NOT") 
    levels(diagnosis$gp54_11) = c("DEFINITELY YES", 
                "YES", "RATHER YES", "NEITHER YES OR NO", "PROBABLY NOT", 
                "NO", "DEFINITELY NOT") 
    levels(diagnosis$gp54_12) = c("DEFINITELY YES", "YES", 
                "RATHER YES", "NEITHER YES OR NO", "PROBABLY NOT", "NO", 
                "DEFINITELY NOT")
    levels(diagnosis$gp54_13) = c("DEFINITELY AGREE", "AGREE", 
                "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE")
    levels(diagnosis$gp54_14) = c("DEFINITELY AGREE", 
                "AGREE", "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_15) = c("DEFINITELY AGREE", 
                "AGREE", "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_16) = c("DEFINITELY YES", 
                "YES", "RATHER YES", "NEITHER YES OR NO", "PROBABLY NOT", 
                "NO", "DEFINITELY NOT")
    levels(diagnosis$gp54_17) = c("DEFINITELY AGREE", 
                "AGREE", "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_18) = c("DEFINITELY AGREE", 
                "AGREE", "RATHER AGREE", "NEITHER AGREE NOR DISAGREE", "RATHER DISAGREE", 
                "DISAGREE", "DEFINITELY DISAGREE") 
    levels(diagnosis$gp54_19) = c("DEFINITELY YES", 
                "YES", "RATHER YES", "NEITHER YES OR NO", "PROBABLY NOT", 
                "NO", "DEFINITELY NOT") 
    levels(diagnosis$gp54_20) = c("DEFINITELY YES", "YES", 
                "RATHER YES", "NEITHER YES NOR NOT", "RATHER NOT", "NO", 
                "DEFINITELY NOT")
    levels(diagnosis$gp54_21) = c("DEFINITELY YES", "YES", "RATHER YES", 
                "NEITHER YES NOR NOT", "RATHER NOT", "NO", "DEFINITELY NOT"
    )
    levels(diagnosis$gp54_22) = c("DEFINITELY YES", "YES", "RATHER YES", "NEITHER YES NOR NOT", 
                "RATHER NOT", "NO", "DEFINITELY NOT")
    
    levels(auta2012$Skrzynia.biegow) = c("", "automatic", "manual")
    levels(auta2012$Pojazd.uszkodzony) = c("", "Yes")
    levels(auta2012$Rodzaj.paliwa) = c("petrol", "petrol+LPG", "ethanol", "hybrid", "electric", "diesel")
    levels(auta2012$Kolor) <- c("", "sand", "sand-metallic", "white", "white-metallic", 
                                "dark red", "dark red-metallic", "brown", "brown-metallic", 
                                "black", "black-metallic", "red", "red-metallic", 
                                "violet", "violet-metallic", "graphite", "graphite-metallic", 
                                "dark blue", "dark blue-metallic", "blue", "blue-metallic", 
                                "orange", "orange-metallic", "pink", "pink-metallic", 
                                "silver", "silver-metallic", "grey", "grey-metallic", "cherry", 
                                "cherry-metallic", "green", "green-metallic", "yellow", 
                                "yellow-metallic", "gold", "gold-metallic")
    auta2012$Kolor <- factor(as.character(auta2012$Kolor))
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="el. lusterka", replacement="electric mirrors")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="klimatyzacja", replacement="air conditioning")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="alufelgi", replacement="alloy wheels")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="centralny zamek", replacement="central locking")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="poduszka powietrzna", replacement="airbag")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="wspomaganie kierownicy", replacement="power steering")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="komputer", replacement="computer")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="przyciemniane szyby", replacement="tinted windows")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="skorzana tapicerka", replacement="leather upholstery")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="tempomat", replacement="cruise control")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="hak", replacement="hook")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="el. szyby", replacement="el. windows")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="welurowa tapicerka", replacement="velor upholstery")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="lwiatla przeciwmglowe", replacement="fog lights
                                           ")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="kierownica wielofunkcyjna", replacement="multifunction steering wheel")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="pod. przednia szyba", replacement="the windshield")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="podgrzewane fotele", replacement="heated seats")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="czujnik parkowania", replacement="parking sensor")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="czujnik deszczu", replacement="rain sensor")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="system nawigacji", replacement="navigation system")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="ksenony", replacement="xeons")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="szyberdach", replacement="sunroof")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="niezalezne ogrzewanie", replacement="independent heating")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="bagaznik na dach", replacement="trunk on the roof")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="blokada skrzyni biegAlw", replacement="gearbox lock")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="reg. wysokole podwozia", replacement="height adjustable chassis")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="blokada dyferencjalu", replacement="differential lock")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="blokada skrzyni biegow", replacement="transmission lock box")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="instalacja gazowa", replacement="gas-fittings")
    auta2012$Wyposazenie.dodatkowe <- gsub(auta2012$Wyposazenie.dodatkowe, pattern="klatka", replacement="cage")
    
    pushGlobal("seriesIMDB", value = seriale_eng)
    pushGlobal("warsaw_names", value = imiona_warszawa_eng)
    pushGlobal("votes2014", value = mandatySejmik2014_eng)
    pushGlobal("cats_birds", value = cats_birds)
    pushGlobal("diagnosis", value = diagnosis)
    pushGlobal("diagnosisDict", value = diagnosisDict)
  } else {
    colnames(auta2012_eng) <- c("Cena", "Waluta", "Cena.w.PLN", "Brutto.netto", "KM", "kW", 
      "Marka", "Model", "Wersja", "Liczba.drzwi", "Pojemnosc.skokowa", 
      "Przebieg.w.km", "Rodzaj.paliwa", "Rok.produkcji", "Kolor", "Kraj.aktualnej.rejestracji", 
      "Kraj.pochodzenia", "Pojazd.uszkodzony", "Skrzynia.biegow", "Status.pojazdu.sprowadzonego", 
      "Wyposazenie.dodatkowe")
    
    colnames(TCGA_BRCA_eng) <- c("TP53", "plec", "czy.zyje", "dni.do.smierci", "czy.nowy.guz")
    
    colnames(pearson_eng) <- c("syn", "ojciec")
    colnames(galton_eng) <- c("syn", "sr_rodzic")
 #   colnames(imiona_warszawa_eng) <- c("imie", "plec", "rok", "miesiac", "liczba")  
    colnames(WIG_eng) <- c("Data", "Nazwa", "Kurs.otwarcia", "Kurs.maksymalny",
                           "Kurs.minimalny", "Kurs.zamkniecia", "Zmiana", "Wartosc.obrotu.w.tys.zl")

#    colnames(mandatySejmik2014_eng) <- c("Wojewodztwo", "PSL", "PiS", "PO", "SLD", "Inne", "ProcentWaznychGlosow", 
#                                         "long", "lat")   
#    colnames(cats_birds) <- c("gatunek", "waga", "dlugosc", "predkosc", "habitat", "zywotnosc", "druzyna" )
#    pushGlobal("koty_ptaki", value = cats_birds)
  }
  
  pushGlobal("pearson", value = pearson_eng)
  pushGlobal("galton", value = galton_eng)
  pushGlobal("TCGA_BRCA", value = TCGA_BRCA_eng)
  pushGlobal("WIG", value = WIG_eng)
  pushGlobal("auta2012", value = auta2012_eng)
  invisible(0)
}
