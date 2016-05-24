context ("Remote API Client: avisApiBusAvanzada")

testSpeciesId <- 480
expectedSpecies <- "Pica pica"

response <- .avisApiBusAvanzada(list(id_especie = testSpeciesId))

test_that(".avisApiBusAvanzada returns the expected header",{ 
    
    expectedNames <- c(
    	"Id..Obs.",
		"Nombre.comun",
		"Especie",
		"Fecha",
		"Numero",
		"Com.Aut",
		"Provincia",
		"UTM",
		"Observador",
		"Periodo",
		"Hora",
		"Edad",
		"Sexo",
		"Interes",
		"Grado.Reprod.",
		"Categ.Fenol.",
		"Habitat",
		"Codigo.de.Habitat",
		"Notas"
	);

    expect_equal (names(response), expectedNames)
})

test_that(".avisApiBusAvanzada returns a dataframe",{ 
    
    expect_true(is.data.frame (response))
})

test_that(".avisApiBusAvanzada returns expected species",{ 
    
    # note: test coupled to remote database
    expect_equal(as.character(response$Especie[1]), expectedSpecies)
})
