Sys.setlocale("LC_COLLATE", "C")

test_that("uniqtag aaaaaab k=3", expect_equal(
	unname(uniqtag(c("aaaaaab", "aaab"), k = 3)),
	c("aaa-1", "aaa-2")))

test_that("uniqtag aaaaaab k=4", expect_equal(
	unname(uniqtag(c("aaaaaab", "aaab"), k = 4)),
	c("aaaa", "aaab")))

states <- sub(" ", "", state.name)

states3 <- setNames(c(
	"aba-1", "las-1", "Ari-1", "Ark-1", "Cal-1", "Col-1", "Con-1", "Del-1",
	"Flo-1", "Geo-1", "Haw-1", "Ida-1", "Ill-1", "Ind-1", "Iow-1", "Kan-1",
	"Ken-1", "Lou-1", "Mai-1", "Mar-1", "Mas-1", "Mic-1", "Min-1", "ipp-1",
	"our-1", "Mon-1", "Neb-1", "Nev-1", "Ham-1", "Jer-1", "Mex-1", "Yor-1",
	"Car-1", "Dak-1", "Ohi-1", "Okl-1", "Ore-1", "Pen-1", "Isl-1", "Car-2",
	"Dak-2", "Ten-1", "Tex-1", "Uta-1", "Ver-1", "Vir-1", "Was-1", "Wes-1",
	"Wis-1", "Wyo-1"), states)
test_that("uniqtag states k=3", expect_equal(uniqtag(states, k = 3), states3))

states4 <- setNames(c(
	"Alab", "Alas", "Ariz", "Arka", "Cali", "Colo", "Conn", "Dela",
	"Flor", "Geor", "Hawa", "Idah", "Illi", "Indi", "Iowa", "Kans",
	"Kent", "Loui", "Main", "Mary", "Mass", "Mich", "Minn", "ippi",
	"isso", "Mont", "Nebr", "Neva", "Hamp", "Jers", "Mexi", "NewY",
	"rthC", "rthD", "Ohio", "Okla", "Oreg", "Penn", "Isla", "uthC",
	"uthD", "Tenn", "Texa", "Utah", "Verm", "Virg", "Wash", "West",
	"Wisc", "Wyom"), states)
test_that("uniqtag states k=4", expect_equal(uniqtag(states, k = 4), states4))
