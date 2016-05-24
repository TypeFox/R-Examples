library(testthat)
library(lettercase)
library(magrittr)

context( "str_snake_case" )

str_snake_case( 'mission_of_burma' ) %>% expect_equal( 'mission_of_burma' ) 
str_snake_case( 'mission-of-burma' ) %>% expect_equal( 'mission_of_burma' )
str_snake_case( 'mission.of.burma' ) %>% expect_equal( 'mission_of_burma' )
str_snake_case( 'Mission of Burma' ) %>% expect_equal( 'mission_of_burma' )
str_snake_case( ' .mission.of.burma. ' ) %>% expect_equal( 'mission_of_burma' )

# str_snake_case( 'mission-of-burma' ) %>% expect_equal( 'mission_of_burma' )
# str_snake_case( 'mission_of_burma' ) %>% expect_equal( 'mission_of_burma' )
# str_snake_case( 'MissionOfBurma' )   %>% expect_equal( 'mission_of_burma' )
