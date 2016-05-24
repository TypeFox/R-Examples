# Copyright PGE ToolBox (Cai et. al) !
# Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
# Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127

codontable <- function(){

PROTEIN <- c("F F L L S S S S C C W * Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -", # standard nuclear Dna
"F F L L S S S S C C W W Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M M T T T T S S * * N N K K -",
"F F L L S S S S C C W W Y Y * * T T T T P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M M T T T T S S R R N N K K -",
"F F L L S S S S C C W W Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -",
"F F L L S S S S C C W W Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M M T T T T S S S S N N K K -",
"F F L L S S S S C C W * Y Y Q Q L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -",
"F F L L S S S S C C W W Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S S S N N K N -",
"F F L L S S S S C C W C Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -",
"F F L L S S S S C C W * Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -",
"F F L L S S S S C C W * Y Y * * L L S L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -",
"F F L L S S S S C C W W Y Y * * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M M T T T T S S G G N N K K -",
"F F L L S S S S C C W W Y Y * Y L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S S S N N K N -",
"F F L L S S S S C C W * Y Y Q * L L L L P P P P R R R R H H Q Q V V V V A A A A G G G G D D E E I I M I T T T T S S R R N N K K -")

Protein <- matrix(unlist(strsplit(PROTEIN," ")),ncol=65,byrow=13)

Triplets <- matrix(
c(
1, 1, 1,
1, 1, 2,
1, 1, 3, 
1, 1, 4, 
1, 2, 1, 
1, 2, 2,
1, 2, 3, 
1, 2, 4, 
1, 3, 1,
1, 3, 2, 
1, 3, 3, 
1, 3, 4, 
1, 4, 1, 
1, 4, 2, 
1, 4, 3, 
1, 4, 4, 
2, 1, 1, 
2, 1, 2, 
2, 1, 3, 
2, 1, 4, 
2, 2, 1, 
2, 2, 2, 
2, 2, 3, 
2, 2, 4, 
2, 3, 1, 
2, 3, 2, 
2, 3, 3, 
2, 3, 4, 
2, 4, 1, 
2, 4, 2, 
2, 4, 3, 
2, 4, 4, 
3, 1, 1, 
3, 1, 2, 
3, 1, 3, 
3, 1, 4, 
3, 2, 1, 
3, 2, 2, 
3, 2, 3, 
3, 2, 4, 
3, 3, 1, 
3, 3, 2, 
3, 3, 3, 
3, 3, 4, 
3, 4, 1, 
3, 4, 2, 
3 ,4, 3, 
3, 4, 4, 
4, 1, 1,
4 ,1 ,2, 
4, 1, 3, 
4, 1, 4, 
4, 2, 1, 
4, 2, 2, 
4, 2, 3,
4, 2, 4, 
4, 3, 1, 
4, 3, 2, 
4, 3, 3, 
4, 3, 4, 
4, 4, 1, 
4, 4, 2, 
4, 4, 3, 
4, 4, 4),ncol=3,byrow=3)

Triplets <- rbind(Triplets,c(5,5,5))

return(list(Protein=Protein,Triplets=Triplets))

}

