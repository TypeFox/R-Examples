
# Example of using the postscript terminal, with special characters and symbols specified as PostScript symbols
Gprun("set terminal postscript size 8cm,6cm eps enhanced color font \"Helvetica,10\"\nset output \"SinCos.eps\"\nset style line 1 linecolor rgb \"blue\" linetype 1 linewidth 5  # blue\nset style line 2 linecolor rgb \"red\" linetype 1 linewidth 5  # red\nset xlabel \"{/Helvetica-Italic x}\"\nset ylabel \"{/Helvetica-Italic y}\"\nset xtics (\"-2{/Symbol p}\" -2*pi, \"-{/Symbol p}\" -pi, 0, \"{/Symbol p}\" pi, \"2{/Symbol p}\" 2*pi)\nplot [-2*pi:2*pi][-1.5:1.5] sin(x) tit \"sin({/Helvetica-Italic x})\" with lines ls 1, cos(x) tit \"cos({/Helvetica-Italic x})\" with lines ls 2")

# Example of using the latex terminal with special characters and symbols specified as Latex symbols
Gprun("set terminal latex size 8cm, 6cm color courier 10\nset output \"SinCosV0.tex\"\nset style line 1 linecolor rgb \"blue\" linetype 1 linewidth 5  # blue\nset style line 2 linecolor rgb \"red\" linetype 1 linewidth 5  # red\nset xlabel \"$x$\"\nset ylabel \"$y$\"\nset xtics (\"$-2\\\\pi$\" -2*pi, \"$-\\\\pi$\" -pi, 0, \"$\\\\pi$\" pi, \"$2\\\\pi$\" 2*pi)\nplot [-2*pi:2*pi][-1.5:1.5] sin(x) tit \"$\\\\sin(x)$\" with lines ls 1, cos(x) tit \"$\\\\cos(x)$\" with lines ls 2")

# Compile the following code with Tex -> DVI -> PDF to display the previous result \begin{lstlisting}[language=R,caption={Latex code to show
# SinCosV0.tex},label={lst:textshowSinCosV0}] \documentclass{article} \begin{document} \begin {figure} \input{SinCosV0} \end {figure} \end{document} \end{lstlisting}






 
