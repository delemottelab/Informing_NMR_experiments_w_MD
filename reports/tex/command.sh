# Clean-up
latexmk -C
# Compile
latexmk -r latexmkrc
# Make diff file only showing addition
git latexdiff HEAD~1 --main aip_text.tex -o diff.pdf --bibtex --allow-spaces --latexmk --preamble latexdiff_onlyblue --graphics-markup=2
