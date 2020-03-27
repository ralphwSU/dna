.PHONY : c
c :
	latex src.tex; dvips -o src.ps src.dvi; ps2pdf src.ps

.PHONY : s
s :
	latex Corona_genomes_22_Feb_2020.tex;
	dvips -o Corona_genomes_22_Feb_2020.ps Corona_genomes_22_Feb_2020.dvi;
	ps2pdf Corona_genomes_22_Feb_2020.ps
	mv Corona_genomes_22_Feb_2020.pdf tmp.pdf
	pdftk A=tmp.pdf B=cv.pdf C=ebola_sample.pdf cat A B C output Corona_genomes_22_Feb_2020.pdf

.PHONY : e
e :
	latex editTest.tex; dvips -o editTest.ps editTest.dvi;
	ps2pdf editTest.ps

.PHONY : clean
clean :
	rm *.log *.dvi *.ps *.aux *.out *~
