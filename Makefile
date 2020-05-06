all:
	maple basisConversion.mpl
	g++ -o trim trim.cpp
	./trim
	rm -rf trim
	mv out.txt output.txt
