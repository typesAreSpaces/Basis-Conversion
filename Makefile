all: trim
	maple basisConversion.mpl
	./trim
	mv out.txt output.txt

trim: trim.cpp
	g++ -o trim trim.cpp

.PHONY: clean
clean:
	rm -rf output.txt trim

