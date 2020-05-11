all: trim_output
	maple basisConversion2.mpl
	./trim_output
	mv out.txt output.txt

trim_output: trim_output.cpp	
	g++ -o $@ $^

.PHONY: clean
clean:
	rm -rf output.txt trim_output

