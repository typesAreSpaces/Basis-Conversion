#include <iostream>
#include <fstream>
#include <string>

int main(){

  std::ifstream in("output.txt");
  std::ofstream out("out.txt");
  std::string line;

  if(in.is_open()){
    while(getline(in, line))
      if(line[0] != '>' && line.size() > 0)
        out << line << std::endl;
    in.close();
  }
  else
    std::cout << "Unable to open file" << std::endl;
  
  return 0;
}
