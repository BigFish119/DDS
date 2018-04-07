# DDS
The algorithm DDS is the abbreviation of Distance Dynamic Synchronization. DDS is used to detect non-overlapping community sturcture in the network. A network is saved in the file with a .txt suffix. Each line consisting of two numbers in the file is an edge of the network. Each number in the file is the id of a vertex.
1.INPUT FILE FORMAT
(1)The input file must be end with .txt suffix;
(2)Each number in the file is the id of a vertex;
(3)Each line consisting of two numbers in the file is an undirected edge in the network. The two numbers are separated by tab.
(4)The first number in each line must be smaller than the second number;
(5)Each edge appear only once in the file.
(6)All numbers must be consecutive integers starting at 1;
  
2.OUTPUT FILE FORMAT
(1)The output file must be end with .txt suffix;
(2)Each line contains two numbers in the output file: the first number is the vertex id, and the second number is its community id;
  
3.STEPS OF RUNNING DDS
(1)Change the value of global variable NetNodeNum(it is number of vertices existing in the network) in line 8;
(2)Compile the source code;
  For Example:on windows 7 OS:Compile with Virtual Studio 2013;
              on Linux:g++ DDS.cpp -o DDS
(3)Execute DDS
  For Example:on Windows 7 OS:use the cd command to the directory of DDS.exe in the command line, and input the commands "DDS.exe karate.txt resultKarate.txt 78" and press the "Enter" key to run the algorithm.
              on Linux:Input the commands "./DDS karate.txt resultKarate.txt 78" in the comand line, and press the "Enter" key to run the algorithm.

NOTE:
the second input is the input file name, the third input is the output file name, the forth input is the number of edges in the network. The Karate network contains 78 edges, so the forth input is 78.
