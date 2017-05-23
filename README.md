# STAG-CNS
An Order-Aware Conserved Non-coding Sequences (CNS) Discovery Tool For Arbitrary Numbers of Species

# Compile the program
g++ -o stagCNS stagCNS.cpp -std=c++11

# Input
The first n ( n = number of sequences) lines contain some additional information about the genes. It contains seven fields separated by space. The fields are gene name, chromosome ID, forward(+)/reverse(-) strand, range of 10kbp upstream, rage of 10kbp downstream respectively. Please check the sample test file (input.fasta).

# Run the program
./stagCNS -file fasta_file -mem minimum_CNS_length  -out prefix_output_file   >  out.txt

Example: ./stagCNS -file  input.fasta  -mem 8   -out TEST  > out.txt

# Outputs
1. TEST_CNS_8.csv (contains CNS information)
2. TEST_MEM_1_8.csv (contains all MEMs)
3. TEST_MEM_2_8.csv (contains all MEMs without serial number)
4. TEST_LPMEM_8_csv (contains all CNSs without actual chromosome location ) 
5. TEST_MEM_8.html  (used for browser visualization )
6. TEST_LPMEM_8.html (used for browser visualization)
7. out.txt (contains all console outputs and error logs)

The first file contains the CNS information. The other files are used for visualization.

# Visualization
Two ways to visualiza the output

1. Using XDAT: 
   Run the collowing command on command line

   java -jar xdat.jar
   
   It will open an window where the TEST_MEM_2_8.csv and TEST_LPMEM_8_csv  files can be uploaded (data) to see the visualizations
   (chart->parallel coordinate set). 
   Don't forget to change the settings (import settings) to accept comma separated files.

2. Using browser (firefox)
 Run the following commands:

     firefox TEST_MEM_8.html
 
     firefox TEST_LPMEM_8.html

3. Gobe Visualization will be added
