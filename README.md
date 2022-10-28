# MCCD-WN
MCCD-WN stands for Maximal Clique based Community Detection algorithm for undirected Weighted Network. This algorithm uses a similarity matrix constructed by incorporating information obtained from maximal-cliques and leveraging both the local and global importance of nodes. This algorithm is developed for a sparse (OTU/ ASV  co-abundance) network. <br/>
In order to run the code, please make sure that the following R pacakges are installed: <br/>
I. data.table <br/>
II. igraph <br/>
III. rARPACK <br/>
The code was tested on Ubuntu 22.04.1 LTS (GNU/Linux 5.15.0-48-generic x86_64) having gcc version 11.2.0 and R version 4.2.1 installed on it. <br/>
You just need to execute the shell script (command.sh) <br/>
Before executing the script, kindly make all necessary changes in the file itself. <br/>
If you face any difficulties to run this code or want to apply the entire workflow as used in [1] for analysing your OTU/ ASV co-abundance networks, please feel free to write an email to anirban.bhar@uni-greifswald.de and we would be happy to assist you or share the rest of the code. <br/>
If you find MCCD-WN and/ or the entire workflow useful for your project, please cite the following publication: <br/> <br/> <br/>
[1] Bhar A, Gierse LC, Meene A, Wang H, Karte C, Schwaiger T, Schröder C, Mettenleiter TC, Urich T, Riedel K and Kaderali L. (2022) Application of a
maximal-clique based community detection algorithm to gut microbiome data reveals driver microbes during influenza A virus infection. Front. Microbiol. 13:979320. doi: 10.3389/fmicb.2022.979320


