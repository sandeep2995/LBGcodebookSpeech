# LBGcodebookSpeech
Generates codebook using LBG algorithm

All the following files and graphs are in "LBG/" directory

config.h---> configuration file to set various things such as delimiter, filename, thresold distortion allowed, precision of the code vectors, etc.

lbg---> implements lbg algorithm on universe of cepstral coefficient vectors.

set the values in config.h and run the project.

distortion will be saved in "logdistortion.txt" (average distortion per cluster)

select all and paste in excel sheet to see the graph of average distortion Vs. number of iterations.(see various graphs lbg1.png, lbg2.png, lbg3.png and lbg4.png attached)

"codebook.txt" contains the final code vectors.

"logclustersize.txt" contains intermediate result indicating number of vectors associated with each cluster after every iteration.



P.S:

Program that we used to compute the universe of cepstral coefficients also attached.
"cep.txt" contains the universe of cepstral coeffients separated by space.
"universe.txt" contains the universe of cepstral coeffients separated by tab.
