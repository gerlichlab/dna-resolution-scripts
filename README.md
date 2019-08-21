# dna-resolution-scripts

## Local_separation_along_linescan
The script needs linescans in csv format. File directory can be choosen in the block under the import block.
Linescans were  generated with the imageJ macro "allProfils.ijm". Linescans need to have at least the column Distance, Hoechst and f-ara-EdU.
The Script reads the csv file and applies a z-standardization.

![Illustration of Programm](https://github.com/gerlichlab/dna-resolution-scripts/blob/master/local%20co%20localization.jpg)


The script recognizes everything over a z-score= -1 as chromosomal mass. 
Since its possible, that there is more than one chromosome along the line, the biggest consecutive block of intensity over -1 is determined as chromosome of interest. If two consecutive blocks havethe same size, the dataset is skipped.
The chromosome dataset is devided into two equal parts. For EdU and Hoechst curve, the integrals for each side are calculated. The side with a bigger EdU integral is said to be the EdU containing sister.
To calculate percentages the Integral for the EdU containing side is divided by the some of both. This is done for EdU and Hoechst respectively.

The Output of the script is a pandas dataframe, that is saved to the choosen directory as a csv file.
Columns are File, Left intergal, right integral, ratio (EdU/non-EdU) and percentage. Data is given for EdU and Hoechst respectively combined in one Dataframe.
