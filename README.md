# Border-Analysis
Spatial scan methods are commonly used to identify and analyze geographic disease patterns. Most often, they are used to identify clusters or areas of high disease incidence. When testing for clusters, there are many proposed ways to construct the list of locations to test i.e., the set of candidate zones. Two of the most popular methods are the circular scan test and the elliptic scan test, which achieve their popularity from their accessibility, efficiency, and accuracy. One limitation of the scan methods is the inability to determine the exact border of a cluster. To combat this border uncertainty, Oliveira et al. proposed the F function, which calculates a value between [0, 1] for each region in the study area, where higher values indicate more evidence that the region belongs to the true cluster. We have created R scripts for the F function for both the circular scan test and the elliptic scan test. 
