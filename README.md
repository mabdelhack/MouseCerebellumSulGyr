# MouseCerebellumSulGyr

These are the codes and data for reproducing results of the manuscript [Comparative Morphology of Dendritic Arbors in Populations of Purkinje Cells in Mouse Sulcus and Apex](http://dx.doi.org/10.1155/2013/948587)

To start, run the code Neighbors_Overlap.m changing the values of the following variables
- s1: to choose the mouse. Possible values are P121, P123 for young and P1501, P1502 for old.
- Area: sulcus or apex

and then for each run of the previous script run slice_overlap2.manuscript to generate the output file


The code depends on the library [Polygon Intersection](https://www.mathworks.com/matlabcentral/fileexchange/18173-polygon_intersection).  It is reuploaded here for convenience.

## Output
Excel sheet is generated with all the overlaps and self-overlaps in the folder of each region and mouse

More details to follow.
