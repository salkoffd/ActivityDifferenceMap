# ActivityDifferenceMap
This function computes the intensity difference and statistics for images of two conditions. It is a permutation test that follows the logic of
this paper: 
[Nonparametric permutation tests for functional neuroimaging: a primer with examples](https://www.ncbi.nlm.nih.gov/pubmed/11747097).
The output adm.map is the algebraic difference between condition 1 and 2. To derive statistics (without assuming a Gaussian distribution),
a permutation test is performed to test if each pixel is significantly larger (or smaller) than expected if the 
two conditions are identical. Finally, clusters of significant pixels are tested to see if they are significantly larger than expected by
chance.

As an example, I tested to see if neural activity was significantly different after a visual stimulus.
This image (automatically generated from the function) shows the dorsal aspect of the left hemisphere of a mouse. Neural activity was 
significantly higher after the stimulus (condition 1) in the anterior (left) and medial region, outlined in black.

![Alt Text](https://i.imgur.com/8otYK9S.png)

**Usage**
-----

You can run this by doing:

    ActivityDifferenceMap(responseHits, responseMisses, nIterations)

Inputs  
- responseHits is a h x w x t1 matrix (3D) of images from condition 1. hxw is the height and width of image. t1 is the number of images.  
- videosMisses is a h x w x t2 matrix (3D) of images from condition 2. t2 is the number of images and does not need to be equal to t1.
- nIterations is the number of iterations for the permutations test (default 1000).

Outputs  

adm is a structure containing the following:
- adm.map is the map of activity difference.  
- adm.pmap is the p-value estimate per pixel.  
- adm.nPixelsPval is the p-value associated with the total number of pixels with significant difference (compared with null distribution)  
- adm.clusterSizePvals is the p-value associated with the size of each cluster (in pixels)  
- adm.clusterMassPvals is the p-value associated with the mass of each cluster  
- adm.clustersSigPos is a logical map showing where positive sign clusters are (Hit trial activity significantly larger than Miss trial activity)  
- adm.clustersSigPos is a logical map showing where negative sign clusters are (Hit trial activity significantly smaller than Miss trial activity)  
