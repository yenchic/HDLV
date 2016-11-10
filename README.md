# HDLV
The level set estimation and visualization. 

- Paper reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Density Level Sets: Asymptotics, Inference, and Visualization." arXiv preprint arXiv:1504.05438 (2015).
- Contact: yenchic@uw.edu

## HDLV_RS.R


### HDLV
- The main function for computing level set estimation and visualization.
- Inputs:
  - data: Input data matrix.
  - h: Smoothing parameter. Default \emph{NULL} will choose from the normal reference rule.
  - lv_seq: A sequence of density levels. Default \emph{NULL} will choose from a sequence of level.
  - lv_min: Minimal density level (to stablized the level sets). Default \emph{NULL} will pick 0.05*p_max.
  - eps: The tolerance. If mean shift moves less than this value, we will consider it done.
  - max.iterations: Maximal number of iteration for mean shift.
  - cut: The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
- Outputs:
  - An S4 object "HDLV" consisting of the following attributes:
    - data: The original data.
    - h: The smoothing parameter.
    - lv.min: The minimal level.
    - density: The density by kernel density estimator.
    - MS.labels: The cluster labels from mode clustering.
    - modes: The local modes corresopnding to each cluster.
    - sig.clu: The significant clusters (local modes).
    - lv_seq: The sequence of level sets.
    - size: The matrix for the sizes of each significant clusters at different density levels.
    - conn: The list of connectivity matrices at each density level.


