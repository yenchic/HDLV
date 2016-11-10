# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@uw.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Density Level Sets: Asymptotic, Inference, and Visualization."
# Date: 04/20/2015
#High dimension level set
library(scatterplot3d)
source("HDLV_RS.R") ## For R-studio
#source("HDLV.R") ## For ordinary R

  ## Five cluster dataset (d = 3+d.add dimensions)
X = five_cluster(d.add=7)
  ## Plotting the first three coordinates
v.lim=c(-0.2,1.2)
S3d = scatterplot3d(X[,1:3], cex.symbols=0.5,tick.marks=FALSE, xlab="", ylab="", zlab="",mar=c(0.1,0.1,0.1,0.1), xlim=v.lim, ylim=v.lim,zlim=v.lim)


	## show that this is indeed a high dimensional dataset
head(X)



  ## HDLV: High dimensional level set
X_HDLV = HDLV(X)
  ## Visualization
plot(X_HDLV)



  ## Customized level sequences.
plot(X_HDLV, lv_seq= (10:50)/100*max(X_HDLV$density))



  ## Summary information for HDLV.
summary(X_HDLV)
names(X_HDLV)




  ## (R-studio only) An interactive level set visualization.
Intplot.HDLV(X_HDLV, r0=1, l_size =30)


