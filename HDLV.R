# Author: Yen-Chi Chen, Christopher R. Genovese, Larry Wasserman
# Maintainer: Yen-Chi Chen <yenchic@andrew.cmu.edu>
# Reference: Chen, Yen-Chi, Christopher R. Genovese, and Larry Wasserman. "Density Level Sets: Asymptotic, Inference, and Visualization."
# Date: 04/20/2015
### High dimensional level set
### libraries
#' @import RANN
library(RANN)	# knn
#' @import TDA
library(TDA)	# kde
#' @import Rcpp
library(Rcpp) # writing c-code
#' @import mvtnorm
library(mvtnorm)  # generating five clusters dataset
#' @import plotrix
library(plotrix)	# for plotting circle

cppFunction('NumericVector MSCpp(NumericMatrix data, NumericMatrix query, double h, int max_iterations , double eps){
	int n = data.nrow(), d = data.ncol(), m = query.nrow();
	double dist_tmp; 
	NumericVector K(n);
	double K_tot;
	
	NumericVector MS(d);
	NumericVector newPt(d);
	NumericVector oldPt(d);
	int iter_now;
	double err_now;
	
	NumericVector result(d*m);
	

	for(int w=0; w<m; w++){
		// for each point
		for(int j = 0; j<d;j++){
			newPt(j) = query(w,j);
		}
		err_now = 1e14;
		iter_now = 0; 
		
		while((iter_now < max_iterations)&&(err_now > eps)){	// ignore those nearly not change
		
			for(int j =0; j<d; j++){
				MS(j) =0;
				oldPt(j) = newPt(j);
			}
			
			K_tot = 0;
			for(int i = 0; i<n; i++){
				dist_tmp = 0;
				for(int j =0; j<d; j++){
					dist_tmp += (data(i,j) - newPt(j))*(data(i,j) - newPt(j))/h/h;
				}
				K(i) = exp(-1*dist_tmp/2);
				K_tot += K(i);
				for(int j =0; j<d; j++){
					MS(j) += data(i,j)*K(i);
				}
			}
				
			// updates & errors
			err_now = 0;
			for(int j =0; j<d; j++){
				newPt(j) = MS(j)/K_tot;
				err_now =+ (newPt(j) - oldPt(j)) * (newPt(j) - oldPt(j));
			}
			err_now = sqrt(err_now);
			iter_now++;
		}
		
		for(int j =0; j<d; j++){
			//result(w*d+j) = newPt(j);
			result(w+m*j) = newPt(j);
		}
	}
	
	return result;

}')


#' Mean shift algorithm using Rcpp
#' 
#' @param data Input data matrix.
#' @param query The mesh points that you want to apply mean shift to.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @return The mesh points after mean shift.
#' @examples
#' x = matrix(rnorm(1000), ncol=2)
#' x_MS = MS(x, x, 0.5)
#'  # mean shift
#' 
#' plot(x, cex=0.5)
#' points(x_MS, col="red", pch=20, cex=2)
#'  # plot the shifted result versus original case
#' 
#' @export
MS = function(data, query, h, max.iterations=200, eps= 1e-15){
	tmp = MSCpp(data=data, query=query, h=h, max_iterations= max.iterations, eps=eps)
	return(matrix(tmp, ncol=ncol(query)))
}

#' Fast mean shift using heirachical clustering.
#' 
#' @param data Input data matrix.
#' @param query The mesh points that you want to apply mean shift to.
#' @param h Smoothing parameter.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @param cut The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
#' @return The mode clustering result; a list consisting of
#'  \item{label}{the cluster labels for query points.}
#'  \item{modes}{the local modes corresponding to each label.}
#'  
#' @export
fms = function(data, query, h, eps=1.0e-8, max.iterations=100, cut = 0.1){

	result = list()

	pt_ms = MS(data, query, h=h, eps=eps, max.iterations= max.iterations)
	D1 = dist(pt_ms)
	D1_hclust = hclust(D1)
	cluster_lab= cutree(D1_hclust, h=0.1*h)
		# find the cluster lable
		
	modes = matrix(NA, nrow=max(cluster_lab), ncol=ncol(X))
	for(i in 1:max(cluster_lab)){
		w_tmp = which(cluster_lab==i)
		if(length(w_tmp)>1){
			modes[i,] = colSums(pt_ms[w_tmp,])/length(w_tmp)
		}
		if(length(w_tmp)==1){
			modes[i,] = pt_ms[w_tmp,]
		}
	}
	result$label = cluster_lab
	result$modes = modes
	
	return(result)
}

#' High dimensional five clusters dataset
#' @param N.c Sample size per each cluster
#' @param N.f Sample size per filaments
#' @param dis.c Randomness around clusters
#' @param dis.f Randomness around filaments
#' @param d.add Added dimensions.
#' @return A five clusters dataset with 3+d.add dimensions.
#' @export
five_cluster = function(N.c=200, N.f=100, dis.c=0.08, dis.f=0.03, d.add=7){
	  ## setting clusters
  C0.1 = c(0,0,0)
	C0.2 = c(1,0,0)
	C0.3 = c(0,1,0)
	C0.4 = c(0,0,1)
	C0.5 = c(0,1,1)
	C0 = cbind(rbind(C0.1,C0.2,C0.3,C0.4,C0.5), matrix(0,nrow=5, ncol=d.add))
	
    ## settting total dimensions
	d= d.add+3
	
    ## adding edges
	C = NULL
	for(i.c in 1:5){
	  C = rbind(C, rmvnorm(N.c,C0[i.c,],sigma = diag(rep(dis.c^2,d))))
	}
	U0.1= runif(N.f)
	E1 = t((C0.2-C0.1)%o%U0.1)
	U0.1= runif(N.f)
	E3 = t((C0.4-C0.1)%o%U0.1)
	U0.1= runif(N.f)
	E4 = t((C0.5-C0.4)%o%U0.1+C0.4)
	E0 = cbind(rbind(E1,E3,E4),matrix(0,nrow=3*N.f, ncol=d.add))	
	E = E0+matrix(rnorm(nrow(E0)*ncol(E0),sd=dis.f), nrow=nrow(E0), ncol=ncol(E0))
	
	X = rbind(C,E)
	return(X)
}


#' Normal reference rule
#' @param data Input data.
#' @return Smoothing parameter by normal reference rule.
#' @export
h_NR = function(data){
		d = ncol(data)
		n = nrow(data)
		h = (4/(d+2))^(1/(d+4))/n^(1/(d+4))*mean(apply(data,2,sd))
		return(h)
}

#' Connection detection.
#' 
#' Checking if a pair of sample points are connected in the upper level set.
#' @param X1 First sample.
#' @param X2 Second sample.
#' @param X1.den Density for first data.
#' @param X2.den Density for second data.
#' @param lv Density level.
#' @param data Total data.
#' @param h Smoothing parameter.
#' @return 0 or 1. 1 indicates that the two sample points are connected.
#' 
pairConn = function(X1, X2, X1.den, X2.den, lv, data, h){
	d = ncol(data)
	if(length(X1)==d){
		X1= rbind(X1,X1)
		X1.den = c(X1.den,X1.den)
	}
	if(length(X2)==d){
		X2= rbind(X2,X2)
		X2.den = c(X2.den, X2.den)
	}
		
	knn12 = nn2(X1,X2, k=1)
	knn21 = nn2(X2,X1, k=1)
	idx1_2 = which(knn12$nn.dist==min(knn12$nn.dist)) #on X2
	idx1_1 = knn12$nn.idx[which(knn12$nn.dist==min(knn12$nn.dist))] #on X1
	
	idx2_1 = which(knn21$nn.dist==min(knn21$nn.dist)) #on X1
	idx2_2 = knn21$nn.idx[which(knn21$nn.dist==min(knn21$nn.dist))] #on X2
	
	q1 = X1[idx1_1,]
	q2 = X2[idx1_2,]
	
	if(min(knn12$nn.dist)> min(knn21$nn.dist)){
		q1 = X1[idx2_1,]
		q2 = X2[idx2_2,]
	}
	min_dist = min(min(knn12$nn.dist), min(knn21$nn.dist))
	
	nl_seq = 20
	nl = sapply(1:d, function(x){
		tmp = seq(from=q1[x], to=q2[x], length.out=nl_seq)
		return(tmp)
	})
	
	nl_den = kde(data,nl,h)
	
	result = 1
	if(sum(nl_den<lv)>0||min_dist>2*h){
		# if distance is too large (2*h is roughly the influence radius for Gaussian kernel), we will not accept it
		result = 0
	}	
	return(result)
}


#' Cluster connectivity at certain level.
#' @param data Input data.
#' @param data.den Density for input data.
#' @param data.label Mode clustering lavel for each data point.
#' @param lv Density level.
#' @param h Smoothing parameter.
#' @param lv_min Minimal density level (to stablized the level sets).
#' @return A matrix of 0 and 1 indicating if two clusters are connected (1: connected).
CluConn = function(data, data.den, data.label, lv, h, lv_min=0.05*max(data.den)){
	cl_sig = intersect(unique(data.label[which(data.den>lv_min)]), which(table(data.label)>2))
	
	Conn_m = matrix(0, nrow=length(cl_sig), ncol=length(cl_sig))
	colnames(Conn_m) = cl_sig
	rownames(Conn_m) = cl_sig	
	for(i in 1:(length(cl_sig)-1)){
		for(j in (i+1):length(cl_sig)){
			idx_i = which(data.label==cl_sig[i])
			idx_j = which(data.label==cl_sig[j])
			Conn_m[i,j] = pairConn(data[idx_i,], data[idx_j,], data.den[idx_i], data.den[idx_j], lv=lv, data=data, h=h)
			Conn_m[j,i] = Conn_m[i,j]
		}
	}
	return(Conn_m)
}


#' Summary information about clusters under various density levels.
#' 
#' @param data Input data.
#' @param data.den Density for input data.
#' @param data.label Mode clustering lavel for each data point.
#' @param lv_seq A sequence of density levels.
#' @param h Smoothing parameter.
#' @param lv_min Minimal density level (to stablized the level sets).
#' @return A list consisting of 
#' \item{Conn}
#' {Matrices of 0 and 1 indicating if two clusters are connected (1: connected).}
#' \item{S_matrix}
#' {A matrix indicating the size (number of data points) of clusters at each level.}

CluInfo = function(data, data.den, data.label, lv_seq, h, lv_min=0.05*max(data.den)){
	result = NULL
	cl_sig = intersect(unique(data.label[which(data.den>lv_min)]), which(table(data.label)>2))
	
	l_conn = lapply(lv_seq, function(x){
		tmp = CluConn(data=X, data.den=data.den, data.label=data.label, lv=x, h=h)
		return(tmp)
		}
	)
	
	S_matrix = matrix(NA, ncol=length(cl_sig), nrow=length(lv_seq))
	for(i_lv in 1:length(lv_seq)){
		w_tmp = which(data.den>lv_seq[i_lv])
		for(i_cl in 1:length(cl_sig)){
			S_matrix[i_lv, i_cl] = length(intersect(w_tmp, which(data.label==cl_sig[i_cl])))
		}
	}

	colnames(S_matrix) = cl_sig
	
	result$Conn = l_conn
	result$Size = S_matrix
	return(result)
}

#' Creating an HDLV object.
#' @param data Input data.
#' @param ... See \emph{HDLV.default}.
#' @return An HDLV object.
#' @export
HDLV = function(data, ...) UseMethod("HDLV")

#' High dimensional density level set.
#' @param data Input data matrix.
#' @param h Smoothing parameter. Default \emph{NULL} will choose from the normal reference rule.
#' @param lv_seq A sequence of density levels. Default \emph{NULL} will choose from a sequence of level.
#' @param lv_min Minimal density level (to stablized the level sets). Default \emph{NULL} will pick 0.05*p_max.
#' @param eps The tolerance. If mean shift moves less than this value, we will consider it done.
#' @param max.iterations Maximal number of iteration for mean shift.
#' @param cut The cut for heirachical clustering (we cut the dedrogram by height = cut*h).
#' @return An S4 object "HDLV" consisting of the following attributes:
#' \item{data}
#' {The original data.}
#' \item{h}
#' {The smoothing parameter.}
#' \item{lv.min}
#' {The minimal level.}
#' \item{density}
#' {The density by kernel density estimator.}
#' \item{MS.labels}
#' {The cluster labels from mode clustering.}
#' \item{modes}
#' {The local modes corresopnding to each cluster.}
#' \item{sig.clu}
#' {The significant clusters (local modes).}
#' \item{lv_seq}
#' {The sequence of level sets.}
#' \item{size}
#' {The matrix for the sizes of each significant clusters at different density levels.}
#' \item{conn}
#' {The list of connectivity matrices at each density level.}
HDLV.default = function(data,h = NULL, lv_seq=NULL, lv_min=NULL, eps=1.0e-8, max.iterations=1000, cut=0.1, ...){
	result = list()
	
	if(is.null(h)){
		h = h_NR(data)
	}
	
	X_MS = fms(data, data, h=h, eps=eps, max.iterations = max.iterations, cut= cut)
	X_kde = TDA:::kde(data, data, h=h)

	if(is.null(lv_min)){
		lv_min = 0.05*max(X_kde)
	}
	if(is.null(lv_seq)){
		lv_seq = (1:20)/20*max(X_kde)
		flag_seq = T
	}
	
	cl_sig = intersect(unique(X_MS$label[which(X_kde>lv_min)]), which(table(X_MS$label)>2))
	C_info = CluInfo(data, X_kde, X_MS$label, lv_seq, h=h, lv_min=lv_min)
	
	result$data = data
	result$h = h
	result$lv.min = lv_min
	
	result$density = X_kde
	result$MS.labels = X_MS$label
	result$modes = X_MS$modes
	result$sig.clu = cl_sig
	if(flag_seq)
		result$lv_seq.type="p_max"
		
	result$lv_seq = lv_seq
	result$size = C_info$Size
	result$conn = C_info$Conn

	
	class(result) = "HDLV"
	return(result)
}
print.HDLV = function(object,...){
	cat("Density Estimate:\n")
	print(object$density)
	cat("\n")
}
summary.HDLV = function(object,...){
	cat("Density Estimate:\n")
	print(object$density)
	cat("\n")
	cat("Mode Clustering:\n")	
	print(object$MS.labels)
	cat("\n")
	cat("Modes:\n")
	print(object$modes)
	cat("\n")
	cat("Significant clusters:\n")
	print(object$sig.clu)
	cat("\n")
	cat("Level sequence type:\n")
	print(object$lv_seq.type)
	cat("\n")
	cat("Cluster size at each level:\n")
	print(object$size)
	cat("\n")
	cat("To see how mode clusters are connected into upper level set cluster, use:\n object$conn")
}


#' Plotting HDLV.
#' @param x An HDLV object.
#' @param lv_seq A sequence of density levels. Default \emph{NULL} will choose from a sequence of level.
#' @param col_code The color for the level set. Default \emph{NULL} will automatically generate.
#' @param r0 Visualization size for each cluster.
#' @param l_size Line width for the connection between clusters.
#' @param legend.default True or False. To deplay the default legend or not.
plot.HDLV = function(x, lv_seq=NULL, col_code=NULL, r0 = 0.8, l_size = 10, legend.default=T, ...){
	flag_seq=T
	flag_max="n"
	if(is.null(lv_seq)){
		lv_seq = x$lv_seq
		flag_seq = F
		flag_max = x$lv_seq.type
	}
	if(is.null(col_code))
		col_code = colorRampPalette(c("lightgreen","limegreen","yellow","orange","red"))(length(lv_seq))

	M_sig = x$modes[x$sig.clu,]
	M_mds = cmdscale(dist(M_sig))
	row.names(M_mds) = x$sig.clu

	## plotting parameters
	S_max = max(x$size)
	minM_dist = min(dist(M_mds))	
	N_size = x$size/S_max
	conn_list = x$conn
	## if exist input level sequence
	if(flag_seq){
		C_info = CluInfo(x$data, x$density, x$MS.label, lv_seq=lv_seq, h=x$h, lv_min=x$lv.min)
		N_size = C_info$Size/max(C_info$Size)
		conn_list = C_info$Conn
	}


	plot(M_mds, xlim=c(-1,1), ylim=c(-1,1), xlab="",ylab="", axes=F, frame.plot=T)
	for(i_size in 1:length(lv_seq)){
		for(i_mode in 1:length(x$sig.clu)){
		draw.circle(x= M_mds[i_mode,1], y= M_mds[i_mode,2], radius = N_size[i_size,i_mode]* minM_dist/2*r0, border=NA, col= col_code[i_size])
			}
			for(i_tmp in 1:(length(x$sig.clu)-1)){
				for(j_tmp in (i_tmp+1):length(x$sig.clu)){
					if(conn_list[[i_size]][i_tmp, j_tmp]==1){
									lines(x=c(M_mds[i_tmp,1], M_mds[j_tmp,1]), y= c(M_mds[i_tmp,2], M_mds[j_tmp,2]), lwd=l_size*(N_size[i_size,i_tmp]+ N_size[i_size,j_tmp]), col=col_code[i_size])
	
						}
					}
				}
		}
	points(M_mds, cex=2, pch=10)
	if(flag_max=="p_max"& legend.default){
		C5 = c(1,5,9,13,17,20)/20
		col5 = colorRampPalette(c("lightgreen","limegreen","yellow","orange","red"))(6)
		legend("topleft", legend=paste(C5, "*p_max", sep=""), col=col5, lwd=rep(15,6) ,  bg="gray", cex=1)
	}
		
}

