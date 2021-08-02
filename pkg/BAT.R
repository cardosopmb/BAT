#####BAT - Biodiversity Assessment Tools
#####Version 2.7.0 (2021-08-02)
#####By Pedro Cardoso, Stefano Mammola, Francois Rigal, Jose Carlos Carvalho
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P., Rigal, F. & Carvalho, J.C. (2015) BAT - Biodiversity Assessment Tools, an R package for the measurement and estimation of alpha and beta taxon, phylogenetic and functional diversity. Methods in Ecology and Evolution, 6: 232-236.
#####Reference: Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#####Changed from v2.6.1:
#####Added functions r2, aic, and ses
#####Added function hyper.build
#####Added function kernel.hotspots
#####Added option cores to kernel.build
#####Added option comp to *.beta
#####Improved function evenness

library("ape")
library("geometry")
library("graphics")
library("hypervolume")
library("MASS")
library("nls2")
library("parallel")
library("raster")
library("stats")
library("utils")
library("vegan")
#' @import geometry
#' @import graphics
#' @import hypervolume
#' @import nls2
#' @import parallel
#' @import stats
#' @import utils
#' @import vegan
#' @importFrom ape as.phylo pcoa
#' @importFrom MASS stepAIC
#' @importFrom raster cellStats
#' @importFrom raster raster
#' @importFrom raster rasterize
#' @importFrom raster rasterToPoints

#####auxiliary functions
prep <- function(comm, xtree, abund = TRUE){
	len <- xtree[[1]] 							## length of each branch
	A <- xtree[[2]]									## matrix species X branches
	minBranch <- min(len[colSums(A)==1]) 	## minimum branch length of terminal branches
	if(is.data.frame(comm))
		comm = as.matrix(comm)
	BA <- comm%*%A 												## matrix samples X branches
	if (!abund)	BA = ifelse(BA >= 1, 1, 0)
	return (list(lenBranch = len, sampleBranch = BA, speciesBranch = A, minBranch = minBranch))
}

clean <- function(comm, tree = NA){
  if(is.vector(comm))
    comm <- matrix(comm, nrow = 1)
  comm <- as.matrix(comm)
  if (!missing(tree)){
    comm = reorderComm(comm, tree)
    tree <- xTree(tree)
  }
  return(list(comm, tree))
}

reorderComm <- function(comm, tree = NULL){
  
  if(is.vector(comm))
    comm = as.matrix(comm, nrow = 1)
  if(class(tree) == "hclust")
    tree = ape::as.phylo(tree)

  if (class(tree) == "phylo"){
    if(!is.null(tree$tip.label) && !is.null(colnames(comm))){ ##if both tree and comm have species names match and reorder species (columns) in comm
      
      #if some species are missing from comm add 0s
      if(length(tree$tip.label) > ncol(comm) && all(colnames(comm) %in% tree$tip.label)){
        miss = tree$tip.label[which(!(tree$tip.label %in% colnames(comm)))]
        addComm = matrix(0, nrow = nrow(comm), ncol = length(miss))
        colnames(addComm) = miss
        comm = cbind(comm, addComm)
      }
      if(length(dim(comm)) == 2)
        comm <- comm[,match(tree$tip.label, colnames(comm)), drop = FALSE]
      else
        comm <- comm[,match(tree$tip.label, colnames(comm)),, drop = FALSE]
      if (any(tree$tip.label != colnames(comm)))
        warning("Species names of comm and tree do not match!")
    }
  } else if (class(tree) == "dist"){
    if(!is.null(colnames(tree)) && !is.null(colnames(comm))){ ##if both tree and comm have species names match and reorder species (columns) in comm
      comm <- comm[,match(colnames(tree), colnames(comm))]
      if (any(colnames(tree) != colnames(comm)))
        warning("Species names of comm and distance do not match!")
    }
  }
  return(comm)
}

nMin <- function(comm){
	n <- sum(comm)
	for (s in 1:nrow(comm))
		n <- min(n, sum(comm[s,]))
	return(n)
}

rss <- function(x, y){
	return (sum((x-y)^2))
}

logit <- function(x){
	return(log(x/(1-x)))
}

revLogit <- function(x){
	return(exp(x)/(1+exp(x)))
}

euclid <- function(x, y){
	return(sqrt(sum((x - y) ^ 2)))
}

#####xTree function partly adapted from http://owenpetchey.staff.shef.ac.uk/Code/Code/calculatingfd_assets/Xtree.r
#####by Jens Schumacher (described in Petchey & Gaston 2002, 2006)
xTree <- function(tree) {
  if (class(tree) == "hclust"){
  	nSpp <- nrow(as.data.frame(tree['order']))
  	sppEdges <- matrix(0, nSpp, 2 * nSpp - 2)
  	lenEdges <- vector("numeric", 2 * nSpp - 2)
  	for(i in 1:(nSpp - 1)) {
  		if(tree$merge[i, 1] < 0) {
  			lenEdges[2 * i - 1] <- tree$height[order(tree$height)[i]]
  			sppEdges[ - tree$merge[i, 1], 2 * i - 1] <- 1
  		} else {
  			lenEdges[2 * i - 1] <- tree$height[order(tree$height)[i]] - tree$height[order(tree$height)[tree$merge[i, 1]]]
  			sppEdges[, 2 * i - 1] <- sppEdges[, 2 * tree$merge[i, 1] - 1] + sppEdges[ , 2 * tree$merge[i, 1]]
  		}
  		if(tree$merge[i, 2] < 0) {
  			lenEdges[2 * i] <- tree$height[order(tree$height)[i]]
  			sppEdges[ - tree$merge[i, 2], 2 * i] <- 1
  		} else {
  			lenEdges[2 * i] <- tree$height[order(tree$height)[i]] - tree$height[order(tree$height)[tree$merge[i, 2]]]
  			sppEdges[, 2 * i] <- sppEdges[, 2 * tree$merge[i, 2] - 1] + sppEdges[, 2 *tree$merge[i, 2]]
  		}
  	}
  	rownames(sppEdges) <- tree$labels
  	list(lenEdges, sppEdges)
  } else if (class(tree) == "phylo"){
    lenEdges <- tree$edge.length
    nSpp <- length(tree$tip.label)
    nEdges <- length(tree$edge.length)
    root <- nSpp + 1
    sppEdges <- matrix(0, nSpp, nEdges)
    for(i in 1:nSpp){
      find = i                                    #start by finding the ith species
      repeat{
        row = which(tree$edge[,2] == find)        #locate in which row of the edge table is our species or edge to be found
        sppEdges[i, row] = 1
        find = tree$edge[row,1]                   #find next edge if any until reaching the root
        if(find == root) break                    #all edges of this species were found, go to next species
      }
    }
    rownames(sppEdges) <- tree$tip.label
    list(lenEdges, sppEdges)
  } else {
    cat("Unrecognized tree object!")
  }
}

#####observed diversity
sobs <- function(comm, xtree){
	if (is.vector(comm))
		comm = matrix(comm, nrow = 1)
	if (missing(xtree)){
		return(length(colSums(comm)[colSums(comm) > 0]))
	} else {
		data <- prep(comm, xtree)
		value <- ifelse (colSums(data$sampleBranch) > 0, 1, 0) # vector of observed branches
		return (sum(value*data$lenBranch))
	}
}

#####observed abundance
nobs <- function(comm, xtree){
  if (is.vector(comm))
    comm = matrix(comm, nrow = 1)
  if (missing(xtree)){
    return(sum(comm))
  } else {
    data <- prep(comm, xtree)
    value <- colSums(data$sampleBranch) # vector of observed branches
    return (sum(value*data$lenBranch))
  }
}

#####hill numbers
hillobs <- function(comm, q = 0){
  #comm must be a vector
  comm = comm[comm > 0]
  comm = comm / sum(comm)    #convert to proportions
  if (q == 1)
    res = exp(-1*sum(comm * log(comm)))
  else
    res = (sum(comm^q))^(1/(1-q))
  return(res)
}

#####rao quadratic entropy
raoobs <- function(comm, distance){
  
  #comm must be a vector
  s = length(comm)
  comm = comm / sum(comm)    #convert to proportions
  distance = as.matrix(distance)
  res = 0
  for(i in 1:(s-1)){
    for(j in (i+1):s){
      res = res + (distance[i,j] * comm[i] * comm[j])
    }
  }
  return(res)
}

#####diversity of rare species for abundance - singletons, doubletons, tripletons, etc
srare <- function(comm, xtree, n = 1){
	if(missing(xtree)){
		return(length(colSums(comm)[colSums(comm) == n]))
	} else {
		data <- prep(comm, xtree)
		value <- ifelse (colSums(data$sampleBranch) == n, 1, 0) # vector of branches with given abundance
		return (sum(value*data$lenBranch))
	}
}

#####diversity of rare species for incidence - uniques, duplicates, triplicates, etc
qrare <- function(comm, xtree, n = 1){
	if(missing(xtree)){
		comm <- ifelse(comm > 0, 1, 0)
		return(length(colSums(comm)[colSums(comm) == n]))
	} else {
		data <- prep(comm, xtree, FALSE)
		value <- ifelse (colSums(data$sampleBranch) == n, 1, 0) # vector of branches with given incidence
		return (sum(value*data$lenBranch))
	}
}

#####minimum terminal branch length, = 1 in case of TD
minBranch <- function(comm, xtree){
	if (missing(xtree)){
		return(1)
	} else {
		data <- prep(comm, xtree)
		return(data$minBranch)
	}
}

#####non-parametric estimators
chao <- function(obs, s1, s2, mb){
	return(obs + (s1*(s1-mb))/(2*(s2+mb)))
}

jack1ab <- function(obs, s1){
	return(obs + s1)
}

jack1in <- function(obs, q1, q){
	return(obs + q1 * ((q-1)/q))
}

jack2ab <- function(obs, s1, s2){
	return(obs + 2*s1 - s2)
}

jack2in <- function(obs, q1, q2, q){
	if (q > 1)	return(obs + (q1*(2*q-3)/q - q2*(q-2)^2/(q*(q-1))))
	else return(obs + 2*q1 - q2)
}

pcorr <- function(obs, s1){
	return(1+(s1/obs)^2)
}

#####observed beta (a = shared species/edges, b/c = species/edges exclusive to either site, comm is a 2sites x species matrix)
betaObs <- function(comm, xtree, func = "jaccard", abund = TRUE, comp = FALSE){
  if(sum(comm) == 0)                                ##if no species on any community return 0
    return(list(Btotal = 0, Brepl = 0, Brich = 0))
	if (!abund || max(comm) == 1) {										##if incidence data
		obs1 <- sobs(comm[1,,drop=FALSE], xtree)
		obs2 <- sobs(comm[2,,drop=FALSE], xtree)
		obsBoth <- sobs(comm, xtree)
		a <- obs1 + obs2 - obsBoth
		b <- obsBoth - obs2
		c <- obsBoth - obs1
	} else if (abund & missing(xtree)){								##if abundance data
		a <- 0
		b <- 0
		c <- 0
		for (i in 1:ncol(comm)){
		  minComm <- min(comm[1,i], comm[2,i])
		  a <- a + minComm
		  b <- b + comm[1,i] - minComm
		  c <- c + comm[2,i] - minComm
		}
	} else {																					##if abundance and tree
		##due to the way Soerensen doubles the weight of the a component, using a tree or not will be the same with abundance data.
		data <- prep(comm, xtree)
		a = sum(data$lenBranch * apply(data$sampleBranch,2,min))
		diff = data$lenBranch * (data$sampleBranch[1,] - data$sampleBranch[2,])
		b = sum(replace(diff, diff < 0, 0))
		c = sum(replace(diff, diff > 0, 0) * -1)
	}
	denominator <- a + b + c
	if(tolower(substr(func, 1, 1)) == "s")
		denominator <- denominator + a
	betaValues = (list(Btotal = (b+c)/denominator, Brepl = 2*min(b,c)/denominator, Brich = abs(b-c)/denominator))
	if(comp){
	  betaValues$Shared = a
	  betaValues$Unique1 = b
	  betaValues$Unique2 = c
	}
  return(betaValues)
}

##create latitude layer
raster.lat <- function(layers){
	lat <- layers[[1]]
	x <- rasterToPoints(lat)[,1:2]
	lat <- rasterize(x[,1:2], lat, x[,2])
	names(lat) <- "latitude"
	return(lat)
}

##create longitude layer
raster.long <- function(layers){
	long <- layers[[1]]
	x <- rasterToPoints(long)[,1:2]
	long <- rasterize(x[,1:2], long, x[,1])
	names(long) <- "longitude"
	return(long)
}

##################################################################################
##################################MAIN FUNCTIONS##################################
##################################################################################

#' Alpha diversity (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Observed richness with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details TD is equivalent to species richness. Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' The rarefaction option is useful to compare communities with much different numbers of individuals sampled, which might bias diversity comparisons (Gotelli & Colwell 2001)
#' @return A matrix of sites x diversity values (either "Richness" OR "Mean, Median, Min, LowerCL, UpperCL and Max").
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology Letters, 4, 379-391.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples comm <- matrix(c(0,0,1,1,0,0,2,1,0,0), nrow = 2, ncol = 5, byrow = TRUE)
#' trait = 1:5
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha(comm)
#' alpha(comm, raref = 0)
#' alpha(comm, trait)
#' alpha(comm, tree)
#' alpha(comm, tree, 2, 100)
#' @export
alpha <- function(comm, tree, raref = 0, runs = 100){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

  #first organize the data
  if(!missing(tree)){
    cleanData = clean(comm, tree)
    comm = cleanData[[1]]
    tree = cleanData[[2]]
  }

  #now let's go for what matters
  nComm <- nrow(comm)
	if(raref < 1){						# no rarefaction if 0 or negative
		results <- matrix(0, nComm, 1)
		for (s in 1:nComm){
			results[s,1] <- sobs(comm[s,, drop=FALSE], tree)
		}
		rownames(results) <- rownames(comm)
		colnames(results) <- "Richness"
		return (results)
	}
	if (raref == 1)
		raref <- nMin(comm)				# rarefy by minimum n among all communities
	results <- matrix(0, nComm, 6)
	for (s in 1:nComm){
		res <- c()
		for (r in 1:runs){
			res <- c(res,sobs(rrarefy(comm[s,], raref), tree))
		}
		results[s,] <- c(mean(res), quantile(res, 0.5), min(res), quantile(res, 0.025), quantile(res, 0.975), max(res))
	}
	rownames(results) <- rownames(comm)
	colnames(results) <- c("Mean", "Median", "Min", "LowerCL", "UpperCL", "Max")
	return (results)
}

#' Alpha diversity accumulation curves (observed and estimated).
#' @description Estimation of alpha diversity of a single site with accumulation of sampling units.
#' @param comm A sampling units x species matrix, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param func The class of estimators to be used:
#' If func is partial match of "curve", TD, PD or FD are based on extrapolating the accumulation curve of observed diversity.
#' If func is partial match of "nonparametric", TD, PD or FD are based on non-parametric estimators.
#' If func is partial match of "completeness", PD or FD estimates are based on the completeness of TD (requires a tree to be used).
#' If not specified, default is "nonparametric.
#' @param target True diversity value to calculate the accuracy of curves (scaled mean squared error). If not specified do not calculate accuracy (default), -1 uses the total observed diversity as true diversity and any other value is the true known diversity.
#' @param runs Number of random permutations to be made to the sampling order. If not specified, default is 100.
#' @param prog Present a text progress bar in the R console.
#' @details Observed diversity often is an underestimation of true diversity. Several approaches have been devised to estimate species richness (TD) from incomplete sampling.
#' These include: (1) fitting asymptotic functions to randomised accumulation curves (Soberon & Llorente 1993; Flather 1996; Cardoso et al. in prep.)
#' (2) the use of non-parametric estimators based on the incidence or abundance of rare species (Heltshe & Forrester 1983; Chao 1984, 1987; Colwell & Coddington 1994).
#' A correction to non-parametric estimators has also been recently proposed, based on the proportion of singleton or unique species
#' (species represented by a single individual or in a single sampling unit respectively; Lopez et al. 2012).
#' Cardoso et al. (2014) have proposed a way of adapting these approaches to estimate PD and FD, also adding a third possible approach for
#' these dimensions of diversity: (3) correct PD and FD values based on the completeness of TD, where completeness equals the proportion of estimated true diversity that was observed.
#' Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of sampling units x diversity values (sampling units, individuals, observed and estimated diversity).
#' The values provided by this function are:
#' @return Sampl - Number of sampling units;
#' @return Ind - Number of individuals;
#' @return Obs - Observed diversity;
#' @return S1 - Singletons;
#' @return S2 - Doubletons;
#' @return Q1 - Uniques;
#' @return Q2 - Duplicates;
#' @return Jack1ab - First order jackknife estimator for abundance data;
#' @return Jack1in - First order jackknife estimator for incidence data;
#' @return Jack2ab - Second order jackknife estimator for abundance data;
#' @return Jack2in - Second order jackknife estimator for incidence data;
#' @return Chao1 - Chao estimator for abundance data;
#' @return Chao2 - Chao estimator for incidence data;
#' @return Clench - Clench or Michaelis-Menten curve;
#' @return Exponential - Exponential curve;
#' @return Rational - Rational function;
#' @return Weibull - Weibull curve;
#' @return The P-corrected version of all non-parametric estimators is also provided.
#' @return Accuracy - if accuracy is to be calculated a list is returned instead, with the second element being the scaled mean squared error of each estimator.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, 5: 452-461.
#' @references Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
#' @references Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' @references Colwell, R.K. & Coddington, J.A. (1994) Estimating terrestrial biodiversity through extrapolation. Phil. Trans. Roy. Soc. London B 345, 101-118.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Flather, C. (1996) Fitting species-accumulation functions and assessing regional land use impacts on avian diversity. Journal of Biogeography, 23, 155-168.
#' @references Heltshe, J. & Forrester, N.E. (1983) Estimating species richness using the jackknife procedure. Biometrics, 39, 1-11.
#' @references Lopez, L.C.S., Fracasso, M.P.A., Mesquita, D.O., Palma, A.R.T. & Riul, P. (2012) The relationship between percentage of singletons and sampling effort: a new approach to reduce the bias of richness estimates. Ecological Indicators, 14, 164-169.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @references Soberon, M.J. & Llorente, J. (1993) The use of species accumulation functions for the prediction of species richness. Conservation Biology, 7, 480-488.
#' @examples comm <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha.accum(comm)
#' alpha.accum(comm, func = "nonparametric")
#' alpha.accum(comm, tree, "completeness")
#' alpha.accum(comm, tree, "curve", runs = 1000)
#' alpha.accum(comm, target = -1)
#' @export
alpha.accum <- function(comm, tree, func = "nonparametric", target = -2, runs = 100, prog = TRUE){

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  #first organize the data
  if(!missing(tree)){
    cleanData = clean(comm, tree)
    comm = cleanData[[1]]
    tree = cleanData[[2]]
  }

  #####function options:
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	#####completeness (PD/FD with TD completeness correction)
	#####curve (TD/PD/FD with curve fitting)
	func <- match.arg(func, c("nonparametric", "completeness", "curve"))

	#####nonparametric (TD/PD/FD with non-parametric estimators)
	switch(func, nonparametric = {
		resultsArray <- array(0, dim = c(nrow(comm), 19, runs))
		if(target > -2){
		  smse <- matrix(0, runs, 19)
		  smsew <-  smse
		}
		if (prog) pb <- txtProgressBar(0, runs, style = 3)
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]			#shuffle rows (sampling units)
			data <- matrix(0,1,ncol(comm))
			runData <- matrix(0,nrow(comm),19)
			colnames(data) = colnames(comm)
			for (q in 1:nrow(comm)){
				data <- rbind(data, comm[q,])
				n <- sum(rowSums(data))
				obs <- sobs(data, tree)
				s1 <- srare(data, tree, 1)
				s2 <- srare(data, tree, 2)
				q1 <- qrare(data, tree, 1)
				q2 <- qrare(data, tree, 2)
				mb <- minBranch(data, tree)
				j1ab <- jack1ab(obs, s1)
				j1abP <- j1ab * pcorr(obs, s1)
				j1in <- jack1in(obs, q1, q)
				j1inP <- j1in * pcorr(obs, q1)
				j2ab <- jack2ab(obs, s1, s2)
				j2abP <- j2ab * pcorr(obs, s1)
				j2in <- jack2in(obs, q1, q2, q)
				j2inP <- j2in * pcorr(obs, q1)
				c1 <- chao(obs, s1, s2, mb)
				c1P <- c1 * pcorr(obs, s1)
				c2 <- chao(obs, q1, q2, mb)
				c2P <- c2 * pcorr(obs, q1)
				runData[q,] <- c(q, n, obs, s1, s2, q1, q2, j1ab, j1abP, j1in, j1inP, j2ab, j2abP, j2in, j2inP, c1, c1P, c2, c2P)
			}
			resultsArray[,,r] <- runData
			if(exists("smse")){					##if accuracy is to be calculated
				if(r == 1){
					if(target == -1){
						truediv <- runData[nrow(runData),3]
					}else{
						truediv <- target
					}
				}
				s <- accuracy(runData, truediv)
				smse[r,3] <- s[1,1]
				smse[r,8:19] <- s[1,-1]
				smsew[r,3] <- s[2,1]
				smsew[r,8:19] <- s[2,-1]
			}
			if (prog) setTxtProgressBar(pb, r)
		}
		if (prog) close(pb)

		#####calculate averages or medians of all runs
		results <- matrix(0,nrow(comm),19)
		v <- array(0, dim = c(runs))
		for (i in 1:nrow(comm)){
			for (j in 1:19){
				for (k in 1:runs){
					v[k] <- resultsArray[i,j,k]
				}
				if (j < 16 || missing(tree))
					results[i,j] <- mean(v)
				else
					results[i,j] <- median(v)
			}
		}
		if(exists("smse")){						##calculate accuracy
			smse <- colMeans(smse)
			smsew <- colMeans(smsew)
		}

		#####completeness (PD/FD with TD completeness correction)
	}, completeness = {
		if (missing(tree))
			stop("Completeness option not available without a tree...")
		results <- alpha.accum(comm, runs = runs)
		obs <- matrix(0,nrow(comm),1)
		for (r in 1:runs){
			comm <- comm[sample(nrow(comm)),, drop=FALSE]			#shuffle rows (sampling units)
			for (s in 1:nrow(comm)){
				obs[s,1] <- obs[s,1] + sobs(comm[1:s,], tree)
			}
		}
		obs <- obs / runs
		for (i in 8:19)
			results[,i] <- obs * (results[,i] / results[,3])
		results[,3] <- obs

		#####curve (TD/PD/FD with curve fitting)
	}, curve = {
    results <- matrix(NA,nrow(comm),7)
    results[,1] <- seq(1,nrow(comm))  ##fill samples column
    results[,2] <- seq(sum(comm)/nrow(comm),sum(comm), sum(comm)/nrow(comm))  ##fill individuals column
    runObs <- rep(0,nrow(comm))
    if (prog) pb <- txtProgressBar(0, runs, style = 3)
    for (r in 1:runs){
 		  comm <- comm[sample(nrow(comm)),, drop=FALSE]		#shuffle rows (sampling units)
 		  for (s in 1:nrow(comm)){
 		    runObs[s] <- runObs[s] + sobs(comm[1:s,,drop=FALSE], tree)
 		  }
 		  if (prog) setTxtProgressBar(pb, r)
    }
    if (prog) close(pb)
		results[,3] <- runObs / runs

    rich <- results[nrow(comm),3]

		for (s in 3:nrow(results)){				##fit curves only with 3 or more sampling units
		  ## curve fitting
      x <- results[1:s,1]
			y <- results[1:s,3]
			##Clench
			stlist <- data.frame(a = rich, b = c(0.1, 0.5, 1))
			form <- y ~ (a*x)/(b+x)
			mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
			if(class(curve) != "try-error"){
				a <- coef(curve)[1]
				results[s,5] <- a
			}
			##Negative exponential
      form <- y ~ a*(1-exp(-b*x))
			mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
			if(class(curve) != "try-error"){
				a <- coef(curve)[1]
				results[s,6] <- a
			}
			##Rational
			stlist <- data.frame(a = rich, b = c(0.1, 0.5, 1, 5, 10), c = c(1, 10, 100, 1000, 10000))
			form <- y ~ (c+(a*x))/(b+x)
			mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
			if(class(curve) != "try-error"){
				a <- coef(curve)[1]
				results[s,4] <- a
			}
			##Weibull
 			stlist <- data.frame(a = rich, b = c(0,1,10), c = c(0,0.1,1))
 			form <- y ~ a*(1-exp(-b*(x^c)))
      mod <- try(nls2(form, start = stlist, algorithm = "random-search"), silent = TRUE)
 			curve <- try(nls2(form, start = mod, algorithm = "default"), silent = TRUE)
  			if(class(curve) != "try-error"){
  				a <- coef(curve)[1]
  				results[s,7] <- a
  			}
		}
		colnames(results) <- c("Sampl", "Ind", "Obs", "Clench", "Exponential", "Rational", "Weibull")
		return (results)
	})
	colnames(results) <- c("Sampl", "Ind", "Obs", "S1", "S2", "Q1", "Q2", "Jack1ab", "Jack1abP", "Jack1in", "Jack1inP", "Jack2ab", "Jack2abP", "Jack2in", "Jack2inP", "Chao1", "Chao1P", "Chao2", "Chao2P")
	if(exists("smse")){
		smse <- rbind(smse, smsew)
		colnames(smse) <- colnames(results)
		rownames(smse) <- c("Raw", "Weighted")
		smse <- smse[,-c(1:2,4:7)]
		return(list(results, smse))
	}	else {
		return(results)
	}
}

#' Alpha diversity estimates.
#' @description Estimation of alpha diversity of multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundances or number of incidences.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param func The class of estimators to be used:
#' If func is partial match of "nonparametric", TD, PD or FD are based on non-parametric estimators.
#' If func is partial match of "completeness", PD or FD estimates are based on the completeness of TD (requires a tree to be used).
#' If not specified, default is "nonparametric".
#' @details Observed diversity often is an underestimation of true diversity.
#' Non-parametric estimators based on the incidence or abundance of rare species have been proposed to overcome the problem of undersampling (Heltshe & Forrester 1983; Chao 1984, 1987; Colwell & Coddington 1994).
#' A correction to non-parametric estimators has also been recently proposed, based on the proportion (P) of singleton or unique species
#' (species represented by a single individual or in a single sampling unit respectively; Lopez et al. 2012).
#' Cardoso et al. (2014) have proposed a way of adapting non-parametric species richness estimators to PD and FD. They have also proposed correcting PD and FD values based on the completeness of TD, where completeness equals the proportion of estimated true diversity that was observed.
#' Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of sites x diversity values (individuals, observed and estimated diversity).
#' The values provided by this function are:
#' @return Ind - Number of individuals;
#' @return Obs - Observed diversity;
#' @return S1 - Singletons;
#' @return S2 - Doubletons;
#' @return Jack1ab - First order jackknife estimator for abundance data;
#' @return Jack2ab - Second order jackknife estimator for abundance data;
#' @return Chao1 - Chao estimator for abundance data.
#' @return The P-corrected version of all estimators is also provided.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, 5: 452-461.
#' @references Chao, A. (1984) Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics, 11, 265-270.
#' @references Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' @references Colwell, R.K. & Coddington, J.A. (1994) Estimating terrestrial biodiversity through extrapolation. Phil. Trans. Roy. Soc. London B 345, 101-118.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Heltshe, J. & Forrester, N.E. (1983) Estimating species richness using the jackknife procedure. Biometrics, 39, 1-11.
#' @references Lopez, L.C.S., Fracasso, M.P.A., Mesquita, D.O., Palma, A.R.T. & Riul, P. (2012) The relationship between percentage of singletons and sampling effort: a new approach to reduce the bias of richness estimates. Ecological Indicators, 14, 164-169.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples comm <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' alpha.estimate(comm)
#' alpha.estimate(comm, tree)
#' alpha.estimate(comm, tree, func = "completeness")
#' @export
alpha.estimate <- function(comm, tree, func = "nonparametric"){

	if (max(comm) == 1)
		stop("No estimates are possible without abundance or incidence frequency data")

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
	#####function options:
	#####nonparametric (TD/PD/FD with non-parametric estimators)
	#####completeness (PD/FD with TD completeness correction)
	func <- match.arg(func, c("nonparametric", "completeness"))

	#####nonparametric (TD/PD/FD with non-parametric estimators)
	switch(func, nonparametric = {
	  #first organize the data
	  if(!missing(tree)){
	    cleanData = clean(comm, tree)
	    comm = cleanData[[1]]
	    tree = cleanData[[2]]
	  }
	  
	  #now let's go for what matters
	  results <- matrix(0,0,10)
		for (s in 1:nrow(comm)){
			data <- comm[s,,drop = FALSE]
			obs <- sobs(data, tree)
			n <- sum(data)
			s1 <- srare(data, tree, 1)
			s2 <- srare(data, tree, 2)
			mb <- minBranch(data, tree)
			j1ab <- jack1ab(obs, s1)
			j1abP <- j1ab * pcorr(obs, s1)
			j2ab <- jack2ab(obs, s1, s2)
			j2abP <- j2ab * pcorr(obs, s1)
			c1 <- chao(obs, s1, s2, mb)
			c1P <- c1 * pcorr(obs, s1)
			results <- rbind(results, c(n, obs, s1, s2, j1ab, j1abP, j2ab, j2abP, c1, c1P))
		}

		#####completeness (PD/FD with TD completeness correction)
	}, completeness = {
	  if(is.vector(comm))
	    comm <- matrix(comm, nrow = 1)
	  comm <- as.matrix(comm)
		if (missing(tree))
			stop("Completeness option not available without a tree...")
		results <- alpha.estimate(comm, tree, "nonparametric")
	  if(!is.null(tree$labels) && !is.null(colnames(comm))) ##if both tree and comm have species names match and reorder species (columns) in comm
	    comm <- comm[,match(tree$labels, colnames(comm))]
	  tree <- xTree(tree)
		obs <- matrix(0,nrow(comm),1)
		for (s in 1:nrow(comm))
			obs[s,1] <- obs[s,1] + sobs(comm[s,], tree)
		for (i in 5:10)
			results[,i] <- obs[,1] * (results[,i] / results[,2])
		results[,2] <- obs[,1]
	})
	rownames(results) <- rownames(comm)
	colnames(results) <- c("Ind", "Obs", "S1", "S2", "Jack1ab", "Jack1abP", "Jack2ab", "Jack2abP", "Chao1", "Chao1P")
	return(results)
}

#' Rao quadratic entropy.
#' @description Rao quadratic entropy for Phylogenetic or Functional richness.
#' @param comm A sites x species matrix, with abundance data.
#' @param tree A phylo or hclust object (used only for PD or FD).
#' @param distance A dist object representing the phylogenetic or functional distance between species or alternatively a species x traits matrix or data.frame to calculate distances.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details Rao quadratic entropy (Rao, 1982) measures diversity based on the abundance of species and the dissimilarity between them.
#' @return A matrix of sites x diversity values (either "Rao" OR "Mean, Median, Min, LowerCL, UpperCL and Max").
#' @references Rao, C.R. (1982). Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21: 24-43.
#' @examples comm <- matrix(c(0,0,1,1,0,0,100,1,2,0), nrow = 2, ncol = 5, byrow = TRUE)
#' distance = dist(1:5)
#' rao(comm)
#' rao(comm, , distance)
#' rao(comm, hclust(distance), raref = 1)
#' @export
rao <- function(comm, tree, distance, raref = 0, runs = 100){
  
  #convert traits to distance if needed
  if(!missing(distance) && (is.matrix(distance) || is.data.frame(distance) || is.vector(distance)))
    distance = gower(distance)

  if(!missing(tree))
    distance = cophenetic(tree)
  else if(missing(distance))
    distance = as.dist(matrix(1, ncol(comm), ncol(comm)))
  
  nComm <- nrow(comm)
  if(raref < 1){						# no rarefaction if 0 or negative
    results <- matrix(0, nComm, 1)
    for (s in 1:nComm){
      results[s,1] <- raoobs(comm[s,], distance)
    }
    rownames(results) <- rownames(comm)
    colnames(results) <- "Rao"
    return (results)
  }
  if (raref == 1)
    raref <- nMin(comm)				# rarefy by minimum n among all communities
  results <- matrix(0, nComm, 6)
  for (s in 1:nComm){
    res <- c()
    for (r in 1:runs){
      rarefiedComm = rrarefy(comm[s,], raref)
      res <- c(res, raoobs(rarefiedComm, distance))
    }
    results[s,] <- c(mean(res), quantile(res, 0.5), min(res), quantile(res, 0.025), quantile(res, 0.975), max(res))
  }
  rownames(results) <- rownames(comm)
  colnames(results) <- c("Mean", "Median", "Min", "LowerCL", "UpperCL", "Max")
  return (results)
}


#' Hill numbers.
#' @description Hill numbers with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with abundance data.
#' @param q Hill number order: q(0) = species richness, q(1) ~ Shannon diversity, q(2) ~ Simpson diversity.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details Hill numbers are based on the number of equally abundant species  that would match the current diversity.
#' Depending on the single parameter they give more or less weight to rare species (Jost 2002).
#' @return A matrix of sites x diversity values (either "Hill q" OR "Mean, Median, Min, LowerCL, UpperCL and Max").
#' @references Hill, M.O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54: 427-432.
#' @examples comm <- matrix(c(0,0,1,1,0,0,100,1,0,0), nrow = 2, ncol = 5, byrow = TRUE)
#' hill(comm)
#' hill(comm, q = 1)
#' hill(comm, q = 2, 1)
#' @export
hill <- function(comm, q = 0, raref = 0, runs = 100){
  
  nComm <- nrow(comm)
  if(raref < 1){						# no rarefaction if 0 or negative
    results <- matrix(0, nComm, 1)
    for (s in 1:nComm){
      results[s,1] <- hillobs(comm[s,], q)
    }
    rownames(results) <- rownames(comm)
    colnames(results) <- paste("Hill", q)
    return (results)
  }
  if (raref == 1)
    raref <- nMin(comm)				# rarefy by minimum n among all communities
  results <- matrix(0, nComm, 6)
  for (s in 1:nComm){
    res <- c()
    for (r in 1:runs){
      rarefiedComm = rrarefy(comm[s,], raref)
      res <- c(res, hillobs(rarefiedComm, q))
    }
    results[s,] <- c(mean(res), quantile(res, 0.5), min(res), quantile(res, 0.025), quantile(res, 0.975), max(res))
  }
  rownames(results) <- rownames(comm)
  colnames(results) <- c("Mean", "Median", "Min", "LowerCL", "UpperCL", "Max")
  return (results)
}

#' Beta diversity (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Beta diversity with possible rarefaction, multiple sites simultaneously.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param func Partial match indicating whether the Jaccard or Soerensen family of beta diversity measures should be used.  If not specified, default is Jaccard.
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @param comp Boolean indicating whether beta diversity components (shared and unique fractions) should be returned
#' @details The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' The rarefaction option is useful to compare communities with much different numbers of individuals sampled, which might bias diversity comparisons (Gotelli & Colwell 2001)
#' @return Three distance matrices between sites, one per each of the three beta diversity measures (either "Obs" OR "Mean, Median, Min, LowerCL, UpperCL and Max").
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology Letters, 4, 379-391.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta(comm)
#' beta(comm, abund = FALSE, comp = TRUE)
#' beta(comm, tree)
#' beta(comm, raref = 1)
#' beta(comm, tree, "s", abund = FALSE, raref = 2)
#' @export
beta <- function(comm, tree, func = "jaccard", abund = TRUE, raref = 0, runs = 100, comp = FALSE){

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

  #first organize the data
  if(!missing(tree)){
    cleanData = clean(comm, tree)
    comm = cleanData[[1]]
    tree = cleanData[[2]]
  }
  
  #now let's go for what matters
  nComm <- nrow(comm)

	if(raref < 1){						# no rarefaction if 0 or negative
		results <- array(0, dim=c(nComm, nComm, 3))
		rownames(results) = colnames(results) = rownames(comm)
		comps <- results
		for (i in 1:(nComm-1)){
			for (j in (i+1):nComm){
				commBoth <- as.matrix(rbind(comm[i,], comm[j,]))
				betaValues <- betaObs(commBoth, tree, func, abund, comp)
				results[j,i,] <- unlist(betaValues)[1:3]
				if(comp)
				  comps[j,i,] <- unlist(betaValues)[4:6]
			}
		}
		results <- list(Btotal = as.dist(results[,,1]),Brepl = as.dist(results[,,2]),Brich = as.dist(results[,,3]))
		if (comp){
		  results$Shared = round(as.dist(comps[,,1]),3)
		  results$Unique1 = round(as.dist(comps[,,2]),3)
		  results$Unique2 = round(as.dist(comps[,,3]),3)
		}
		return(results)
	}
	if (raref == 1)
		raref <- nMin(comm)				# rarefy by minimum n among all communities
	results <- array(0, dim=c(nComm, nComm, 3, 6))

	for (i in 1:(nComm-1)){
		for (j in (i+1):nComm){
			run <- matrix(0, runs, 3)
			for (r in 1:runs){
				commBoth <- as.matrix(rbind(rrarefy(comm[i,], raref), rrarefy(comm[j,], raref)))
				betaValues <- betaObs(commBoth, tree, func, abund)
				run[r,1] <- betaValues$Btotal
				run[r,2] <- betaValues$Brepl
				run[r,3] <- betaValues$Brich
			}
			for (b in 1:3){
			  results[j,i,b,1] <- mean(run[,b])
			  results[j,i,b,2] <- quantile(run[,b], 0.5)
				results[j,i,b,3] <- min(run[,b])
				results[j,i,b,4] <- quantile(run[,b], 0.025)
				results[j,i,b,5] <- quantile(run[,b], 0.975)
				results[j,i,b,6] <- max(run[,b])
			}
		}
	}
	results.total <- list(Btotal.mean = as.dist(results[,,1,1]), Btotal.median = as.dist(results[,,1,2]), Btotal.min = as.dist(results[,,1,3]), Btotal.lowCL = as.dist(results[,,1,4]), Btotal.upCL = as.dist(results[,,1,5]), Btotal.max = as.dist(results[,,1,6]))
	results.repl <- list(Brepl.mean = as.dist(results[,,2,1]), Brepl.median = as.dist(results[,,2,2]), Brepl.min = as.dist(results[,,2,3]), Brepl.lowCL = as.dist(results[,,2,4]), Brepl.upCL = as.dist(results[,,2,5]), Brepl.max = as.dist(results[,,2,6]))
	results.rich <- list(Brich.mean = as.dist(results[,,3,1]), Brich.median = as.dist(results[,,3,2]), Brich.min = as.dist(results[,,3,3]), Brich.lowCL = as.dist(results[,,3,4]), Brich.upCL = as.dist(results[,,3,5]), Brich.max = as.dist(results[,,3,6]))
	results <- c(results.total, results.repl, results.rich)
	return (results)
}

#' Beta diversity accumulation curves.
#' @description Beta diversity between two sites with accumulation of sampling units.
#' @param comm1 A sampling units x species matrix for the first site, with either abundance or incidence data.
#' @param comm2 A sampling units x species matrix for the second site, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param func Partial match indicating whether the Jaccard or Soerensen family of beta diversity measures should be used. If not specified, default is jaccard.
#' @param abund A boolean (T/F) indicating whether abundance data should be used or converted to incidence before analysis.
#' @param runs Number of random permutations to be made to the sampling order. If not specified, default is 100.
#' @param prog Present a text progress bar in the R console.
#' @details As widely recognized for species richness, beta diversity is also biased when communities are undersampled.
#' Beta diversity accumulation curves have been proposed by Cardoso et al. (2009) to test if beta diversity has approached an asymptote when comparing two undersampled sites.
#' The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone;
#' Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm1 and comm2 must be the same as in tree. Also, the number of sampling units should be similar in both sites.
#' @return Three matrices of sampling units x diversity values, one per each of the three beta diversity measures (sampling units, individuals and observed diversity).
#' @references Cardoso, P., Borges, P.A.V. & Veech, J.A. (2009) Testing the performance of beta diversity measures based on incidence data: the robustness to undersampling. Diversity and Distributions, 15, 1081-1090.
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta.accum(comm1, comm2)
#' beta.accum(comm1, comm2, func = "Soerensen")
#' beta.accum(comm1, comm2, tree)
#' beta.accum(comm1, comm2, abund = FALSE)
#' beta.accum(comm1, comm2, tree,, FALSE)
#' @export
beta.accum <- function(comm1, comm2, tree, func = "jaccard", abund = TRUE, runs = 100, prog = TRUE){

  if(nrow(comm1) < 2 || nrow(comm1) != nrow(comm2))
		stop("Both communities should have multiple and the same number of sampling units")
  comm1 <- as.matrix(comm1)
  comm2 <- as.matrix(comm2)
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

  #first organize the data
  if(!missing(tree)){
    cleanData = clean(comm1, tree)
    comm1 = cleanData[[1]]
    cleanData = clean(comm2, tree)
    comm2 = cleanData[[1]]
    tree = cleanData[[2]]
  }
  
  #now let's go for what matters
	nSamples <- nrow(comm1)
	results <- matrix(0,nSamples, 4)
	colnames(results) <- c("Sampl", "Btotal", "Brepl", "Brich")
	if (prog) pb <- txtProgressBar(0, runs, style = 3)
	for (r in 1:runs){
		comm1 <- comm1[sample(nSamples),, drop=FALSE]			#shuffle sampling units of first community
		comm2 <- comm2[sample(nSamples),, drop=FALSE]			#shuffle sampling units of second community
		for (q in 1:nSamples){
			commBoth <- as.matrix(rbind(colSums(comm1[1:q,,drop=FALSE]),colSums(comm2[1:q,,drop=FALSE])))
			results[q,1] <- results[q,1] + q
			betaValues <- betaObs(commBoth, tree, func, abund)
			results[q,2] <- results[q,2] + betaValues$Btotal
			results[q,3] <- results[q,3] + betaValues$Brepl
			results[q,4] <- results[q,4] + betaValues$Brich
		}
		if (prog) setTxtProgressBar(pb, r)
	}
	if (prog) close(pb)
	results <- results/runs
	return(results)
}

#' Beta diversity among multiple communities.
#' @description Beta diversity with possible rarefaction - multiple sites measure calculated as the average or variance of all pairwise values.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param func Indicates whether the Jaccard or Soerensen family of beta diversity measures should be used. If not specified, default is jaccard.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details Beta diversity of multiple sites simultaneously is calculated as either the average or the variance among all pairwise comparisons (Legendre, 2014).
#' The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone;
#' Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in comm must be the same as in tree.
#' @return A matrix of beta measures x diversity values (average and variance).
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Legendre, P. (2014) Interpreting the replacement and richness difference components of beta diversity. Global Ecology and Biogeography, 23: 1324-1334.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples comm <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' beta.multi(comm)
#' beta.multi(comm, func = "Soerensen")
#' beta.multi(comm, tree)
#' beta.multi(comm, raref = 1)
#' beta.multi(comm, tree, "s", FALSE, raref = 2)
#' @export
beta.multi <- function(comm, tree, func = "jaccard", abund = TRUE, raref = 0, runs = 100){
	pairwise <- beta(comm, tree, func, abund, raref, runs)
	Btotal.avg <- mean(pairwise$Btotal)
	Brepl.avg <- mean(pairwise$Brepl)
	Brich.avg <- mean(pairwise$Brich)
	Btotal.var <- sum(pairwise$Btotal)/(ncol(comm)*(ncol(comm)-1))
	Brepl.var <- sum(pairwise$Brepl)/(ncol(comm)*(ncol(comm)-1))
	Brich.var <- sum(pairwise$Brich)/(ncol(comm)*(ncol(comm)-1))
	results <- matrix(c(Btotal.avg, Brepl.avg, Brich.avg, Btotal.var, Brepl.var, Brich.var), nrow = 3, ncol = 2)
	colnames(results) <- c("Average", "Variance")
	rownames(results) <- c("Btotal", "Brepl", "Brich")
	return(results)
}

#' Beta diversity evenness (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Difference of evenness between pairs of sites.
#' @param comm A sites x species matrix, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist or matrix object representing the phylogenetic or functional distance between species. If both tree and distance are missing, taxonomic evenness is calculated.
#' @param method Calculate evenness using "expected" values (default) or values based on "contribution" of species to the tree.
#' @param func Calculate evenness using "Camargo" (default) or "Bulla" index.
#' @param abund A boolean (T/F) indicating whether evenness should be calculated using abundance data.
#' @details This measure is simply the pairwise difference of evenness calculated based on the index of Camargo (1993) or Bulla (1994) using the values of both species abundances and edge lengths in the tree (if PD/FD).
#' @details If no tree or distance is provided the result is the original index.
#' @return Distance matrix between sites.
#' @references Bulla, L. (1994) An index of evenness and its associated diversity measure. Oikos, 70: 167-171.
#' @references Camargo, J.A. (1993) Must dominance increase with the number of subordinate species in competitive interactions? Journal of Theoretical Biology, 161: 537-542.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,1,1,1,1,100), nrow = 4, byrow = TRUE)
#' distance <- dist(c(1:5), method = "euclidean")
#' tree <- hclust(distance, method = "average")
#' beta.evenness(comm)
#' beta.evenness(comm, tree)
#' beta.evenness(comm, tree, method = "contribution")
#' beta.evenness(comm, tree, abund = FALSE)
#' @export
beta.evenness <- function(comm, tree, distance, method = "expected", func = "camargo", abund = TRUE){
  return(dist(evenness(comm, tree, distance, method, func, abund)))
}

#' Phylogenetic/functional originality of species or individuals.
#' @description Average dissimilarity between a species or individual and all others in a community.
#' @param comm A sites x species matrix, with either abundance or incidence data. If missing, the originality using the full tree or distance matrix is calculated.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist object representing the phylogenetic or functional distance between species.
#' @param abund A boolean (T/F) indicating whether originality should be calculated per individual (T) or species (F).
#' @param relative A boolean (T/F) indicating whether originality should be relative to the maximum distance between any two species in the tree or distance matrix.
#' @details This is the originality measure of Pavoine et al. (2005) without replacement.
#' @return A matrix of sites x species values.
#' @references Pavoine, S., Ollier, S. & Dufour, A.-B. (2005) Is the originality of a species measurable? Ecology Letters, 8: 579-586.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,1,1,1), nrow = 4, byrow = TRUE)
#' distance <- dist(c(1:5), method="euclidean")
#' tree <- tree.build(1:5)
#' originality(tree = tree)
#' originality(distance = distance)
#' originality(comm, tree)
#' originality(comm, tree, abund = FALSE)
#' originality(comm, tree, abund = FALSE, relative = FALSE)
#' @export
originality <- function(comm, tree, distance, abund = TRUE, relative = TRUE){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  if(missing(comm)){
		if(!missing(distance))
			comm = rep(1, attributes(distance)$Size)
		if(!missing(tree))
			comm = rep(1, length(tree$order))
	}
	if(is.vector(comm))
		comm <- matrix(comm, nrow = 1)
	
	if(!missing(tree)){
	  comm = reorderComm(comm, tree)
	  distance <- cophenetic(tree)			     								#cophenetic distances of species
	}else if(missing(distance)){
		return(warning("Need one of tree or distance!"))
	}
	distance <- as.matrix(distance)													#convert distance to matrix
  if(!abund){
  	comm <- ifelse(comm > 0, 1, 0)
  	for(i in 1:nrow(distance))
  		distance[i,i] = NA
  }

  original <- matrix(NA,nrow(comm),ncol(comm))
  for (r in 1:nrow(comm)){    					                  #cycle through all sites/samples
    present <- which(comm[r,]>0)                          #which species exist in this site
    nSpp <- length(present)                               #how many species are present in this site
    proportion <- comm[r,present]/sum(comm[r,present])                  #proportion incidence/abundance of species in this site
    for (c in present){
    	original[r,c] <- sum(distance[present,c] * proportion, na.rm=T)
    	if(!abund)
    		original[r,c] <- original[r,c] * nSpp / (nSpp-1) #correct not to take distance to self into account
    }
  }
  if(relative)
    original <- original / max(distance, na.rm=T)
	return(original)
}

#' Phylogenetic/functional uniqueness of species.
#' @description Dissimilarity between each species and the single closest in a community.
#' @param comm A sites x species matrix, with either abundance or incidence data. If missing, the uniqueness using the full tree or distance matrix is calculated.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist object representing the phylogenetic or functional distance between species.
#' @param relative A boolean (T/F) indicating whether uniqueness should be relative to the maximum distance between any two species in the tree or distance matrix.
#' @details This is equivalent to the originality measure of Mouillot et al. (2013).
#' @return A matrix of sites x species values.
#' @references Mouillot, D., Graham, N.A., Villeger, S., Mason, N.W. & Bellwood, D.R. (2013) A functional approach reveals community responses to disturbances. Trends in Ecology and Evolution, 28: 167-177.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,1,1,1), nrow = 4, byrow = TRUE)
#' distance <- dist(c(1:5), method="euclidean")
#' tree <- hclust(distance, method="average")
#' uniqueness(tree = tree)
#' uniqueness(distance = distance)
#' uniqueness(comm, tree)
#' @export
uniqueness <- function(comm, tree, distance, relative = TRUE){
	
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  if(missing(comm)){
		if(!missing(distance))
			comm = rep(1, attributes(distance)[1])
		if(!missing(tree))
			comm = rep(1, length(tree$order))
	}
	if(is.vector(comm))
		comm <- matrix(comm, nrow = 1)
	
	if(!missing(tree)){
	  comm = reorderComm(comm, tree)
	  distance <- cophenetic(tree)			     								#cophenetic distances of species
	}else if(missing(distance)){
		return(warning("Need one of tree or distance!"))
	}
	distance <- as.matrix(distance)													#convert distance to matrix

	comm <- ifelse(comm > 0, 1, 0)
	for(i in 1:nrow(distance))
		distance[i,i] = NA
	
	unique <- matrix(NA,nrow(comm),ncol(comm))
	
	for (r in 1:nrow(comm)){    					                  #cycle through all sites/samples
		present <- which(comm[r,]>0)                          #which species exist in this site
		for (c in present){
			unique[r,c] <- min(distance[present,c], na.rm=T)
		}
	}
	if(relative){
		unique <- unique / max(distance, na.rm = T)
	}
	return(unique)
}

#' Contribution of species or individuals to total phylogenetic/functional diversity.
#' @description Contribution of each species or individual to the total PD or FD of a number of communities.
#' @param comm A sites x species matrix, with either abundance or incidence data. If missing, the contribution of all species to the full tree is calculated.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param abund A boolean (T/F) indicating whether contribution should be weighted by abundance of each species.
#' @param relative A boolean (T/F) indicating whether contribution should be relative to total PD or FD (proportional contribution per individual or species). If FALSE, the sum of contributions for each site is equal to total PD/FD, if TRUE it is 1.
#' @details Contribution is equivalent to the evolutionary distinctiveness index (ED) of Isaac et al. (2007) if done by species and to the abundance weighted evolutionary distinctiveness (AED) of Cadotte et al. (2010) if done by individual.
#' @return A matrix of sites x species values (or values per species if no comm is given).
#' @references Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. & Baillie, J.E.M. (2007) Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS One, 2: e296.
#' @references Cadotte, M.W., Davies, T.J., Regetz, J., Kembel, S.W., Cleland, E. & Oakley, T.H. (2010) Phylogenetic diversity metrics for ecological communities: integrating species richness, abundance and evolutionary history. Ecology Letters, 13: 96-105.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,1,1,1), nrow = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' contribution(tree = tree)
#' contribution(comm, tree)
#' contribution(comm, tree, FALSE)
#' contribution(comm, tree, abund = FALSE, relative = FALSE)
#' @export
contribution <- function(comm, tree, abund = TRUE, relative = TRUE){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  if(missing(comm))
		comm = rep(1, length(tree$order))
	if(is.vector(comm))
		comm <- matrix(comm, nrow = 1)
	comm <- as.matrix(comm)
	
	if(!abund)
		comm <- ifelse(comm > 0, 1, 0)
	if (!missing(tree)){
		if (class(tree) == "phylo"){
			nEdges <- length(tree$edge.length)
			if(!is.null(tree$tip.label) && !is.null(colnames(comm))) ##if both tree and comm have species names match and reorder species (columns) in comm
				comm <- comm[,match(tree$tip.label, colnames(comm))]
		} else {
			nEdges <- length(tree$merge)
			if(!is.null(tree$labels) && !is.null(colnames(comm))) ##if both tree and comm have species names match and reorder species (columns) in comm
				comm <- comm[,match(tree$labels, colnames(comm))]
		}
	} else {
		tree = hclust(as.dist(matrix(1,ncol(comm),ncol(comm))))
	}
	
	contrib <- matrix(0,nrow(comm),ncol(comm))
	for (i in 1:nrow(comm)){											#cycle through all sites/samples
		dataSample <- prep(comm[i,], xTree(tree), TRUE)
		valueBranch <- dataSample$lenBranch / dataSample$sampleBranch
		valueBranch <- ifelse(valueBranch == Inf, 0, valueBranch)
		valueBranch <- ifelse(is.na(valueBranch), 0, valueBranch)
		for (j in 1:ncol(comm)){										#cycle through all species
			for (k in 1:nEdges){	            				#cycle through all branches
				contrib[i,j] <- contrib[i,j] + (dataSample$speciesBranch[j,k] * valueBranch[k] * comm[i,j])
			}
		}
	}
	if(relative){
	  for(r in 1:nrow(comm))
	    contrib[r,] <- contrib[r,] / c(alpha(comm[r,], tree))
	}
	if(abund){                      #contribution weighted by abundance
	  for (r in 1:nrow(comm)){											#cycle through all sites/samples
	    relAbund = comm[r,] / sum(comm[r,])
	    contrib[r,] = (contrib[r,] * relAbund) / sum(contrib[r,] * relAbund)
	  }
	}
	
	contrib[comm[] == 0] = NA
	return(contrib)
}

#' Phylogenetic/functional dispersion of species or individuals.
#' @description Average dissimilarity between any two species or individuals randomly chosen in a community.
#' @param comm A sites x species matrix, with either abundance or incidence data. If missing, the dispersion using the full tree or distance matrix is calculated.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist object representing the phylogenetic or functional distance between species.
#' @param func Calculate dispersion using originality (default), uniqueness or contribution.
#' @param abund A boolean (T/F) indicating whether dispersion should be calculated using individuals (T) or species (F).
#' @param relative A boolean (T/F) indicating whether dispersion should be relative to the maximum distance between any two species in the tree or distance matrix.
#' @details If abundance data is used and a tree is given, dispersion is the quadratic entropy of Rao (1982).
#' If abundance data is not used but a tree is given, dispersion is the phylogenetic dispersion measure of Webb et al. (2002).
#' @return A vector of values per site (or a single value if no comm is given).
#' @references Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21: 24-43.
#' @references Webb, C.O., Ackerly, D.D., McPeek, M.A. & Donoghue, M.J. (2002) Phylogenies and community ecology. Annual Review of Ecology and Systematics, 33: 475-505.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,1,1,1), nrow = 4, byrow = TRUE)
#' distance <- dist(c(1:5), method="euclidean")
#' tree <- hclust(distance, method="average")
#' dispersion(tree = tree)
#' dispersion(distance = distance)
#' dispersion(comm, tree)
#' dispersion(comm, tree, abund = FALSE)
#' dispersion(comm, tree, abund = FALSE, relative = FALSE)
#' @export
dispersion <- function(comm, tree, distance, func = "originality", abund = TRUE, relative = TRUE){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  if(missing(comm)){
		if(!missing(distance))
			comm = rep(1, attributes(distance)[1])
		if(!missing(tree))
			comm = rep(1, length(tree$order))
	}
	if(is.vector(comm))
		comm <- matrix(comm, nrow = 1)
	if(func == "originality")
		funcValue <- originality(comm, tree, distance, abund, relative)
	else if (func == "uniqueness")
		funcValue <- uniqueness(comm, tree, distance, relative)
	else if (func == "contribution")
		funcValue <- contribution(comm, tree, abund, relative)
	else
    stop(sprintf("Function %s not recognized.", func))

	disp <- rep(0,nrow(comm))
	
	for (r in 1:nrow(comm)){  						                         #cycle through all sites/samples
		present <- which(comm[r,]>0)                                 #which species exist in this site
		proportion <- comm[r,present]/sum(comm[r,present])           #proportion incidence/abundance of species in this site
		disp[r] <- sum(funcValue[r,present]*proportion)
	}
	
	disp = matrix(disp, ncol = 1)
	rownames(disp) <- rownames(comm)
	colnames(disp) <- "Dispersion"
	return(disp)
}

#' Taxonomic/phylogenetic/functional evenness of species or individuals.
#' @description Regularity of abundances and distances (if PD/FD) between species in a community.
#' @param comm A sites x species matrix, with either abundance or incidence data. If missing, the evenness using the full tree or distance matrix is calculated.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist or matrix object representing the phylogenetic or functional distance between species. If both tree and distance are missing, taxonomic evenness is calculated.
#' @param method Calculate evenness using "expected" values (default) or values based on "contribution" of species to the tree.
#' @param func Calculate evenness using "Camargo" (default) or "Bulla" index.
#' @param abund A boolean (T/F) indicating whether evenness should be calculated using abundance data.
#' @details Evenness is calculated based on the index of Camargo (1993) or Bulla (1994) using the values of both species abundances and edge lengths in the tree (if PD/FD).
#' @details If no tree or distance is provided the result is the original index.
#' @details If any site has < 2 species its value will be NA.
#' @return A vector of values per site (or a single value if no comm is given).
#' @references Bulla, L. (1994) An index of evenness and its associated diversity measure. Oikos, 70: 167-171.
#' @references Camargo, J.A. (1993) Must dominance increase with the number of subordinate species in competitive interactions? Journal of Theoretical Biology, 161: 537-542.
#' @examples comm <- matrix(c(1,2,0,0,0,1,1,0,0,0,0,2,2,0,0,1,1,1,1,100), nrow = 4, byrow = TRUE)
#' distance <- dist(c(1:5), method = "euclidean")
#' tree <- hclust(distance, method = "average")
#' evenness(comm)
#' evenness(tree = tree, func = "bulla")
#' evenness(comm, tree)
#' evenness(comm, tree, method = "contribution")
#' evenness(comm, tree, abund = FALSE)
#' @export
evenness <- function(comm, tree, distance, method = "expected", func = "camargo", abund = TRUE){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

	if(missing(comm))
		comm = rep(1, length(tree$order))
	if(is.vector(comm))
		comm <- matrix(comm, nrow = 1)
	if(!abund)
		comm <- ifelse(comm > 0, 1, 0)
	comm[is.na(comm)] = 0
	
	if(!missing(tree)){
	  comm = reorderComm(comm, tree)
	} else if (!missing(distance)){
			tree = hclust(distance, method = "average")
	} else {
		tree = hclust(as.dist(matrix(1,ncol(comm),ncol(comm))))
		tree$labels = colnames(comm)
	}
	
	evenness <- rep(0, nrow(comm))
	for (i in 1:nrow(comm)){								#cycle through all sites/samples
	  thisComm = comm[i,comm[i,] > 0]			  #redo this comm
	  if(length(thisComm) < 2){             #if comm has less than 2 species evenness is NA
	    evenness[i] = NA
	    next
	  }              
	  thisTree = as.matrix(cophenetic(tree))[comm[i,] > 0, comm[i,] > 0]    #redo this tree
	  thisTree = hclust(as.dist(thisTree))
	  if(method == "expected"){               #if expected
	    thisTree = prep(thisComm, xTree(thisTree), abund)
	    thisEdges = which(thisTree$lenBranch > 0 & thisTree$sampleBranch > 0)
		  thisObs = c()
		  for(j in thisEdges){											    			#cycle through all edges of this site/sample
		    #calculate the observed values as avg abundance per species of edge / length of edge 
		    thisObs = c(thisObs, (thisTree$sampleBranch[j] / sum(thisTree$speciesBranch[,j]) / thisTree$lenBranch[j]))
		  }
		  thisObs = thisObs / sum(thisObs)
		  if(func == "bulla"){
		    ##calculate the expected values as avg length of tree edges
		    thisExp = 1 / length(thisEdges)
		    #calculate evenness as the sum of minimum values between observed and expected with correction from Bulla, 1994
		    evenness[i] = (sum(apply(cbind(thisObs, rep(thisExp, length(thisObs))), 1, min)) - (1/length(thisEdges))) / (1-1/length(thisEdges))
		  } else if(func == "camargo"){      #if Camargo
		    nEdges = length(thisObs)
		    for(j in 1:(nEdges-1)){
		      for(k in (j+1):nEdges){
		        evenness[i] = evenness[i] + abs(thisObs[j] - thisObs[k])
		      } 
		    }
		    evenness[i] = 1 - (evenness[i] / (nEdges*(nEdges-1)/2))
		  } else {
		    stop(sprintf("Function %s not recognized.", func))
		  }
		} else if (method == "contribution") {          #if using the contribution of species
		  contrib = contribution(thisComm, thisTree, abund = abund)
		  nSp = length(thisComm)
		  if(func == "bulla"){
		    ##calculate the expected contribution as 1/nSp
		    thisExp = 1 / nSp
		    #calculate evenness as the sum of minimum values between observed and expected with correction from Bulla, 1994
		    evenness[i] = (sum(apply(cbind(contrib, rep(thisExp, nSp)), 1, min)) - (1/nSp)) / (1-1/nSp)
		  } else if(func == "camargo"){      #if Camargo
		    for(j in 1:(nSp-1)){
		      for(k in (j+1):nSp){
		        evenness[i] = evenness[i] + abs(contrib[j] - contrib[k])
		      } 
		    }
		    evenness[i] = 1 - (evenness[i] / (nSp*(nSp-1)/2))
		  } else {
		    stop(sprintf("Function %s not recognized.", func))
		  }
		} else {
		  stop(sprintf("Method %s not recognized.", method))
		}
	}
	
	evenness = matrix(evenness, ncol = 1)
	rownames(evenness) <- rownames(comm)
	colnames(evenness) <- "Evenness"
	return(evenness)
}

#' Contribution of each species or individual to the total taxonomic/phylogenetic/functional evenness.
#' @description Contribution of each observation to the regularity of abundances and distances (if PD/FD) between species in a community (or individuals in a species).
#' @param comm A sites x species matrix, with either abundance or incidence data. If missing, the evenness using the full tree or distance matrix is calculated.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist or matrix object representing the phylogenetic or functional distance between species. If both tree and distance are missing, taxonomic evenness is calculated.
#' @param method Calculate evenness using "expected" values (default) or values based on "contribution" of species to the tree.
#' @param func Calculate evenness using "Camargo" (default) or "Bulla" index.
#' @param abund A boolean (T/F) indicating whether evenness should be calculated using abundance data.
#' @details Contribution to evenness is calculated using a leave-one-out approach, whereby the contribution of a single observation is the total evenness minus the evenness calculated without that observation. Evenness is based on the index of Camargo (1993) or Bulla (1994) using the values of both species abundances and edge lengths in the tree (if PD/FD).
#' Note that the contribution of a species or individual can be negative, if the removal of an observation increases the total evenness.  
#' @details If no tree or distance is provided the result is calculated for taxonomic evenness using the original index.
#' @return A matrix of sites x species (or a vector if no comm is given).
#' @references Bulla, L. (1994) An index of evenness and its associated diversity measure. Oikos, 70: 167-171.
#' @references Camargo, J.A. (1993) Must dominance increase with the number of subordinate species in competitive interactions? Journal of Theoretical Biology, 161: 537-542.
#' @examples comm <- matrix(c(1,2,0,5,5,1,1,0,0,0,0,2,2,0,0,1,1,1,1,100), nrow = 4, byrow = TRUE)
#' distance <- dist(c(1:5), method = "euclidean")
#' tree <- hclust(distance, method = "average")
#' evenness.contribution(comm)
#' evenness.contribution(tree = tree, func = "bulla")
#' evenness.contribution(comm, tree)
#' evenness.contribution(comm, tree, method = "contribution")
#' evenness.contribution(comm, tree, abund = FALSE)
#' @export
evenness.contribution <- function(comm, tree, distance, method = "expected", func = "camargo", abund = TRUE){
  
  #check if right data is provided
  if(missing(comm))
    comm <- rep(1, length(tree$order))
  if(is.vector(comm))
    comm <- matrix(comm, nrow = 1)
  if(!abund)
    comm <- ifelse(comm > 0, 1, 0)
  comm[is.na(comm)] = 0
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  if(!missing(tree)){
    comm = reorderComm(comm, tree)
  } else if (!missing(distance)){
    tree = hclust(distance, method = "average")
  } else {
    tree = hclust(as.dist(matrix(1,ncol(comm),ncol(comm))))
    tree$labels = colnames(comm)
  }
  
  #extract total evenness
  tot_evenness <- evenness(comm = comm, tree = tree, distance = distance, method = method, func = func, abund = abund)

  #leave-one-out
  evenness.contrib <- comm
  evenness.contrib[] <- NA
  
  for(i in 1:nrow(evenness.contrib))  {
    if(length(comm[i, comm[i,]>0]) < 3) {
      warning(paste("Community ", as.character(i), " contains less than 3 species. Cannot evaluate contribution to evenness." ) )
    } else  {
      for(k in 1:length(evenness.contrib[i,])) {
        comm2 <- comm[i,]
        comm2[k] <- 0 #for each iteration, assign one species to 0 
        evenness.contrib[i,k] <- (tot_evenness[i] - evenness(comm = comm2, tree = tree, distance = distance, method = method, func = func, abund = abund))
      }
    }
  }
  
  evenness.contrib[comm[] == 0] <- NA
  
  rownames(evenness.contrib) <- rownames(comm)
  colnames(evenness.contrib) <- colnames(comm)
  return(evenness.contrib)
  
}

#' Alpha diversity using convex hull hypervolumes.
#' @description Estimation of functional richness of one or multiple sites, based on convex hull hypervolumes.
#' @param comm A 'convhulln' object or list, preferably built with function hull.build.
#' @details Estimates the functional richness (alpha FD) of one or more communities using convex hull hypervolumes.
#' Functional richness is expressed as the total volume of the convex hull.
#' @return One value or a vector of alpha diversity values for each site.
#' @examples comm = rbind(c(1,3,0,5,3), c(3,2,5,0,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = hull.build(comm[1,], trait)
#' hull.alpha(hv)
#' hvlist = hull.build(comm, trait, axes = 2)
#' hull.alpha(hvlist)
#' @export
hull.alpha <- function(comm){

  #check if right data is provided
  if (!(class(comm) %in% c("list", "convhulln")))
    stop("A convhulln or list is needed as input data.")
  
  #if single comm
  if (class(comm) == "convhulln")
    return(comm$vol)

  #if multiple comm
  alphaValues <- c()
  for (i in 1:length(comm))
    alphaValues <- c(alphaValues, comm[[i]]$vol)
  
  #finalize
  names(alphaValues) <- names(comm)
  return(alphaValues)
}

#' Beta diversity partitioning using convex hull hypervolumes.
#' @description Pairwise beta diversity partitioning into replacement and net difference in amplitude components of convex hulls.
#' @param comm A list of 'convhulln' objects, preferably built with function hull.build.
#' @param func Partial match indicating whether the Jaccard (default) or Soerensen family of beta diversity measures should be used.
#' @param comp Boolean indicating whether beta diversity components (shared and unique fractions) should be returned.
#' @details Computes a pairwise decomposition of the overall differentiation among kernel hypervolumes into two components: the replacement (shifts) of space between hypervolumes and net differences between the amount of space enclosed by each hypervolume.
#' The beta diversity measures used here follow the FD partitioning framework where Btotal = Breplacement + Brichness. Beta diversity ranges from 0 (when hypervolumes are identical) to 1 (when hypervolumes are fully dissimilar).
#' See Carvalho & Cardoso (2020) and Mammola & Cardoso (2020) for the full formulas of beta diversity used here.
#' @return Three pairwise distance matrices, one per each of the three beta diversity components.
#' @references Carvalho, J.C. & Cardoso, P. (2020) Decomposing the causes for niche differentiation between species using hypervolumes. Frontiers in Ecology and Evolution. https://doi.org/10.3389/fevo.2020.00243
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13424
#' @examples comm <- rbind(c(1,1,1,1,1), c(1,1,1,1,1), c(0,0,1,1,1),c(0,0,1,1,1))
#' colnames(comm) = c("SpA","SpB","SpC","SpD", "SpE")
#' rownames(comm) = c("Site 1","Site 2","Site 3","Site 4")
#'
#' trait <- cbind(c(2.2,4.4,6.1,8.3,3),c(0.5,1,0.5,0.4,4))
#' colnames(trait) = c("Trait 1","Trait 2")
#' rownames(trait) = colnames(comm)
#'
#' hvlist = hull.build(comm, trait)
#' hull.beta(hvlist)
#' hvlist = hull.build(comm, trait, axes = 2)
#' hull.beta(hvlist, comp = TRUE)
#' @export
hull.beta <- function(comm, func = "jaccard", comp = FALSE){
  
  #check if right data are provided
  if (!(class(comm) %in% c("list")))
    stop("A list of convhulln is needed as input data.")
  
  #create matrices to store results
  nComm <- length(comm)
  Btotal <- Brepl <- Brich <- compA <- compB <- compC <- matrix(NA, nrow = nComm, ncol = nComm)
  
  #calculate beta values
  for (i in 1:(nComm-1)){
    for(j in (i+1):nComm){
      intersection <- geometry::intersectn(comm[[i]]$p, comm[[j]]$p, options = "FA")$ch$vol
      unique1 <- comm[[i]]$vol - intersection
      unique2 <- comm[[j]]$vol - intersection
      union <- unique1 + unique2 + intersection
      if(comp){
        compA[j,i] = union - unique1 - unique2
        compB[j,i] = unique1
        compC[j,i] = unique2
      }
      if(tolower(substr(func, 1, 1)) == "s")
        union <- 2 * union - unique1 - unique2
      Btotal[j,i] <- (unique1 + unique2) / union
      Brepl[j,i]  <- 2 * min(unique1, unique2) / union 
      Brich[j,i]  <- abs(unique1 - unique2) / union 
    }
  }
  
  #finalize
  rownames(Btotal) <- colnames(Btotal) <- rownames(Brepl) <- colnames(Brepl) <- rownames(Brich) <- colnames(Brich) <- names(comm)
  betaValues <- list(Btotal = round(as.dist(Btotal),3), Brepl = round(as.dist(Brepl),3), Brich = round(as.dist(Brich),3))
  if (comp){
    rownames(compA) <- colnames(compA) <- rownames(compB) <- colnames(compB) <- rownames(compC) <- colnames(compC) <- names(comm)
    betaValues$Shared = round(as.dist(compA),3)
    betaValues$Unique1 = round(as.dist(compB),3)
    betaValues$Unique2 = round(as.dist(compC),3)
  }
  return(betaValues)
}

#' Contribution of each observation to a convex hull hypervolume.
#' @description Contribution of each species or individual to the total volume of one or more convex hulls.
#' @param comm A 'convhulln' object or list, preferably built with function hull.build.
#' @param relative A boolean (T/F) indicating whether contribution should be relative to total PD or FD (proportional contribution per individual or species). If FALSE, the sum of contributions for each site is equal to total PD/FD, if TRUE it is 1.
#' @details The contribution of each observation (species or individual) to the total volume of a convex hull, calculated as the difference in volume between the total convex hull and a second hypervolume lacking this specific observation (i.e., leave-one-out approach; Mammola & Cardoso, 2020). 
#' @return A vector or matrix with the contribution values of each species or individual for each site.
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.13424
#' @examples comm = rbind(c(1,3,0,5,3), c(3,2,5,1,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = hull.build(comm[1,], trait)
#' hull.contribution(hv)
#' hvlist = hull.build(comm, trait, axes = 2)
#' hull.contribution(hvlist, relative = TRUE)
#' @export
hull.contribution = function(comm, relative = FALSE){
  
  #check if right data is provided
  if (!(class(comm) %in% c("list", "convhulln")))
    stop("A convhulln or list is needed as input data.")
  
  #if a convex hull is provided go for it.
  if (class(comm) == "convhulln"){
    contrib <- c()
    for (i in 1:nrow(comm$p))
      contrib <- c(contrib, (comm$vol - geometry::convhulln(comm$p[-i,], options = "FA")$vol))
    names(contrib) = rownames(comm$p)
    if(relative)
      contrib = contrib/sum(contrib)
    
  #if a comm matrix is provided just call this same function using hypervolumes.
  } else {
    
    #get species names from hypervolumes and order them alphabetically
    spp = c()
    for(i in 1:length(comm)){
      spp = c(spp, rownames(comm[[i]]$p))
    }
    spp = unique(spp)
    spp = spp[order(spp)]
    
    #get values calling this same function one hv at a time
    contrib <- matrix(NA, nrow = length(comm), ncol = length(spp))
    rownames(contrib) = names(comm)
    colnames(contrib) = spp
    for (i in 1:length(comm)){
      contr = hull.contribution(comm[[i]], relative)
      contrib[i, which(spp %in% names(contr))] = contr
    }
  }
  
  return(contrib)
}

#' Alpha diversity using kernel density hypervolumes.
#' @description Estimation of functional richness of one or multiple sites, based on n-dimensional hypervolumes.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @details Estimates the functional richness (alpha FD) of one or more communities using kernel density hypervolumes, as implemented in Blonder et al. (2014, 2018).
#' Functional richness is expressed as the total volume of the n-dimensional hypervolume (Mammola & Cardoso, 2020). Note that the hypervolume is dimensionless, and that only hypervolumes with the same number of dimensions can be compared in terms of functional richness.
#' Given that the density and positions of stochastic points in the hypervolume are probabilistic, the functional richness of the trait space will intimately depend on the quality of input hypervolumes (details in Mammola & Cardoso, 2020).
#' @return A value or vector of alpha diversity values for each site.
#' @references Blonder, B., Lamanna, C., Violle, C. & Enquist, B.J. (2014) The n-dimensional hypervolume. Global Ecology and Biogeography, 23: 595-609.
#' @references Blonder, B., Morrow, C.B., Maitner, B., Harris, D.J., Lamanna, C., Violle, C., ... & Kerkhoff, A.J. (2018) New approaches for delineating n-dimensional hypervolumes. Methods in Ecology and Evolution, 9: 305-319.
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @examples \dontrun{
#' comm = rbind(c(1,3,0,5,3), c(3,2,5,0,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' kernel.alpha(hv)
#' hvlist = kernel.build(comm, trait, axes = 0.8)
#' kernel.alpha(hvlist)
#' }
#' @export
kernel.alpha <- function(comm){

  #check if right data is provided
  if (!(class(comm) %in% c("HypervolumeList", "Hypervolume")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")

  #convert data if needed
  if (class(comm) == "Hypervolume")
    return(get_volume(comm))

  #calculate alpha values and give them a name
  alphaValues <- c()
  for (i in 1:length(comm@HVList)){
    alphaValues <- c(alphaValues, get_volume(comm@HVList[[i]]))
    names(alphaValues[i]) <- comm@HVList[[i]]@Name
  }
  
  #return alpha values
  return(alphaValues)
}

#' Beta diversity partitioning using kernel density hypervolumes.
#' @description Pairwise beta diversity partitioning into replacement and net difference in amplitude components of n-dimensional hypervolumes.
#' @param comm A 'HypervolumesList' object, preferably built using function kernel.build.
#' @param func Partial match indicating whether the Jaccard or Soerensen family of beta diversity measures should be used.  If not specified, default is Jaccard.
#' @param comp Boolean indicating whether beta diversity components (shared and unique fractions) should be returned
#' @details Computes a pairwise decomposition of the overall differentiation among kernel hypervolumes into two components: the replacement (shifts) of space between hypervolumes and net differences between the amount of space enclosed by each hypervolume.
#' The beta diversity measures used here follow the FD partitioning framework developed by Carvalho & Cardoso (2020), where Btotal = Breplacement + Brichness. Beta diversity ranges from 0 (when hypervolumes are identical) to 1 (when hypervolumes are fully dissimilar).
#' See Carvalho & Cardoso (2020) and Mammola & Cardoso (2020) for the full formulas of beta diversity used here.
#' @return Three pairwise distance matrices, one per each of the three beta diversity components. If comp = TRUE also three distance matrices with beta diversity components.
#' @references Carvalho, J.C. & Cardoso, P. (2020) Decomposing the causes for niche differentiation between species using hypervolumes. Frontiers in Ecology and Evolution, 8: 243.
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @examples \dontrun{
#' comm <- rbind(c(1,1,1,1,1), c(1,1,1,1,1), c(0,0,1,1,1),c(0,0,1,1,1))
#' colnames(comm) = c("SpA","SpB","SpC","SpD", "SpE")
#' rownames(comm) = c("Site 1","Site 2","Site 3","Site 4")
#'
#' trait <- cbind(c(2.2,4.4,6.1,8.3,3),c(0.5,1,0.5,0.4,4),c(0.7,1.2,0.5,0.4,5),c(0.7,2.2,0.5,0.3,6))
#' colnames(trait) = c("Trait 1","Trait 2","Trait 3","Trait 4")
#' rownames(trait) = colnames(comm)
#'
#' hvlist = kernel.build(comm, trait)
#' kernel.beta(hvlist)
#' hvlist = kernel.build(comm, trait, axes = 0.9)
#' kernel.beta(hvlist, comp = TRUE)
#' }
#' @export
kernel.beta = function(comm, func = "jaccard", comp = FALSE){
  
  #check if right data is provided
  if (!(class(comm) %in% c("HypervolumeList")))
    stop("A HypervolumeList is needed as input data.")
  
  #create matrices to store results
  nComm <- length(comm@HVList)
  Btotal <- Brepl <- Brich <- compA <- compB <- compC <- matrix(NA, nrow = nComm, ncol = nComm)

  #calculate beta values and give them a name
  commNames <- c()
  commNames[nComm] <- comm@HVList[[nComm]]@Name
  for (i in 1:(nComm-1)){
    for(j in (i+1):(nComm)){
      hyperSet <- hypervolume_set(comm@HVList[[i]], comm@HVList[[j]], check.memory = FALSE, verbose = FALSE, num.points.max = 10000)
      union   <- hyperSet[[4]]@Volume
      unique1 <- hyperSet[[5]]@Volume
      unique2 <- hyperSet[[6]]@Volume
      if(comp){
        compA[j,i] = union - unique1 - unique2
        compB[j,i] = unique1
        compC[j,i] = unique2
      }
      if(tolower(substr(func, 1, 1)) == "s")
        union <- 2 * union - unique1 - unique2
      Btotal[j,i] <- (unique1 + unique2) / union
      Brepl[j,i]  <- 2 * min(unique1, unique2) / union 
      Brich[j,i]  <- abs(unique1 - unique2) / union
    }
    commNames[i] <- comm@HVList[[i]]@Name
  }
  
  #finalize
  rownames(Btotal) <- colnames(Btotal) <- rownames(Brepl) <- colnames(Brepl) <- rownames(Brich) <- colnames(Brich) <- commNames
  betaValues <- list(Btotal = round(as.dist(Btotal),3), Brepl = round(as.dist(Brepl),3), Brich = round(as.dist(Brich),3))
  if (comp){
    rownames(compA) <- colnames(compA) <- rownames(compB) <- colnames(compB) <- rownames(compC) <- colnames(compC) <- commNames
    betaValues$Shared = round(as.dist(compA),3)
    betaValues$Unique1 = round(as.dist(compB),3)
    betaValues$Unique2 = round(as.dist(compC),3)
  }
  return(betaValues)
}

#' Functional beta diversity evenness using kernel density hypervolumes.
#' @description Difference of evenness between pairs of sites, measuring the regularity of stochastic points distribution within the total functional space.
#' @param comm A 'HypervolumesList' object, preferably built using function kernel.build.
#' @details This measure is simply the pairwise difference of evenness calculated based on the functional evenness (Mason et al., 2005) of a n-dimensional hypervolume, namely the regularity of stochastic points distribution within the total trait space (Mammola & Cardoso, 2020).
#' Evenness is calculated as the overlap between the observed hypervolume and a theoretical hypervolume where traits and abundances are evenly distributed within the range of their values (Carmona et al., 2016, 2019).
#' @return Distance matrix between sites.
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2016) Traits without borders: integrating functional diversity across scales. Trends in Ecology and Evolution, 31: 382-394.
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2019) Trait probability density (TPD): measuring functional diversity across scales based on TPD with R. Ecology, 100: e02876.
#' @references Mason, N.W.H., Mouillot, D., Lee, W.G. & Wilson, J.B. (2005) Functional richness, functional evenness and functional divergence: the primary components of functional diversity. Oikos, 111: 112-118.
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @examples \dontrun{
#' comm <- rbind(c(1,1,1,1,1), c(1,1,1,1,1), c(0,0,1,1,1),c(0,0,1,1,1))
#' colnames(comm) = c("SpA","SpB","SpC","SpD", "SpE")
#' rownames(comm) = c("Site 1","Site 2","Site 3","Site 4")
#'
#' trait <- cbind(c(2.2,4.4,6.1,8.3,3),c(0.5,1,0.5,0.4,4),c(0.7,1.2,0.5,0.4,5),c(0.7,2.2,0.5,0.3,6))
#' colnames(trait) = c("Trait 1","Trait 2","Trait 3","Trait 4")
#' rownames(trait) = colnames(comm)
#' 
#' hvlist = kernel.build(comm, trait)
#' kernel.beta.evenness(hvlist)
#' hvlist = kernel.build(comm, trait, axes = 0.9)
#' kernel.beta.evenness(hvlist)
#' }
#' @export
kernel.beta.evenness <- function(comm){
    if (!(class(comm) %in% c("HypervolumeList")))
      stop("A HypervolumeList is needed as input data.")
    return(dist(kernel.evenness(comm)))
}

#' Functional originality of observations in kernel density hypervolumes.
#' @description Average dissimilarity between a species or individual and a sample of random points within the boundaries of the n-dimensional hypervolume.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @param frac A value between 0.01 and 1, indicating the fraction of random points to be used in the estimation of originality. Default is 0.1.
#' @param relative A boolean (T/F) indicating whether originality should be relative to the most original species in the community.
#' @details A measure of the originality (sensu Pavoine et al., 2005) of each observation (species or individuals) used to construct the n-dimensional hypervolume. In a probabilistic hypervolume, originality is calculated as the average distance between each observation to a sample of stochastic points within the boundaries of the n-dimensional hypervolume (Mammola & Cardoso, 2020).
#' Originality is a measure of functional rarity (sensu Violle et al., 2017; Carmona et al., 2017) that allows to map the contribution of each observation to the divergence components of FD (Mammola & Cardoso, 2020).
#' The number of sample points to be used in the estimation of the originality is controlled by the frac parameter. Increase frac for less deviation in the estimation, but mind that computation time also increases.
#' @return A vector or matrix with the originality values of each species or individual in each site.
#' @references Carmona, C.P., de Bello, F., Sasaki, T., Uchida, K. & Partel, M. (2017) Towards a common toolbox for rarity: A response to Violle et al. Trends in Ecology and Evolution, 32: 889-891. 
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @references Pavoine, S., Ollier, S. & Dufour, A.-B. (2005) Is the originality of a species measurable? Ecology Letters, 8: 579-586.
#' @references Violle, C., Thuiller, W., Mouquet, N., Munoz, F., Kraft, N.J.B., Cadotte, M.W., ... & Mouillot, D. (2017) Functional rarity: the ecology of outliers. Trends in Ecology and Evolution, 32: 356-367.
#' @examples \dontrun{
#' comm = rbind(c(1,3,0,5,3), c(3,2,5,1,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' kernel.originality(hv)
#' hvlist = kernel.build(comm, trait)
#' kernel.originality(hvlist)
#' kernel.originality(hvlist, relative = TRUE)
#' }
#' @export
kernel.originality = function(comm, frac = 0.1, relative = FALSE){

  #check if right data is provided
  if (!(class(comm) %in% c("Hypervolume", "HypervolumeList")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")
  
  #check if right frac parameter is provided
  if (frac < 0.01 | frac > 1)
    stop("The frac parameter should be a number between 0.01 and 1.")

  if (class(comm) == "Hypervolume") {    
    
    sample.points <- comm@RandomPoints[sample(1:nrow(comm@RandomPoints), nrow(comm@RandomPoints)*frac), ]
    
    originality <- c()
    for (i in 1:length(unique(rownames(comm@Data)))){
      originality_run <- c()
      subHvData <- comm@Data[rownames(comm@Data)[i], ]
      for (r in 1:nrow(sample.points))
        originality_run <- c(originality_run, dist(c(subHvData, sample.points[r,1:ncol(sample.points)])))
      originality <- c(originality, mean(originality_run))
    }
    names(originality) = rownames(comm@Data)
    
  } else if (class(comm) == "HypervolumeList") {
    if (all(rownames(comm@HVList[[1]]@Data) == 1:nrow(comm@HVList[[1]]@Data)))
      warning("Species names are probably missing from hypervolumes, consider renaming them.")
    spp = c()
    commNames = c()
    for(i in 1:length(comm@HVList)){
      spp = c(spp, rownames(comm@HVList[[i]]@Data))
      commNames = c(commNames, comm@HVList[[i]]@Name)
    }
    spp = unique(spp)
    spp = spp[order(spp)]
    
    originality = matrix(NA, nrow = length(comm@HVList), ncol = length(spp))
    rownames(originality) = commNames
    colnames(originality) = spp
    for(i in 1:length(comm@HVList)){
      originality[i, which(spp %in% rownames(comm@HVList[[i]]@Data))] <- kernel.originality(comm = comm@HVList[[i]], frac = frac)
    }
  }
  
  #finalize
  if(relative){
    if(is.matrix(originality))
      originality = t(apply(originality, 1, function(x) x / max(x, na.rm = TRUE)))
    else
      originality <- originality/max(originality, na.rm = TRUE)
  }
  return(originality)
}

#' Contribution of each observation to the kernel density hypervolume.
#' @description Contribution of each species or individual to the total volume of one or more kernel hypervolumes.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @param func Calculate contribution using either closest "neighbor" or leave "one out" approach. 
#' @param relative A boolean (T/F) indicating whether contribution should be relative to total FD (proportional contribution per individual or species). If FALSE, the sum of contributions for each site is equal to total FD, if TRUE it is 1.
#' @details Contribution is a measure of functional rarity (sensu Violle et al., 2017; Carmona et al., 2017) that allows to map the contribution of each observation to the richness components of FD (Mammola & Cardoso, 2020).
#' If using func = "neighbor", each random point will be attributed to the closest species. The contribution of each species will be proportional to the number of its points. The sum of contributions of all species is equal to total richness.
#' Note that the contribution of a species or individual can be negative if leave-one-out approach is taken, if the removal of an observation increases the total volume (see Figure 2d in Mammola & Cardoso 2020).  
#' This might happen, although not always, in cases when the presence of a given species decreases the average distance between all the species in the community, i.e., when a given species is close to the "average" species of that community, making that community less diverse in some sense (Mammola & Cardoso, 2020).
#' @return A matrix with the contribution values of each species or individual for each site.
#' @references Carmona, C.P., de Bello, F., Sasaki, T., Uchida, K. & Partel, M. (2017) Towards a common toolbox for rarity: A response to Violle et al. Trends in Ecology and Evolution, 32(12): 889-891. 
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @references Violle, C., Thuiller, W., Mouquet, N., Munoz, F., Kraft, N.J.B., Cadotte, M.W., ... & Mouillot, D. (2017) Functional rarity: The ecology of outliers. Trends in Ecology and Evolution, 32: 356-367.
#' @examples \dontrun{
#' comm = rbind(c(1,3,0,5,3), c(3,2,5,1,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' kernel.contribution(hv)
#' hvlist = kernel.build(comm, trait, axes = 2)
#' kernel.contribution(hvlist)
#' kernel.contribution(hvlist, relative = TRUE)
#' }
#' @export
kernel.contribution = function(comm, func = "neighbor", relative = FALSE){
  
  #check if right data is provided
  if (!(class(comm) %in% c("Hypervolume", "HypervolumeList")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")
  
  #if hypervolume is provided go for it.
  if (class(comm) == "Hypervolume"){
    contrib  <- c()
    if(substring(func, 1, 3) == "nei"){
      
      #check if any data points are similar and create vector with similar
      simPoints = c()
      for(i in 1:nrow(comm@Data)){
        simPoints[i] = which(paste(comm@Data[i, ], collapse = "") == apply(comm@Data, 1, paste, collapse = ""))[1]
      }
      
      #calculate distance between each random point and each data point
      distPoints = matrix(NA, ncol = nrow(comm@Data), nrow = nrow(comm@RandomPoints))
      for(d in 1:nrow(comm@Data)){
        for(r in 1:nrow(comm@RandomPoints)){
          distPoints[r, d] = dist(rbind(comm@RandomPoints[r,], comm@Data[d,]))
        }
      }
      
      #attribute each random point to each data point, repeated points will have 0
      distPoints = apply(distPoints, 1, which.min)
      
      ##calculate contribution of each data point
      for(d in 1:nrow(comm@Data)){
        contrib[d] = sum(distPoints == d)
      }
      contrib = contrib / length(distPoints) * kernel.alpha(comm)
      
      ##divide contribution of each data point by similar species
      temp = c()
      for(i in 1:nrow(comm@Data)){
        dups = which(simPoints == simPoints[i])
        temp[i] = contrib[dups[1]] / length(dups)
      }
      contrib = temp
      
    } else if (substring(func, 1, 3) == "one"){
      for (i in 1:nrow(comm@Data)){
        if (comm@Method == "Box kernel density estimate")
          contrib <- c(contrib, comm@Volume - hypervolume_box(comm@Data[-i,], verbose = FALSE)@Volume)
        else if (comm@Method == "Gaussian kernel density estimate")
          contrib <- c(contrib, comm@Volume - hypervolume_gaussian(comm@Data[-i,], verbose = FALSE)@Volume)
        else if (comm@Method == "One-class support vector machine")
          contrib <- c(contrib, comm@Volume - hypervolume_svm(comm@Data[-i,], verbose = FALSE)@Volume)
      }
    } else {
      stop("func not recognized, must be 'neighbor' or 'one out'")
    }
    
    names(contrib) = rownames(comm@Data)
    
  } else if (class(comm) == "HypervolumeList") {
    if (all(rownames(comm@HVList[[1]]@Data) == 1:nrow(comm@HVList[[1]]@Data)))
      warning("Species names are probably missing from hypervolumes, consider renaming them.")
    spp = c()
    commNames = c()
    for(i in 1:length(comm@HVList)){
      spp = c(spp, rownames(comm@HVList[[i]]@Data))
      commNames = c(commNames, comm@HVList[[i]]@Name)
    }
    spp = unique(spp)
    spp = spp[order(spp)]
    
    contrib = matrix(NA, nrow = length(comm@HVList), ncol = length(spp))
    rownames(contrib) = commNames
    colnames(contrib) = spp
    for(i in 1:length(comm@HVList)){
      contrib[i, which(spp %in% rownames(comm@HVList[[i]]@Data))] <- kernel.contribution(comm = comm@HVList[[i]], func = func)
    }
  }

  #finalize
  if(relative){
    if(is.matrix(contrib))
      contrib = t(apply(contrib, 1, function(x) x / sum(x, na.rm = TRUE)))
    else
      contrib <- contrib/sum(contrib, na.rm = TRUE)
  }
  return(contrib)
}

#' Functional dispersion of kernel density hypervolumes.
#' @description Average distance to centroid or dissimilarity between random points within the boundaries of the kernel density hypervolume.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @param func Function for calculating dispersion. One of 'divergence', 'dissimilarity' or 'regression'.
#' @param frac A value between 0.01 and 1, indicating the fraction of random points to be used. Default is 0.1.
#' @details This function calculates dispersion either: i) as the average distance between stochastic points within the kernel density hypervolume and the centroid of these points (divergence; Laliberte & Legendre, 2010; see also Carmona et al., 2019); ii) as the average distance between all points (dissimilarity, see also function BAT::dispersion); or iii) as the average distance between stochastic points within the kernel density hypervolume and a regression line fitted through the points.
#' The number of stochastic points is controlled by the 'frac' parameter (increase this number for less deviation in the estimation).
#' @return A value or vector of dispersion values for each site.
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2019) Trait probability density (TPD): measuring functional diversity across scales based on TPD with R. Ecology, 100: e02876.
#' @references Laliberte, E. & Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. Ecology 91: 299-305.
#' @examples \dontrun{
#' comm = rbind(c(1,3,0,5,3), c(3,2,5,1,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' kernel.dispersion(hv)
#' hvlist = kernel.build(comm, trait, axes = 2)
#' kernel.dispersion(hvlist)
#' kernel.dispersion(hvlist, func = "divergence")
#' }
#' @export
kernel.dispersion = function(comm, func = 'dissimilarity', frac = 0.1) {
  
  #check if right data is provided
  if (!(class(comm) %in% c("HypervolumeList", "Hypervolume")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")
  
  #check if right frac parameter is provided
  if (frac < 0.01 | frac > 1)
    stop("Frac parameter should be a number between 0.01 and 1.")

  if (class(comm) == "Hypervolume"){
    
    random_points = comm@RandomPoints[sample(1:nrow(comm@RandomPoints), nrow(comm@RandomPoints)*frac), ]
    
    if (func == "dissimilarity"){
      disp <- mean(dist(random_points))
    } else if (func == "divergence"){
      cent <- get_centroid(comm)
      disp <- c()
      for (k in 1:comm@Dimensionality){
        disp <- cbind(disp, (cent[k] - random_points[, k])^2)
      }
      disp <- mean(rowSums(disp)^0.5)
    } else if (func == "regression"){
      disp <- c()
      for(m in 1:(ncol(random_points)-1)){           # build all bivariate predictor combinations
        for(n in (m+1):ncol(random_points)){
          disp <- append(disp,summary(lm(random_points[,m]~random_points[,n]))$r.squared) #get the R^2
        }
      }
      disp <- mean(disp)
    } else {
      stop(sprintf("Function %s not recognized.", func))
    }
    names(disp) <- comm@Name
    
  } else if (class(comm) == "HypervolumeList"){
    disp <- c()
    for (i in 1:length(comm@HVList))
      disp <- c(disp, kernel.dispersion(comm@HVList[[i]], func = func, frac = frac))
  }
  return(disp)
}

#' Functional evenness of kernel density hypervolumes.
#' @description Functional evenness of a community, measuring the regularity of stochastic points distribution within the total functional space.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @details This function measures the functional evenness (Mason et al., 2005) of a n-dimensional hypervolume, namely the regularity of stochastic points distribution within the total trait space (Mammola & Cardoso, 2020).
#' Evenness is calculated as the overlap between the observed hypervolume and a theoretical hypervolume where traits and abundances are evenly distributed within the range of their values (Carmona et al., 2016, 2019).
#' @return A value or vector of evenness values for each site.
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2016) Traits without borders: integrating functional diversity across scales. Trends in Ecology and Evolution, 31: 382-394.
#' @references Carmona, C.P., de Bello, F., Mason, N.W.H. & Leps, J. (2019) Trait probability density (TPD): measuring functional diversity across scales based on TPD with R. Ecology, 100: e02876.
#' @references Mason, N.W.H., Mouillot, D., Lee, W.G. & Wilson, J.B. (2005) Functional richness, functional evenness and functional divergence: the primary components of functional diversity. Oikos, 111: 112-118.
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @examples \dontrun{
#' comm = rbind(c(100,3,0,5,3), c(3,2,5,1,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' kernel.evenness(hv)
#' hv = kernel.build(comm[1,], trait, abund = FALSE)
#' kernel.evenness(hv)
#' hvlist = kernel.build(comm, trait, axes = 2)
#' kernel.evenness(hvlist)
#' }
#' @export
kernel.evenness = function(comm) {
  
  #check if right data is provided
  if (!(class(comm) %in% c("HypervolumeList", "Hypervolume")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")

  #if hypervolume is provided go for it
  if (class(comm) == "Hypervolume") {
    
    #creating a perfectly even hypervolume within the distribution of traits
    ref_even <- comm@Data
    ref_even[] <- NA
    for(j in 1:ncol(comm@Data)){
      space <- (range(comm@Data[, j])[2] - range(comm@Data[, j])[1]) / (length(comm@Data[, j]) - 1)
      ref_even[, j] <- seq(from = range(comm@Data[,j])[1], to = range(comm@Data[,j])[2], by = space)
    }
    
    if (comm@Method == "Box kernel density estimate")
      hv_ref <- hypervolume_box(ref_even, verbose = FALSE)
    else if (comm@Method == "Gaussian kernel density estimate")
      hv_ref <- hypervolume_gaussian(ref_even, verbose = FALSE)
    else if (comm@Method == "One-class support vector machine")
      hv_ref <- hypervolume_svm(ref_even, verbose = FALSE)
    
    #Checking the overlap between the hypervolume and the even hypervolume
    set <- hypervolume_set(hv_ref, comm, check.memory = FALSE, distance.factor = 1, verbose = FALSE)
    even <- hypervolume_overlap_statistics(set)[1]
    names(even) <- comm@Name
    
  } else if (class(comm) == "HypervolumeList"){
    even <- c()
    for(i in 1:length(comm@HVList)){
      even <- c(even, kernel.evenness(comm@HVList[[i]]))
      message(paste("Evenness of hypervolume", as.character(i), "out of", as.character(length(comm@HVList)), "has been estimated.\n"))
    }
  }
  return(even)
}

#' Contribution of each observation to the evenness of a kernel density hypervolume.
#' @description Contribution of each species or individual to the evenness of one or more kernel hypervolumes.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @details The contribution of each observation (species or individual) to the total evenness of a kernel hypervolume. Contribution to evenness is calculated as the difference in evenness between the total hypervolume and a second hypervolume lacking this specific observation (i.e., leave-one-out approach; Mammola & Cardoso, 2020). 
#' Note that the contribution of a species or individual can be negative, if the removal of an observation increases the total evenness.  
#' @return A vector or matrix with the contribution values of each species or individual for each community or species respectively.
#' @references Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.
#' @examples \dontrun{
#' comm = rbind(c(100,3,0,5,3), c(3,2,5,1,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' kernel.evenness.contribution(hv)
#' hvlist = kernel.build(comm, trait)
#' kernel.evenness.contribution(hvlist)
#' hvlist = kernel.build(comm, trait, axes = 0.8)
#' kernel.evenness.contribution(hvlist)
#' }
#' @export
kernel.evenness.contribution = function(comm){
  
  #check if right data is provided
  if (!(class(comm) %in% c("Hypervolume", "HypervolumeList")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")

  #if hypervolume is provided go for it.
  if (class(comm) == "Hypervolume"){

    #extract total evenness:
    comm.evenness <- kernel.evenness(comm)
    
    #leave-one-out:
    contrib <- c()
    for (i in 1:nrow(comm@Data)){
      if (comm@Method == "Box kernel density estimate")
        contrib <- c(contrib, comm.evenness - kernel.evenness(hypervolume_box(comm@Data[-i, ], verbose = FALSE)))
      else if (comm@Method == "Gaussian kernel density estimate")
        contrib <- c(contrib, comm.evenness - kernel.evenness(hypervolume_gaussian(comm@Data[-i, ], verbose = FALSE)))
      else if (comm@Method == "One-class support vector machine")
        contrib <- c(contrib, comm.evenness - kernel.evenness(hypervolume_svm(comm@Data[-i, ], verbose = FALSE)))
    }
    
    names(contrib) = rownames(comm@Data)

  }	else if (class(comm) == "HypervolumeList"){
    if (all(rownames(comm@HVList[[1]]@Data) == 1:nrow(comm@HVList[[1]]@Data)))
      warning("Species names are probably missing from hypervolumes, consider renaming them.")
    spp = c()
    commNames = c()
    for(i in 1:length(comm@HVList)){
      spp = c(spp, rownames(comm@HVList[[i]]@Data))
      commNames = c(commNames, comm@HVList[[i]]@Name)
    }
    spp = unique(spp)
    spp = spp[order(spp)]
    
    contrib = matrix(NA, nrow = length(comm@HVList), ncol = length(spp))
    rownames(contrib) = commNames
    colnames(contrib) = spp
    for(i in 1:length(comm@HVList)){
      contrib[i, which(spp %in% rownames(comm@HVList[[i]]@Data))] <- kernel.evenness.contribution(comm = comm@HVList[[i]])
    }
  }

  return(contrib)
}

#' Pairwise similarity among kernel density hypervolumes.
#' @description Calculate pairwise distance metrics (centroid and minimum distance) and similarity indices (Intersection, Jaccard, Soerensen-Dice) among n-dimensional hypervolumes.
#' @param comm A 'HypervolumesList' object, preferably built using function kernel.build.
#' @details Computes a pairwise comparison between kernel density hypervolumes of multiple species or communities, based on the distance and similarity metrics implemented in hypervolume R package (Blonder et al., 2014, 2018).
#' See Mammola (2019) for a description of the different indices, and a comparison between their performance. Note that computation time largely depends on the number of 'Hypervolume' objects in the list, and scales almost exponentially with the number of hypervolume axes.
#' @return Five pairwise distance matrices, one per each of the distance and similarity indices (in order: distance between centroids, minimum distance, Jaccard overlap, Soerensen-Dice overlap, and Intersection among hypervolumes).
#' @references Blonder, B., Lamanna, C., Violle, C. & Enquist, B.J. (2014) The n-dimensional hypervolume. Global Ecology and Biogeography, 23: 595-609.
#' @references Blonder, B., Morrow, C.B., Maitner, B., Harris, D.J., Lamanna, C., Violle, C., ... & Kerkhoff, A.J. (2018) New approaches for delineating n-dimensional hypervolumes. Methods in Ecology and Evolution, 9: 305-319.
#' @references Mammola, S. (2019) Assessing similarity of n-dimensional hypervolumes: Which metric to use?. Journal of Biogeography, 46: 2012-2023.
#' @examples \dontrun{
#' comm <- rbind(c(1,1,1,1,1), c(1,1,1,1,1), c(0,0,1,1,1),c(0,0,1,1,1))
#' colnames(comm) = c("SpA","SpB","SpC","SpD", "SpE")
#' rownames(comm) = c("Site 1","Site 2","Site 3","Site 4")
#'
#' trait <- cbind(c(2.2,4.4,6.1,8.3,3),c(0.5,1,0.5,0.4,4),c(0.7,1.2,0.5,0.4,5),c(0.7,2.2,0.5,0.3,6))
#' colnames(trait) = c("Trait 1","Trait 2","Trait 3","Trait 4")
#' rownames(trait) = colnames(comm)
#'
#' hvlist = kernel.build(comm, trait)
#' kernel.similarity(hvlist)
#' hvlist = kernel.build(comm, trait, axes = 0.9)
#' kernel.similarity(hvlist)
#' }
#' @export
kernel.similarity <- function(comm) {
  
  #check if right data is provided
  if (!(class(comm) %in% c("HypervolumeList")))
    stop("A HypervolumeList is needed as input data.")
  
  #create matrices to store results
  nComm  <- length(comm@HVList)
  dist_c <- matrix(nrow = nComm, ncol = nComm, NA)
  dist_m <- matrix(nrow = nComm, ncol = nComm, NA)
  int    <- matrix(nrow = nComm, ncol = nComm, NA)
  jac    <- matrix(nrow = nComm, ncol = nComm, NA)
  sor    <- matrix(nrow = nComm, ncol = nComm, NA)
  
  #calculate similarity values and give them a name
  commNames <- c()
  for (i in 1:nComm){
    for(j in i:nComm){
      dst_cent  <- hypervolume_distance(comm@HVList[[i]], comm@HVList[[j]], type = "centroid", num.points.max = 1000, check.memory = TRUE)
      dst_min   <- hypervolume_distance(comm@HVList[[i]], comm@HVList[[j]], type = "minimum", check.memory = FALSE)
      set       <- hypervolume_set(comm@HVList[[i]], comm@HVList[[j]], check.memory = FALSE, verbose = FALSE)
      dist_c[j,i] <- dst_cent
      dist_m[j,i] <- dst_min
      int[j,i]   <- get_volume(set)[[3]]
      jac[j,i]   <- hypervolume_overlap_statistics(set)[1]
      sor[j,i]   <- hypervolume_overlap_statistics(set)[2]
    }
    commNames <- c(commNames, comm@HVList[[i]]@Name)
    message(paste("Similarity values of the hypervolume", as.character(i), "out of", as.character(nComm), "have been calculated.\n"))
  }
  
  #finalize
  rownames(dist_c) <- colnames(dist_c) <- rownames(dist_m) <- colnames(dist_m) <- rownames(jac) <- colnames(jac) <- rownames(sor) <- colnames(sor) <- rownames(int) <- colnames(int) <- commNames
  similarity <- list(Distance_centroids = as.dist(dist_c), Minimum_distance = as.dist(dist_m), Intersection = as.dist(int), Jaccard = as.dist(jac), Sorensen = as.dist(sor))
  return(similarity)
}

#' Hotspots in hypervolumes.
#' @description Identify hotspots in kernel density hypervolumes based on minimum volume needed to cover a given proportion of random points.
#' @param comm A 'Hypervolume' or 'HypervolumesList' object, preferably built using function kernel.build.
#' @param prop Proportion of random points to be included.
#' @details Estimates the hotspots of one or more communities using kernel density hypervolumes as in Carmona et al. (2021).
#' @return A 'Hypervolume' or 'HypervolumeList' with the hotspots of each site.
#' @references Carmona, C.P., et al. (2021) Erosion of global functional diversity across the tree of life. Science Advances, 7: eabf2675. DOI: 10.1126/sciadv.abf2675
#' @examples \dontrun{
#' comm = rbind(c(1,3,0,5,3), c(3,2,5,0,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' plot(hv)
#' kernel.alpha(hv)
#' 
#' hot = kernel.hotspots(hv, 0.5)
#' plot(hot)
#' kernel.alpha(hot)
#' 
#' hvlist = kernel.build(comm, trait)
#' hot = kernel.hotspots(hvlist, 0.1)
#' kernel.alpha(hot)
#' }
#' @export
kernel.hotspots <- function(comm, prop = 0.5){
  
  #check if right data is provided
  if (!(class(comm) %in% c("HypervolumeList", "Hypervolume")))
    stop("A Hypervolume or HypervolumeList is needed as input data.")
  
  #single volume
  if (class(comm) == "Hypervolume"){
    
    hot = hypervolume_threshold(comm, num.thresholds = 1000, quantile.requested = prop, quantile.requested.type = "probability", plot = FALSE, verbose = FALSE)
    hot = hot$HypervolumesThresholded
    hot@Name = comm@Name

  } else if (class(comm) == "HypervolumeList") {
    #if list call this same function
    for (i in 1:length(comm@HVList)){
      newHot <- kernel.hotspots(comm@HVList[[i]], prop)
      cat(paste("Hotspot", as.character(i), "out of", as.character(length(comm@HVList)), "has been calculated.\n"))
      if(i == 1){
        hot <- hypervolume_join(newHot)
      } else {
        hot <- hypervolume_join(hot,newHot)
      }
    }
  }

  #return hot
  return(hot)
}

#' Community Weighted Mean.
#' @description Average value of each of a series of traits in multiple communities.
#' @param comm A sites x species matrix, with incidence or abundance data about the species in the community.
#' @param trait A species x traits matrix, with trait values for each species in comm.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis. If not specified, default is TRUE.
#' @param na.rm Remove NA values before calculating cwm.
#' @details Community weighted mean is used to compare communities in terms of their "typical" trait values.
#' @return A sites x trait matrix with mean value per site and trait.
#' @examples comm <- matrix(c(2,5,0,0,0,1,1,0,0,0,0,1,2,0,0,0,0,0,10,1), nrow = 4, ncol = 5, byrow = TRUE)
#' rownames(comm) = c("Site1","Site2","Site3","Site4")
#' colnames(comm) = c("Sp1","Sp2","Sp3","Sp4","Sp5")
#' trait <- data.frame(Trait1 = c(1,0,0,2,0), Trait2 = c(rep("A",2), rep("B",3)))
#' rownames(trait) = colnames(comm)
#' cwm(comm, trait)
#' cwm(comm, trait, FALSE)
#' @export
cwm <- function(comm, trait, abund = TRUE, na.rm = FALSE){
  trait = dummy(trait)
  if(!abund)
      comm[comm > 1] = 1
  nSites = nrow(comm)
  nTraits = ncol(trait)
  nSp = rowSums(comm)
  results = matrix(NA, nrow = nSites, ncol = nTraits)
  rownames(results) = rownames(comm)
  colnames(results) = colnames(trait)
  for (s in 1:nSites)
    for (t in 1:nTraits)
      results[s, t] = sum(comm[s,] * trait[,t], na.rm = na.rm) / nSp[s]
  return(results)
}

#' Community Weighted Dispersion.
#' @description Standard deviation value of each of a series of traits in multiple communities.
#' @param comm A sites x species matrix, with incidence or abundance data about the species in the community.
#' @param trait A species x traits matrix, with trait values for each species in comm.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis. If not specified, default is TRUE.
#' @param na.rm Remove NA values before calculating cwd.
#' @details Community weighted dispersion is used to compare communities in terms of their dispersion of trait values around a mean, reflecting individual trait variability or diversity.
#' @return A sites x trait matrix with sd value per site and trait.
#' @examples comm <- matrix(c(2,5,0,0,0,1,1,0,0,0,0,1,2,0,0,0,0,0,10,1), nrow = 4, ncol = 5, byrow = TRUE)
#' rownames(comm) = c("Site1","Site2","Site3","Site4")
#' colnames(comm) = c("Sp1","Sp2","Sp3","Sp4","Sp5")
#' trait <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 5, ncol = 4, byrow = TRUE)
#' rownames(trait) = colnames(comm)
#' colnames(trait) = c("Trait1","Trait2","Trait3","Trait4")
#' cwd(comm, trait)
#' cwd(comm, trait, FALSE)
#' @export
cwd <- function(comm, trait, abund = TRUE, na.rm = FALSE){
  trait = dummy(trait)
  if(!abund)
    comm[comm > 1] = 1
  nSites = nrow(comm)
  nTraits = ncol(trait)
  nSp = rowSums(comm)     
  results = matrix(NA, nrow = nSites, ncol = nTraits)
  rownames(results) = rownames(comm)
  colnames(results) = colnames(trait)
  cwmean = cwm(comm, trait, abund, na.rm)
  for (s in 1:nSites)
    for (t in 1:nTraits)
      results[s, t] = (sum(comm[s,] * (trait[,t] - cwmean[s,t])^2, na.rm = na.rm) / nSp[s])^0.5
  return(results)
}

#' Community Weighted Evenness.
#' @description Evenness value of each of a series of traits in multiple communities.
#' @param comm A sites x species matrix, with incidence or abundance data about the species in the community.
#' @param trait A species x traits matrix, with trait values for each species in comm.
#' @param func Calculate evenness using Camargo (1993, default) or Bulla (1994) index.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis. If not specified, default is TRUE.
#' @param na.rm Remove NA values before calculating cwe.
#' @details Community weighted evenness is used to compare communities in terms of their evenness of trait values, reflecting trait abundance and distances between values.
#' @return A sites x trait matrix with evenness value per site and trait.
#' @references Bulla, L. (1994) An index of evenness and its associated diversity measure. Oikos, 70: 167-171.
#' @references Camargo, J.A. (1993) Must dominance increase with the number of subordinate species in competitive interactions? Journal of Theoretical Biology, 161: 537-542.
#' @examples comm <- matrix(c(1,1,1,1,0,1,1,0,0,0,0,1,2,0,0,0,0,0,10,1), nrow = 4, ncol = 5, byrow = TRUE)
#' rownames(comm) = c("Site1","Site2","Site3","Site4")
#' colnames(comm) = c("Sp1","Sp2","Sp3","Sp4","Sp5")
#' trait <- matrix(c(4,1,3,4,2,2,2,1,3,3,2,0,1,4,0,0,5,5,2,1), nrow = 5, ncol = 4, byrow = TRUE)
#' rownames(trait) = colnames(comm)
#' colnames(trait) = c("Trait1","Trait2","Trait3","Trait4")
#' cwe(comm, trait)
#' cwe(comm, trait, abund = FALSE)
#' cwe(comm, trait, "bulla")
#' @export
cwe <- function(comm, trait, func = "camargo", abund = TRUE, na.rm = FALSE){
  trait = dummy(trait)
  if(!abund)
    comm[comm > 1] = 1
  nSites = nrow(comm)
  nTraits = ncol(trait)
  results = matrix(NA, nrow = nSites, ncol = nTraits)
  rownames(results) = rownames(comm)
  colnames(results) = colnames(trait)
  
  for (s in 1:nSites){
    for (t in 1:nTraits){

      #clean stuff for this run
      thisComm = comm[s, comm[s,] > 0]		                   #filter comm
      thisTrait = trait[comm[s,] > 0, t]                     #filter trait values
      thisComm = thisComm[order(thisTrait)]                  #order comm by trait values
      thisTrait = thisTrait[order(thisTrait)]                #order trait by trait values
      
      #if any trait values are similar, merge in same functional unit
      i = 1
      while(i < length(thisTrait)){
        if(thisTrait[i] == thisTrait[i+1]){
          thisComm[i+1] = thisComm[i] + thisComm[i+1]
          thisComm = thisComm[-i]
          thisTrait = thisTrait[-i]
        } else {
          i = i + 1
        }
      }
      
      #if only 1 functional category skip, as evenness does not make sense
      nDist = length(thisComm) - 1                           #number of links
      if(nDist == 0) next

      #if only 2 categories use regular evenness without the functional part
      if(nDist == 1){
        #calculate the observed values as proportional abundance per species
        thisObs = thisComm / sum(thisComm, na.rm = na.rm)
        if(func == "bulla"){
          thisExp = 1 / length(thisComm)
          results[s,t] = as.numeric((sum(apply(cbind(thisObs, rep(thisExp, length(thisObs))), 1, min), na.rm = na.rm) - thisExp) / (1 - thisExp))
        } else if(func == "camargo"){
          results[s,t] = as.numeric(1 - (abs(thisObs[1] - thisObs[2])))
        }
        next
      }

      #if more than 2 categories proceed with the functional part

      #calculate distances between trait values
      disTraits = c()                                        
      for(i in 1:nDist)
        disTraits[i] = thisTrait[i+1] - thisTrait[i]

      #calculate the observed values as proportional abundance per species / distance
      thisObs = c()
      for(i in 1:nDist)											   #cycle through all distances of this site/sample
        thisObs[i] = mean(as.numeric(thisComm[c(i, i+1)]), na.rm = na.rm) / disTraits[i]
      thisObs = thisObs / sum(thisObs, na.rm = na.rm)         #sum all observations to 1
      
      if(func == "bulla"){
        ##calculate the expected values as average length of distances between observations
        thisExp = 1 / nDist
        #calculate evenness as the sum of minimum values between observed and expected with correction from Bulla, 1994
        results[s,t] = (sum(apply(cbind(thisObs, rep(thisExp, length(thisObs))), 1, min), na.rm = na.rm) - thisExp) / (1 - thisExp)
      } else if(func == "camargo"){
        results[s,t] = 0
        for(j in 1:(nDist - 1)){
          for(k in (j + 1):nDist){
            results[s,t] = results[s,t] + abs(thisObs[j] - thisObs[k])
          }
        }
        results[s,t] = 1 - (results[s,t] / (nDist * (nDist - 1) / 2))

      } else {
        stop(sprintf("Function %s not recognized.", func))
      }
    }
  }

  return(results)
}

#' Scaled mean squared error of accumulation curves.
#' @description Accuracy (scaled mean squared error) of accumulation curves compared with a known true diversity value (target).
#' @param accum A matrix resulting from the alpha.accum or beta.accum functions (sampling units x diversity values).
#' @param target The true known diversity value, with which the curve will be compared. If not specified, default is the diversity observed with all sampling units.
#' @details Among multiple measures of accuracy (Walther & Moore 2005) the SMSE presents several advantages, as it is (Cardoso et al. 2014):
#' (i) scaled to true diversity, so that similar absolute differences are weighted according to how much they represent of the real value;
#' (ii) scaled to the number of sampling units, so that values are independent of sample size;
#' (iii) squared, so that small, mostly meaningless fluctuations around the true value are down-weighted; and
#' (iv) independent of positive or negative deviation from the real value, as such differentiation is usually not necessary.
#' For alpha diversity accuracy may also be weighted according to how good the data is predicted to be. The weight of each point in the curve is proportional to its sampling intensity (i.e. n/Sobs).
#' @return Accuracy values (both raw and weighted) for all observed and estimated curves.
#' @references Cardoso, P., Rigal, F., Borges, P.A.V. & Carvalho, J.C. (2014) A new frontier in biodiversity inventory: a proposal for estimators of phylogenetic and functional diversity. Methods in Ecology and Evolution, 5: 452-461.
#' @references Walther, B.A. & Moore, J.L. (2005) The concepts of bias, precision and accuracy, and their use in testing the performance of species richness estimators, with a literature reviewof estimator performance. Ecography, 28, 815-829.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' acc.alpha = alpha.accum(comm1)
#' accuracy(acc.alpha)
#' accuracy(acc.alpha, 10)
#' acc.beta = beta.accum(comm1, comm2, tree)
#' accuracy(acc.beta)
#' accuracy(acc.beta, c(1,1,0))
#' @export
accuracy <- function(accum, target = -1){
  if(ncol(accum) > 5 || accum[nrow(accum), 3] > 1){		#if alpha
		if (target == -1)
			target <- accum[nrow(accum), 3]
		intensTotal = accum[nrow(accum), 2] / accum[nrow(accum), 3]	#sampling intensity = final n / final S
		if(ncol(accum) > 10){              #if non-parametric
      smse <- matrix(0, 13, nrow = 2)
      for (i in 1:nrow(accum)){
      	intensity = accum[i, 2] / accum[i, 3] / intensTotal
      	error = (accum[i,3] - target)^2 / (target^2 * nrow(accum))
      	smse[1,1] <- smse[1,1] + error
      	smse[2,1] <- smse[2,1] + error * intensity
      	for (j in 2:13){
      		error = (accum[i,j+6] - target)^2 / (target^2 * nrow(accum))
      		smse[1,j] <- smse[1,j] + error
      		smse[2,j] <- smse[2,j] + error * intensity
      	}
      }
      rownames(smse) <- c("Raw", "Weighted")
     	colnames(smse) <- c("Obs", "Jack1ab", "Jack1abP", "Jack1in", "Jack1inP", "Jack2ab", "Jack2abP", "Jack2in", "Jack2inP", "Chao1", "Chao1P", "Chao2", "Chao2P")
    }
    else{                              #if curve
      smse <- matrix(0, 5, nrow = 2)
      for (i in 3:nrow(accum)){
      	intensity = accum[i, 2] / accum[i, 3] / intensTotal
      	for (j in 1:5){
          if (!is.na(accum[i,j+2])){
          	error = (accum[i,j+2] - target)^2 / (target^2 * nrow(accum))
          	smse[1,j] <- smse[1,j] + error
          	smse[2,j] <- smse[2,j] + error * intensity
          }
      	}
      }
      rownames(smse) <- c("Raw", "Weighted")
      colnames(smse) <- c("Obs", "Clench", "Exponential", "Rational", "Weibull")
    }
	} else {																						#if beta
		if (target[1] == -1)
			target <- accum[nrow(accum), 2:4]
		smse <- rep(0, 3)
		for (i in 1:nrow(accum)){
			for (j in 1:3)
				smse[j] <- smse[j] + (accum[i,j+1] - target[j])^2
		}
		smse <- smse / nrow(accum)
		smse <- list(Btotal=smse[1], Brepl=smse[2], Brich=smse[3])
		smse <- c(unlist(smse))
	}
	return(smse)
}

#' Slope of accumulation curves.
#' @description This is similar to the first derivative of the curves at each of its points.
#' @param accum A matrix resulting from the alpha.accum or beta.accum functions (sampling units x diversity values).
#' @details Slope is the expected gain in diversity when sampling a new individual. The slope of an accumulation curve, of either observed or estimated diversity, allows verifying if the asymptote has been reached (Cardoso et al. 2011).
#' This is an indication of either the completeness of the inventory (low final slopes of the observed curve indicate high completeness) or reliability of the estimators (stability of the slope around a value of 0 along the curve indicates reliability).
#' @return A matrix of sampling units x slope values.
#' @references Cardoso, P., Pekar, S., Jocque, R. & Coddington, J.A. (2011) Global patterns of guild composition and functional diversity of spiders. PLoS One, 6, e21710.
#' @examples comm1 <- matrix(c(2,2,0,0,0,1,1,0,0,0,0,2,2,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' comm2 <- matrix(c(1,1,0,0,0,0,2,1,0,0,0,0,2,1,0,0,0,0,2,1), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' acc.alpha = alpha.accum(comm1)
#' slope(acc.alpha)
#' acc.beta = beta.accum(comm1, comm2, tree)
#' slope(acc.beta)
#' @export
slope <- function(accum){
	if(ncol(accum) > 5 || accum[nrow(accum), 3] > 1){			#if alpha
		sl <- accum[,-2]
		accum <- rbind(rep(0,ncol(accum)), accum)
		for (i in 1:nrow(sl)){
			sl[i,1] <- i
			for (j in 2:ncol(sl)){
				sl[i,j] <- (accum[i+1,j+1]-accum[i,j+1])/(accum[i+1,2]-accum[i,2])
			}
		}
	} else {																							#if beta
		sl <- accum
		sl[1,] <- 0
		sl[1,1] <- 1
		for (i in 2:nrow(sl)){
			for (j in 2:ncol(sl)){
				sl[i,j] <- (accum[i,j]-accum[i-1,j])
			}
		}
	}
	return(sl)
}

#' Coverage of datasets.
#' @description Coverage is a measure of completeness of a dataset.
#' @param comm A matrix of sites x species with abundance values.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @details Calculated as the estimated proportion of individuals that belong to the species (or phylogenetic, or functional diversity) already collected (Chao and Jost 2012).
#' @return A vector with coverage values per site.
#' @references Chao, A. & Jost, L. (2012). Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size. Ecology, 93: 2533-2547.
#' @examples comm <- matrix(c(2,1,0,0,100,1,2,0,0,3,1,2,4,0,0,0,0,0,2,2), nrow = 4, ncol = 5, byrow = TRUE)
#' tree <- hclust(dist(c(1:5), method="euclidean"), method="average")
#' coverage(comm)
#' coverage(comm, tree)
#' @export
coverage <- function(comm, tree){

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
  if(!missing(tree))
    tree = xTree(tree)
  cover = comm[,1,drop = FALSE]

  for(r in 1:nrow(comm)){
    data <- comm[r,,drop = FALSE]
    n <- nobs(data, tree)
    s1 <- srare(data, tree, 1)
    s2 <- srare(data, tree, 2)
    cover[r,1] = 1 - (s1/n)*(((n-1)*s1)/((n-1)*s1+2*s2))
  }
  return(cover)
}

#' Optimization of alpha diversity sampling protocols.
#' @description Optimization of alpha diversity sampling protocols when different methods and multiple samples per method are available.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param base A vector defining a base protocol from which to build upon (complementarity analysis) (length must be equal to number of methods).
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @param prog Present a text progress bar in the R console.
#' @details Often a combination of methods allows sampling maximum plot diversity with minimum effort, as it allows sampling different sub-communities, contrary to using single methods.
#' Cardoso (2009) proposed a way to optimize the number of samples per method when the target is to maximize sampled alpha diversity. It is applied here for TD, PD and FD, and for one or multiple sites simultaneously.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix of samples x methods (values being optimum number of samples per method). The last column is the average alpha diversity value, rescaled to 0-1 if made for several sites, where 1 is the true diversity of each site.
#' @references Cardoso, P. (2009) Standardization and optimization of arthropod inventories - the case of Iberian spiders. Biodiversity and Conservation, 18, 3949-3962.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2), c(4,3,2))
#' colnames(comm) <- c("Sp1","Sp2","Sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels <- colnames(comm)
#' optim.alpha(comm,,methods)
#' optim.alpha(comm, tree, methods)
#' optim.alpha(comm,, methods = methods, base = c(0,0,1), runs = 100)
#' @export
optim.alpha <- function(comm, tree, methods, base, runs = 0, prog = TRUE){

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

	##preliminary stats
	methods <- as.vector(t(methods))
	nSamples <- length(methods)							          ##number of samples
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)							          ##number of methods
	if (missing(base))										            ##if no samples to start with for complementarity analysis
		samples <- rep(0,metNum)
	else
		samples <- base
	nMiss <- nSamples - sum(samples)          				##number of samples missing
	nSamplesMet <- rep(0,metNum)						          ##samples per method
	for (m in 1:metNum)
		nSamplesMet[m] <- sum(methods == metUnique[m])

	##accumulation process
	if (prog)
	  pb <- txtProgressBar(max = nMiss + 1, style = 3)
	div <- rep(0, nMiss + 1)									      	##diversity along the optimal accumulation curve
	if (sum(samples) > 0)
		div[1] <- optim.alpha.stats(comm, tree, methods, samples, runs)
	if (prog)
	  setTxtProgressBar(pb, 1)
	for (s in 2:(nMiss+1)){
		samples <- rbind (samples, rep(0,metNum))
		samples[s,] <- samples[s-1,]
		metValue <- rep(0, metNum)										  #diversity when adding each method
		for (m in 1:metNum){
			if (samples[s,m] < nSamplesMet[m]){
				samples[s,m] <- samples[s,m] + 1
				metValue[m] <- optim.alpha.stats(comm, tree, methods, samples[s,], runs)
				samples[s,m] <- samples[s,m] - 1
			}
		}
		div[s] <- max(metValue)
		best <- which(metValue == div[s])
		if (length(best) > 1)
			best = best[sample(1:length(best),1)]					#if tie, choose one of the best methods randomly
		samples[s, best] <- samples[s, best] + 1
		if (prog) setTxtProgressBar(pb, s)
	}
	if (prog) close(pb)
	colnames(samples) <- metUnique
	rownames(samples) <- (0:nMiss+sum(samples[1,]))
	samples <- cbind(samples, div)
	return(samples)
}

#' Efficiency statistics for alpha-sampling.
#' @description Average alpha diversity observed with a given number of samples per method.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param samples A vector defining the number of samples per method to be evaluated (length must be equal to number of methods).
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @details Different combinations of samples per method allow sampling different sub-communities.
#' This function allows knowing the average TD, PD or FD values for a given combination, for one or multiple sites simultaneously.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A single average alpha diversity value. Rescaled to 0-1 if made for several sites, where 1 is the true diversity of each site.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2), c(4,3,2))
#' colnames(comm) <- c("Sp1","Sp2","Sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels <- colnames(comm)
#' optim.alpha.stats(comm,,methods, c(1,1,1))
#' optim.alpha.stats(comm, tree, methods = methods, samples = c(0,0,1), runs = 100)
#' @export
optim.alpha.stats <- function(comm, tree, methods, samples, runs = 0){

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

	##preliminary stats
	if (!missing(tree)){
	  comm = reorderComm(comm, tree)
	  tree <- xTree(tree)
	}
	if(length(dim(comm)) == 3)					                          ##number of sites
		nSites <- dim(comm)[3]
	else
		nSites <- 1
	methods <- as.vector(t(methods))
	metUnique <- as.vector(t(unique(methods)))				            ##list of methods
	metNum <- length(metUnique)					                          ##number of methods
	div <- 0														                          ##average diversity obtained using this particular combination of samples per method

	for (i in 1:nSites){
		if (nSites > 1){
			site <- as.matrix(comm[,,i])
			true <- sobs(site, tree) 			                          	##true diversity of each site
		} else {
			site <- as.matrix(comm)
			true <- 1
		}
	  
    for (r in 1:runs){
		 	addSample <- rep(0, ncol(comm))
		 	for (m in 1:metNum){
		 		if (samples[m] > 0){
		 			filterList <- site[which(methods == metUnique[m]),,drop=F]						##filter by method m
		 			filterList <- filterList[sample(nrow(filterList),samples[m]),,drop=F]	##randomly select rows
		 			addSample <- rbind(addSample, filterList)															##add random samples
		 		}
		 	}
		 	div <- div + sobs(addSample, tree) / runs / nSites / true
		}
	}
	return(div)
}

#' Optimization of beta diversity sampling protocols.
#' @description Optimization of beta diversity sampling protocols when different methods and multiple samples per method are available.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param base Allows defining a base mandatory protocol from which to build upon (complementarity analysis). It should be a vector with length = number of methods.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis.
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @param prog Present a text progress bar in the R console.
#' @details Often, comparing differences between sites or the same site along time (i.e. measure beta diversity) it is not necessary to sample exhaustively. A minimum combination of samples targeting different sub-communities (that may behave differently) may be enough to perceive such differences, for example, for monitoring purposes.
#' Cardoso et al. (in prep.) introduce and differentiate the concepts of alpha-sampling and beta-sampling. While alpha-sampling optimization implies maximizing local diversity sampled (Cardoso 2009), beta-sampling optimization implies minimizing differences in beta diversity values between partially and completely sampled communities.
#' This function uses as beta diversity measures the Btotal, Brepl and Brich partitioning framework (Carvalho et al. 2012) and respective generalizations to PD and FD (Cardoso et al. 2014).
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix of samples x methods (values being optimum number of samples per method). The last column is the average absolute difference from real beta.
#' @references Cardoso, P. (2009) Standardization and optimization of arthropod inventories - the case of Iberian spiders. Biodiversity and Conservation, 18, 3949-3962.
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Cardoso, P., et al. (in prep.) Optimal inventorying and monitoring of taxon, phylogenetic and functional diversity.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm3 <- matrix(c(2,0,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2, comm3), c(4,3,3))
#' colnames(comm) <- c("sp1","sp2","sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels <- colnames(comm)
#' optim.beta(comm, methods = methods, runs = 100)
#' optim.beta(comm, tree, methods = methods, abund = FALSE, base = c(0,0,1), runs = 100)
#' @export
optim.beta <- function(comm, tree, methods, base, abund = TRUE, runs = 0, prog = TRUE){

	##preliminary stats
	methods <- as.vector(t(methods))
	nSamples <- length(methods)							          ##number of samples
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)						          	##number of methods

	if (missing(base))										            ##if no samples to start with
		samples <- rep(0,metNum)
	else
		samples <- base
	nMiss <- nSamples - sum(samples)				          ##number of samples missing
	nSamplesMet <- rep (0, metNum)					          ##samples per method
	for (m in 1:metNum)
		nSamplesMet[m] <- sum(methods == metUnique[m])

	##accumulation process
	if (prog) pb <- txtProgressBar(max = nMiss+1, style = 3)
	diff <- rep(0,nMiss+1)														#absolute difference along the optimal accumulation curve
	diff[1] <- optim.beta.stats(comm, tree, methods, samples, abund, runs)
	if (prog) setTxtProgressBar(pb, 1)
	if (diff[1] == "NaN")
		diff[1] = 1
	for (s in 2:(nMiss+1)){
	  samples <- rbind (samples, rep(0,metNum))
		samples[s,] <- samples[s-1,]
		metValue <- rep(1, metNum)										  #absolute difference when adding each method
		for (m in 1:metNum){
			if (samples[s,m] < nSamplesMet[m]){
				samples[s,m] <- samples[s,m] + 1
				metValue[m] <- optim.beta.stats(comm, tree, methods, samples[s,], abund, runs)
				samples[s,m] <- samples[s,m] - 1
			}
		}
		diff[s] <- min(metValue)
		best <- which(metValue == diff[s])
		if (length(best) > 1)
			best = best[sample(1:length(best),1)]					#if tie, choose one of the best methods randomly
		samples[s, best] <- samples[s, best] + 1
		if (prog) setTxtProgressBar(pb, s)
	}
	if (prog) close(pb)
	colnames(samples) <- metUnique
	rownames(samples) <- (0:nMiss+sum(samples[1,]))
	samples <- cbind(samples, diff)
	return(samples)
}

#' Efficiency statistics for beta-sampling.
#' @description Average absolute difference between sampled and real beta diversity when using a given number of samples per method.
#' @param comm A samples x species x sites array, with either abundance or incidence data.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param methods A vector specifying the method of each sample (length must be equal to nrow(comm))
#' @param samples The combination of samples per method we want to test. It should be a vector with length = number of methods.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis.
#' @param runs Number of random permutations to be made to the sample order. Default is 1000.
#' @details Different combinations of samples per method allow sampling different sub-communities.
#' This function allows knowing the average absolute difference between sampled and real beta diversity for a given combination, for one or multiple sites simultaneously.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A single average absolute beta diversity difference value.
#' @examples comm1 <- matrix(c(1,1,0,2,4,0,0,1,2,0,0,3), nrow = 4, ncol = 3, byrow = TRUE)
#' comm2 <- matrix(c(2,2,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm3 <- matrix(c(2,0,0,3,1,0,0,0,5,0,0,2), nrow = 4, ncol = 3, byrow = TRUE)
#' comm <- array(c(comm1, comm2, comm3), c(4,3,3))
#' colnames(comm) <- c("sp1","sp2","sp3")
#' methods <- c("Met1","Met2","Met2","Met3")
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels <- colnames(comm)
#' optim.beta.stats(comm,,methods, c(1,1,1))
#' optim.beta.stats(comm, tree, methods = methods, samples = c(0,0,1), runs = 100)
#' @export
optim.beta.stats <- function(comm, tree, methods, samples, abund = TRUE, runs = 0){

  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)

	##preliminary stats
	if(length(dim(comm)) == 3){					              ##number of sites
		nSites <- dim(comm)[3]
	}else{
		return(message("Need sample data from at least two sites to perform analyses."))
	}
	methods <- as.vector(t(methods))
	metUnique <- as.vector(t(unique(methods)))				##list of methods
	metNum <- length(metUnique)					              ##number of methods
	diff <- 0													              	##average absolute difference between observed and true diversity obtained using this particular combination of samples per method
  
	if(!missing(tree))
  	comm = reorderComm(comm, tree)
	
	##calculate true beta values
	sumComm <- matrix(0, nrow = nSites, ncol = ncol(comm))
	for (i in 1:nSites){
		sumComm[i,] <- colSums(comm[,,i])
	}
	true <- beta(sumComm, tree, abund)

	##calculate absolute difference between sampled and true beta values
	for (r in 1:runs){
		sumComm <- matrix(0, nrow = nSites, ncol = ncol(comm))
		for (m in 1:metNum){
			if (samples[m] > 0){
				filterList <- comm[which(methods == metUnique[m]),,,drop=F] 							##filter by method m
				filterList <- filterList[sample(nrow(filterList),samples[m]),,,drop=F]		##randomly select rows
				for (i in 1:nSites){
					sumComm[i,] <- sumComm[i,] + colSums(filterList[,,i,drop=F])
				}
			}
		}
		sampleBeta <- beta(sumComm, tree, abund)
		for(i in 1:3){
			diff <- diff + mean(abs(sampleBeta[[i]] - true[[i]])) / 3 / runs
		}
	}
	return(diff)
}

#' Optimization of spatial sampling.
#' @description Optimization of sampling site distribution in space based on environmental (or other) variables.
#' @param layers A Raster* object (typically a multi-layer type: RasterStack or RasterBrick).
#' @param n The number of intended sampling sites (clusters).
#' @param latlong Boolean indicating whether latitude and longitude should be taken into account when clustering.
#' @param clusterMap Boolean indicating whether to build a new raster with clusters.
#' @details Optimizing the selection of sampling sites often requires maximizing the environmental diversity covered by them.
#' One possible solution to this problem, here adopted, is performing a k-means clustering using environmental data and choosing the sites closest to the multidimensional environmental centroid of each cluster for sampling (Jimenez-Valverde & Lobo 2004)
#' @return Either a matrix of cells x clusters (also indicating distance to centroid, longitude and latitude of each cell) or a list with such matrix plus the clusterMap.
#' @references Jimenez-Valverde, A., & Lobo, J. M. (2004) Un metodo sencillo para seleccionar puntos de muestreo con el objetivo de inventariar taxones hiperdiversos: el caso practico de las familias Araneidae y Thomisidae (Araneae) en la comunidad de Madrid, Espana. Ecologia, 18: 297-305.
#' @export
optim.spatial <- function(layers, n, latlong = TRUE, clusterMap = TRUE){
  for(i in 1:length(layers))              ##transform all layers to a scale [0,1]
    layers[[i]] <- (layers[[i]]-cellStats(layers[[i]], min))/(cellStats(layers[[i]], max)-cellStats(layers[[i]], min))
  dataMat <- as.matrix(layers)
	dataMat <- dataMat[complete.cases(dataMat),]
	dataMat <- cbind(dataMat, rasterToPoints(layers[[1]])[,1:2])					##add latlong
	if (latlong)
		res <- kmeans(dataMat, n)        ##do k-means
	else
		res <- kmeans(dataMat[,-c((ncol(dataMat)-1),ncol(dataMat))], n)        ##do k-means

	cl = c()
	for(c in 1:n){
		cData <- dataMat[res$cluster==c,]           #filter to cluster c
		cCenter <- res$centers[c,]
		dist2centroid <- c()
		for(r in 1:nrow(cData))
			dist2centroid[r] = dist(rbind(cData[r,], cCenter))
		cData <- cbind(rep(c, nrow(cData)), dist2centroid, cData[,ncol(cData)], cData[,(ncol(cData)-1)])
		colnames(cData) <- c("cluster", "dist2centroid", "lat", "long")
		cData <- cData[sort.list(cData[,2]), ]
		cl <- rbind(cl, cData)
	}

	#output raster with clusters
	if(clusterMap){
		map <- rasterize(rasterToPoints(layers[[1]])[,1:2], layers[[1]], res$cluster)
		names(map) <- "clusters"
		cl <- list(cl, map)
	}
	return(cl)
}

#' Maps of alpha diversity (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Observed alpha diversity using rasters of species distributions (presence/absence).
#' @param layers A Raster* object of species distributions (typically a multi-layer type: RasterStack or RasterBrick).
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @details TD is equivalent to species richness. Calculations of PD and FD are based on Faith (1992) and Petchey & Gaston (2002, 2006), which measure PD and FD of a community as the total branch length of a tree linking all species represented in such community.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in layers must be the same as in tree.
#' @return A raster object representing richness in space.
#' @references Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @references Petchey, O.L. & Gaston, K.J. (2002) Functional diversity (FD), species richness and community composition. Ecology Letters, 5, 402-411.
#' @references Petchey, O.L. & Gaston, K.J. (2006) Functional diversity: back to basics and looking forward. Ecology Letters, 9, 741-758.
#' @examples sp1 <- raster::raster(matrix(c(NA,1,1,1,1,0,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' sp2 <- raster::raster(matrix(c(0,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3, byrow = TRUE))
#' sp3 <- raster::raster(matrix(c(0,0,0,1,1,1,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' spp <- raster::stack(sp1, sp2, sp3)
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels = c("Sp1", "Sp2", "Sp3")
#' names(spp) = tree$labels
#' raster.alpha(spp)
#' raster.alpha(spp, tree)
#' @export
raster.alpha <- function(layers, tree){
	res = raster::raster(matrix(NA, nrow = nrow(layers), ncol = ncol(layers)))
	for(r in 1:nrow(layers)){
		for(c in 1:ncol(layers)){
			if(is.na(sum(layers[r,c])))
				res[r,c] = NA
			else
				res[r,c] = alpha(layers[r,c], tree)
		}
	}
	return(res)
}

#' Maps of beta diversity (Taxon, Phylogenetic or Functional Diversity - TD, PD, FD).
#' @description Observed beta diversity using rasters of species distributions (presence/absence or abundance).
#' @param layers A Raster* object of species distributions (typically a multi-layer type: RasterStack or RasterBrick).
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param func Partial match indicating whether the Jaccard or Soerensen family of beta diversity measures should be used. If not specified, default is Jaccard.
#' @param neighbour Either 8 (default) or 4 cells considered to calculate beta diversiy of each focal cell.
#' @param abund A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis.
#' @details The beta diversity measures used here follow the partitioning framework independently developed by Podani & Schmera (2011) and Carvalho et al. (2012)
#' and later expanded to PD and FD by Cardoso et al. (2014), where Btotal = Brepl + Brich.
#' Btotal = total beta diversity, reflecting both species replacement and loss/gain;
#' Brepl = beta diversity explained by replacement of species alone; Brich = beta diversity explained by species loss/gain (richness differences) alone.
#' PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric). The path to the root of the tree is always included in calculations of PD and FD.
#' The number and order of species in layers must be the same as in tree.
#' @return A raster.stack object representing Btotal, Brepl and Brich in space.
#' @references Cardoso, P., Rigal, F., Carvalho, J.C., Fortelius, M., Borges, P.A.V., Podani, J. & Schmera, D. (2014) Partitioning taxon, phylogenetic and functional beta diversity into replacement and richness difference components. Journal of Biogeography, 41, 749-761.
#' @references Carvalho, J.C., Cardoso, P. & Gomes, P. (2012) Determining the relative roles of species replacement and species richness differences in generating beta-diversity patterns. Global Ecology and Biogeography, 21, 760-771.
#' @references Gotelli, N.J. & Colwell, R.K. (2001) Quantifying biodiversity: procedures and pitfalls in the measurement and comparison of species richness. Ecology Letters, 4, 379-391.
#' @references Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. Oikos, 120, 1625-1638.
#' @examples sp1 <- raster::raster(matrix(c(NA,1,1,1,1,0,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' sp2 <- raster::raster(matrix(c(0,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3, byrow = TRUE))
#' sp3 <- raster::raster(matrix(c(0,0,0,1,1,1,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' spp <- raster::stack(sp1, sp2, sp3)
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels = c("Sp1", "Sp2", "Sp3")
#' names(spp) = tree$labels
#' raster.beta(spp)
#' raster.beta(spp, tree)
#' @export
raster.beta <- function(layers, tree, func = "jaccard", neighbour = 8, abund = FALSE){
	resTotal = raster::raster(matrix(NA, nrow = nrow(layers), ncol = ncol(layers)))
	resRepl = resTotal
	resRich = resTotal
	for(c in 1:(raster::ncell(layers))){
		if(is.na(sum(layers[c]))){
			resTotal[c] = NA
			resRepl[c] = NA
			resRich[c] = NA
		} else {
			betaValue = matrix(ncol=3)
			adj = raster::adjacent(layers, c, neighbour)[,2]
			for(a in adj)
				if(!is.na(sum(layers[a])))
					betaValue = rbind(betaValue, beta(rbind(layers[c], layers[a]), tree, func = func, abund = abund))
			betaValue = betaValue[-1,]
			resTotal[c] = mean(unlist(betaValue[,1]))
			resRepl[c] = mean(unlist(betaValue[,2]))
			resRich[c] = mean(unlist(betaValue[,3]))
		}
	}
	res = raster::stack(resTotal, resRepl, resRich)
	names(res) = c("Btotal", "Brepl", "Brich")
	return(res)
}

#' Maps of phylogenetic/functional dispersion of species or individuals.
#' @description Average dissimilarity between any two species or individuals randomly chosen in a community using rasters of species distributions (presence/absence or abundance).
#' @param layers A Raster* object of species distributions (typically a multi-layer type: RasterStack or RasterBrick).
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist object representing the phylogenetic or functional distance between species.
#' @param func Calculate dispersion using originality (default), uniqueness or contribution.
#' @param abund A boolean (T/F) indicating whether dispersion should be calculated using individuals (T) or species (F).
#' @param relative A boolean (T/F) indicating whether dispersion should be relative to the maximum distance between any two species in the tree or distance matrix.
#' @details If abundance data is used and a tree is given, dispersion is the quadratic entropy of Rao (1982).
#' If abundance data is not used but a tree is given, dispersion is the phylogenetic dispersion measure of Webb et al. (2002).
#' @return A raster object representing dispersion in space.
#' @references Rao, C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21: 24-43.
#' @references Webb, C.O., Ackerly, D.D., McPeek, M.A. & Donoghue, M.J. (2002) Phylogenies and community ecology. Annual Review of Ecology and Systematics, 33: 475-505.
#' @examples sp1 <- raster::raster(matrix(c(NA,1,1,1,1,0,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' sp2 <- raster::raster(matrix(c(0,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3, byrow = TRUE))
#' sp3 <- raster::raster(matrix(c(0,0,0,1,1,1,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' spp <- raster::stack(sp1, sp2, sp3)
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels = c("Sp1", "Sp2", "Sp3")
#' names(spp) = tree$labels
#' raster.dispersion(spp, tree)
#' @export
raster.dispersion <- function(layers, tree, distance, func = "originality", abund = FALSE, relative = FALSE){
	res = raster::raster(matrix(NA, nrow = nrow(layers), ncol = ncol(layers)))
	for(r in 1:nrow(layers)){
		for(c in 1:ncol(layers)){
			if(is.na(sum(layers[r,c])))
				res[r,c] = NA
			else
				res[r,c] = dispersion(layers[r,c], tree, distance, func, abund, relative)
		}
	}
	return(res)
}

#' Maps of phylogenetic/functional evenness of species or individuals.
#' @description Regularity of distance and abundance between any two species in a community using rasters of species distributions (presence/absence or abundance).
#' @param layers A Raster* object of species distributions (typically a multi-layer type: RasterStack or RasterBrick).
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree.
#' @param distance A dist object representing the phylogenetic or functional distance between species.
#' @param method Calculate dispersion using "expected" values (default) or values based on "contribution" of species to the tree.
#' @param func Calculate dispersion using "Camargo" (default) or "Bulla" index.
#' @param abund A boolean (T/F) indicating whether evenness should be calculated using abundance data.
#' @details Evenness is calculated based on the index of Bulla (1994) using the values of both edge lengths in the tree and their abundance.
#' @details If no tree or distance is provided the result is the original index of Bulla with correction.
#' @return A raster object representing evenness in space.
#' @references Bulla, L. (1994) An index of evenness and its associated diversity measure. Oikos, 70: 167-171.
#' @examples sp1 <- raster::raster(matrix(c(NA,1,1,1,1,0,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' sp2 <- raster::raster(matrix(c(0,0,0,0,1,1,1,1,1), nrow = 3, ncol = 3, byrow = TRUE))
#' sp3 <- raster::raster(matrix(c(0,0,0,1,1,1,0,0,0), nrow = 3, ncol = 3, byrow = TRUE))
#' spp <- raster::stack(sp1, sp2, sp3)
#' tree <- hclust(dist(c(1:3), method="euclidean"), method="average")
#' tree$labels = c("Sp1", "Sp2", "Sp3")
#' names(spp) = tree$labels
#' raster.evenness(spp, tree)
#' @export
raster.evenness <- function(layers, tree, distance, method = "expected", func = "camargo", abund = TRUE){
	res = raster::raster(matrix(NA, nrow = nrow(layers), ncol = ncol(layers)))
	for(r in 1:nrow(layers)){
		for(c in 1:ncol(layers)){
			if(is.na(sum(layers[r,c])) || sum(ifelse(layers[r,c] > 0, 1, 0)) < 2)
				res[r,c] = NA
			else
				res[r,c] = evenness(layers[r,c], tree, distance, method, func, abund)
		}
	}
	return(res)
}

#' Species-abundance distribution (SAD).
#' @description Fits the SAD to community abundance data with possible rarefaction.
#' @param comm Either a vector with the abundance per species, or a sites x species matrix.
#' @param raref An integer specifying the number of individuals for rarefaction (individual based).
#' If raref < 1 no rarefaction is made.
#' If raref = 1 rarefaction is made by the minimum abundance among all sites.
#' If raref > 1 rarefaction is made by the abundance indicated.
#' If not specified, default is 0.
#' @param runs Number of resampling runs for rarefaction. If not specified, default is 100.
#' @details Classes defined as n = 1, 2-3, 4-7, 8-15, .... Rarefaction allows comparison of sites with different total abundances.
#' @return A vector or matrix with the different values per class per community.
#' @examples comm1 <- c(20,1,3,100,30)
#' comm2 <- c(1,2,12,0,45)
#' comm <- rbind(comm1, comm2)
#' sad(comm1)
#' sad(comm)
#' sad(comm, raref = 1)
#' @export
sad <- function(comm, raref = 0, runs = 100){
	if(is.vector(comm))
		comm <- matrix(comm, nrow = 1)
	if(raref == 1)
		raref = min(rowSums(comm))

	nClasses = as.integer(max(log(comm, 2))) + 1
	res = matrix(0, nrow = nrow(comm), ncol = nClasses)
	
	for(i in 1:nrow(comm)){
		#if rarefying data
		if(raref > 0){
			r = matrix(NA, nrow = runs, ncol = nClasses)
			for(j in 1:runs){
				samp = rrarefy(comm[i,], sample = raref)
				samp = samp[samp > 0]
				newSad = sad(samp)
				r[j,1:length(newSad)] = newSad
			}
			res[i,] = apply(r, 2, median, na.rm = TRUE)
			res[i, is.na(res[i,])] = 0
		} else {
			r = c()
			thisComm = comm[i,]
			thisComm = thisComm[thisComm > 0]
			thisComm = as.integer(log(thisComm, 2)) + 1
			for(j in 1:max(thisComm))
				r[j] = sum(thisComm == j)
			res[i,1:length(r)] = r
		}
	}
	return(res)
}

#' Species-area relationship (SAR).
#' @description Fits and compares several of the most supported models for the species (or PD, or FD) -area relationship.
#' @param comm Either a vector with the diversity values per site, or a sites x species matrix.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree (used only to fit the PD or FD-area relationships, requires comm to be a sites x species matrix).
#' @param area A vector with the area per site.
#' @details Larger areas (often islands) usually carry more species. Several formulas were proposed in the past to describe this relationship (Arrhenius 1920, 1921; Gleason 1922).
#' Recently, the same approach began to be used for other measures of diversity, namely phylogenetic (PD) and functional (FD) diversity (Whittaker et al. 2014).
#' The function compares some of the most commonly used and theoretically or empirically suported models.
#' The relationships for PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix with the different model parameters and explanatory power.
#' @references Arrhenius, O. (1920) Distribution of the species over the area. Meddelanden fran Vetenskapsakadmiens Nobelinstitut, 4: 1-6.
#' @references Arrhenius, O. (1921) Species and area. Journal of Ecology, 9: 95-99.
#' @references Gleason, H.A. (1922) On the relation between species and area. Ecology, 3: 158-162.
#' @references Whittaker, R.J., Rigal, F., Borges, P.A.V., Cardoso, P., Terzopoulou, S., Casanoves, F., Pla, L., Guilhaumon, F., Ladle, R. & Triantis, K.A. (2014) Functional biogeography of oceanic islands and the scaling of functional diversity in the Azores. Proceedings of the National Academy of Sciences USA, 111: 13709-13714.
#' @examples div <- c(1,2,3,4,4)
#' comm <- matrix(c(2,0,0,0,3,1,0,0,2,4,5,0,1,3,2,5,1,1,1,1), nrow = 5, ncol = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:4), method="euclidean"), method="average")
#' area <- c(10,40,80,160,160)
#' sar(div,,area)
#' sar(comm,,area)
#' sar(comm,tree,area)
#' @export
sar <- function(comm, tree, area){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
	if(is.vector(comm)){
		div = comm
	} else if (missing(tree)){
		div = alpha(comm)
	} else {
		div = alpha(comm, tree)
	}
  div = as.vector(div)
	if (!missing(tree)){
	  comm = reorderComm(comm, tree)
	  tree <- xTree(tree)
	}
	results <- matrix(NA, 6, 7)
	colnames(results) <- c("c", "z", "r2", "AIC", "\U0394 AIC", "AICc", "\U0394 AICc")
	rownames(results) <- c("Linear", "Linear (origin)", "Exponential", "Exponential (origin)", "Power", "Power (origin)")
	k <- c(3,2,3,2,3,2)
	model <- list()
	model[[1]] <- try(nls(div ~ c + z*area, start = data.frame(c = 0, z = 1)))
	model[[2]] <- try(nls(div ~ z*area, start = data.frame(z = 1)))
	model[[3]] <- try(nls(div ~ c + z*log(area), start = data.frame(c = 0, z = 1)))
	model[[4]] <- try(nls(div ~ z*log(area), start = data.frame(z = 1)))
	model[[5]] <- try(nls(div ~ c + area^z, start = data.frame(c = 0, z = 1)))
	model[[6]] <- try(nls(div ~ area^z, start = data.frame(z = 1)))
	for(m in 1:length(model)){
		if(k[m] == 3){
			results[m,1] <- coef(summary(model[[m]]))[1,1]
			results[m,2] <- coef(summary(model[[m]]))[2,1]
		} else {
			results[m,2] <- coef(summary(model[[m]]))[1,1]
		}
		est <- predict(model[[m]], area=area)
		results[m,3] <- r2(div, est)
		results[m,4] <- aic(div, est, k[m])
		results[m,6] <- aic(div, est, k[m], correct = TRUE)
	}
	for(m in 1:length(model)){
		results[m,5] <- results[m,4] - min(results[,4])
		results[m,7] <- results[m,6] - min(results[,6])
	}
	return(results)
}

#' General dynamic model of oceanic island biogeography (GDM).
#' @description Fits and compares several of the most supported models for the GDM (using TD, PD or FD).
#' @param comm Either a vector with the diversity values per island, or an island x species matrix.
#' @param tree A phylo or hclust object (used only for PD or FD) or alternatively a species x traits matrix or data.frame to build a functional tree (used only to fit the PD or FD GDM, requires comm to be a sites x species matrix).
#' @param area A vector with the area of islands.
#' @param time A vector with the age of islands. If not given, the species-area relationship is returned instead.
#' @details The general dynamic model of oceanic island biogeography was proposed to account for diversity patterns within and across oceanic archipelagos as a function of area and age of the islands (Whittaker et al. 2008).
#' Several different equations have been found to describe the GDM, extending the different SAR models with the addition of a polynomial term using island age and its square (TT2), depicting the island ontogeny.
#' The first to be proposed was an extension of the exponential model (Whittaker et al. 2008), the power model extensions following shortly after (Fattorini 2009; Steinbauer et al. 2013), as was the linear model (Cardoso et al. 2020).
#' The relationships for PD and FD are calculated based on a tree (hclust or phylo object, no need to be ultrametric).
#' @return A matrix with the different model parameters and explanatory power.
#' @references Cardoso, P., Branco, V.V., Borges, P.A.V., Carvalho, J.C., Rigal, F., Gabriel, R., Mammola, S., Cascalho, J. & Correia, L. (2020) Automated discovery of relationships, models and principles in ecology. Frontiers in Ecology and Evolution, 8: 530135.
#' @references Fattorini, S. (2009) On the general dynamic model of oceanic island biogeography. Journal of Biogeography, 36: 1100-1110.
#' @references Steinbauer, M.J, Klara, D., Field, R., Reineking, B. & Beierkuhnlein, C. (2013) Re-evaluating the general dynamic theory of oceanic island biogeography. Frontiers of Biogeography, 5: 185-194.
#' @references Whittaker, R.J., Triantis, K.A. & Ladle, R.J. (2008) A general dynamic theory of oceanic island biogeography. Journal of Biogeography, 35: 977-994.
#' @examples div <- c(1,3,5,8,10)
#' comm <- matrix(c(2,0,0,0,3,1,0,0,2,4,5,0,1,3,2,5,1,1,1,1), nrow = 5, ncol = 4, byrow = TRUE)
#' tree <- hclust(dist(c(1:4), method="euclidean"), method="average")
#' area <- c(10,40,80,160,160)
#' time <- c(1,2,3,4,5)
#' gdm(div,,area,time)
#' gdm(comm,tree,area,time)
#' gdm(div,,area)
#' @export
gdm <- function(comm, tree, area, time){
  
  #convert traits to a tree if needed
  if(!missing(tree) && (is.matrix(tree) || is.data.frame(tree) || is.vector(tree)))
    tree = tree.build(tree)
  
	if(missing(time))
		return(sar(comm,tree,area))
	if(is.vector(comm)){
		div = comm
	} else if (missing(tree)){
		div = alpha(comm)
	} else {
		div = alpha(comm, tree)
	}
  div = as.vector(div)
	if (!missing(tree)){
	  comm = reorderComm(comm, tree)
	  tree <- xTree(tree)
	}

	results <- matrix(NA, 4, 9)
	colnames(results) <- c("c", "z", "x", "y", "r2", "AIC", "\U0394 AIC", "AICc", "\U0394 AICc")
	rownames(results) <- c("Linear", "Exponential", "Power (area)", "Power (area, time)")
	k <- 5
	model <- list()
	model[[1]] <- try(nls(div ~ c + z*area + x*time + y*time^2, start = data.frame(c=1, z=1, x=1, y=0)))
	model[[2]] <- try(nls(div ~ c + z*log(area) + x*time + y*time^2, start = data.frame(c=1, z=1, x=1, y=0)))
	model[[3]] <- try(nls(div ~ exp(c + z*log(area) + x*time + y*time^2), start = data.frame(c=1, z=1, x=1, y=0)))
	model[[4]] <- try(nls(div ~ exp(c + z*log(area) + x*log(time) + y*log(time)^2), start = data.frame(c=1, z=1, x=1, y=0)))
	for(m in 1:length(model)){
		results[m,1] <- coef(summary(model[[m]]))[1,1]
		results[m,2] <- coef(summary(model[[m]]))[2,1]
		results[m,3] <- coef(summary(model[[m]]))[3,1]
		results[m,4] <- coef(summary(model[[m]]))[4,1]
		est <- predict(model[[m]], area=area, time=time)
		results[m,5] <- r2(div, est)
		results[m,6] <- aic(div, est, k)
		results[m,8] <- aic(div, est, k, correct = TRUE)
	}
	for(m in 1:length(model)){
		results[m,7] <- results[m,6] - min(results[,6])
		results[m,9] <- results[m,8] - min(results[,8])
	}
	return(results)
}

#' Interspecific abundance-occupancy relationship (IAOR).
#' @description Fits and compares several of the most supported models for the IAOR.
#' @param comm A sites x species matrix with abundance values.
#' @details Locally abundant species tend to be widespread while locally rare species tend to be narrowly distributed.
#' That is, for a given species assemblage, there is a positive interspecific abundance-occupancy relationship (Brown 1984).
#' This function compares some of the most commonly used and theoretically or empirically suported models (Nachman 1981; He & Gaston 2000; Cardoso et al. 2020).
#' @return A matrix with the different model parameters and explanatory power.
#' @references Brown, J.H. (1984) On the relationship between abundance and distribution of species. American Naturalist, 124: 255-279.
#' @references Cardoso, P., Branco, V.V., Borges, P.A.V., Carvalho, J.C., Rigal, F., Gabriel, R., Mammola, S., Cascalho, J. & Correia, L. (2020) Automated discovery of relationships, models and principles in ecology. Frontiers in Ecology and Evolution, 8: 530135.
#' @references He, F.L. & Gaston, K.J. (2000) Estimating species abundance from occurrence. American Naturalist, 156: 553-559.
#' @references Nachman, G. (1981) A mathematical model of the functional relationship between density and spatial distribution of a population. Journal of Animal Ecology, 50: 453-460.
#' @examples comm <- matrix(c(4,3,2,1,5,4,3,2,3,2,1,0,6,3,0,0,0,0,0,0), nrow = 5, ncol = 4, byrow = TRUE)
#' iaor(comm)
#' @export
iaor <- function(comm){
  results <- matrix(NA, 4, 7)
  colnames(results) <- c("a", "b", "r2", "AIC", "\U0394 AIC", "AICc", "\U0394 AICc")
  rownames(results) <- c("Linear", "Exponential", "Negative Binomial", "SR")
  k <- c(3,3,2,2)
  abund <- colMeans(comm)                   #mean abundance per species (including sites with 0 individuals)
  occup <- colMeans(ifelse(comm>0,1,0))     #proportion occupancy per species

  model <- list()
  model[[1]] <- try(nls(logit(occup) ~ a+b*log(abund), start = data.frame(a = 1, b = 1))) #linear
  model[[2]] <- try(nls(occup ~ 1-exp(a*abund^b), start = data.frame(a = -1, b = 1))) #exponential
  model[[3]] <- try(nls(occup ~ 1-(1+(abund/a))^(0-a), start = data.frame(a = 0))) #negative binomial
  model[[4]] <- try(nls(occup ~ abund/(a+abund), start = data.frame(a = 0))) #SR = Clench with asymptote 1
  for(m in 1:length(model)){
    if(m < 3){
      results[m,1] <- coef(summary(model[[m]]))[1,1]
      results[m,2] <- coef(summary(model[[m]]))[2,1]
    } else {
      results[m,1] <- coef(summary(model[[m]]))[1,1]
    }
    est <- predict(model[[m]], abund=abund)
    if(m==1) est = revLogit(est)
    results[m,3] <- r2(occup, est)
    results[m,4] <- aic(occup, est, k[m])
    results[m,6] <- aic(occup, est, k[m], correct = TRUE)
  }
  for(m in 1:length(model)){
    results[m,5] <- results[m,4] - min(results[,4])
    results[m,7] <- results[m,6] - min(results[,6])
  }
  return(results)
}

#' Create Linnean tree.
#' @description Creates a Linnean tree from taxonomic hierarchy.
#' @param taxa A taxonomic matrix with columns ordered according to linnean hierarchy starting with the highest.
#' @param distance A vector with distances between levels starting with the highest. If not provided distances will be evenly distributed from 1 to 0.
#' @return An hclust with all species.
#' @examples family <- c("Nemesiidae", "Nemesiidae", "Zodariidae", "Zodariidae")
#' genus <- c("Iberesia", "Nemesia", "Zodarion", "Zodarion")
#' species <- c("Imachadoi", "Nungoliant", "Zatlanticum", "Zlusitanicum")
#' taxa <- cbind(family, genus, species)
#' par(mfrow = c(1, 2))
#' plot(linnean(taxa))
#' plot(linnean(taxa, c(2, 0.5, 0.3)))
#' @export
linnean <- function(taxa, distance = NULL){
	if(is.null(distance))
		distance = seq(from = 1, to = 1/ncol(taxa), by = -1*1/ncol(taxa))
	nspp = nrow(taxa)
	distTable = matrix(NA, nrow = nspp, ncol = nspp)
	colnames(distTable) = rownames(distTable) = taxa[,ncol(taxa)]
	for(i in 1:nspp){
		for(j in 1:nspp){
			level = 0
			for(k in 1:ncol(taxa))
				if(taxa[i,k] != taxa[j,k])
					level = level + 1
			if(level == 0)
				distTable[i,j] = 0
			else
				distTable[i,j] = distance[length(distance) - level + 1]
		}
	}
	tree = hclust(as.dist(distTable))
	return(tree)
}

#' Dummify variables.
#' @description Convert factor variables to dummy variables.
#' @param trait A species x traits matrix or data.frame.
#' @param convert A vector of column numbers, usually categorical variables, to be converted to dummy variables.
#' @param weight Indicates whether weights of variables should be returned (TRUE/FALSE) or a vector with weights per variable.
#' @details If convert is given the algorithm will convert these column numbers to dummy variables. Otherwise it will convert all columns with factors or characters.
#' @return A matrix with variables converted or, if weight == TRUE or a vector, a list also with weights.
#' @examples trait = data.frame(length = c(2,4,6,3,1), wing = c("A", "B", "A", "A", "B"))
#' dummy(trait)
#' dummy(trait, weight = TRUE)
#' dummy(trait, convert = 2, weight = c(0.9, 0.1))
#' @export
dummy <- function(trait, convert = NULL, weight = FALSE){
  
  traitNames = colnames(trait)
  if (is.matrix(trait))
    trait <- as.data.frame(trait)
  if (length(weight) > 1){
    if(length(weight) != ncol(trait))
      stop("Length of 'weight' must be equal to the number of columns in 'trait'.")
    vecWeight = weight
  } else {
    vecWeight = rep(1, ncol(trait))
  }
  
  newTrait = c()
  newWeight = c()
  for(i in 1:length(traitNames)){
    if(is.null(convert) || (i %in% convert)){
      if(is.numeric(trait[,i]) && !(i %in% convert)){
        newTrait = cbind(newTrait, trait[,i])
        colnames(newTrait)[i] = traitNames[i]
        newWeight = c(newWeight, vecWeight[i])
      }
      if(is.logical(trait[,i]) && !(i %in% convert)){
        newTrait = cbind(newTrait, ifelse(trait[,i] == T, 1, 0))
        colnames(newTrait)[i] = traitNames[i]
        newWeight = c(newWeight, vecWeight[i])
      }
      if(is.factor(trait[,i]) || is.character(trait[,i]) || (i %in% convert)){
        colTrait <- factor(trait[,i])
        nValues <- length(colTrait)
        nLevels <- length(levels(colTrait))
        y <- matrix(0L, nrow = nValues, ncol = nLevels)
        colnames(y) <- paste(rep(colnames(trait)[i], nLevels), levels(colTrait), sep = ".")
        y[cbind(seq_len(nValues), as.integer(colTrait))] <- 1L
        colTrait = y
        newTrait = cbind(newTrait, colTrait)
        newWeight = c(newWeight, rep((1/nLevels), nLevels) *  vecWeight[i])
      }
    } else {
      newTrait = cbind(newTrait, trait[,i])
      colnames(newTrait)[ncol(newTrait)] = colnames(trait)[i]
      newWeight = c(newWeight, vecWeight[i])
    }
  }
  if(length(weight) > 1 || weight == TRUE)
    return(list(trait = newTrait, weight = newWeight/sum(newWeight)))
  else
    return(newTrait)
}

#' Filling missing data.
#' @description Estimation of missing trait values (NA) based on different methods.
#' @param trait A species x traits matrix (a species or individual for each row and traits as columns).
#' @param method Method for imputing missing data. One of "mean" (mean value of the trait), "median" (median value of the trait), "similar" (input from closest species), "regression" (linear regression), "w_regression" (regression weighted by species distance), or "PCA" (Principal Component Analysis).
#' @param group A vector (string of characters, factorial, etc.) whose values indicate which species belong to the same group as the missing and should be used in the estimation of missing data. If NULL all species will be used.
#' @param weight A hclust, phylo or dist object to calculate the distance between species and use as weights.
#' Note that the order of tip labels in trees or of species in the distance matrix should be the same as the order of species in trait.
#' @param step A boolean (T/F) indicating if a stepwise regression model based on AIC should be performed. Ignored is regression is not used.
#' @details Inputs missing data in the trait matrix based on different methods (see Taugourdeau et al. 2014; Johnson et al. 2021 for comparisons among the performance of different methods). 
#' The simplest approach is the average imputation ("mean" or "median"), calculating the mean/median of the values for that trait based on all the observations that are non-missing. It has the advantage of keeping the same mean and the same sample size, but many disadvantages. 
#' The "similar" method inputs a systematically chosen value from the closest species who has similar values on other variables.
#' The default method is linear regression ("regression"), where the predicted value is obtained by regressing the missing variable on other variables. This preserves relationships among variables involved in the imputation model, but not variability around predicted values (i.e., may lead to extrapolations).
#' The "w_regression" takes into account the relative distance among species in the imputation of missing traits, based on the phylogenetic or functional distance between missing and non-missing species.
#' The "PCA" method performs PCA with incomplete data sensu Podani et al. (2021).
#' Note that for PCA and regressions methods the performance of the prediction increases as the number of collinear traits increase.
#' @return A trait matrix with missing data (NA) filled with predicted values.
#' If method = "PCA" the function returns the standard output of a principal component analysis as a list with:
#' Eigenvalues
#' Positive eigenvalues
#' Positive eigenvalues as percent
#' Square root of eigenvalues
#' Eigenvectors
#' Component scores
#' Variable scores
#' Object scores in a biplot
#' Variable scores in a biplot
#' @references Johnson, T.F., Isaac, N.J., Paviolo, A. & Gonzalez-Suarez, M. (2021). Handling missing values in trait data. Global Ecology and Biogeography, 30: 51-62.
#' @references Podani, J., Kalapos, T., Barta, B. & Schmera, D. (2021). Principal component analysis of incomplete data. A simple solution to an old problem. Ecological Informatics, 101235.
#' @references Taugourdeau, S., Villerd, J., Plantureux, S., Huguenin-Elie, O. & Amiaud, B. (2014). Filling the gap in functional trait databases: use of ecological hypotheses to replace missing data. Ecology and Evolution, 4: 944-958.
#' @examples \dontrun{
#' trait <- iris[,-5]
#' group <- iris[,5]
#' 
#' #Generating some random missing data
#' for (i in 1:10)
#' trait[sample(nrow(trait), 1), sample(ncol(trait), 1)] <- NA
#' 
#' #Estimating the missing data with different methods
#' fill(trait, "mean")
#' fill(trait, "mean", group)
#' fill(trait, "median")
#' fill(trait, "median", group)
#' fill(trait, "similar")
#' fill(trait, "similar", group)
#' fill(trait, "regression", step = FALSE)
#' fill(trait, "regression", group, step = TRUE)
#' fill(trait, "w_regression", step = TRUE)
#' fill(trait, "w_regression", weight = dist(trait), step = TRUE)
#' fill(trait, "PCA")
#' }
#' @export
fill <- function(trait, method = "regression", group = NULL, weight = NULL, step = TRUE) {
  
  n_sp <- nrow(trait) 
  
  #Initial checking of group
  if(is.null(group) == TRUE){
    group <- as.factor(rep("group", n_sp))
  } else {
    if (!(class(group) %in% c("character", "factor", "vector", "numeric")))
      stop("Group should be a vector, string of characters, or a factor.")
    if(length(group) != n_sp)
      stop("The number of elements in the trait matrix should be equal to the number of elements in the grouping vector.")
    if(any(is.na(group) == TRUE))
      stop("The group cannot contain missing data.")
    group <- as.factor(group) 
  }
  
  #Convert characters to factors if needed
  trait[sapply(trait, is.character)] <- lapply(trait[sapply(trait, is.character)], as.factor)
  
  #Convert <NA> to NA if needed
  trait[] <- lapply(trait, function(x) {
    is.na(levels(x)) <- levels(x) == "NA"
    x
  })

  #Check if there are rows with only missing data
  trait <- trait[rowSums(is.na(trait)) != ncol(trait),]  
  if(nrow(trait) != n_sp)
    stop("There are rows in the dataset that contain exclusively missing data. Please omit those and re-try.")
  
  #Provide % of missing data
  for(i in 1:ncol(trait))
    message(paste("-- Column ", colnames(trait)[i], " contains ", round((sum(is.na(trait[,i]) == TRUE) / nrow(trait) )*100,2), "% of missing data.", sep=''))
  
  #Missing Data Imputation with mean
  if(method == "mean") {               
    for(i in 1:nlevels(group))
      for(j in 1:ncol(trait))
        trait[is.na(trait[group==levels(group)[i], j]), j] <- mean(trait[group==levels(group)[i], j], na.rm = TRUE)
    
    #Missing Data Imputation with median
  } else if(method == "median") {      
    for(i in 1:nlevels(group))
      for(j in 1:ncol(trait))
        trait[is.na(trait[group==levels(group)[i], j]), j] <- median(trait[group==levels(group)[i], j], na.rm = TRUE)
    
    #Missing Data Imputation based on the closest species
  } else if(method == "similar") {
    trait_total <- data.frame()
    for (m in 1 : nlevels(group)) {
      trait_m <- trait[group == levels(group)[m], ]
      for (i in 1 : ncol(trait)){
        NA_df <- subset(trait_m, is.na(trait_m[,i]) == TRUE)
        if(nrow(NA_df) > 0) {
          missing_value <- c()
          for (j in 1 : nrow(NA_df)){
            close <- which.min(gower(rbind(NA_df[j,], trait_m[-j,]))[1:nrow(trait_m)])
            missing_value <- append(missing_value, trait_m[ close, i ] )
          } 
          trait_m[c(rownames(NA_df)), i] <- missing_value
        } 
      } 
      trait_total <- rbind(trait_total,trait_m)
    }
    trait <- trait_total
    
    #Missing Data Imputation with linear regression
  } else if(method == "regression") {
    
    #Initial checking
    if(ncol(trait) < 2)
      stop("A minimum of two traits are needed for predicting missing values using regression.")
    
    #With group
    if(nlevels(group) > 1) { #If a grouping factor is provided
      trait2 <- data.frame(trait, group)
      for (i in 1:ncol(trait)){
        NA_df <- subset(trait2, is.na(trait2[,i]) == TRUE)
        #if the column has missing data, make prediction for NA:
        if(nrow(NA_df) > 0) {
          missing_value <- c()
          for (j in 1:nrow(NA_df)) {
            data_j <- NA_df[j,]
            #column without missing data
            column_names <-  names(which(colSums(is.na(NA_df[j,])) == 0)) 
            
            formula.lm <- as.formula(paste(colnames(trait2)[i], " ~ ", paste(column_names, collapse=" + "), sep=""))
            #model fit  
            model <- lm(formula.lm, data = data.frame(na.omit(trait2)))
            #stepwise selection
            if(step)
              model <- MASS::stepAIC(model, direction = "both",  trace = FALSE)
            #make the prediction
            pred_data <- data.frame(NA_df[j,column_names]) ; colnames(pred_data) <- column_names
            missing_value <- append( missing_value , predict(model, newdata = pred_data))
          }
          trait[c(rownames(NA_df)),i] <- missing_value
        }
      }
      
      #Without group
    } else { 
      for (i in 1:ncol(trait)){
        NA_df <- subset(trait, is.na(trait[,i]) == TRUE)
        if(nrow(NA_df) > 0) {
          #if the column has missing data, make prediction for NA:
          missing_value <- c()
          for (j in 1:nrow(NA_df)) {
            #column without missing data
            column_names <- names(which(colSums(is.na(NA_df[j,])) == 0))
            #model formula 
            formula.lm <- as.formula(paste(colnames(trait)[i], " ~ ", paste(column_names, collapse=" + "), sep=""))
            #model fit
            model <- lm(formula.lm, data = data.frame(na.omit(trait)))
            #stepwise selection
            if(step)
              model <- MASS::stepAIC(model, direction = "both",  trace = FALSE)
            #make the prediction
            pred_data <- data.frame(NA_df[j,column_names]) ; colnames(pred_data) <- column_names
            missing_value <- append(missing_value , predict(model, newdata = pred_data))
          }
          trait[c(rownames(NA_df)),i] <- missing_value
        }
      }
    }
    
    #Missing Data Imputation with weighted linear regression
  } else if(method == "w_regression") {  
    
    # Initial checking
    if(ncol(trait) < 2)
      stop("A minimum of two traits are needed for predicting missing values using weighted regression.")
    
    trait <- data.frame(trait)
    
    for (i in 1:ncol(trait)){
      NA_df <- subset(trait, is.na(trait[,i]) == TRUE) 
      #if the column has missing data, make prediction for NA:
      if(nrow(NA_df) > 0) {
        missing_value <- c()
        for(j in rownames(NA_df)){
          j <- as.numeric(j)
          #column without missing data
          column_names <-  names(which(colSums(is.na(trait[j,])) == 0))
          #model formula 
          formula.lm <- as.formula(paste(colnames(trait)[i], " ~ ", paste(column_names, collapse=" + "), sep=""))
          #preparing the weight
          if(class(weight) == "dist")
            weight <- (1 / as.matrix(weight)[j,])
          else if(class(weight) == "phylo" || class(weight) == "hclust")
            weight <- (1 / as.matrix(cophenetic(weight))[j,])
          else
            weight <- (1 - as.matrix(gower(trait))[j,])
          
          weight[is.infinite(weight) == TRUE] <- max(weight[which(!is.infinite(weight))])
          trait2 <- na.omit(data.frame(trait,weight))
          #model fit  
          model <- lm(formula.lm, data = trait2, weights = weight)
          #stepwise selection
          if(step)
            model <- MASS::stepAIC(model, direction = "both", trace = FALSE)
          #make the prediction
          pred_data <- data.frame(trait[j,column_names]) ; colnames(pred_data) <- column_names
          missing_value <- append(missing_value , predict(model, newdata =  pred_data))
        } 
        trait[c(rownames(NA_df)),i] <- missing_value
      }
    } 
    
    #Missing Data Imputation based on a PCA
  } else if(method == "PCA") {
    
    # Note that this chunk of code for PCA is adapted from Podani et al. (2021)
    # Initial checking
    if(ncol(trait) < 2)
      stop("A minimum of two traits are needed for predicting missing values using PCA")
    
    #scaling
    X <- scale(trait, center = TRUE, scale = TRUE)
    #correlation
    C <- cor(X , use="pairwise.complete.obs")
    #Eigenvalue
    Eigenvalues<-eigen(C)$values
    Eigenvalues.pos<-Eigenvalues[Eigenvalues>0]
    Eigenvalues.pos.as.percent<-100*Eigenvalues.pos/sum(Eigenvalues.pos)
    #Eigenvectors
    V <- eigen(C)$vectors
    #Principal components
    X2<-X
    X2[is.na(X2)] <- 0
    PC <- as.matrix(X2) %*% V
    #object.standardized
    PCstand1 <- PC[,Eigenvalues>0]/sqrt(Eigenvalues.pos)[col(PC[,Eigenvalues>0])]
    PCstand2 <- PCstand1 / sqrt(nrow(PC) - 1)
    #loadings
    loadings<-cor(X,PC,use="pairwise.complete.obs")
    #arrows for biplot
    arrows<-cor(X,PC,use="pairwise.complete.obs")*sqrt(nrow(X) - 1)
    #output
    PCA <- list()
    PCA$Correlation.matrix<-C
    PCA$Eigenvalues<-Eigenvalues
    PCA$Positive.Eigenvalues<-Eigenvalues.pos
    PCA$Positive.Eigenvalues.as.percent<-100*Eigenvalues.pos/sum(Eigenvalues.pos)
    PCA$Square.root.of.eigenvalues <- sqrt(Eigenvalues.pos)
    PCA$Eigenvectors<-V
    PCA$Component.scores<-PC
    PCA$Variable.scores<-loadings
    PCA$Biplot.objects<-PCstand2
    PCA$Biplot.variables<-arrows
    #Final output
    trait <- PCA
  } else {
    stop(sprintf("Method %s not recognized.", method))
  }
  return(trait)
}

#' Standardize variables.
#' @description Standardize (or normalize) variables in different ways.
#' @param trait A species x traits matrix or data.frame.
#' @param method One of "standard" (standardize to mean = 0 and sd = 1, i.e., use z-score), "range" (rescale with range 0-1), or "rank" (rescale with range 0-1 after ranking).
#' @param convert A vector of column numbers to be standardized. If NULL all will be standardized.
#' @details Standardizing values allows to directly compare variables of interest with inherently different ranges, avoiding artificial distortions of distances between observations.
#' @return A matrix with variables standardized.
#' @examples trait = data.frame(body = c(20,40,60,30,10), beak = c(NA,4,6,3,1))
#' standard(trait)
#' standard(trait, method = "range")
#' standard(trait, method = "rank")
#' @export
standard <- function(trait, method = "standard", convert = NULL){

  if(is.vector(trait))
    trait = matrix(trait, ncol = 1)
  if(is.null(convert))
    convert = 1:ncol(trait)
  
  for(i in convert){
    
    #if standardization with mean 0 and sd 1 
    if(method == "standard"){
      trait[,i] = (trait[,i] - mean(trait[,i], na.rm = TRUE)) / sd(trait[,i], na.rm = TRUE)
    
    #if standardization with range 0-1  
    } else if(method == "range"){
      trait[,i] = (trait[,i] - min(trait[,i], na.rm = TRUE)) / (max(trait[,i], na.rm = TRUE) - min(trait[,i], na.rm = TRUE))
    
    #if standardization with ranking of values
    } else if(method == "rank"){
      trait[,i] = rank(trait[,i], ties.method = "average", na.last = "keep")
      trait[,i] = standard(trait[,i], method = "range")
      
    } else {
      stop("Method not recognized, must be one of 'standard', 'range' or 'rank'.")
    }
  }
  return(trait)
}

#' Gower distance.
#' @description Calculates Gower distances between observations.
#' @param trait A species x traits matrix or data.frame.
#' @param convert A vector of column numbers, usually categorical variables, to be converted to dummy variables.
#' @param weight A vector of column numbers with weights for each variable. Its length must be equal to the number of columns in trait.
#' @details The Gower distance allows continuous, ordinal, categorical or binary variables, with possible weighting (Pavoine et al. 2009).
#' NAs are allowed as long as each pair of species has at least one trait value in common.
#' If convert is given the algorithm will convert these column numbers to dummy variables. Otherwise it will convert all columns with factors or characters as values.
#' @return A dist object with pairwise distances between species.
#' @references Pavoine et al. (2009) On the challenge of treating various types of variables: application for improving the measurement of functional diversity. Oikos, 118: 391-402.
#' @examples trait = data.frame(body = c(NA,2,3,4,4), beak = c(1,1,1,1,2))
#' gower(trait)
#' gower(trait, weight = c(1, 0))
#' @export
gower <- function(trait, convert = NULL, weight = NULL){
  
  #prepare data
  if(is.vector(trait))
    trait = matrix(trait, ncol = 1)
  if(is.null(weight))
    weight = rep(1, ncol(trait))
  spNames = rownames(trait)
  
  #dummify, standardize, and get weights
  trait = dummy(trait, convert, weight = weight)
  weight = trait$weight
  trait = standard(trait$trait, method = "range")

  #calculate gower
  nSp = nrow(trait)
  res = matrix(0, nrow = nSp, ncol = nSp)
  for(i in 1:(nSp-1)){
    for(j in (i+1):nSp){
      sumWeights = 0
      for(t in 1:ncol(trait)){
        value1 = trait[i,t]
        value2 = trait[j,t]
        if(!any(is.na(c(value1, value2)))){
          res[j,i] = res[j,i] + ((trait[i,t] - trait[j,t])^2 * weight[t])
          sumWeights = sumWeights + weight[t]
        }
      }
      res[j,i] = (res[j,i] / sumWeights)^0.5
    }
  }
  rownames(res) = colnames(res) = spNames
  return(as.dist(res))
}

#' Model R2.
#' @description Calculates R2 from the summed squared differences between observed and estimated values.
#' @param obs Either a model or a vector with observed values.
#' @param est A vector with estimated values. Only used if obs is not a model.
#' @param param Number of parameters in the model to calculate the adjusted R2 if > 0. If obs is a model param will be ignored and the number of parameters will be calculated from the model.
#' @details Useful for models or functions that do not provide r2 values.
#' @return The r2 value.
#' @examples obs = c(1,4,5,6)
#' est = c(0,1,4,7)
#' 
#' #example using values
#' r2(obs, est)
#' r2(obs, est, param = 1)
#' 
#' #example using model
#' mod = lm(obs ~ est)
#' r2(mod)
#' summary(mod)$r.squared
#' r2(mod, param = 1)
#' summary(mod)$adj.r.squared
#' 
#' @export
r2 <- function(obs, est = NULL, param = 0){
  if(!is.vector(obs)){
    if(param > 0)
      param = length(obs$coefficients) - 1
    est = predict(obs)
    obs = obs$model[,1]
  }
  SSn <- rss(obs, est)
  SSd <- sum((obs-mean(obs))^2)
  res <- 1-(SSn/SSd)
  if(param > 0){
    n = length(obs)
    res <- 1-(1-res)*((n-1)/(n-param-1))
  }
  return(res)
}

#' Akaike Information Criterion.
#' @description Calculates the Akaike Information Criterion (AIC) of any model based on observed and estimated values.
#' @param obs Either a model or a vector with observed values.
#' @param est A vector with estimated values. Only used if obs is not a model.
#' @param param Number of parameters in the model. If obs is a model param will be ignored and the number of parameters will be calculated from the model.
#' @param correct Boolean indicating whether the corrected version of AIC (AICc) should be calculated, mostly for models with few observations.
#' @details Useful for models or functions that do not provide logLik values.
#' @return The AIC or AICc value.
#' @examples obs = c(1,4,5,6)
#' est = c(0,1,4,7)
#' 
#' #example using values
#' aic(obs, est)
#' aic(obs, est, param = 1)
#' aic(obs, est, param = 1, correct = TRUE)
#' 
#' #example using model
#' mod = lm(obs ~ est)
#' aic(mod)
#' extractAIC(mod)[2]
#' aic(mod, correct = TRUE)
#' 
#' @export
aic <- function(obs, est = NULL, param = 0, correct = FALSE){
  if(!is.vector(obs)){
    param = length(obs$coefficients) - 1
    est = predict(obs)
    obs = obs$model[,1]
  }
  k = param + 1
  n = length(obs)
  res = n * log(rss(obs, est)/n) + 2 * k
  
  if(correct)
    res = res + (2*k*(k+1))/(n-k-1)
  return(res)
}

#' Standard Effect Size.
#' @description Calculates the standard effect size from observed and estimated values.
#' @param obs A single observed value.
#' @param est A vector with estimated values.
#' @param param Value is calculated with parametric or non-parametric method. The later is preferable when distribution of estimated values is not normally distributed.
#' @param p Boolean indicating whether the p-value should be calculated based on the rank of 'obs' in 'est'. As the test is two-sided, usual significance is indicated by p < 0.025 or p > 0.975.
#' @return The ses value or a vector with ses and p-value.
#' @examples est = rnorm(1000, 500, 100)
#' 
#' ses(100, est)
#' ses(100, est, param = FALSE)
#' ses(500, est)
#' ses(500, est, param = FALSE)
#' ses(900, est, p = TRUE)
#' ses(900, est, param = FALSE, p = TRUE)
#' @export
ses <- function(obs, est, param = TRUE, p = TRUE){
  if(param){
    res = (obs - mean(est)) / (sd(est, na.rm = TRUE))
    pval = pnorm(res)   #converts ses to p
  } else {
    pval = (sum(est < obs) + sum(est == obs)/2) / (length(est) + 1)
    if(pval == 0)
      pval = pval + (1 / (length(est)))
    res = qnorm(pval)   #converts p to ses
  }

  if(p){
    res = c(res, pval)
    names(res) = c("ses", "p-value")
  }
  
  return(res)
}

#' Build functional tree.
#' @description Builds a functional tree from trait data.
#' @param trait A species x traits matrix or data.frame.
#' @param convert A vector of column numbers, usually categorical variables, to be converted to dummy variables.
#' @param weight A vector of column numbers with weights for each variable. Its length must be equal to the number of columns in trait.
#' @details The tree will be built using UPGMA over gower distances (Pavoine et al. 2009) after traits are dummyfied (if needed) and standardized (always).
#' Gower distance allows continuous, ordinal, categorical or binary variables, with possible weighting.
#' NAs are allowed as long as each pair of species has at least one trait value in common.
#' If convert is given the algorithm will convert these column numbers to dummy variables. Otherwise it will convert all columns with factors or characters as values.
#' @return An hclust object representing a functional tree.
#' @references Pavoine et al. (2009) On the challenge of treating various types of variables: application for improving the measurement of functional diversity. Oikos, 118: 391-402.
#' @examples trait = data.frame(body = c(NA,2,3,4,4), beak = c(1,1,1,1,2))
#' plot(tree.build(trait))
#' plot(tree.build(trait, weight = c(0, 1)))
#' @export
tree.build <- function(trait, convert = NULL, weight = NULL){
  trait = gower(trait, convert, weight)
  res = hclust(trait, method = "average")
  return(res)
}

#' Build hyperspace.
#' @description Builds hyperspace by transforming trait data to use with either hull.build or kernel.build.
#' @param trait A species x traits matrix or data.frame.
#' @param weight A vector of column numbers with weights for each variable. Its length must be equal to the number of columns in trait. Only used if axes > 0.
#' @param axes If 0, no transformation of data is done.
#' If 0 < axes <= 1 a PCoA is done with Gower distances and as many axes as needed to achieve this proportion of variance explained are selected.
#' If axes > 1 these many axes are selected.
#' @param convert A vector of column numbers, usually categorical variables, to be converted to dummy variables. Only used if axes > 0.
#' @details The hyperspace can be constructed with the given data or data can be transformed using PCoA over Gower distances (Pavoine et al. 2009) after traits are dummyfied (if needed) and standardized (always).
#' Gower distance allows continuous, ordinal, categorical or binary variables, with possible weighting.
#' NAs are allowed as long as each pair of species has at least one trait value in common.
#' If convert is given the algorithm will convert these column numbers to dummy variables. Otherwise it will convert all columns with factors or characters as values.
#' Note that each community should have at least 3 species and more species than traits or axes (if axes > 0) to build convex hull hypervolumes.
#' Transformation of traits is recommended if (Carvalho & Cardoso, 2020):
#' 1) Some traits are not continuous;
#' 2) Some traits are correlated; or
#' 3) There are less species than traits + 1, in which case the number of axes should be smaller.
#' @return A matrix with trait data.
#' @references Carvalho, J.C. & Cardoso, P. (2020) Decomposing the causes for niche differentiation between species using hypervolumes. Frontiers in Ecology and Evolution, 8: 243.
#' @references Pavoine et al. (2009) On the challenge of treating various types of variables: application for improving the measurement of functional diversity. Oikos, 118: 391-402.
#' @examples
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' 
#' hs = hyper.build(trait, weight = c(1,2), axes = 2)
#' plot(hs)
#' @export
hyper.build <- function(trait, weight = NULL, axes = 1, convert = NULL){
  
  #do gower and pcoa if axes is larger than 0
  if(axes > 0){
    trait = gower(trait, convert, weight)
    trait = ape::pcoa(trait)

    if(axes <= 1){
      axes = which(trait$values$Cumul_eig >= axes)[1]
    } else {
      axes = min(axes, ncol(trait$vectors))
    }
    
    trait = trait$vectors[,(1:axes)]
  }
  
  if(is.null(colnames(trait)))
    colnames(trait) = paste("Trait", 1:ncol(trait), sep = "")
  if(is.null(rownames(trait)))
    rownames(trait) = paste("Sp", 1:nrow(trait), sep = "")
  
  return(trait)
}

#' Build convex hull hypervolumes.
#' @description Builds convex hull hypervolumes for each community from incidence and trait data.
#' @param comm A sites x species matrix, data.frame or vector, with incidence data about the species in the community.
#' @param trait A species x traits matrix or data.frame.
#' @param weight A vector of column numbers with weights for each variable. Its length must be equal to the number of columns in trait. Only used if axes > 0.
#' @param axes If 0, no transformation of data is done.
#' If 0 < axes <= 1 a PCoA is done with Gower distances and as many axes as needed to achieve this proportion of variance explained are selected.
#' If axes > 1 these many axes are selected.
#' @param convert A vector of column numbers, usually categorical variables, to be converted to dummy variables. Only used if axes > 0.
#' @details The hypervolumes can be constructed with the given data or data can be transformed using PCoA over Gower distances (Pavoine et al. 2009) after traits are dummyfied (if needed) and standardized (always).
#' Beware that if transformations are required, all communities to be compared should be built simultaneously to guarantee comparability. In such case, one might want to first run hyper.build and use the resulting data in different runs of hull.build.
#' See function hyper.build for more details.
#' @return A 'convhulln' object or a list, representing the hypervolumes of each community.
#' @examples comm = rbind(c(1,3,0,5,3), c(3,2,5,0,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site 1", "Site 2")
#' 
#' trait = data.frame(body = c(1,2,3,4,4), beak = c(1,5,4,1,2))
#' rownames(trait) = colnames(comm)
#' 
#' hv = hull.build(comm[1,], trait)
#' plot(hv)
#' hvlist = hull.build(comm, trait)
#' plot(hvlist[[2]])
#' hvlist = hull.build(comm, trait, axes = 2, weight = c(1,2))
#' plot(hvlist[[1]])
#' @export
hull.build <- function(comm, trait, weight = NULL, axes = 0, convert = NULL){
  
  trait = hyper.build(trait, weight, axes, convert)
  if(is.vector(comm))
    comm = matrix(comm, nrow = 1)

  #check if there are communities with less species than traits + 1 and remove them
  fewSpp = as.vector(which(rowSums(ifelse(comm > 0, 1, 0)) <= ncol(trait)))
  if(length(fewSpp) > 0){
    if(nrow(comm) == length(fewSpp))
      stop(paste("There are no communities with enough species for convex hull delineation.\n"))
    comm2 = comm[-(fewSpp), , drop = FALSE]
    for(i in fewSpp)
      warning(paste("Site", fewSpp, "does not contain enough species for convex hull delineation and has been removed prior to convex hull estimation.\n"))
    comm <- comm2
  }
  
  if(is.null(colnames(comm)))
    colnames(comm) = rownames(trait)
  if(is.null(rownames(comm)))
    rownames(comm) = paste("Comm", 1:nrow(comm), sep = "")

  #build the convex hulls for each community
  hull_list <- list()
  for (s in 1:nrow(comm)) {
    subComm <- trait[comm[s,] > 0, ]
    hull_list[[s]] <- geometry::convhulln(subComm, options = "FA")
    rownames(hull_list[[s]]$p) = rownames(subComm)
  }
  names(hull_list) = rownames(comm)

  #return
  if(length(hull_list) > 1){
    return(hull_list)
  } else {
    return(hull_list[[1]])
  }
}

#' Build kernel hypervolumes.
#' @description Builds kernel density hypervolumes from trait data.
#' @param comm A sites x species matrix, data.frame or vector, with incidence or abundance data about the species in the community.
#' @param trait A species x traits matrix or data.frame.
#' @param method Method for constructing the 'Hypervolume' object. One of "gaussian" (Gaussian kernel density estimation, default), "box" (box kernel density estimation), or "svm" (one-class support vector machine). See respective functions of the hypervolume R package for details.
#' @param abund A boolean (T/F) indicating whether abundance data should be used as weights in hypervolume construction. Only works if method = "gaussian".
#' @param weight A vector of column numbers with weights for each variable. Its length must be equal to the number of columns in trait. Only used if axes > 0.
#' @param axes If 0, no transformation of data is done.
#' @param ... further arguments to be passed to hypervolume::hypervolume
#' If 0 < axes <= 1 a PCoA is done with Gower distances and as many axes as needed to achieve this proportion of variance explained are selected.
#' If axes > 1 these many axes are selected.
#' @param convert A vector of column numbers, usually categorical variables, to be converted to dummy variables. Only used if axes > 0.
#' @param cores Number of cores to be used in parallel processing. If = 0 all available cores are used. Beware that multicore for Windows is not optimized yet and it often takes longer than single core.
#' @details The hypervolumes can be constructed with the given data or data can be transformed using PCoA over Gower distances (Pavoine et al. 2009) after traits are dummyfied (if needed) and standardized (always).
#' Beware that if transformations are required, all communities to be compared should be built simultaneously to guarantee comparability. In such case, one might want to first run hyper.build and use the resulting data in different runs of kernel.build.
#' See function hyper.build for more details.
#' @return A 'Hypervolume' or 'HypervolumeList', representing the hypervolumes of each community.
#' @examples \dontrun{
#' comm = rbind(c(1,1,0,5,1), c(3,2,5,0,0))
#' colnames(comm) = c("SpA", "SpB", "SpC", "SpD", "SpE")
#' rownames(comm) = c("Site1", "Site2")
#' 
#' trait = data.frame(body = c(1,2,3,1,2), beak = c(1,2,3,2,1))
#' rownames(trait) = colnames(comm)
#' 
#' hv = kernel.build(comm[1,], trait)
#' plot(hv)
#' hvlist = kernel.build(comm, trait, abund = FALSE, cores = 0)
#' plot(hvlist)
#' hvlist = kernel.build(comm, trait, method = "box", weight = c(1,2), axes = 2)
#' plot(hvlist)
#' }
#' @export
kernel.build <- function(comm, trait, method = "gaussian", abund = TRUE, weight = NULL, axes = 0, convert = NULL, cores = 1, ... ){

  trait = hyper.build(trait, weight, axes, convert)
  if(is.vector(comm))
    comm = matrix(comm, nrow = 1)

  #check if there are communities with no species
  fewSpp = as.vector(which(rowSums(comm) == 0))
  if(length(fewSpp) > 0){
    if(nrow(comm) == length(fewSpp))
      stop(paste("There are no communities with species.\n"))
    comm2 = comm[-(fewSpp), , drop = FALSE]
    for(i in fewSpp)
      warning(paste("Site", fewSpp, "does not contain any species and has been removed prior to hypervolume estimation.\n"))
    comm <- comm2
  }
  
  #general function for lapply (serial), mcapply (Mac/Linux) or parLapply (Win)
  parbuild <- function(subComm, trait, method, abund, ... ){
    subTrait <- trait[subComm > 0, ] ##Select traits
    subComm <- subComm[subComm > 0]
     
    if (method == "box"){
      newHv <- hypervolume_box(subTrait, verbose = FALSE, ... )
    } else if (method == "svm"){
      newHv <- hypervolume_svm(subTrait, verbose = FALSE, ... )
    } else if (method == "gaussian"){
      if(abund){
        weight = subComm / sum(subComm)
        newHv <- hypervolume_gaussian(subTrait, weight = weight, verbose = FALSE, ... )
      } else {
        newHv <- hypervolume_gaussian(subTrait, verbose = FALSE, ... )
      }
    } else {
      stop(sprintf("Method %s not recognized.", method))
    }
    return(newHv)
  }
  
  #get number of cores, if cores = 0 uses all available
  if(cores == 0){
    cores = detectCores()
  } else {
    cores = min(cores, detectCores())
  }
  
  #make list to go for parallel if required
  commList = as.list(as.data.frame(t(comm)))
  if(cores == 1){
    hv = lapply(commList, parbuild, trait = trait, method = method, abund = abund, ... )
  } else if (Sys.info()[['sysname']] == 'Windows'){
    cl = makeCluster(cores)
    func = paste("hypervolume", method, sep = "_")
    clusterExport(cl, c("trait", "method", "abund", func))
    hv = parLapply(cl, commList, parbuild, trait = trait, method = method, abund = abund, ...)
    stopCluster(cl)
  } else {
    hv = mclapply(commList, parbuild, trait = trait, method = method, abund = abund, mc.cores = cores, ... )
  }
  
  #name hypervolumes and convert list to HypervolumeList
  if(length(hv) > 1){
    for(i in 1:length(hv))
      hv[[i]]@Name = rownames(comm)[i]
  }
  hv = hypervolume_join(hv)
  
  #return
  if(length(hv@HVList) > 1){
    return(hv)
  } else {
    return(hv[[1]])
  }
}

#' Simulation of species abundance distributions (SAD).
#' @description Creates artificial communities following given SADs.
#' @param n total number of individuals.
#' @param s number of species.
#' @param sad The SAD distribution type (lognormal, uniform, broken stick or geometric). Default is lognormal.
#' @param sd The standard deviation of lognormal distributions. Default is 1.
#' @details Species Abundance Distributions may take a number of forms. A lognormal SAD probably is the most supported by empirical data, but we include other common types useful for testing multiple algorithms including several of the functions in BAT.
#' @return A matrix of species x abundance per species.
#' @examples comm1 <- sim.sad(10000, 100)
#' comm2 <- sim.sad(10000, 100, sd = 2)
#' comm3 <- sim.sad(10000, 100, sad = "uniform")
#' par(mfrow=c(1,3))
#' hist(log(comm1$Freq))
#' hist(log(comm2$Freq))
#' hist(log(comm3$Freq))
#' @export
sim.sad <- function(n, s, sad = "lognormal", sd = 1) {
	if (s > n)
		stop("Number of species can't be larger than number of individuals")
	sppnames = paste("Sp", 1:s, sep="") ##species names
	sad <- match.arg(sad, c("lognormal", "uniform", "broken", "geometric"))

	##lognormal distribution
	switch(sad, lognormal = {
		comm = sample(sppnames, size = n, replace = T, prob = c(rlnorm(s, sdlog = sd)))

		##uniform distribution
	}, uniform = {
		comm = sample(sppnames, size = n, replace = T)

		##broken stick distribution
	}, broken = {
		broken.stick <- function(p){
			result = NULL
			for(j in 1:p) {
				E = 0
				for(x in j:p)
					E = E+(1/x)
				result[j] = E/p
			}
			return(result)
		}
		broken.prob = broken.stick(s)
		comm = sample(sppnames, size = n, replace = TRUE, prob = c(broken.prob))
	}, geometric = {
		geo.ser <- function(s, k = 0.3){
			result = NULL
			for (x in 1:s) {
				result[x] = k*(1-k)^(x-1)/(1-(1-k)^s)
			}
			return(result)
		}
		geo.prob = geo.ser(s)
		comm = sample(sppnames, size = n, replace = TRUE, prob = c(geo.prob))
	})
	return(as.data.frame(table(comm)))
}

#' Simulation of species spatial distributions.
#' @description Creates artificial communities with given SAD and spatial clustering.
#' @param n total number of individuals.
#' @param s number of species.
#' @param sad The SAD distribution type (lognormal, uniform, broken stick or geometric). Default is lognormal.
#' @param sd The standard deviation of lognormal distributions. Default is 1.
#' @param distribution The spatial distribution of individual species populations (aggregated, random, uniform or gradient). Default is aggregated.
#' @param clust The clustering parameter if distribution is either aggregated or gradient (higher values create more clustered populations). Default is 1.
#' @details The spatial distribution of individuals of given species may take a number of forms.
#' Competitive exclusion may cause overdispersion, specific habitat needs or cooperation may cause aggregation and environmental gradients may cause abundance gradients.
#' @return A matrix of individuals x (species, x coords and y coords).
#' @examples par(mfrow = c(3 ,3))
#' comm = sim.spatial(100, 9, distribution = "uniform")
#' for(i in 1:9){
#' 	sp <- comm[comm[1] == paste("Sp", i, sep = ""), ]
#' 	plot(sp$x, sp$y, main = paste("Sp", i), xlim = c(0,1), ylim = c(0,1))
#' }
#' comm = sim.spatial(1000, 9, sad = "lognormal", sd = 0.5, distribution = "aggregated", clust = 2)
#' for(i in 1:9){
#' 	sp <- comm[comm[1] == paste("Sp", i, sep=""), ]
#' 	plot(sp$x, sp$y, main = paste("Sp", i), xlim = c(0,1), ylim = c(0,1))
#' }
#' @export
sim.spatial <- function(n, s, sad = "lognormal", sd = 1, distribution = "aggregated", clust = 1){
	repeat{
		simsad <- sim.sad(n, s, sad, sd)
		coords <- matrix(ncol = 2)
		for (j in 1:nrow(simsad)){	#species by species
			spCoords <- matrix(c(runif(1),runif(1)), ncol=2)
			if(simsad[j,2] > 1){
				for(i in 2:simsad[j,2]){
					repeat{
						newcoords <- c(runif(1),runif(1))

						##aggregated distribution
						if(distribution == "aggregated"){
							mindist = 1
							for(r in 1:nrow(spCoords)){
								mindist <- min(mindist, euclid(newcoords, spCoords[r,]))
							}
							thres = abs(rnorm(1, sd = 1/clust))/10
							if(mindist < thres){
								spCoords <- rbind(spCoords, newcoords)
								break
							}
						
						##random distribution
						} else if (distribution == "random"){
							spCoords <- rbind(spCoords, newcoords)
							break
						
						##uniform distribution
						} else if (distribution == "uniform"){
							mindist = 1
							for(r in 1:nrow(spCoords)){
								mindist <- min(mindist, euclid(newcoords, spCoords[r,]))
							}
							thres = runif(1)
							if(mindist > thres){
								spCoords <- rbind(spCoords, newcoords)
								break
							}
						
						##gradient distribution
						} else if (distribution == "gradient") {
							thres = runif(1)^(1/clust)
							if(newcoords[2] > thres){
								spCoords <- rbind(spCoords, newcoords)
								break
							}
							

							
						} else {
						  return(message("distribution not recognized"))
						} 
					}
				}
			}
			coords <- rbind(coords, spCoords)
		}
		spp <- rep(as.character(simsad[,1]), simsad[,2])
		comm <- data.frame(Spp = spp, x = coords[-1,1], y = coords[-1,2])
		if(nrow(comm) == n && length(which(is.na(comm$x))) == 0){
			break
		}
	}
	return(comm)
}

#' Plots of simulated species spatial distributions.
#' @description Plots individuals from artificial communities with given SAD and spatial clustering.
#' @param comm artificial community data from function sim.spatial.
#' @param sad boolean indicating if the SAD plot should also be shown. Default is FALSE.
#' @param s number of species to plot simultaneously. Default is the number of species in comm.
#' @details Function useful for visualizing the results of sim.spatial.
#' @examples comm <- sim.spatial(1000, 24)
#' sim.plot(comm)
#' sim.plot(comm, sad = TRUE)
#' sim.plot(comm, s = 9)
#' @export
sim.plot <- function(comm, sad = FALSE, s = 0){
	spp <- length(unique(comm$Spp))
	if(s < 1){
		if(!sad)
			side = ceiling(sqrt(spp))
		else
			side = ceiling(sqrt(spp+1))
	}	else{
		side = ceiling(sqrt(s))
	}
	par(mfrow = c(side,side), mar = c(1,1,2,1))
	if(sad)
		hist(log(table(comm[,1])), main = "All species", xlab = "Abundance (log)", xaxt = "n")
	for(i in 1:spp){
		sp <- comm[comm[1] == paste("Sp", i, sep=""), ]
		plot(sp$x, sp$y, main = paste("Sp", i), xlim = c(0,1), ylim = c(0,1), xlab="", ylab="", xaxt="n", yaxt="n")
	}
}

#' Simulation of sampling from artificial communities.
#' @description Simulates a sampling process from artificial communities.
#' @param comm simulated community data from function sim.spatial.
#' @param cells number of cells to divide the simulated space into. Default is 100.
#' @param samples number of samples (cells) to randomly extract. Default is the number of cells (the entire community).
#' @details The space will be divided in both dimensions by sqrt(cells).
#' @details Function useful for simulating sampling processes from the results of sim.spatial.
#' @details May be used as direct input to other functions (e.g. alpha, alpha.accum, beta, beta.accum) to test the behavior of multiple descriptors and estimators.
#' @return A matrix of samples x species (values are abundance per species per sample).
#' @examples comm <- sim.spatial(1000, 10)
#' sim.sample(comm)
#' sim.sample(comm, cells = 10, samples = 5)
#' @export
sim.sample <- function(comm, cells = 100, samples = 0){
	side <- round(sqrt(cells),0)
	cells = side^2

	comm$ind <- 0
	xv <- cut(comm$x, seq(0, 1, 1/side))
	yv <- cut(comm$y, seq(0, 1, 1/side))
	grid1 <- data.frame(table(xv, yv))
	grid1 <- grid1[,-3]

	s <- 1:cells
	for (i in 1:cells){
		id <- NULL
		id <- which (xv == grid1$xv[s[i]] & yv == grid1$yv[s[i]])
		comm$ind[id] <- paste("Sample", i, sep="")
	}
	comm <- table(comm$ind, comm$Spp) ##entire community
	comm <- comm[rownames(comm) != "0", ]
	if (samples < 1 || samples > nrow(comm))
		samples = nrow(comm)

	##number of samples to take
	samp <- comm[sample(nrow(comm), samples, replace = FALSE),] ## sampled community

	return(samp)
}

#' Simulation of phylogenetic or functional tree.
#' @description Simulates a random tree.
#' @param s number of species.
#' @param m a structural parameter defining the average difference between species. Default is 100. Lower numbers create trees dominated by increasingly similar species, higher numbers by increasingly dissimilar species.
#' @details A very simple tree based on random genes/traits.
#' @return An hclust object.
#' @examples tree <- sim.tree(10)
#' plot(as.dendrogram(tree))
#' tree <- sim.tree(100,10)
#' plot(as.dendrogram(tree))
#' tree <- sim.tree(100,1000)
#' plot(as.dendrogram(tree))
#' @export
sim.tree <- function(s, m = 100){
	sim.matrix <- matrix(sample(0:m, ceiling(s*m/50), replace = TRUE), nrow = s, ncol = m)
	tree <- hclust(dist(sim.matrix), method = "average")
	tree$height <- tree$height / max(tree$height)
	return(tree)
}

#' Sample data of spiders in Arrabida (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 sampling units. Details are described in:
#' Cardoso, P., Gaspar, C., Pereira, L.C., Silva, I., Henriques, S.S., Silva, R.R. & Sousa, P. (2008) Assessing spider species richness and composition in Mediterranean cork oak forests. Acta Oecologica, 33: 114-127.
#'
#' @docType data
#' @keywords datasets
#' @name arrabida
#' @usage data(arrabida)
#' @format A data frame with 320 sampling units (rows) and 338 species (variables).
NULL

#' Sample data of spiders in Geres (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 sampling units. Details are described in:
#' Cardoso, P., Scharff, N., Gaspar, C., Henriques, S.S., Carvalho, R., Castro, P.H., Schmidt, J.B., Silva, I., Szuts, T., Castro, A. & Crespo, L.C. (2008) Rapid biodiversity assessment of spiders (Araneae) using semi-quantitative sampling: a case study in a Mediterranean forest. Insect Conservation and Diversity, 1: 71-84.
#'
#' @docType data
#' @keywords datasets
#' @name geres
#' @usage data(geres)
#' @format A data frame with 320 sampling untis (rows) and 338 species (variables).
NULL

#' Sample data of spiders in Guadiana (Portugal)
#'
#' A dataset containing the abundance of 338 spider species in each of 320 sampling units. Details are described in:
#' Cardoso, P., Henriques, S.S., Gaspar, C., Crespo, L.C., Carvalho, R., Schmidt, J.B., Sousa, P. & Szuts, T. (2009) Species richness and composition assessment of spiders in a Mediterranean scrubland. Journal of Insect Conservation, 13: 45-55.
#'
#' @docType data
#' @keywords datasets
#' @name guadiana
#' @usage data(guadiana)
#' @format A data frame with 192 sampling units (rows) and 338 species (variables).
NULL

#' Functional tree for 338 species of spiders
#'
#' A dataset representing the functional tree for 338 species of spiders captured in Portugal.
#' For each species were recorded: average size, type of web, type of hunting, stenophagy, vertical stratification in vegetation and circadial activity. Details are described in:
#' Cardoso, P., Pekar, S., Jocque, R. & Coddington, J.A. (2011) Global patterns of guild composition and functional diversity of spiders. PLoS One, 6: e21710.
#'
#' @docType data
#' @keywords datasets
#' @name functree
#' @usage data(functree)
#' @format An hclust object with 338 species.
NULL

#' Taxonomic tree for 338 species of spiders (surrogate for phylogeny)
#'
#' A dataset representing an approximation to the phylogenetic tree for 338 species of spiders captured in Portugal.
#' The tree is based on the linnean hierarchy, with different suborders separated by 1 unit, families by 0.75, genera by 0.5 and species by 0.25.
#'
#' @docType data
#' @keywords datasets
#' @name phylotree
#' @usage data(phylotree)
#' @format An hclust object with 338 species.
NULL