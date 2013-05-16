### Scuba Helper Scripts
###
## load with: source('~/Documents/r_scripts/scuba_scripts.r', chdir = TRUE)

normalit <- function (m) (m - min(m))/(max(m)-min(m))


surfing_score <- function(m, k) {
	require(cluster)

	# distanz matrix berechnen
	d <- as.matrix(daisy(m, metric = c("euclidean"), stand = FALSE, type = list()))

	# distanzen normieren
	d <- t(t(d) / max(d))
	
	d_sorted <- t(apply(d,1,sort))
	# d_sorted <- d_sorted[,-1]
	# d_sorted <- d_sorted[,-(k+1:ncol(d_sorted))]
  # max_Distances <- d_sorted[, ncol(d_sorted)]

	# getting only the max Distances (per definition: kNN = max(d))
	max_Distances <- d_sorted[, k]

	# calculating the mean of kNN distances
	mean_of_maxDistances <- mean(max_Distances)
	
	meanDistances <- sapply(max_Distances, function(x) abs(mean_of_maxDistances - x))
	
	diff_mean <- 0.5 * sum(meanDistances)
	below_s <- length(which(max_Distances < mean_of_maxDistances))
	if (below_s > 0) {
		quality_s <- diff_mean / (below_s * mean_of_maxDistances)
	} else {
		quality_s <- 0
	}

	return(quality_s)
}

calc_entropy <- function(m, num_xbins) {
	require(hexbin)
	require(entropy)
	
	if (ncol(m) > 2) {
		stop("Matrix dimensions > 2 -- is Subspace projected?")
	}
	
	bin<-hexbin(m[,1],m[,2], xbins = num_xbins, shape = 1, xbnds = range(m[,1]), ybnds = range(m[,2]),xlab = NULL, ylab = NULL, IDs = TRUE)
	vect   <- c(bin@count, rep(0, length(m[,1]) - length(bin@count)))
	counts <- as.vector(table(vect))
	entr <- entropy.empirical(counts, unit = c("log2"))
	return (entr)
}

data_topology <- function(m, n, sample_size, number_of_neighbours) {
	require(cluster)
	int_sample_size      = as.integer(nrow(m) * sample_size)
	int_num_neighbours   = as.integer(nrow(m) * number_of_neighbours)

	# calculate distance matrix of given subspace (matrix m)
	dMat1 <- as.matrix(daisy(m, metric = c("euclidean"), stand = FALSE, type = list()))
	dMat2 <- as.matrix(daisy(n, metric = c("euclidean"), stand = FALSE, type = list()))

	# select random points
	random_points <- sample(1:nrow(m), int_sample_size, replace=F)

	avg <- vector()

	for(i in random_points) {
		first_neighbours  <- sort(abs(dMat1[-i,i] - dMat1[i,i]), index.return=TRUE)$ix[1:int_num_neighbours] + 1
		second_neighbours <- sort(abs(dMat2[-i,i] - dMat2[i,i]), index.return=TRUE)$ix[1:int_num_neighbours] + 1

		int <- intersect(first_neighbours, second_neighbours)
		res <- length(int) / length(first_neighbours)

		avg <- append(avg, res)
	}

	return (1-mean(avg))
}

vis_data_topology <- function(m, n, sample_size, number_of_neighbours) {
	require(cluster)
	int_sample_size      = as.integer(nrow(m) * sample_size)
	int_num_neighbours   = as.integer(nrow(m) * number_of_neighbours)

	# calculate mds projections
	mds1 <- cmdscale(daisy(m, metric = c("euclidean"), stand = FALSE, type = list()), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
	mds2 <- cmdscale(daisy(n, metric = c("euclidean"), stand = FALSE, type = list()), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

	# calculate distance matrix of the projections
	dMat1 <- as.matrix(daisy(mds1, metric = c("euclidean"), stand = FALSE, type = list()))
	dMat2 <- as.matrix(daisy(mds2, metric = c("euclidean"), stand = FALSE, type = list()))

	# select random points
	random_points <- sample(1:nrow(m), int_sample_size, replace=F)

	avg <- vector()

	for(i in random_points) {
		first_neighbours  <- sort(abs(dMat1[-i,i] - dMat1[i,i]), index.return=TRUE)$ix[1:int_num_neighbours] + 1
		second_neighbours <- sort(abs(dMat2[-i,i] - dMat2[i,i]), index.return=TRUE)$ix[1:int_num_neighbours] + 1

		int <- intersect(first_neighbours, second_neighbours)
		res <- length(int) / length(first_neighbours)

		avg <- append(avg, res)
	}

	return (1-mean(avg))
}

factorial <- function (x) gamma(1 + x)

subspace_from_indices <- function(sub) {
	m <- NULL
	for (j in seq(1, length(sub))) {
		m <- cbind(m, pl.data[, sub[j]])
	}
	
	return (m)
}

create_sampled_subspaces <- function(num_of_subspaces_per_cardinality, cardinalities, total_attributes) {	
	z <- vector("list", length(cardinalities) * num_of_subspaces_per_cardinality)
	
	for (i in seq(1, length(cardinalities))) {
		m <- subspace_sampling(cardinalities[i], total_attributes, num_of_subspaces_per_cardinality)

		for (j in seq(1, ncol(m))) {
			index <- ((i-1) * num_of_subspaces_per_cardinality) + j
			z[[index]] <- m[,j]
		}
	}
	
	return (z)
}


subspace_sampling <- function(cardinality, total_attributes, num_subspaces) {
	max_subs <- factorial(total_attributes) / (factorial(total_attributes-cardinality) * factorial(cardinality))

	if (num_subspaces > max_subs) {
		stop("Requested number of subspaces exceeds the permutation limit!")
	}

	m <- matrix(cbind(sort(sample(total_attributes, cardinality))))
	
	while (ncol(m) < num_subspaces) {
		vec <- sort(sample(total_attributes, cardinality))
		
		val <- FALSE
		
		for (i in ncol(m)) {
			m_vec <- m[,i]
			
			if (identical(vec, m_vec)) {
				break;
			}
		}
		
		if (val == FALSE) {
			m <- cbind(m, vec)
		}
	}
	return (m)
}

indices_of_subspaces_containing_dimensions <- function(dimensions, subspaces) {
	indices <- vector()
	for (i in seq(1, length(subspaces))) {
		if(all(dimensions %in% subspaces[[i]])) {
			indices <- append(indices, i)
		}
	}
	return (indices)
}

assign_class_to_indices <-function(subspaces, indices, class) {
	ind_temp <- subspaces$z
	
	for (i in seq(1, length(indices))) {
		ind_temp[indices[i]] <- class
	}

	return (ind_temp)
}

dimension_count_for_cluster_indices <- function(subspaces, cluster_indices) {
	dimensions <- vector()
	for (i in seq(1, length(cluster_indices))) {
		dimensions <- append(dimensions, subspaces[[cluster_indices[i]]])
	}
	return (dimensions)
}

fix_entropy_values <- function(x) {
	val <- ifelse (x < 0.5,   x+0.5, x )
	return(val)
}

color_vector <- function(length, red_indices) {
	vec <- rep("black", length)
	for (i in seq(1, length(red_indices))) {
		vec[red_indices[i]] <- 'red'
	}
	return (vec)
}