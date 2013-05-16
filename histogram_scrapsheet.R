print (paste("Max cluster ID:", max(d$cluster)))

# extract cluster indices: x1-n 
for (i in seq(1, max(d$cluster))) {
 	indices <- which(d$cluster == i)
	assign(paste('x', i, sep=''), indices)
}

# assign classes to cluster indices
for (i in seq(1, max(d$cluster))) {
	indices <- get(paste('x', i, sep=''))
	for (j in seq(1, length(indices))) {
		dataf$z[indices[j]] <- i
	}
}

dataf$z <- as.character(dataf$z)

ggplot(dataf, aes(x,y,colour=z)) + coord_fixed() + geom_point(alpha=.8, size=5) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank()) + labs(title = paste(filename_base))

# ==== Histogram erstellen =============

dimensions <- vector()
for (i in seq(1, length(x3))) {
	dimensions <- append(dimensions, z[[x3[i]]])
}

h <-hist(dimensions, breaks=seq(0,ncol(pl.data)), col=color_vector(ncol(pl.data), z[[153]]))
sort(h$counts, index.return=TRUE, decreasing=TRUE)$ix

# sort subspace by max-dimensions
histsort <- sort(h$counts, index.return=TRUE, decreasing=TRUE)$ix
intersect(histsort, z[[150]])

# sort dimension names by sorted histogram values
names(pl.data)[sort(histcount, index.return=TRUE, decreasing=TRUE)$ix]
names(pl.data)[31]
# ==== ============= =============

ggplot(dataf, aes(x,y,colour=z)) + coord_fixed() + geom_point(alpha=.8, size=5) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.title = element_blank()) + labs(title = paste(filename_base)) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank())