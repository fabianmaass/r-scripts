# ////// Start:
# in Z die "subspaces" laden -> hier card: 10,20,30 mit jeweils 50 subspaces
	z    <- create_sampled_subspaces(50, c(10,20,30), ncol(pl.data))

# dann das script einmal laufen lassen um die simmilarity matrix zu berechnen

mvec <- vector("list", length(z))
subspaces <- vector("list", length(z))

################################################################
### MDS Projection der Subspaces berechnen

print("calculate MDS projection of subspaces")

for (i in seq(1, length(z))) {
	print (i)
	subspace <- as.matrix(subspace_from_indices(z[[i]]))
	
	mdstemp <- cmdscale(daisy(subspace, metric = c("euclidean"), stand = FALSE, type = list()), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
	mvec[[i]] <- mdstemp
}

################################################################
### Store all subspaces in a vector

for (i in seq(1, length(z))) {
	print (i)
	subspace <- as.matrix(subspace_from_indices(z[[i]]))
	subspaces[[i]] <- subspace
}
		
###################################################################################################################
### MDS mit data topology similarity

dt_sim <- matrix(0, nrow=length(z), ncol=length(z))

surfing_scores <- vector()
entropy_scores <- vector()

for (i in seq(1, length(z) - 1)) {
	print (paste("i: " , i))

	sub_p <- mvec[[i]]

	surfing        <- surfing_score(sub_p, 20)
	entropy        <- calc_entropy(sub_p, 30)
	surfing_scores <- append(surfing_scores, surfing)
	entropy_scores <- append(entropy_scores, entropy)

	for (j in seq(i+1, length(z))) {
		print (paste("j: " , j))

		sub_q <- mvec[[j]]
		
		# use VISUAL data topology measure (working directly on the MDS projection)
		dt <- data_topology(sub_p, sub_q, 0.1, 0.1)

		dt_sim[i,j] <- dt
		dt_sim[j,i] <- dt
	}
}
sub_p <- mvec[[length(z)]]

surfing        <- surfing_score(sub_p, 20)
entropy        <- calc_entropy(sub_p, 30)
surfing_scores <- append(surfing_scores, surfing)
entropy_scores <- append(entropy_scores, entropy)

mds <- cmdscale(as.dist(dt_sim), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)

dataf = as.data.frame(cbind(normalit(mds[,1]),normalit(mds[,2]), normalit(surfing_scores)))
names(dataf) = c("x","y","z")

#######################################################################################
#######################################################################################

##### jetzt hast du in:
### dt_sim: die simmilarity matrix
### mds:    die mds koordinaten
### dataf:  ebenfalls die mds koordinaten, nur als datenstruktur eine datatable, mit einer extra spalte

## dann lässt du auf dem mds ein dbscan laufen zum erkennen der cluster 
## durch showplot kannst du ja mit den epsilon-werten etwas rumspielen

d <- dbscan(mds, 0.05, MinPts = 3, showplot = 2)

## jetz hast du die dbscan ergebnisse in der variable d
## anschließend holen wir uns die cluster indizes und speichern sie in extra vectoren
## ist vielleicht überflüssig .. ich will aber nich in der d variable rummachen, das ist mir zu umständlich

# extract cluster indices: x1-n 
for (i in seq(1, max(d$cluster))) {
 	indices <- which(d$cluster == i)
	assign(paste('x', i, sep=''), indices)
}

## jetz weisen wir den indizes in der struktur dataf, die cluster indizes aus den x vektoren zu
## funktion ist in scuba_scripts.r enthalten

# assign classes to cluster indices
dataf$z <- assign_class_to_indices(dataf, x1, 1)
dataf$z <- assign_class_to_indices(dataf, x2, 2)
dataf$z <- assign_class_to_indices(dataf, x3, 3)
dataf$z <- assign_class_to_indices(dataf, x4, 4)

dataf$z <- as.character(dataf$z)

## so und jetz kannst du die dataf struktur plotten, und kriegst dazu die farb-codes aus dem dbscan

ggplot(dataf, aes(x,y,colour=z)) + coord_fixed() + geom_point(alpha=.8, size=5) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(),  axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.title = element_blank()) + labs(title = "random cluster injection") + theme(panel.background=element_rect(fill="transparent",colour="black"))

## und wenn du das alles hast, dann kannst du neue subspaces hinzufügen indem du dir die generieren lässt und an die z struktur unten anfügst:

z_temp    <- create_sampled_subspaces(50, c(10,20,30), ncol(pl.data))

z_length <- length(z)
for (i in seq(1, length(z_temp))) {
  z[[z_length+i]] <- z_temp[[i]]
}

## und dann kannst du entweder das ganze script hier ab zeile 7-72 nochmal laufen lassen ... 
## oder die simmilarity matrix nur um die zeilen und spalten der hinzugefügten subspaces erweitern ... 

## jetz noch die cluster indizes in dataf: zeile 101-106 (zeile ) und die neue indizes 
## ich bin hier jetzt mal davon ausgegangen das vorher 150 subspaces in dataf waren und die neuen von 151-200 sind...
dataf$z <- assign_class_to_indices(dataf,  seq(151,200), 5)

## und dann kannst du das wieder plotten...
ggplot(dataf, aes(x,y,colour=z)) + coord_fixed() + geom_point(alpha=.8, size=5) + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(),  axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), legend.title = element_blank()) + labs(title = "random cluster injection") + theme(panel.background=element_rect(fill="transparent",colour="black"))

### und so weiter ;)

## wenn du dann die indizes der gruppen willst, kannst du die funktion benutzen:
## damit bekommst du die un-aggregierten dimensions-indizes ... dh. du kannst da dann ein hist drauf laufen lassen

dimensions_x1 <- dimension_count_for_cluster_indices(z, x1)
h <-hist(dimensions_x1, breaks=seq(0,ncol(pl.data)))