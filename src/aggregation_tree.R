# to build an aggregation tree based on the similarity matrix
library("phangorn")
setwd("/Users/shaopeng/Desktop/mcibs_project")

### read raw data
df_otu <- read.csv("otutable.txt", sep="\t", row.names=1)
df_meta <- read.csv("task.txt", sep="\t")
# filter unmatched samples (some baseline control)
keep_otu <- intersect(colnames(df_otu), df_meta$X.SampleID)
df_otu <- df_otu[ ,keep_otu]
df_meta <- df_meta[df_meta$X.SampleID %in% keep_otu, c(1,2)]
rm(keep_otu)

### check numbers
species =  strsplit(row.names(df_otu), split="_")
species = sapply(species,'[[',4)
length(unique(species)) #752
genus = strsplit(row.names(df_otu), split="_")
genus = sapply(genus, '[[', 3)
length(unique(genus)) #384

###clean data (prepare into analysis-ready format)
row.names(df_meta) <- df_meta[,1]
df_meta <- df_meta[, -1, drop=F] 
colnames(df_meta) <- "label"
original_otu <- row.names(df_otu) 
fuso_otu <- original_otu[grep("Fusobacterium", original_otu)]
fuso_otu_id <- unname(sapply(fuso_otu, function(x){paste(unlist(strsplit(x, split="_"))[1:4], collapse = "_")}))
row.names(df_otu) <- sapply(row.names(df_otu), function(x){paste(unlist(strsplit(x, split="_"))[1:4], collapse = "_")})
otu_count <- t(df_otu)
otu_count <- merge(df_meta, otu_count, by=0)
row.names(otu_count) <- otu_count$Row.names
otu_count <- otu_count[, -1]
f.otu.count <- otu_count[rowSums(otu_count[, -1])>2000, ]


write.csv(f.otu.count, file="raw_input.csv")

### sample/feature exclusion
# 3.1 sample size (rowSum): mean 4501, median 4297, max 11339, min(237)
# indicating some samples are extremely sparse, can't be used
# arbitrary cutoff: >2000 (129 / 156 remaining)
par(mfrow=c(1,2))
plot(density(rowSums(otu_count[, -1])), main="density of library size", col="orange")
abline(v=2000, col="red")
boxplot(rowSums(otu_count[, -1]), main="Boxplot of library size", col="grey")
abline(h=2000, col="red")
f.otu.count <- otu_count[rowSums(otu_count[, -1])>2000, ]

# 3.2 remove low abundance species (biological reasoning)
par(mfrow=c(1,1))
plot(density(colSums(f.otu.count[,-1])))
quantile(colSums(f.otu.count[,-1]))
### so many variables appear very few times
# 0~1: not reliable, itu is not sensitive at low-abundance level
# 1~20: the appearance is too low to be biological informative (0.5% total popu)
#f.otu.count <- f.otu.count[, c(TRUE, colSums(f.otu.count[,-1])>20)]
f.otu.count <- f.otu.count[, c(TRUE, colSums(f.otu.count[,-1]>5)>5)]

write.csv(f.otu.count, file="cleaned_input.csv")

# corrplot
library(corrplot)
corrplot(cor(f.otu.count[,-1]), method="circle", title="Corrplot of features", tl.pos = "n", mar=c(0,0,1,0))

### check numbers
species =  strsplit(colnames(f.otu.count)[-1], split="_")
species = sapply(species,'[[',4)
length(unique(species)) #182
filtered_species = species
genus = strsplit(colnames(f.otu.count)[-1], split="_")
genus = sapply(genus, '[[', 3)
length(unique(genus)) #110


### generate similarity clustering
k_list= seq(10,60,5)
cluster_num = append(seq(5,40,5), seq(50,100,10))
#k = 30
matching_file = read.table("GCA_genus_species_matching.txt", col.names = c("GCA", "genus", "species"))
row.names(matching_file) = matching_file$GCA

plot_phylo <- function(phylo_obj, name='upgma.pdf'){
  pdf(name)
  opar=par(no.readonly=TRUE)
  par(mfrow=c(2, 2), col.main="red", family="serif")
  par(mai=c(0.2, 0.2, 0.2, 0.2))
  plot(phylo_obj, type="phylogram", main="phylogram", show.tip.label=F) 
  plot(phylo_obj, type="fan", main="fan", show.tip.label=F)
  plot(phylo_obj, type="unrooted", main="unrooted", show.tip.label=F)
  plot(phylo_obj, type="radial", main="radial", show.tip.label=F)
  par(opar)
  dev.off()
}

get_phylo_obj <- function(k){
  file=paste0("JI_matrix_of_input_k",k,".csv")
  df_temp <- read.csv(file)
  #change names to species name
  ids <- sapply(colnames(df_temp), function(x){gsub("X.data.sml6467.github.MCIBS556_20FA_project.src.NCBI_database.genomes.","",x)})
  ids =  strsplit(ids, split="_")
  ids = sapply(ids,'[[',2)
  df_temp_species= rep('none', length(ids))
  for ( i in seq(1:length(ids))) {
    match_gca = paste0("GCA_", ids[i])
    df_temp_species[i] = as.character(matching_file[match_gca, "species"])
  }
  if ('none' %in% df_temp_species) {
    print("There is unclassified species")
    print(paste0("The first apperance is at ",which(df_temp_species == "none")))
  }
  # species-level dedup: some parallel species still pass bash filter
  df_temp <- df_temp[!duplicated(df_temp_species),!duplicated(df_temp_species)]
  colnames(df_temp) <- df_temp_species[!duplicated(df_temp_species)]
  row.names(df_temp) <- df_temp_species[!duplicated(df_temp_species)]
  up = upgma(df_temp)
  plot_phylo(up, name=paste0("Phylo_plot_k",k,".pdf"))
  out_list <- list('up' =up, 'ji'=df_temp)
  return(out_list)
}

get_clade <- function(k, m_list){
  # species list in the otu table
  temp_species =  strsplit(colnames(f.otu.count)[-1], split="_")
  otu_genus = sapply(temp_species,'[[',3)
  otu_species = sapply(temp_species,'[[',4)
  phylo_out = get_phylo_obj(k)
  phylo_obj = phylo_out$up

  for (m in m_list) {
    temp_hclust = as.hclust.phylo(phylo_obj)
    temp_clusters = cutree(temp_hclust, k=m)
    names_clusters = names(temp_clusters)
    # match cluster info into otu table
    otu_group=rep('none', length(otu_species))
    
    for (i in seq(1:length(otu_species))) {
      if ( otu_species[i] %in% names_clusters ) {
        otu_group[i] <- as.character(temp_clusters[which(names_clusters == otu_species[i])])
      } 
      else {
        otu_group[i] <- otu_species[i]
      }
    }
    
    modified_otu = f.otu.count[, -1]
    out_otu = data.frame(t(rowsum(t(modified_otu), otu_group)))
    out_otu$label <- f.otu.count$label
    temp_size = dim(out_otu)[2]
    out_otu <- out_otu[,c(temp_size, 1:temp_size-1)]
    write.csv(out_otu, file=paste0("aggregated_k",k,"_m",m,".csv"))
    id_otu <- cbind(colnames(modified_otu), otu_group)
    write.csv(id_otu, file=paste0("group_infor_aggregated_k",k,"_m",m,".csv"))
    jpeg(file=paste0("Group_size_k",k, "_m",m, ".png"))
    barplot(table(temp_clusters), main=paste0("Group_size_k",k, "_m",m), xlab = "group_id", ylab="number of species")
    dev.off()
  }
}

for (k in k_list) {
  get_clade(k, m_list=cluster_num)
}



### generate a genus-level aggregation file
genus_otu <- f.otu.count[,-1]
temp_genus =  strsplit(colnames(genus_otu), split="_")
temp_genus = sapply(temp_genus,'[[',3)
out_otu = data.frame(t(rowsum(t(genus_otu), temp_genus)))  #110 features
out_otu$label <- f.otu.count$label
temp_size = dim(out_otu)[2]
out_otu <- out_otu[,c(temp_size, 1:temp_size-1)]
write.csv(out_otu, file="genus_level_otu_table.csv")
jpeg(file="genus_level_aggregation.png")
barplot(table(temp_genus), main="genus level aggregation", xlab = "genus", ylab="number of species")
dev.off()



### because RF is the best model, there might be no high order relationships
# within group correlation
cor_within_cluster <- function(group_infor, name="input data") {
  corr_vec=vector()
  input_data=f.otu.count[,-1]
  for (sub_grp in group_infor) {
    sub_input <- input_data[ , group_infor == sub_grp, drop=FALSE]
    if (dim(sub_input)[2] > 1) {
      cor_max <- cor(sub_input)
      corr_vec <- append(corr_vec, cor_max[upper.tri(cor_max)])
    }
  }
  
  #plot the results
  jpeg(file=paste0("Density_plot_of_CORR_from_", name, ".jpeg"))
  d <- density(corr_vec)
  plot(d, main=paste0("Density plot of CORR from ", name), xlab = "Within group pairwise correlation")
  polygon(d, col="red", border="blue")
  dev.off()
  # return obj 
  return(corr_vec)
}
# genus cor
cor1 <- cor_within_cluster(temp_genus, name="genus_aggregation")
# km combination
cor_for_km <- function(k, m) {
  phylo_out = get_phylo_obj(k)
  phylo_obj = phylo_out$up
  # aggregation cluster
  temp_hclust = as.hclust.phylo(phylo_obj)
  temp_clusters = cutree(temp_hclust, k=m)
  names_clusters = names(temp_clusters)
  # match to colname
  temp_species =  strsplit(colnames(f.otu.count)[-1], split="_")
  otu_species = sapply(temp_species,'[[',4)
  otu_group=rep('none', length(otu_species))
  for (i in seq(1:length(otu_species))) {
    if ( otu_species[i] %in% names_clusters ) {
      otu_group[i] <- as.character(temp_clusters[which(names_clusters == otu_species[i])])
    } 
    else {
      otu_group[i] <- otu_species[i]
    }
  }
  # get cor plot
  out_cor <- cor_within_cluster(otu_group, name=paste0("similarity_k", k, "_m", m,"_aggregation"))
  return(out_cor)
}
cor2 <- cor_for_km(20, 30)
cor3 <- cor_for_km(40, 35)
cor4 <- cor_for_km(60, 20)

# plot 2 lines together
compare_2_cors <- function(c1, c2, name1, name2) {
  jpeg(file=paste0("Compare_2_corrplot_from_", name1, "_", name2, ".jpeg"))
  d1 <- density(c1)
  d2 <- density(c2)
  plot(d1, col=2, ylim=c(0,max(max(d1$y), max(d2$y))+1), main="Compare 2 corrplot", xlab="Within group pairwise correlation")
  lines(d2, col=4)
  legend("topright", legend=c(paste0(name1, ", N=", length(c1)), paste0(name2, ", N=", length(c2))), fill=c(2,4))
  dev.off()
}
compare_2_cors(cor1, cor2, "genus", "k20m30")
compare_2_cors(cor1, cor3, "genus", "k40m35")
compare_2_cors(cor1, cor4, "genus", "k60m20")




### compare similarity tree with phylo tree
compare_2_trees <- function(k, phylo_tree_file="species_in_ji_matrix.nwk") {
  phylo_out = get_phylo_obj(k)
  simi_tree = as.phylo(phylo_out$up)
  # load phylo tree
  phylo_tree = read.tree(file=phylo_tree_file)
  # prune simi tree for non existing one in phylo tree
  to_prune=vector()
  phylo_species = strsplit(phylo_tree$tip.label, split="_")
  phylo_species = sapply(phylo_species,'[[',2)
  phylo_tree$tip.label <- phylo_species
  for (tree_label in simi_tree$tip.label) {
    if (!tree_label %in% phylo_species) {
      to_prune = append(to_prune, tree_label)
    }
  }
  simi_tree <- drop.tip(simi_tree, to_prune)
  # compare
  compare_out <- comparePhylo(simi_tree, phylo_tree, plot=T, use.edge.length = T)
}



pdf(name)
opar=par(no.readonly=TRUE)
par(mfrow=c(2, 2), col.main="red", family="serif")
par(mai=c(0.2, 0.2, 0.2, 0.2))
plot(phylo_obj, type="phylogram", main="phylogram", show.tip.label=F) 
plot(phylo_obj, type="fan", main="fan", show.tip.label=F)
plot(phylo_obj, type="unrooted", main="unrooted", show.tip.label=F)
plot(phylo_obj, type="radial", main="radial", show.tip.label=F)
par(opar)
dev.off()




################ manual code below
# check # of species in filtered data with complete genomes
length(unique(filtered_species)) #182
all_found = simi_tree$tip.label
length(intersect(filtered_species, all_found))

# look at group 1 in k35m40 output
k=5
m=40
phylo_out = get_phylo_obj(k)
phylo_obj = phylo_out$up
temp_hclust = as.hclust.phylo(phylo_obj)
temp_clusters = cutree(temp_hclust, k=m)
gp1_species <- names(temp_clusters[temp_clusters==1])
# JI matrix
refine_ji = phylo_out$ji
gp1_ji <- refine_ji[gp1_species, gp1_species]
library(corrplot)
corrplot(data.matrix(gp1_ji), method="circle", title="Corrplot of features", tl.pos = "n", mar=c(0,0,1,0))
# conclusion: all group1 should remain untouched

gp2_species <- names(temp_clusters[temp_clusters==2])
gp2_ji <- refine_ji[gp2_species, gp2_species]
corrplot(data.matrix(gp2_ji), method="circle", title="Corrplot of features", tl.pos = "n", mar=c(0,0,1,0))

plot_group <- function(grp_num, grp_vector=temp_clusters, input_ji=refine_ji) {
  gp_species <- names(grp_vector[grp_vector==grp_num])
  gp_ji <- input_ji[gp_species, gp_species]
  corrplot(data.matrix(gp_ji), method="circle", title=paste0("Corrplot of features of group ", grp_num), tl.pos = "n", mar=c(0,0,1,0))
}

plot_group(1)



