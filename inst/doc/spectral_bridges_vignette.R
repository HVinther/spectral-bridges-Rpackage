## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(spectralBridges)
library(factoextra)

X<-iris[,1:4]
True_classes<-iris$Species

res<-spectral_bridges(X,n_classes=3,n_cells=12)
fviz_cluster(res)             
knitr::kable(table(res$cluster,True_classes))

## -----------------------------------------------------------------------------
# Load necessary libraries
library(ClusterR)
library(factoextra)

# Sample data
set.seed(123)
X <- iris[,1:4]
True_classes=iris$Species
# Perform K-means clustering
n_cells <- 12
kmeans_result <- KMeans_rcpp(X, clusters = n_cells, num_init = 3, max_iters = 30, initializer = 'kmeans++')

# Extract cluster centers, labels, and sizes
kmeans_centers <- as.matrix(kmeans_result$centroids)
kmeans_labels <- kmeans_result$clusters
kmeans_size <- kmeans_result$obs_per_cluster

## -----------------------------------------------------------------------------
# Center the data points within each cluster
X.centered <- as.matrix(do.call(rbind, lapply(1:nrow(X), function(i) {
  X[i, ] - kmeans_centers[kmeans_labels[i], ]
})))

# Pre-compute distances between cluster centers
dist_centers <- as.matrix(dist(kmeans_centers))

# Define a function to compute affinity for one center
compute_affinity <- function(k) {
  affinity_row <- numeric(n_cells)
  for (l in 1:n_cells) {
    if (k != l) {
      distkl2 <- dist_centers[k, l]^2
      centered_k <- X.centered[kmeans_labels == k, ]
      centered_l <- X.centered[kmeans_labels == l, ]

      alpha_kl <- pmax(0, (kmeans_centers[l, ] - kmeans_centers[k, ]) %*% t(centered_k)) / distkl2
      alpha_lk <- pmax(0, (kmeans_centers[k, ] - kmeans_centers[l, ]) %*% t(centered_l)) / distkl2

      alphai <- c(alpha_kl, alpha_lk)
      affinity_row[l] <- sqrt(sum(alphai^2) / (kmeans_size[k] + kmeans_size[l]))
    }
  }
  return(affinity_row)
}

# Compute affinity for all centers using a for loop
affinity <- matrix(0, n_cells, n_cells)
for (k in 1:n_cells) {
  affinity[k, ] <- compute_affinity(k)
}

## -----------------------------------------------------------------------------
transform <- "exp"
M <- 1e3

if (transform == "exp") {
  gamma <- log(M) / diff(quantile(affinity, c(0.1, 0.9)))
  affinity <- exp(gamma * (affinity - 0.5 * max(affinity)))
}

## -----------------------------------------------------------------------------
# Normalized Laplacian matrix
D_inv_sqrt <- 1 / sqrt(rowSums(affinity))
L <- diag(n_cells) - t(affinity * D_inv_sqrt) * D_inv_sqrt
eigen.res <- eigen(-L, symmetric = TRUE)

# Determine the number of classes using the kneedle method
library(kneedle)
n_classes<-3
if (is.null(n_classes)) {
  n_classes <- kneedle(x = 1:length(eigen.res$values), y = eigen.res$values)[1] - 1
}
plot(eigen.res$values)
eigvecs <- eigen.res$vectors[, 1:n_classes]
eigvecs <- eigvecs / sqrt(rowSums(eigvecs ^ 2))
labels <- kmeans(eigvecs, nstart = 20, centers = n_classes)$cluster

## -----------------------------------------------------------------------------
# Assign labels based on clustering results
clusters <- labels[kmeans_labels]

# Return result
result <- list(clustering = clusters, data = X, class = "partition")
knitr::kable(table(Est_classes=result$clustering,True_classes=iris$Species))

