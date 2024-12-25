#' @import ClusterR
library(ClusterR)

#' Spectral Bridges Clustering
#'
#' Performs spectral bridges clustering on the given dataset.
#'
#' @param X A numeric matrix of data to be clustered.
#' @param n_classes Number of classes to cluster into. If NULL, the number of classes will be determined automatically.
#' @param n_cells Number of cells of Voronoi tessellation. If NULL, a heuristic is used.
#' @param p Power parameter for affinity computation.
#' @param M A parameter for the transformation.
#'
#' @return A list containing the clustering result and the original data.
#' @export
#'
#' @examples
#' \dontrun{
#' library(spectralBridges)
#' X <- iris[, 1:4]
#' True_classes <- iris$Species
#' res <- spectral_bridges(X, n_cells = 12, n_classes = 3)
#' table(True_classes, Est_classes = res$clustering)
#' }
#'
spectral_bridges <- function(X, n_classes = NULL, n_cells = NULL, p = 2, M = 1e4) {
  # Input validation
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a numeric matrix or data frame.")
  }

  if (!is.numeric(n_classes) || length(n_classes) != 1 || n_classes <= 0) {
    stop("n_classes must be a positive integer or NULL.")
  }

  if (!is.numeric(n_cells) || length(n_cells) != 1 || n_cells <= 0) {
    stop("n_cells must be a positive integer or NULL.")
  }

  if (!is.null(n_cells) && (n_cells <= n_classes)) {
    stop("Number of cells should be greater than the number of clusters.")
  }

  if (!is.numeric(p) || length(p) != 1 || p <= 0) {
    stop("p must be a strictly positive number.")
  }

  if (!is.numeric(M) || length(M) != 1 || M <= 0) {
    stop("M must be a positive number.")
  }

  # Convert data frames to matrices
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  # 1. Vector Quantization (K-Means)
  ###################################
  kmeans_result <- KMeans_rcpp(X, clusters = n_cells, num_init = 3, max_iters = 30, initializer = 'kmeans++')

  kmeans_centers <- as.matrix(kmeans_result$centroids)
  kmeans_labels <- kmeans_result$clusters

  # Center data by cluster centers
  X_centered <- do.call(rbind, lapply(1:n_cells, function(k) {
    sweep(X[kmeans_labels == k, ], 2, kmeans_centers[k, ])
  }))

  # 2. Affinity Computation
  ###################################
  affinity <- matrix(0, n_cells, n_cells)

  for (k in 1:n_cells) {
    segments <- sweep(kmeans_centers, 2, kmeans_centers[k, ])
    dist_kl2 <- rowSums(segments^2)
    dist_kl2[k] <- 1  # Avoid division by zero

    centered_k <- X_centered[kmeans_labels == k, , drop = FALSE]
    projs <- centered_k %*% t(segments)
    projs <- array(pmax(0, projs) / dist_kl2, dim = dim(projs))
    projs <- projs^p

    affinity[k, ] <- colSums(projs)
  }

  # Normalize affinity matrix
  kmeans_sizes <- table(kmeans_labels)
  counts <- outer(kmeans_sizes, kmeans_sizes, "+")
  affinity <- ((affinity + t(affinity)) / counts)^(1 / p)

  gamma <- log(M) / diff(quantile(affinity, c(0.1, 0.9)))
  affinity <- exp(gamma * (affinity - 0.5 * max(affinity)))

  # 3. Spectral Clustering
  ###################################
  D_inv_sqrt <- 1 / sqrt(rowSums(affinity))
  L <- diag(n_cells) - t(affinity * D_inv_sqrt) * D_inv_sqrt
  eigen_res <- eigen(-L, symmetric = TRUE)

  if (is.null(n_classes)) {
    n_classes <- kneedle(x = 1:length(eigen_res$values), y = eigen_res$values)[1] - 1
  }

  eigvecs <- eigen_res$vectors[, 1:n_classes]
  eigvecs <- eigvecs / sqrt(rowSums(eigvecs^2))
  labels <- kmeans(eigvecs, nstart = 20, centers = n_classes)$cluster

  # Map back to original data
  clusters <- labels[kmeans_labels]

  return(list(clustering = clusters, data = X, class = "partition"))
}
