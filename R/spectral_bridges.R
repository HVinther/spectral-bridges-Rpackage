#' @importFrom stats dist kmeans quantile
#' @import ggplot2
#' @import kneedle
#' @import ClusterR
#' @import factoextra
NULL


#' Spectral Bridges Clustering
#'
#' Performs spectral bridges clustering on the given dataset.
#'
#' @param X A numeric matrix of data to be clustered.
#' @param n_classes Number of classes to cluster into. If NULL, the number of classes will be determined automatically.
#' @param n_cells Number of cells of Voronoi tesselation. If NULL, heuristic is used
#' @param M A parameter for the transformation.
#' @param transform Type of transformation to apply.
#'
#' @return A list containing the clustering result and the original data.
#' @export
#'
#' @examples
#' \dontrun{
#' library(spectralBridges)
#' X<-iris[,1:4]
#' True_classes<-iris$Species
#' res<-spectral_bridges(X,n_cells=12,n_classes=3)
#' table(True_classes,Est_classes=res$clustering)
#' }
spectral_bridges <-  function(X,
                              n_classes = NULL,
                              n_cells = NULL,
                              M = 1e3,transform="exp"){
  # Input validation
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a numeric matrix or data frame.")
  }

  if (!is.null(n_classes) && (!is.numeric(n_classes) || length(n_classes) != 1 || n_classes <= 0)) {
    stop("n_classes must be a positive integer or NULL.")
  }

  if (!is.null(n_cells) && (!is.numeric(n_cells) || length(n_cells) != 1 || n_cells <= 0)) {
    stop("n_cells must be a positive integer or NULL.")
  }

  if (!is.null(n_cells) && (n_cells<=n_classes)) {
    stop("Number of cells should be greater than number of clusters.")
  }

  if (!is.numeric(M) || length(M) != 1 || M <= 0) {
    stop("M must be a positive number.")
  }

  if (!transform %in% c("exp", "none")) {
    stop("transform must be one of 'exp' or 'none'.")
  }

  # Convert data frames to matrices
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  # 1. Vector quantization
  ###################################
  n <- nrow(X)
  if (is.null(n_cells)) n_cells<-n_classes*7

  kmeans_result  = KMeans_rcpp(X, clusters = n_cells,
                               num_init = 3, max_iters =30,
                               initializer = 'kmeans++')

  kmeans_centers <- as.matrix(kmeans_result$centroids)
  kmeans_labels <-  kmeans_result$clusters
  kmeans_size <-    kmeans_result$obs_per_cluster
  kmeans_Iw<-kmeans_result$WCSS_per_cluster


  # 2. Affinity computation
  ###################################
  # Centering of X
  X.centered <- as.matrix(do.call(rbind, lapply(1:n, function(i) {
    X[i, ] - kmeans_centers[kmeans_labels[i], ]
  })))

  # Pre-computation of distances between centers
  dist_centers <- as.matrix(dist(kmeans_centers))

  # Affinity
  affinity<-matrix(0,n_cells,n_cells)
  for (l in 1:(n_cells-1))
    for (k in (l+1):n_cells){
      distkl2 <- dist_centers[k, l]^2
      centered_k <- X.centered[kmeans_labels == k, ]
      centered_l <- X.centered[kmeans_labels == l, ]
      alpha_kl <- pmax(0, (kmeans_centers[l, ] - kmeans_centers[k, ]) %*% t(centered_k)) / distkl2
      alpha_lk <- pmax(0, (kmeans_centers[k, ] - kmeans_centers[l, ]) %*% t(centered_l)) / distkl2
      alphai <- c(alpha_kl, alpha_lk)
      affinity[l,k] <- sqrt(sum(alphai^2) / (kmeans_size[k] + kmeans_size[l]))
      affinity[k,l] <- affinity[l,k]
    }


    if (transform=="exp"){
    gamma<- log(M)/diff(quantile(affinity,c(0.1,0.9)))
    affinity<- exp(gamma*(affinity - 0.5*max(affinity)))}

  # 3. Spectral Clustering of the coding vectors
  ###################################
  # Normalized Laplacian matrix
  D_inv_sqrt <- 1 / sqrt(rowSums(affinity))
  L <- diag(n_cells) - t(affinity * D_inv_sqrt) * D_inv_sqrt
  eigen.res<-eigen(-L, symmetric = TRUE)
  ifelse(is.null(n_classes),
         n_classes <- kneedle(x=1:length(eigen.res$values),y=eigen.res$values)[1]-1,
         n_classes <-n_classes)
  eigvecs <- eigen.res$vectors[,1:n_classes]
  eigvecs <- eigvecs / sqrt(rowSums(eigvecs ^ 2))
  labels <-
    kmeans(eigvecs, nstart = 20, centers = n_classes)$cluster

  # 4. Contaminate
  ###################################
  clusters <- labels[kmeans_labels]

  return(list(clustering = clusters,data =X,class="partition"))
}
