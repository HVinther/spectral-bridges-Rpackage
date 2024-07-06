# spectralBridges Rpackage
 Spectral Bridges Clustering R package  implements a novel clustering algorithm combining k-means and spectral clustering techniques. It leverages efficient affinity matrix computation and merges clusters based on a connectivity measure inspired by SVM's margin concept. This package is designed to provide robust clustering solutions, particularly suited for large datasets. 

## Features

-   **Spectral Bridges Algorithm**: Integrates k-means and spectral clustering with efficient affinity matrix calculation for improved clustering results.
-   **Scalability**: Designed to handle large datasets by optimizing cluster formation through advanced affinity matrix computations.
-   **Customizable**: Parameters such as number of clusters, iterations, and random state allow flexibility in clustering configurations.

## Installation

To install the package 

```{r}
devtools::install_github("https://github.com/cambroise/spectral-bridges-Rpackage/")
library(spectralBridges)
```

## Usage

```{r}
X<-iris[,1:4] # Data to cluster
True_classes<-iris$Species
res<-spectral_bridges(X,n_cells=12,n_classes=3) # Partition in 3 classes
table(True_classes,Est_classes=res$clustering)  
```



## Documentation

For detailed documentation, please refer to the package vignette.

## Python version

Spectral Bridges Clustering algorithm is also available in python <https://pypi.org/project/spectral-bridges/>

