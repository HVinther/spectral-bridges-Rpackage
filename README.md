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



## Documentation and help

For detailed documentation, please refer to the package vignette:
```{r}
vignette("spectral_bridges_vignette")
help(spectral_bridges)
```

## Python version

Spectral Bridges Clustering algorithm is also available in python <https://pypi.org/project/spectral-bridges/>


## Principles

The proposed algorithm uses k-means centroids for vector quantization defining Voronoi region, and a strategy is proposed to link these regions, with an "affinity" gauged in terms of minimal margin between pairs of classes. These affinities are considered as weight of edges defining a completely connected graph whose vertices are the regions. Spectral clustering on the region provide a partition of the input space. The sole parameters of the algorithm are the number of Voronoi region and the number of final cluster. 

### Bridge affinity 

The basic idea involves calculating the difference in inertia achieved by projecting onto a segment connecting two centroids, rather than using the two centroids separately (See Figure @fig-balls-bridge). If the difference is small, it suggests a low density between the classes. Conversely, if this diffrence is large, it indicates that the two classes may reside within the same densely populated region.


<figure>
    <figure class="left" style="float:left">
    <img  src="vignettes/figures/balls.png" width="45%" />
    </figure>
    <figure class="right" style="float:right">
    <img src="vignettes/figures/bridge.png" width="45%" />
    </figure>
 <p>Balls (left) versus Bridge (right). The inertia of each structure is the sum of the squared distances represented by grey lines.</p>
</figure>




Let us consider a sample $X=(\boldsymbol x_i)_{i \in \{1,\cdots,n\}}$ of vectors $\boldsymbol x_i \in \mathbb R^p$ and a set of $m$ coding vectors $(\boldsymbol \mu_k)_{k \in \{1,\cdots,m\}}$ defining a partition $P=\{\mathcal{V}_1,\cdots,\mathcal{V}_m \}$ of $\mathbb R^p$ into $m$ Voronoi regions:

$$
\mathcal{V}_k = \left\{ \boldsymbol{x} \in \mathbb{R}^n \mid \|\boldsymbol{x} - \boldsymbol{\mu}_k\| \leq \|\boldsymbol{x} - \boldsymbol{\mu}_j\| \text{ for all } j \neq k \right\}.
$$

In the following a ball denotes the subset of $X$ in a Voronoi region.  The inertia of two balls $\mathcal{V}_k$ and $\mathcal{V}_l$ is 

$$
I_{kl} = \sum_{\boldsymbol x_i\in \mathcal{V}_k} \|\boldsymbol x_i - \boldsymbol \mu_k\|^2  + \sum_{\boldsymbol x_i\in \mathcal{V}_l} \|\boldsymbol x_i - \boldsymbol \mu_l\|^2.
$$ 
We define a bridge as a structure defined by a segment connecting two centroids $\boldsymbol \mu_k$ and $\boldsymbol \mu_l$. The inertia of a bridge between $\mathcal{V}_k$ and $\mathcal{V}_l$ is defined as 

$$
B_{kl} = \sum_{\boldsymbol x_i\in \mathcal{V}_k \cup \mathcal{V}_l} \|\boldsymbol x_i - \boldsymbol p_{kl}(\boldsymbol x_i)\|^2,
$$ where 
$$
\boldsymbol p_{kl}(\boldsymbol x_i) = \boldsymbol \mu_{k} + t_i(\boldsymbol \mu_{l} - \boldsymbol \mu_{k}),

$$ with $$
t_i  = \min\left(1, \max\left(0, \frac{\langle \boldsymbol x_i - \boldsymbol \mu_k | \boldsymbol \mu_l - \boldsymbol \mu_k\rangle}{\|  \boldsymbol \mu_l - \boldsymbol \mu_k \|^2}\right)\right). 
$$ 

Considering two centroïds, the normalized average of the difference betweenn Bridge and balls inertia   constitutes the basis of our affinity measure between to regions:

$$
\begin{aligned}
\frac{B_{kl}- I_{kl}}{(n_k+n_l)\|\boldsymbol \mu_k - \boldsymbol \mu_l\|^2} &=& \frac{\sum_{\boldsymbol{x_i} \in \mathcal V_k} \langle \boldsymbol{x_i} - \boldsymbol{\mu}_k \vert \boldsymbol{\mu}_l - \boldsymbol{\mu}_k \rangle_+^2  \sum_{\boldsymbol{x_i} \in \mathcal V_l} \langle \boldsymbol{x_i} - \boldsymbol{\mu}_l \vert \boldsymbol{\mu}_k - \boldsymbol{\mu}_l\rangle_+^2}{(n_k+n_l)\|\boldsymbol \mu_k - \boldsymbol \mu_l\|^4},\\
&=&  \frac{\sum_{\boldsymbol{x_i} \in \mathcal V_k \cup \mathcal V_l} \alpha_i^2}{n_k+n_l},
\end{aligned}
$$

where 

$$
\alpha_i=
\begin{cases}
t_i, &  \text{ if } t_i\in[0,1/2],\\
1-t_i, & \text{ if } t_i\in]1/2,1].
\end{cases}
$$ 

The basic intuition behind this affinity is that $t_i$ represents the
relative position of the projection of $\boldsymbol x_i$ on the segment
$[\boldsymbol \mu_k,\boldsymbol \mu_l]$. $\alpha_i$  represents the relative position on the segment, 
with the centroid of the class to which $\boldsymbol x_i$ belongs as the starting point.

The boundary that separates the two clusters defined by centroids
$\boldsymbol \mu_k$ and $\boldsymbol \mu_l$ is a hyperplane. This
hyperplane is orthogonal to the line segment connecting the centroids
and intersects this segment at its midpoint.

If we consider all points $\boldsymbol x_i \in \mathcal V_k \cup \mathcal V_l$ which are not projected on centroids but somewhere on the segment, the distance from
a point to the hyperplane is 
$$
\|\boldsymbol p_{kl}(\boldsymbol x_i) - \boldsymbol \mu_{kl}\| = (1/2-\alpha_i) \| \boldsymbol \mu_k-\boldsymbol \mu_l \|.
$$

This distance is similar to the concept of margin in Support Vector Machine [@Cortes1995]. When the $\alpha_i$ values are small (close to zero since $\alpha_i\in [0,1/2]$), the margins to the hyperplane
are large, indicating a low density between the classes. Conversely, if
the margins are small, it suggests that the two classes may reside
within the same densely populated region. Consequently, the sum of the
$\alpha_i$ or $\alpha_i^2$ increases with the density of the region
between the classes. 

Note that the criterion is local and indicates the relative difference in densities between the balls and the bridge, rather than evaluating a global score for the densities of the structures.

Eventually, we define the bridge affinity between centroids $k$
and $l$ as: 

$$
a_{kl}=
\begin{cases}
0, & \text{ if } k=l,\\
 \frac{\sum_{\boldsymbol{x_i} \in \mathcal V_k \cup \mathcal V_l} \alpha_i^2}{n_k+n_l}, & \text{otherwise}.
\end{cases}
$$ 
To allow points with large margin to dominate and make the algorithm more robust to noise and outliers we consider the following exponential transformation:
$$
\tilde{a}_{kl} = g(a_{kl})=\exp(\gamma\sqrt{a_{kl}}).
$$

where $\gamma$ is a scaling factor. This factor is set to ensure a large enough separation between the final coefficients. This factor is determined by the equation:
$$
\gamma = \frac{log(𝑀)}{\sqrt{q_{90}} - \sqrt{q_{10}}}
$$

where $q_{10}$ and $q_{90}$ are respectively the 10th and 90th percentiles of the original affinity matrix and $M > 0$. Thus, since the transformation is order-preserving, the 90th percentile of the newly constructed matrix is $M$ times greater than the 10th percentile. By default, $M$ is arbitrarily set to a large value of $10^4$.

