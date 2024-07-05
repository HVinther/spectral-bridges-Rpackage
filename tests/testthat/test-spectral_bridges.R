test_that("Testing Spectral Bridges Clustering algo", {
  data(iris)
  X<-iris[,1:4]
  spectralBridges::spectral_bridges(X)
})
