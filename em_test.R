n_genes <- 90
n_cells <- 10000
n_classes <- 10

true_means <- exp( matrix( rnorm( n_genes*n_classes, 1, 2 ), nrow=n_genes ) )
true_sizes <- exp( rnorm( n_cells, mean=2.7, sd=2 ) )
true_class_probs <- rje::rdirichlet( 1, rep( 1, n_classes ) )[1,]
true_classes <- sample( 1:n_classes, n_cells, replace = TRUE, prob = true_class_probs )

counts <- 
  sapply( 1:n_cells, function(j)
     rpois( n_genes, true_sizes[j] * true_means[ , true_classes[j] ] ) )

totals <- colSums(counts)

# init
n_classes <- 20
class_probs <- rje::rdirichlet( 1, rep( 100, n_classes ) )[1,]
#class_probs <- rep( 1/n_classes, n_classes )
class_means <- matrix( exp( rnorm( n_classes*n_genes, 1, 2 ) ), nrow=n_genes ) * mean( 1/totals )
#class_means <- true_means * median( true_sizes / totals ) * exp(rnorm(n_classes,sd=1))


# E step
membership_probs_old <- membership_probs
membership_probs <-
  sapply( 1:n_cells, function(cell) {
    a <- colSums( dpois( counts[,cell], totals[cell] * class_means, log=TRUE ) ) + log(class_probs) 
    a <- exp( a - max( a, na.rm=TRUE ) )
    a[is.nan(a)] <- 0
    a <- a / sum( a ) } )
all(is.finite(membership_probs))

table( true_classes, apply( membership_probs, 2, which.max ) )
weighted.mean( log(membership_probs), membership_probs )

# M step
class_probs <- rowMeans( membership_probs )
class_means <- t( t( counts %*% t(membership_probs) ) / ( totals %*% t(membership_probs) )[1,] )

