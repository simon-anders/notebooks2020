n_genes <- 90
n_cells <- 3000
n_classes <- 10

true_means <- exp( matrix( rnorm( n_genes*n_classes, 1, 2 ), nrow=n_genes ) )
true_sizes <- rep( 1, n_cells ) #exp( rnorm( n_cells, mean=.7, sd=2 ) )
true_classes <- sample( 1:n_classes, n_cells, replace = TRUE )#, prob = c( .1, .2, .3, .39, .01 ) )

counts <- 
  sapply( 1:n_cells, function(j)
     rpois( n_genes, true_sizes[j] * true_means[ , true_classes[j] ] ) )

# init
n_classes <- 20
class_probs <- rje::rdirichlet( 1, rep( 100, n_classes ) )[1,]
class_means <- matrix( exp( rnorm( n_classes*n_genes, 1, 2 ) ), nrow=n_genes )

# E step
membership_probs_old <- membership_probs
membership_probs <-
  sapply( 1:n_cells, function(cell) {
    a <- colSums( dpois( counts[,cell], class_means, log=TRUE ) )
    a <- exp( a - max( a, na.rm=TRUE ) )
    a[is.nan(a)] <- 0
    a <- a / sum( a ) } )
all(is.finite(membership_probs))

table( true_classes, apply( membership_probs, 2, which.max ) )
weighted.mean( log(membership_probs), membership_probs )

# M step
class_probs <- rowMeans( membership_probs )
class_means <- t( t( counts %*% t(membership_probs) ) / rowSums(membership_probs) )

