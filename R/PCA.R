#' Centres and optionally standardizes a given numeric matrix X
#'
#' @param X a matrix
#' @param standardize a boolean declaring whether to standardize
#'
#' @return X a centred matrix
#' @export
#'
#' @examples X <- matrix(seq(1:9))
#' centre(X)
centre <- function(X, standardize = F) {
  if (standardize) {
    X <- apply(X, MARGIN  = 2, FUN = function(x) {x <- x - mean(x)
                                          return(x/sd(x)) })
  }

  else {
    X <- apply(X, MARGIN  = 2, FUN = function(x) x - mean(x))
  }
  return(X)
}


#' Title
#'
#' @param X a data matrix X
#' @param standardize boolean declaring whether to standardize data
#' @param threshold_var the desired percentage of variance that needs to be explained by the principal components
#'
#' @return An object of class "principal components" containing transformed data, principal components and relevant information
#' @export
#'
#' @examples data(iris)
#' pca_iris <- pca(iris[,-5], FALSE, 0.95)
pca <- function(X, standardize,  threshold_var) {
  Z <- centre(as.matrix(X), standardize)
  n <- ncol(X)

  S <- svd(Z)

  P <- S$v
  D <- S$d**2
  sds <- S$d/sqrt(max(1, nrow(X)- 1))

  Zstar = Z %*% P

  row.names(P) <- colnames(X)
  colnames(P) <- paste0("PC",1:n)

  var_explained <- cumsum(D)/sum(D)
  upto <- which(var_explained > threshold_var)[1]

  pca_dat <- data.frame(Zstar[,1:upto, drop = F])
  colnames(pca_dat) <- paste0("PC", 1:upto)

  component_weighting <- data.frame("Component" = 1:n, "Weight" = D/sum(D))

  pca <- list(components = pca_dat, initial_data = X, pca_directions = P, weights = component_weighting, sd = sds)

  class(pca) = "principal_components"

  return(pca)

}

#' Print for principal component object
#'
#' @param x An object of class "principal_components"
#' @param ... generic arguments
#'
#' @return Prints basic information regarding the principal components
#' @rdname print.principal_components
#' @export
#'
#' @examples pca_iris <- pca(iris[,-5], FALSE, 0.95)
#' print(pca_iris)
print.principal_components <- function(x, ...) {
  cat(sprintf("Standard deviations:\n"))
  print(x$sd)
  cat(sprintf("\n PC rotations \n"))
  print(x$pca_directions)
}


#' Summary for principal component object
#'
#' @param object An object of class "principal_components"
#' @param ... generic arguments
#'
#' @return A summary providing the weights of each of the components
#' @rdname summary.principal_components
#' @export
#'
#' @examples pca_iris <- pca(iris[,-5], FALSE, 0.95)
#' summary(pca_iris)
summary.principal_components <- function(object, ...) {
  ob <- object
  ob$weights[,"Cumulative Weight"] <- cumsum(ob$weights[,"Weight"])
  print(ob$weights)
}

#' Biplot for principal components
#'
#' @param x An object of class "principal_components"
#' @param ... generic arguments
#' @param choices Indices declaring which of the PC's you want to plot
#' @param vectors A boolean declaring whether to include vectors representing the constituents of the principal components
#' @param vector_names A boolean declaring whether to label the vectors
#'
#' @return A biplot of the principal components equivalent to that described by Gower and Hand (1996)
#' @rdname plot.principal_components
#' @export
#'
#' @import graphics
#' @examples pca_iris <- pca(iris[,-5], FALSE, 0.95)
#' plot(pca_iris, choices = c(1, 2),  vectors = TRUE, vector_names = TRUE)
plot.principal_components <- function(x, ..., choices = c(1, 2), vectors = TRUE, vector_names = TRUE) {
  pca_ob <- x
  i <- choices[1]
  j <- choices[2]

  # Force a square plot so that angles between vectors are better visually interpretable
  opar <- par(no.readonly = TRUE)
  par(pty = "s")
  # Plot PCs against each other
  plot(pca_ob$components[,choices],
       xlim = c(-max(abs(pca_ob$components[,i])), max(abs(pca_ob$components[,i]))),
       ylim = c(-max(abs(pca_ob$components[,j])), max(abs(pca_ob$components[,j]))),
       xlab = paste0("Principal component ", i, " (", round(pca_ob$weights[i,2], 4)*100, "%)"),
       ylab = paste0("Principal component ", j, " (", round(pca_ob$weights[j,2], 4)*100, "%)"),
  )

  title(main = "Biplot for Principal Components", line = 3)

  if (vectors) {
    par(new = TRUE, las = 1)

    plot.window(xlim = c(-1, 1),
                ylim = c(-1, 1),
                asp = 1)

    # Add new axes for the vector plot

    axis(side = 3,
         at = c(-1, 0.5, 0, -0.5, 1),
         labels = TRUE,
         col = "navy",
         col.ticks = NULL,
         lwd = 2,
         col.axis = "navy")

    axis(side = 4,
         at = c(-1, 0.5, 0, -0.5, 1),
         labels = TRUE,
         col = "navy",
         col.ticks = NULL,
         lwd = 2,
         col.axis = "navy")

    # Adding labels for second axis

    mtext((text = paste("PC" ,i, "rotations")),
          side = 3,
          cex = 1,
          font = 2,
          family = "sans",
          col = "gray10",
          line = 2)

    mtext((text = paste("PC" ,j, "rotations")),
          side = 4,
          cex = 1,
          font = 2,
          family = "sans",
          col = "gray10",
          line = 2,
          las = 3)

    # Add straight lines to mark quadrants

    abline(v = 0, h = 0, lty = 2, col = "grey25")

    # Add variable vectors or arrows to the plot

    arrows(x0 = 0, x1 = pca_ob$pca_directions[,i],
           y0 = 0, y1 = pca_ob$pca_directions[,j],
           col = "navy",
           length = 0.1,
           lwd = 2,
           angle = 30)
    if (vector_names) {
      text(x = pca_ob$pca_directions[,i], y = pca_ob$pca_directions[,j],
           labels = row.names(pca_ob$pca_directions),
           cex = 0.8,
           font = 2,
           pos = c(rep(1, nrow(pca_ob$pca_directions))))
    }
  }
  par(opar)
}

#' Screeplot for principal_components object
#'
#' @param x An object of class "principal_components"
#' @param ... generic arguments
#' @param threshold A threshold that you desire the components to cross, will be plotted if value given
#'
#' @return A screeplot showing the weighting and cumulative weighting of the components
#' @rdname screeplot.principal_components
#' @export
#'
#' @import  stats
#'
#' @examples pca_iris <- pca(iris[,-5], FALSE, 0.95)
#' screeplot(pca_iris)
screeplot.principal_components <- function(x, ..., threshold = -1) {
  pca <- x
  x <- pca$weights[,1]
  y <- pca$weights[,2]
  opar <- par(no.readonly = T) #stores current par settings
  par(mar = c(5, 4, 4, 10),                                  # Specify par parameters
      xpd = TRUE) # changes par settings
  plot(x, y, ylim = c(0, 1), xlim = c(1, max(x)), xlab = "Principal Component", ylab = "Weighting",
       pch = 16, col = "red", xaxt = "n")
  axis(side = 1, at = pca$weights[,1], labels = T) # custom x-axis for dealing with integer values
  lines(y~x, col = 'red', lwd = 3)
  points(x, cumsum(y), pch = 16)
  lines(cumsum(y) ~ x, lwd = 3)
  lines(x = c(0.9, length((x) + 0.1)), y = c(threshold, threshold), lty = "dashed", lwd = 1)
  title(main = "Principal component weighting")
  legend("topright", inset=c(-0.4,0), legend=c("Weights","Cumulative Weights"), fill = c("red", "black"), xpd = TRUE)
  par(opar) # restores par settings to original
}



